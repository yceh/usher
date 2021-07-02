#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include "src/new_tree_rearrangements/tree_rearrangement_internal.hpp"
#include "zlib.h"
#include "tbb/concurrent_queue.h"
#include "tbb/flow_graph.h"
#include <atomic>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <mutex>
#include <string>
#include <tbb/concurrent_vector.h>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "import_vcf.hpp"
#include "tbb/parallel_for.h"
#define ZLIB_BUFSIZ 0x10000
typedef tbb::flow::multifunction_node<char*,tbb::flow::tuple<char*>> decompressor_node_t;
typedef tbb::flow::multifunction_node<char*,tbb::flow::tuple<char*,Parsed_VCF_Line*>> line_parser_t;
typedef tbb::flow::multifunction_node<char*,tbb::flow::tuple<char*>> buffer_only_out;
//Decouple parsing (slow) and decompression, segment file into blocks for parallelized parsing
struct Decompressor{
    gzFile* fd;
    size_t cont_read_size;
    size_t init_read_size;
    void operator()(char* buf,decompressor_node_t::output_ports_type& out) const{
        if (gzeof(*fd)) {
            free(buf);
            return;
        }
        int read_size=gzread(*fd, buf, init_read_size);
        if (!read_size) {
            free(buf);
            return;
        }
        //Make sure the last line is complete in the block.
        if(!gzgets(*fd, buf+read_size, cont_read_size)){
            *(buf+read_size)=0;
        }
        std::get<0>(out).try_put(buf);
    }
};
typedef std::vector<tbb::concurrent_vector<MAT::Valid_Mutation>> new_sample_mut_t;
//Parse a block of lines, assuming there is a complete line in the line_in buffer
struct produce_exist_sample_vcf_line{
    typedef line_parser_t::output_ports_type output_type;
    Parsed_VCF_Line * parsed_line;
    produce_exist_sample_vcf_line(const std::string& chromosome,int pos,nuc_one_hot ref){
        parsed_line=new Parsed_VCF_Line{MAT::Mutation(chromosome,pos, 0, 0, 1,
                      ref)};
        parsed_line->mutated.emplace_back(0,0xf);
    }
    void new_sample_mut(long idx,nuc_one_hot mut_nuc){
        parsed_line->mutated.emplace_back(idx,mut_nuc);
    }
    void output(output_type& out){
        std::sort(parsed_line->mutated.begin(),parsed_line->mutated.end(),mutated_t_comparator());
        std::get<1>(out).try_put(parsed_line);
    }
};
struct ignore_existing_sample{
    typedef buffer_only_out::output_ports_type output_type;
    ignore_existing_sample(const std::string& chromosome,int pos,nuc_one_hot ref){}
    void new_sample_mut(long idx,nuc_one_hot mut_nuc){}
    void output(output_type& out){}
};
template<typename existing_sample_ctl>
struct line_parser{
    const std::vector<long>& header;
    new_sample_mut_t& new_samples;
    void put_sample(existing_sample_ctl& existing_sample_out,long idx,char allele,int pos,uint8_t chrom_idx,nuc_one_hot ref)const{
        if (idx>0) {
            existing_sample_out.new_sample_mut(idx,allele);
        }else{
            new_samples[-idx].emplace_back(pos,chrom_idx,ref,allele);
        }
    }
    void operator()(char* line_in, typename existing_sample_ctl::output_type& out)const{
        char* start=line_in;
        while (*line_in!=0) {
            std::vector<nuc_one_hot> allele_translated;
            std::string chromosome;
            int pos=0;
            //Chromosome
            while (*line_in!='\t') {
                chromosome.push_back(*line_in);
                line_in++;
            }
            line_in++;
            //Position
            while (*line_in!='\t') {
                pos=pos*10+(*line_in-'0');
                line_in++;
            }
            line_in++;
            //ID don't care
            while (*line_in!='\t') {
                line_in++;
            }
            line_in++;
            //REF
            auto ref_nuc = MAT::get_nuc_id(*line_in);
            existing_sample_ctl parsed_line(chromosome, pos, ref_nuc);
            line_in++;
            auto chrom_idx=MAT::Mutation::chromosome_map[chromosome];
            //assert(*line_in=='\t');
            line_in++;
            //ALT
            while (*line_in!='\t') {
                allele_translated.push_back(MAT::get_nuc_id(*line_in));
                line_in++;
                if(*line_in==','){
                    line_in++;
                }else{
                    //assert(*line_in=='\t');
                }
            }
            line_in++;
            unsigned int field_idx=5;
            for (; field_idx < 9; field_idx++) {
              while (*line_in != '\t') {
                line_in++;
              }
              line_in++;
            }
            //samples
            bool is_last=false;
            while (!is_last) {
                unsigned int allele_idx=(*line_in-'0');
                line_in++;
                while (std::isdigit(*line_in)) {
                    allele_idx*=10;
                    allele_idx+=(*line_in-'0');
                    line_in++;
                }
                while (*line_in!='\t') {
                    if (*line_in=='\n'||*line_in==0) {
                        is_last=true;
                        break;
                    }
                    line_in++;
                }
                //output prototype of mutation, and a map from sample to non-ref allele
                if (allele_idx>=(allele_translated.size()+1)) {
                    put_sample(parsed_line,header[field_idx],0xf,pos,chrom_idx,ref_nuc);
                }else if (allele_idx) {
                    put_sample(parsed_line,header[field_idx],allele_translated[allele_idx-1],pos,chrom_idx,ref_nuc);
                }
                field_idx++;
                line_in++;
            }
            //assert(field_idx==header.size());
            parsed_line.output(out);
        }
        std::get<0>(out).try_put(start);
    }
};
//tokenize header, get sample name
static int read_header(gzFile* fd,std::vector<std::string>& out){
    int header_len=0;
    char in=gzgetc(*fd);
    in=gzgetc(*fd);
    bool second_char_pong=(in=='#');

    while (second_char_pong) {
        while (in!='\n') {
            in=gzgetc(*fd);
        }
        in=gzgetc(*fd);
        in=gzgetc(*fd);
        second_char_pong=(in=='#');
    }

    bool eol=false;
    while (!eol) {
        std::string field;
        while (in!='\t') {
            if (in=='\n') {
                eol=true;
                break;
            }
            field.push_back(in);
            in=gzgetc(*fd);
            header_len++;
        }
        in=gzgetc(*fd);
        out.push_back(field);
    }
    return header_len;
}
std::atomic<size_t> assigned_count;
struct Assign_State{
    const std::vector<backward_pass_range>& child_idx_range;
    const std::vector<forward_pass_range>& parent_idx;std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>> &output;
    void operator()(const Parsed_VCF_Line* vcf_line)const{
        Fitch_Sankoff_Whole_Tree(child_idx_range,parent_idx,vcf_line->mutation,vcf_line->mutated,output);
        assigned_count.fetch_add(1,std::memory_order_relaxed);
        delete vcf_line;
    }
};
void print_progress(std::atomic<bool>* done,std::mutex* done_mutex){
    while (true) {
        {
            std::unique_lock<std::mutex> lk(*done_mutex);
            if (timed_print_progress) {
                progress_bar_cv.wait_for(lk,std::chrono::seconds(1));

            }else {
                progress_bar_cv.wait(lk);
            }
            if (done->load()) {
                return;
            }
        }
        fprintf(stderr,"\rAssigned %zu locus",assigned_count.load(std::memory_order_relaxed));
    }
}
std::unordered_set<std::string> get_existing_samples(const MAT::Tree& tree) {
    std::unordered_set<std::string> existing_samples;
    for(const auto& node:tree.all_nodes){
        if (node.second->is_leaf()) {
            auto iter=tree.condensed_nodes.find(node.first);
            if (iter==tree.condensed_nodes.end()) {
                existing_samples.insert(node.first);
            }else{
                existing_samples.insert(iter->second.begin(),iter->second.end());
            }
        }
    }
    return existing_samples;
}
#define CHUNK_SIZ 10
void init(size_t header_size,decompressor_node_t& decompressor,tbb::flow::graph& input_graph,gzFile& fd,std::atomic<bool>& done){
        for (int i=0; i<10; i++) {
        auto chunk=(char*)malloc((CHUNK_SIZ+2)*header_size);
        if (!chunk) {
            fputs("VCF parser chunk size too large \n", stderr);
            exit(1);
        }
        decompressor.try_put(chunk);
    }
    input_graph.wait_for_all();
    gzclose(fd);
    done=true;
    progress_bar_cv.notify_all();

}
void VCF_input(const char * name,MAT::Tree& tree,std::vector<New_Sample_t>& new_samples,bool reassign_state){
    assigned_count=0;
    std::vector<std::string> fields;
    //open file set increase buffer size
    gzFile fd=gzopen(name, "r");
    if (!fd) {
        fprintf(stderr, "cannnot open vcf file : %s, exiting.\n",name);
        exit(EXIT_FAILURE);
    }
    gzbuffer(fd,ZLIB_BUFSIZ);
    unsigned int header_size=read_header(&fd, fields);
    tbb::flow::graph input_graph;
    std::atomic<bool> done(false);
    std::mutex done_mutex;
    std::thread progress_meter(print_progress,&done,&done_mutex);
    decompressor_node_t decompressor(input_graph,1,Decompressor{&fd,CHUNK_SIZ*header_size,2*header_size});
    std::vector<long> idx_map(9);

    std::vector<std::string> new_sample_names;
    new_sample_mut_t new_sample_mut;
    if(reassign_state){
    std::vector<MAT::Node*> bfs_ordered_nodes=tree.breadth_first_expansion();
    for (size_t idx=9; idx<fields.size(); idx++) {
        auto iter=tree.all_nodes.find(fields[idx]);
        if (iter==tree.all_nodes.end()) {
            idx_map.push_back(-new_sample_names.size());
            new_sample_names.push_back(std::move(fields[idx]));
        }else {
            idx_map.push_back(iter->second->bfs_index);
        }
    }
    new_sample_mut.resize(new_sample_names.size());
    line_parser_t parser(input_graph,tbb::flow::unlimited,line_parser<produce_exist_sample_vcf_line>{idx_map,new_sample_mut});
    tbb::flow::make_edge(tbb::flow::output_port<0>(parser),decompressor);
    //feed used buffer back to decompressor
    tbb::flow::make_edge(tbb::flow::output_port<0>(decompressor),parser);

    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>> output(bfs_ordered_nodes.size());
    std::vector<backward_pass_range> child_idx_range;
    std::vector<forward_pass_range> parent_idx;
    Fitch_Sankoff_prep(bfs_ordered_nodes,child_idx_range, parent_idx);
    tbb::flow::function_node<Parsed_VCF_Line*> assign_state(input_graph,tbb::flow::unlimited,Assign_State{child_idx_range,parent_idx,output});
    tbb::flow::make_edge(tbb::flow::output_port<1>(parser),assign_state);
    init(header_size, decompressor, input_graph, fd, done);
    //Filling mutation vector
    tbb::affinity_partitioner ap;
        tbb::parallel_for(
        tbb::blocked_range<size_t>(0, bfs_ordered_nodes.size()),
        [&bfs_ordered_nodes, &output](tbb::blocked_range<size_t> r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                const auto &to_refill = output[i];
                bfs_ordered_nodes[i]->refill(to_refill.begin(), to_refill.end(),
                                             to_refill.size());
            }
        },
    ap);
    } else {
        auto existing_samples = get_existing_samples(tree);
        for (size_t idx = 9; idx < fields.size(); idx++) {
            auto iter = existing_samples.find(fields[idx]);
            if (iter == existing_samples.end()) {
                idx_map.push_back(-new_sample_names.size());
                new_sample_names.push_back(std::move(fields[idx]));
            } else {
                idx_map.push_back(1);
            }
        }
        new_sample_mut.resize(new_sample_names.size());
        buffer_only_out parser(input_graph,tbb::flow::unlimited,line_parser<ignore_existing_sample>{idx_map,new_sample_mut});
        tbb::flow::make_edge(tbb::flow::output_port<0>(parser),decompressor);
        //feed used buffer back to decompressor
        tbb::flow::make_edge(tbb::flow::output_port<0>(decompressor),parser);
        init(header_size, decompressor, input_graph, fd, done);
    }
    new_samples.resize(new_sample_names.size());
    for(size_t new_sample_idx=0;new_sample_idx<new_sample_names.size();new_sample_idx++){
        new_samples[new_sample_idx].name=std::move(new_sample_names[new_sample_idx]);
        new_samples[new_sample_idx].mutations.insert(new_samples[new_sample_idx].mutations.end(),new_sample_mut[new_sample_idx].begin(),new_sample_mut[new_sample_idx].end());
        std::sort(new_samples[new_sample_idx].mutations.begin(),new_samples[new_sample_idx].mutations.end());
    }
    progress_meter.join();
}
