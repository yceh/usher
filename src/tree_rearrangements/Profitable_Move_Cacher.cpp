#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstring>
#include <sys/mman.h>
#include <sys/types.h>
#include <tbb/concurrent_vector.h>
#include <tbb/queuing_rw_mutex.h>
#include <unistd.h>
#define POS_OFFSET 0
#define MUT_NUC_OFFSET 4
#define REF_NUC_OFFSET 5
#define PAR_NUC_OFFSET 6
#define LCA_PARENT_STATE_OFFSET 7
#define META_SIZE 8

#define SRC_OFFSET 0
#define DST_OFFSET 1
#define LCA_OFFSET 2
#define RANGE_FIRST_OFFSET 3
#define RANGE_SECOND_OFFSET 4
#define MOVE_META_SIZE 5

static void serialize_FS_Results(const std::vector<Fitch_Sankoff_Result_Final>& to_serialize,int fd,size_t& offset){
    char temp[META_SIZE];
    size_t n_states=to_serialize.size();
    offset+=write(fd, &n_states, sizeof(n_states));
    for(const Fitch_Sankoff_Result_Final& result:to_serialize){
        //0-3 bytes position
        *((int*)temp+POS_OFFSET)=result.mutation.position;
        //4-6 byte mut_nuc
        temp[MUT_NUC_OFFSET]=result.mutation.mut_nuc;
        temp[REF_NUC_OFFSET]=result.mutation.is_missing?(0x80|result.mutation.ref_nuc):result.mutation.ref_nuc;
        temp[PAR_NUC_OFFSET]=result.mutation.par_nuc;
        //7 byte LCA
        temp[LCA_PARENT_STATE_OFFSET]=result.LCA_parent_state;
        offset+=write(fd, temp, META_SIZE);
        size_t dist=result.scores.size();
        assert(result.scores.size()==dist);
        offset+=write(fd, result.scores.data(), dist*sizeof(Fitch_Sankoff::Score_Type));
    }
}

static void deserialize_FS_Results(std::vector<Fitch_Sankoff_Result_Final>& out,size_t dist,const char* in){
    size_t n_states=*((size_t*) in);
    in+=sizeof(size_t);
    out.reserve(n_states);
    for(;n_states>0;n_states--){
        out.emplace_back();
        //0-3 bytes position
        MAT::Mutation& mut=out.back().mutation;
        mut.position=*((int*)in+POS_OFFSET);
        //4-6 byte mut_nuc
        mut.mut_nuc=in[MUT_NUC_OFFSET];
        mut.is_missing=in[REF_NUC_OFFSET]&0x80;
        mut.ref_nuc=in[REF_NUC_OFFSET]&0x7f;
        mut.par_nuc=in[PAR_NUC_OFFSET];
        out.back().LCA_parent_state=in[LCA_PARENT_STATE_OFFSET];
        in+=META_SIZE;
        out.back().scores=Fitch_Sankoff::Scores_Type(dist);
        memcpy(out.back().scores.data(), in,dist*sizeof(Fitch_Sankoff::Score_Type));
    }
}


struct Move_Sorter{
    bool operator()(const std::pair<int, size_t>& first, const std::pair<int, size_t>& second) const{
        return first.first<second.first;
    }
};

static void serialize_profitable_moves(tbb::concurrent_vector<Profitable_Move>& to_serialize,size_t& offset,std::vector<std::pair<int, size_t>>& file_offsets,int fd){
    size_t temp[MOVE_META_SIZE];
    for(Profitable_Move& ind:to_serialize){
        //save offset
        file_offsets.emplace_back(ind.score_change,offset);
        //save path first
        size_t n_hope=ind.path.size();
        offset+=write(fd, &n_hope, sizeof(size_t));
        offset+=write(fd, ind.path.data(), n_hope*(sizeof(MAT::Node*)));
        //save src, dst and LCA
        temp[SRC_OFFSET]=(size_t)ind.src;
        temp[DST_OFFSET]=(size_t)ind.dst;
        temp[LCA_OFFSET]=(size_t)ind.LCA;
        temp[RANGE_FIRST_OFFSET]=ind.range.first;
        temp[RANGE_SECOND_OFFSET]=ind.range.second;
        offset+=write(fd, temp, MOVE_META_SIZE*(sizeof(MAT::Node*)));
        //save Fitch Sankoff results
        serialize_FS_Results(ind.states,fd,offset);
    }
}

void Profitable_Moves_Cacher::get_path(std::vector<MAT::Node*>& out){
    
}
void Profitable_Moves_Cacher::operator()(){
    size_t offset=lseek(fd, 0, SEEK_SET);
    tbb::concurrent_vector<Profitable_Move> to_swap;
    while (true) {
        std::unique_lock<std::mutex> finished_lock(finish_mutex);
        finish_cv.wait_for(finished_lock,std::chrono::seconds(1));
        if(finished){
            finished_lock.unlock();
            break;
        }
        finished_lock.unlock();
        if (to_monitor.size()>20) {
            {
                tbb::queuing_rw_mutex::scoped_lock lock(swap_lock,true);
                to_swap.swap(to_monitor);
                
                lock.release();
            }
            serialize_profitable_moves(to_monitor, offset, file_offsets, fd);
        }
    }
    serialize_profitable_moves(to_monitor, offset, file_offsets, fd);
    std::sort(file_offsets.begin(),file_offsets.end(),Move_Sorter());
    mapped_address=mmap(nullptr, offset, PROT_READ, MAP_SHARED, fd, 0);
}

void Profitable_Moves_Cacher::run(){
        this_thread=std::thread(std::ref(*this));
}
Profitable_Moves_Cacher::Profitable_Moves_Cacher(tbb::concurrent_vector<Profitable_Move>& to_monitor,tbb::queuing_rw_mutex& rw_mutex):to_monitor(to_monitor),swap_lock(rw_mutex),finished(false){
        filename=(char*)malloc(10);
        filename="XXXXXX";
        fd=mkstemp(filename);
    }
void Profitable_Moves_Cacher::finish(){
        std::lock_guard<std::mutex> finish_lock(finish_mutex);
        finished=true;
        this_thread.join();
}
Profitable_Moves_Cacher::~Profitable_Moves_Cacher(){
        munmap(mapped_address, length);
        free(filename);
}