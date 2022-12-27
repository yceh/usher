#include <csignal>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <thread>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>

struct Seq_Descriptor{
    int fd;
    off64_t start;
    off64_t end;
    Seq_Descriptor(){}
    Seq_Descriptor(int fd, off64_t start):fd(fd),start(start){}
    void set_length(off64_t end_offset){
        end=end_offset;
    }
};
struct loader_out {
    std::vector<std::string> file_name;
    std::unordered_map<std::string, Seq_Descriptor> pos_map;
};
std::fstream open_fifo(char* path){
/*    unlink(path);
    if(mkfifo(path, S_IRWXU|S_IRWXG|S_IRWXO)){
        perror("Error making fifo");
    }*/
    std::fstream f(path,std::ios::in);
    if (!f) {
        perror("error open");
    }
    return f;
}
typedef std::unordered_map<std::string, Seq_Descriptor> Loader_Out_T ;
typedef std::vector<std::pair<std::string, Seq_Descriptor>> Each_Input_T ;
void load_file(std::string command, int out_fd, Each_Input_T* offsets_out){
    auto input_f=popen(command.c_str(), "r");
    char* line=NULL;
    size_t line_size=0;
    auto out_f=fdopen(out_fd, "w");
    auto last_read_len=getline(&line, &line_size, input_f);
    while (last_read_len!=-1) {
        if (line[0]=='>') {
            std::string seq_name(line+1);
            seq_name.pop_back();
            auto cur_offset=ftell(out_f);
            if(!offsets_out->empty()){
                offsets_out->back().second.set_length(cur_offset);
            }
            offsets_out->emplace_back(seq_name,Seq_Descriptor(out_fd,cur_offset));
        }else {
            fwrite_unlocked(line, 1, last_read_len, out_f);
        }
        last_read_len=getline(&line, &line_size, input_f);
    }
    if(!offsets_out->empty()){
        offsets_out->back().second.set_length(ftell(out_f));
    }
    fflush(out_f);
}
Loader_Out_T loader(char* input_fifo_path,char* rename_fifo_path ){
    Loader_Out_T loader_out;
    auto input_fifo=open_fifo(input_fifo_path);
    auto rename_fifo=open_fifo(rename_fifo_path);
    std::vector<std::string> temp_paths;
    std::vector<int> temp_fds;
    std::vector<Each_Input_T*> loaded;
    std::vector<std::thread> loading_threads;
    while (input_fifo) {
        std::string temp;        
        std::getline(input_fifo,temp);
        if (temp=="END") {
            break;
        }
        temp_paths.push_back("XXXXXXX");
        temp_fds.push_back(mkstemp(const_cast<char*>(temp_paths.back().c_str())));
        fprintf(stderr, "%s", temp_paths.back().c_str());
        auto single_loaded_out=new Each_Input_T;
        loaded.push_back(single_loaded_out);
        loading_threads.emplace_back(load_file,temp,temp_fds.back(),single_loaded_out);
    }
    std::unordered_map<std::string, std::string> rename_map;
    while (rename_fifo) {
        std::string temp;        
        std::getline(rename_fifo,temp);
        if (temp=="END") {
            break;
        }
        auto iter=temp.find('\t');
        auto old_name=temp.substr(0,iter);
        auto new_name=temp.substr(iter+1,temp.size()-1-(iter+1));
        rename_map.emplace(old_name,new_name);
    }
    for (auto& t : loading_threads) {
        t.join();
    }
    for (auto ptr : loaded) {
        for(auto& seq:*ptr){
            auto iter=rename_map.find(seq.first);
            if(iter!=rename_map.end()){
                loader_out.emplace(iter->second,seq.second);
            }
        }
        delete ptr;
    }
    return loader_out;
}
typedef std::vector<std::vector<std::string>> To_Write_T;
void gather_write(char* fifo_path,To_Write_T& out){
    auto extract_fifo_file=open_fifo(fifo_path);
    std::string temp;        
    std::getline(extract_fifo_file,temp);

    out.resize(std::atoi(temp.c_str()));  
    while (true) {
        std::getline(extract_fifo_file,temp);
        if (temp=="END") {
            break;
        }
        auto iter=temp.find('\t');
        auto fid=atoi(temp.substr(0,iter).c_str());
        auto name=temp.substr(iter+1,temp.size()-1-(iter+1));
        out[fid].push_back(name);
    }
}
void append_file_range(off64_t in_offset,off64_t in_end,int in_fd, int out_fd){
        off_t out_offset=lseek(out_fd, 0, SEEK_CUR);
        while (in_offset<in_end) {
            if(copy_file_range(in_fd, &in_offset,out_fd,&out_offset, in_end-in_offset,0)==-1){
                perror("sendfile");
                raise(SIGTRAP);
            }
        }
        lseek(out_fd, 0, SEEK_END);
}
void write_file(int idx, const std::vector<std::string>& sequence_to_write,char* output_format_string,int ref_fd,size_t ref_size,const Loader_Out_T& seq_idx){
    char out_buf_name[BUFSIZ];
    snprintf(out_buf_name, BUFSIZ, output_format_string,idx);
    auto out_fd=open(out_buf_name, O_WRONLY|O_CREAT|O_TRUNC,S_IRUSR|S_IWUSR|S_IROTH|S_IRGRP);
    append_file_range(0, ref_size, ref_fd, out_fd);
    for (const auto& name : sequence_to_write) {
        auto iter=seq_idx.find(name);
        if (iter==seq_idx.end()) {
            fprintf(stderr, "%s not found\n", name.c_str());
            continue;
        }
        auto to_write=">"+name+"\n";
        write(out_fd, to_write.c_str(), to_write.size());
        append_file_range(iter->second.start, iter->second.end, iter->second.fd, out_fd);
    }
    close(out_fd);
}
int main(int argc, char** argv){
    // 1:input fifo path, 2:rename fifo path, 3:seq_name_output, 4:extract fifo path, 5:output file name template, 6:reference file path
    To_Write_T to_write;
    std::thread gather_write_thread(gather_write,argv[4],std::ref(to_write));
    auto loader_out=loader(argv[1], argv[2]);
    {
        std::fstream seq_name_f(argv[3]);
        for (const auto& ele : loader_out) {
            seq_name_f<<ele.first<<"\n";
        }
    }
    auto ref_fd=open(argv[6],O_RDONLY);
    if (ref_fd==-1) {
        perror("failed to open reference file");
        exit(EXIT_FAILURE);
    }
    struct stat stat_data;
    fstat(ref_fd,&stat_data);
    gather_write_thread.join();
    auto ref_size=stat_data.st_size;
    tbb::parallel_for(tbb::blocked_range<size_t>(0,to_write.size()),[ref_size,&loader_out,&to_write,argv,ref_fd](tbb::blocked_range<size_t>range){
        for (int idx=range.begin(); idx<range.end(); idx++) {
            write_file(idx, to_write[idx], argv[5],ref_fd ,ref_size, loader_out);
        }
    });
}
