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
    unlink(path);
    if(mkfifo(path, S_IRWXU|S_IRWXG|S_IRWXO)){
        perror("Error making fifo");
    }
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
            auto cur_offset=ftell(input_f);
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
    while (true) {
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
    while (true) {
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
To_Write_T gather_write(char* fifo_path){
    auto extract_fifo_file=open_fifo(fifo_path);
    std::string temp;        
    std::getline(extract_fifo_file,temp);
    To_Write_T out(std::atoi(temp.c_str()));  
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
    return out;
}
void copy_file_range(off64_t in_offset,off64_t in_end){
        off_t out_offset=lseek(out_fd, 0, SEEK_CUR);
        while (in_offset<stat_data.st_size) {
            if(copy_file_range(in_fd, &offset,out_fd,&out_offset, stat_data.st_size-offset,0)==-1){
                perror("sendfile");
                exit(1);
            }
        }
}
void write_file(int idx, const std::vector<std::string>& sequence_to_write,char* output_format_string,int ref_fd,size_t ref_size){
    char out_buf_name[BUFSIZ];
    snprintf(out_buf_name, BUFSIZ, output_format_string,idx);
    auto out_fd=open(out_buf_name, O_WRONLY|O_CREAT|O_TRUNC);
    for (const auto& name : sequence_to_write) {
        
    }
}
int main(int argc, char** argv){
    // input fifo path, rename fifo path, seq_name_output, extract fifo path, output file name template, reference file path
    auto loader_out=loader(argv[1], argv[2]);
    {
        std::fstream seq_name_f(argv[3]);
        for (const auto& ele : loader_out) {
            seq_name_f<<ele.first<<"\n";
        }
    }
    auto to_write=gather_write(argv[3]);



    auto out_fd=open(argv[1], O_CREAT|O_WRONLY|O_TRUNC);
    if (out_fd==-1) {
        perror("out fail");
        exit(EXIT_FAILURE);
    }
    auto in_fd=open(argv[2],O_RDONLY);
    if (in_fd==-1) {
        perror("in fail");
        exit(EXIT_FAILURE);
    }
    struct stat stat_data;
    fstat(in_fd,&stat_data);
    auto times=atoi(argv[3]);
    std::string interperse=argv[4];
    interperse+="\n";
    for (int i=0; i<times; i++) {
        write(out_fd,interperse.data(),interperse.size());
        off_t offset=0;
        off_t out_offset=lseek(out_fd, 0, SEEK_CUR);
        while (offset<stat_data.st_size) {
            if(copy_file_range(in_fd, &offset,out_fd,&out_offset, stat_data.st_size-offset,0)==-1){
                perror("sendfile");
                exit(1);
            }
        }
        lseek(out_fd, 0, SEEK_END);
    }
    
}
