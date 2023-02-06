#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
typedef std::unordered_map<std::string, std::unordered_set<std::string>> descendant_map_t;
static void trim(std::string& in){
    if(in[0]==' '){
        in=in.substr(1);
    }
}
int main(int argc, char** argv){
    std::fstream f1(argv[1]);
    descendant_map_t desc_map;
    while (f1) {
        std::string in;
        std::getline(f1,in);
        std::stringstream tab_dlim(in);
        std::string desc;
        std::getline(tab_dlim,desc,'\t');
        trim(desc);
        std::string desc_all;
        std::getline(tab_dlim,desc_all,'\t');
        std::stringstream desc_eles(desc_all);
        std::unordered_set<std::string> descs;
        while (desc_eles) {
            std::string desc_name;
            std::getline(desc_eles,desc_name,',');
            trim(desc_name);
            descs.insert(desc_name);
        }
        desc_map.emplace(desc,descs);
    }
    std::fstream f2(argv[2]);
    while (f2) {
        std::string in;
        std::getline(f2,in);
        std::stringstream tab_dlim(in);
        std::string desc;
        std::getline(tab_dlim,desc,'\t');
        trim(desc);
        auto iter=desc_map.find(desc);
        if (iter==desc_map.end()) {
            fprintf(stderr, "%s not in first file\n",desc.c_str());
            continue;
        }
        std::string desc_all;
        std::getline(tab_dlim,desc_all,'\t');
        std::stringstream desc_eles(desc_all);
        while (desc_eles) {
            std::string desc_name;
            std::getline(desc_eles,desc_name,',');
            trim(desc_name);
            auto desc_iter=iter->second.find(desc_name);
            if (desc_iter==iter->second.end()) {
                fprintf(stderr, "%s not found as descendat of %s in first file\n", desc_name.c_str(),desc.c_str());
            }else {
                iter->second.erase(desc_iter);
            }
        }
        for (const auto& ele : iter->second) {
                fprintf(stderr, "%s not found as descendat of %s in second file\n", ele.c_str(),desc.c_str());
        }
        desc_map.erase(iter);
    }
    for (const auto& ele : desc_map) {
        fprintf(stderr, "%s not in second file\n",ele.first.c_str());
    }
}