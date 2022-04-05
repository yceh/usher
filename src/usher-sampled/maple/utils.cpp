#include "../usher.hpp"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <vector>
void distancesFromRefPunishNs(std::vector<Sample_Muts>& to_sort){
    tbb::parallel_for(tbb::blocked_range<size_t>(0,to_sort.size()),[&to_sort](tbb::blocked_range<size_t> range){
        for (int idx=range.begin(); idx<range.end(); idx++) {
            int count=1000*to_sort[idx].muts.size();
            for (const auto& mut : to_sort[idx].muts) {
                if (mut.mut_nuc==0xf) {
                    count+=mut.range;
                }
            }
            to_sort[idx].sorting_key1=count;
        }
    });
    tbb::parallel_sort(to_sort.begin(),to_sort.end(),[](const Sample_Muts& first,const Sample_Muts& second){
        return first.sorting_key1<second.sorting_key1;
    });
}

