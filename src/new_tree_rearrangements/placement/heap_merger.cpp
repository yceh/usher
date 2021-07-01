#include "placement.hpp"
#include "src/new_tree_rearrangements/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <vector>
#include "../apply_move/allele_count.hpp"
#define HEAP_MERGE_THESHOLD 4
template<typename T>
struct Heap_Merger {
    struct comparator {
        bool operator()(range<T> &op1,
                        range<T> &op2) {
            return op1->get_position() > op2->get_position();
        }
    };
    std::vector<range<T>> heap;
    Heap_Merger(std::vector<range<T>>& in):heap(std::move(in)) {
        std::make_heap(heap.begin(), heap.end(), comparator());
    }
    Allele_Count_t get_one() {
        Allele_Count_t output=*heap.front();
        std::pop_heap(heap.begin(), heap.end(), comparator());
        ++heap.back();
        if (!heap.back()) {
            heap.pop_back();
        } else {
            std::push_heap(heap.begin(), heap.end(), comparator());
        }
        while ((!heap.empty()) && heap.front()->get_position() == output.get_position()) {
            output+=*heap.front();
            std::pop_heap(heap.begin(), heap.end(), comparator());
            ++heap.back();
            if (!heap.back()) {
                heap.pop_back();
            } else {
                std::push_heap(heap.begin(), heap.end(), comparator());
            }
        }
        return output;
    }
    operator bool() const { return !heap.empty(); }
};