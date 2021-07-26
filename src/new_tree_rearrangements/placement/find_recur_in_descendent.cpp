#include "../mutation_annotated_tree.hpp"
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdio>
#include <mutex>
#include <tbb/concurrent_vector.h>
#include <tbb/task.h>
#include <vector>
extern unsigned int search_radius;
extern std::vector<std::atomic<size_t>> early_stop_down_saving;
extern std::vector<std::atomic<size_t>> early_stop_up_saving;
extern std::vector<std::atomic<size_t>> early_stop_sibling_saving;
template <typename value_type> class range {
    const value_type* curr;
    const value_type* end;

  public:
    template<typename T>
    range(const T &container) : curr(container.data()), end(container.data()+container.size()) {}
    operator bool() const { return curr != end; }
    const value_type* operator->() const { return curr; }
    const value_type &operator*() const { return *curr; }
    void operator++() { curr++; }
    void operator++(int) { curr++; }
    size_t size()const{ return end-curr;}
};
namespace MAT = Mutation_Annotated_Tree;
struct Ancestral_Limit{
    MAT::Mutation* mut;
    unsigned int distance;
    unsigned int last_merge;
    Ancestral_Limit(MAT::Mutation* mut):mut(mut),distance(0),last_merge(0){}
    bool increment(bool merged){
        distance++;
        if (merged) {
            last_merge=distance;
        }
        if (distance>(search_radius+1)) {
            mut->set_mut_tip(last_merge);
        }
        return false;
    }
};
struct Mut_Related :public MAT::Mutation{
    MAT::Mutation* last_level_mutation;
    std::array<std::vector<Ancestral_Limit>,4> last_tip_with_state;
    Mut_Related (const MAT::Mutation& self, MAT::Mutation* ptr):MAT::Mutation(self){
        last_level_mutation=ptr;
    }
    void inc_distance(bool is_merge){
        for (char nuc_idx=0; nuc_idx<4; nuc_idx++) {
            auto& nuc_vect=last_tip_with_state[nuc_idx];
            if (nuc_vect.empty()) {
                continue;
            }
            auto begin=nuc_vect.begin();
            auto valid_end=nuc_vect.end()-1;
            while (begin!=valid_end) {
                if (begin->increment(is_merge)) {
                    std::swap(begin,valid_end);
                    valid_end--;
                }else{
                    begin++;
                }
            }
            assert(begin==valid_end);
            if (begin->increment(is_merge)) {
                valid_end--;
            }
            nuc_vect.erase(valid_end+1,nuc_vect.end());
        }
    }
};
struct Heap_Merger {
    MAT::Node* this_node;
    struct comparator {
        bool operator()(const range<Mut_Related> &op1,
                        const range<Mut_Related> &op2) {
            return op1->get_position() > op2->get_position();
        }
    };
    std::vector<range<Mut_Related>> heap;
    size_t total_count;
    Heap_Merger(const std::vector<std::vector<Mut_Related>> &to_merge,MAT::Node* this_node):this_node(this_node) {
        total_count=to_merge.size();
        heap.reserve(total_count);
        for (const auto &m : to_merge) {
            if (!m.empty()) {
                heap.emplace_back(m);
            }
        }
        std::make_heap(heap.begin(), heap.end(), comparator());
    }
    Mut_Related get_one() {
        Mut_Related output(*heap.front());
        std::vector<Mut_Related> to_set;
        char not_set_stat=0;
        size_t visited=0;
        char all_dec_mut=0;
        bool is_first=true;
        bool is_merged=false;
        while ((!heap.empty()) && heap.front()->get_position() == output.get_position()) {
            visited++;
            all_dec_mut|=heap.front()->get_descendent_mut();
            for(int i=0;i<4;i++){
                if (heap.front()->get_descendent_mut()&(1<<i)) {                
                    assert(!heap.front()->last_tip_with_state[i].empty()||(heap.front()->get_par_one_hot()&(1<<i)));
                }
            }
            if (is_first) {
                is_first=false;
            }else{
                for(int i=0;i<4;i++){
                    is_merged|=(!heap.front()->last_tip_with_state[i].empty());
                    output.last_tip_with_state[i].insert(output.last_tip_with_state[i].end(),heap.front()->last_tip_with_state[i].begin(),heap.front()->last_tip_with_state[i].end());
                }
            }
            if (heap.front()->last_level_mutation) {
                to_set.push_back(*heap.front());
            }else {
                not_set_stat|=heap.front()->get_descendent_mut();
            }
            std::pop_heap(heap.begin(), heap.end(), comparator());
            heap.back()++;
            if (!heap.back()) {
                heap.pop_back();
            } else {
                std::push_heap(heap.begin(), heap.end(), comparator());
            }
        }
        for(int i=0;i<4;i++){
            if (all_dec_mut&(1<<i)) {                
            assert(!output.last_tip_with_state[i].empty()||(output.get_par_one_hot()&(1<<i)));
            }
        }
        if (visited<total_count) {
            not_set_stat|=output.get_par_one_hot();
            all_dec_mut|=output.get_par_one_hot();
        }
        output.set_descendent(all_dec_mut);
        if (!to_set.empty()) {
            int to_set_size=to_set.size();
            std::vector<char> backward(to_set_size);
            backward[to_set_size-1]=to_set[to_set_size-1].get_descendent_mut();
            for (int idx=to_set_size-2; idx>=0; idx--) {
                backward[idx]=to_set[idx].get_descendent_mut()|backward[idx+1];
            }
            char forward_acc=0;
            for (int idx=0; idx<to_set_size-1; idx++) {
                to_set[idx].last_level_mutation->set_sibling(backward[idx+1]|not_set_stat|forward_acc);
                forward_acc|=to_set[idx].get_descendent_mut();
            }
            to_set[to_set_size-1].last_level_mutation->set_sibling(not_set_stat|forward_acc);
        }
        output.inc_distance(is_merged);
        return output;
    }
    operator bool() const { return !heap.empty(); }
};
void set_descendent_alleles(range<Mut_Related> descendent_alleles,
                            MAT::Node *node,
                            std::vector<Mut_Related> &out) {
    for (auto &mut : node->mutations) {
        while (descendent_alleles &&
               (descendent_alleles->get_position() < mut.get_position())) {
            out.push_back(std::move(*descendent_alleles));
            out.back().last_level_mutation=nullptr;
            descendent_alleles++;
        }
        char this_descendent_alleles= mut.get_mut_one_hot();
        if (descendent_alleles &&
            (descendent_alleles->get_position() == mut.get_position())) {
            this_descendent_alleles |=
                descendent_alleles->get_descendent_mut();
        }
        mut.set_descendent(this_descendent_alleles);
        out.emplace_back(mut,&mut);
        if (descendent_alleles &&
            (descendent_alleles->get_position() == mut.get_position())) {
            for(int idx=0;idx<4;idx++){
                if (mut.get_mut_one_hot()&(1<<idx)) {
                    auto& corresponding_entry=out.back().last_tip_with_state[idx];
                    //default distance is max radius, so no need to set here
                    corresponding_entry.clear();
                    corresponding_entry.push_back(&mut);
                }else if (mut.get_par_one_hot()&(1<<idx)) {
                    out.back().last_tip_with_state[idx].clear();
                }
                else {
                    out.back().last_tip_with_state[idx]=std::move(descendent_alleles->last_tip_with_state[idx]);
                }
            }
            descendent_alleles++;
        }else {
            char bin_idx=one_hot_to_two_bit(mut.get_mut_one_hot());
            out.back().last_tip_with_state[bin_idx].push_back(&mut);
        }
        for(int i=0;i<4;i++){
            if (out.back().get_descendent_mut()&(1<<i)) {                
                assert(!out.back().last_tip_with_state[i].empty()||(out.back().get_par_one_hot()&(1<<i)));
            }
        }
    }
    while (descendent_alleles) {
        out.push_back(std::move(*descendent_alleles));
        out.back().last_level_mutation=nullptr;
        descendent_alleles++;
    }
}
struct find_descendent_alleles_continuation : public tbb::task {
    MAT::Node *node;
    // All possible state of all loci that had mutation in subtree rooted at
    // each child Only the allele_present_among_descendents field is useful
    std::vector<std::vector<Mut_Related>> child_result;
    // output, including the mutations at node
    std::vector<Mut_Related> &out;
    find_descendent_alleles_continuation(MAT::Node *node,
                                         std::vector<Mut_Related> &out)
        : node(node), child_result(node->children.size()), out(out) {}
    task *execute() override {
        std::vector<Mut_Related> descendent_mutations;
        Heap_Merger merger(child_result,node);
        while (merger) {
            descendent_mutations.emplace_back(merger.get_one());
        }
        set_descendent_alleles(descendent_mutations, node, out);
        return nullptr;
    }
};
struct find_descendent_alleles : public tbb::task {
    MAT::Node *root;
    std::vector<Mut_Related> &out;
    find_descendent_alleles(MAT::Node *root,
                            std::vector<Mut_Related> &out)
        : root(root), out(out) {}
    tbb::task *execute() override {
        find_descendent_alleles_continuation *cont =
            new (allocate_continuation())
                find_descendent_alleles_continuation(root, out);
        std::vector<tbb::task *> tasks;
        for(auto& mut:root->mutations){
            mut.clean_auxilary();
        }
        for (size_t child_idx = 0; child_idx < root->children.size();
             child_idx++) {
            auto c = root->children[child_idx];
            if (c->is_leaf()) {
                auto& child_mut_out=cont->child_result[child_idx];
                child_mut_out.size();
                for (auto &mut : c->mutations) {
                    mut.clean_auxilary();
                    mut.set_descendent(mut.get_mut_one_hot());
                    child_mut_out.emplace_back(mut,&mut);
                    auto bin_idx = one_hot_to_two_bit(mut.get_mut_one_hot());
                    child_mut_out.back().last_tip_with_state[bin_idx].push_back(&mut);
                }
            } else {
                tasks.push_back(new (cont->allocate_child())
                                    find_descendent_alleles(
                                        c, cont->child_result[child_idx]));
            }
        }
        cont->set_ref_count(tasks.size());
        for (auto t : tasks) {
            cont->spawn(*t);
        }
        return (tasks.empty())?cont:nullptr;
    }
};
void placement_prep(MAT::Tree *t) {
    std::vector<Mut_Related> ignored;
    tbb::task::spawn_root_and_wait(*new (tbb::task::allocate_root()) find_descendent_alleles(t->root, ignored));
    for (auto& pos : ignored) {
        for (char base_idx=0; base_idx<4; base_idx++) {
            for(auto& mut:pos.last_tip_with_state[base_idx]){
                mut.mut->set_mut_tip(mut.last_merge);
            }
        }
    }
}
