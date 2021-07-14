#include "../mutation_annotated_tree.hpp"
#include <algorithm>
#include <bits/types/FILE.h>
#include <cstdio>
#include <mutex>
#include <tbb/concurrent_vector.h>
#include <tbb/task.h>
#include <vector>
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
struct Mut_Related :public MAT::Mutation{
    MAT::Mutation* last_level_mutation;
    std::array<std::vector<MAT::Mutation*>,4> last_tip_with_state;
    std::array<MAT::Node*,4> last_tip_updated_node;  
    Mut_Related (const MAT::Mutation& self, MAT::Mutation* ptr):MAT::Mutation(self){
        last_level_mutation=ptr;
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
        int visited=0;
        char all_dec_mut=0;
        bool is_first=true;
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
                    if (!heap.front()->last_tip_with_state[i].empty()) {
                        if (!output.last_tip_with_state[i].empty()) {
                            output.last_tip_updated_node[i]=this_node;
                        }else {
                            output.last_tip_updated_node[i]=heap.front()->last_tip_updated_node[i];
                        }
                    }
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
                    corresponding_entry.clear();
                    corresponding_entry.push_back(&mut);
                    out.back().last_tip_updated_node[idx]=node;
                }else {
                    if (this_descendent_alleles&(1<<idx)) {                        
                        assert(!descendent_alleles->last_tip_with_state[idx].empty());
                    }
                    out.back().last_tip_with_state[idx]=std::move(descendent_alleles->last_tip_with_state[idx]);
                    out.back().last_tip_updated_node[idx]=descendent_alleles->last_tip_updated_node[idx];
                }
            }
            descendent_alleles++;
        }else {
            char bin_idx=one_hot_to_two_bit(mut.get_mut_one_hot());
            out.back().last_tip_with_state[bin_idx].push_back(&mut);
            out.back().last_tip_updated_node[bin_idx]=node;
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
                    child_mut_out.back().last_tip_updated_node[bin_idx]=root;
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
    size_t tip_added=0;
    size_t tip_missed=0;
    std::vector<int> have_more_than_one_mut(t->all_nodes.size(),0);
    for (auto& pos : ignored) {
        if (pos.get_position()==26144) {
            std::fputc('a',stderr);
        }
        for(int i=0;i<4;i++){
            if((1<<i)&pos.get_ref_one_hot()){
                continue;
            }
            if (pos.get_descendent_mut()&(1<<i)) {                
            assert(!pos.last_tip_with_state[i].empty()||(pos.get_ref_one_hot()&(1<<i)));
            }
            if (pos.last_tip_with_state[i].size()==1) {
                pos.last_tip_with_state[i][0]->set_mut_tip();
                tip_added++;
            }else if(pos.last_tip_with_state[i].size()>1){
                have_more_than_one_mut[pos.last_tip_updated_node[i]->bfs_index]++;
                tip_missed++;
            }
        }
    }
    FILE* freq_out=fopen("merge_stat", "w");
    for(auto temp:have_more_than_one_mut){
        if (temp!=0) {
            std::fprintf(freq_out,"%d\n",temp);
        }
    }
    std::fclose(freq_out);
    std::fprintf(stderr,"Tip added: %zu, tip missed %zu\n",tip_added,tip_missed);
}
