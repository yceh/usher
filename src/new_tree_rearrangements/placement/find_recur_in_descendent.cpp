#include "../Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "../mutation_annotated_tree.hpp"
#include <algorithm>
#include <mutex>
#include <tbb/concurrent_vector.h>
#include <tbb/task.h>
#include <vector>
namespace MAT = Mutation_Annotated_Tree;

struct Heap_Merger {
    struct comparator {
        bool operator()(const range<MAT::Valid_Mutation> &op1,
                        const range<MAT::Valid_Mutation> &op2) {
            return op1->position > op2->position;
        }
    };
    std::vector<range<MAT::Valid_Mutation>> heap;
    Heap_Merger(const std::vector<std::vector<MAT::Valid_Mutation>> &to_merge) {
        heap.reserve(to_merge.size());
        for (const auto &m : to_merge) {
            if (!m.empty()) {
                heap.emplace_back(m);
            }
        }
        std::make_heap(heap.begin(), heap.end(), comparator());
    }
    MAT::Valid_Mutation get_one() {
        MAT::Valid_Mutation output(*heap.front());
        while ((!heap.empty()) && heap.front()->position == output.position) {
            output.allele_present_among_descendents |=
                heap.front()->allele_present_among_descendents;
            std::pop_heap(heap.begin(), heap.end(), comparator());
            heap.back()++;
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
void set_descendent_alleles(range<MAT::Valid_Mutation> descendent_alleles,
                            MAT::Node *node,
                            std::vector<MAT::Valid_Mutation> &out) {
    auto &this_node_valid_mut = node->valid_mutations;
    this_node_valid_mut.clear();
    const auto &this_node_all_mut = node->mutations;
    for (const auto &mut : this_node_all_mut) {
        if (!mut.is_valid()) {
            continue;
        }
        while (descendent_alleles &&
               (descendent_alleles->position < mut.get_position())) {
            out.push_back(*descendent_alleles);
            descendent_alleles++;
        }
        out.push_back(mut);
        auto &out_pushed = out.back();
        out_pushed.allele_present_among_descendents |= mut.get_all_major_allele();
        if (descendent_alleles &&
            (descendent_alleles->position == mut.get_position())) {
            out_pushed.allele_present_among_descendents |=
                descendent_alleles->allele_present_among_descendents;
            descendent_alleles++;
        }
        this_node_valid_mut.push_back(out_pushed);
    }
    while (descendent_alleles) {
        out.push_back(*descendent_alleles);
        descendent_alleles++;
    }
}
struct find_descendent_alleles_continuation : public tbb::task {
    MAT::Node *node;
    // All possible state of all loci that had mutation in subtree rooted at
    // each child Only the allele_present_among_descendents field is useful
    std::vector<std::vector<MAT::Valid_Mutation>> child_result;
    // output, including the mutations at node
    std::vector<MAT::Valid_Mutation> &out;
    find_descendent_alleles_continuation(MAT::Node *node,
                                         std::vector<MAT::Valid_Mutation> &out)
        : node(node), child_result(node->children.size()), out(out) {}
    task *execute() override {
        std::vector<MAT::Valid_Mutation> descendent_mutations;
        Heap_Merger merger(child_result);
        while (merger) {
            descendent_mutations.emplace_back(merger.get_one());
        }
        set_descendent_alleles(descendent_mutations, node, out);
        return nullptr;
    }
};
struct find_descendent_alleles : public tbb::task {
    MAT::Node *root;
    std::vector<MAT::Valid_Mutation> &out;
    find_descendent_alleles(MAT::Node *root,
                            std::vector<MAT::Valid_Mutation> &out)
        : root(root), out(out) {}
    tbb::task *execute() override {
        find_descendent_alleles_continuation *cont =
            new (allocate_continuation())
                find_descendent_alleles_continuation(root, out);
        std::vector<tbb::task *> tasks;
        for (size_t child_idx = 0; child_idx < root->children.size();
             child_idx++) {
            auto c = root->children[child_idx];
            if (c->is_leaf()) {
                for (const auto &mut : c->mutations) {
                    if (mut.is_valid()) {
                        c->valid_mutations.push_back(mut);
                        c->valid_mutations.back()
                            .allele_present_among_descendents =
                            mut.get_all_major_allele();
                        c->valid_mutations.back().set_mut_one_hot(
                            mut.get_all_major_allele());
                    }
                }
                cont->child_result[child_idx] = c->valid_mutations;
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
    std::vector<MAT::Valid_Mutation> ignored;
    tbb::task::spawn_root_and_wait(*new (tbb::task::allocate_root()) find_descendent_alleles(t->root, ignored));
}
