#include "placement.hpp"
#include "src/new_tree_rearrangements/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <vector>
#include "../apply_move/allele_count.hpp"
#include "../Fitch_Sankoff.hpp"
void process_leaf_node(MAT::Node *node) {
    for (const MAT::Valid_Mutation &mut : node->valid_mutations) {
        node->mutations.push_back(MAT::Mutation(mut.chrom_idx, mut.position, 0,
                                                mut.get_mut_one_hot()));
        node->mutations.back().set_auxillary(mut.get_mut_one_hot(),
                                             0xf ^ (mut.get_mut_one_hot()));
    }
}
struct Heap_Merger {
    struct comparator {
        bool operator()(range<MAT::Mutation> &op1,
                        range<MAT::Mutation> &op2) {
            return op1->get_position() > op2->get_position();
        }
    };
    std::vector<range<MAT::Mutation>> heap;
    Heap_Merger(MAT::Node* node) {
        heap.reserve(node->children.size());
        for(auto c:node->children){
            heap.push_back(c->mutations);
        }
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
        while ((!heap.empty()) && heap.front().get_position() == output.get_position()) {
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
bool get_major_allele_polytomy(MAT::Node *node) {
    MAT::Mutations_Collection new_majoralleles_out;
    std::vector<Allele_Count_t> allele_count;
    Heap_Merger merger(node);
    while (merger) {
        auto allele_cnt=merger.get_one();
        auto par_nuc = allele_cnt.base.get_par_one_hot();
        int par_idx = one_hot_to_two_bit(par_nuc);
        allele_cnt.count[par_idx] += (node->children.size() - allele_cnt.node_cnt);
        uint8_t boundary1_major_allele;
        set_state_from_cnt(allele_cnt.count, boundary1_major_allele);

    }

    std::vector<int> backward_mut;
    for (const auto &allele : allele_count) {
        while (iter != end && iter->get_position() < allele.get_position()) {
            //original mutation shared by all children
            rewind_ori_mut_ploytomy(new_major_alleles_out, iter, changed);
            iter++;
        }
        nuc_one_hot major_alleles = boundary1_major_allele & 0xf;
        //Have a matching original mutation, set par_nuc from it
        if (iter != end && iter->get_position() == allele.get_position()) {
            MAT::Mutation altered = *iter;
            if (major_alleles != iter->get_all_major_allele()) {
                changed = true;
                nuc_one_hot ori_mut_nuc = altered.get_mut_one_hot();
                nuc_one_hot par_nuc = altered.get_par_one_hot();
                //set mut_nuc
                if (major_alleles & par_nuc) {
                    //follow parent if it can
                    altered.set_mut_one_hot(major_alleles & par_nuc);
                } else if (!(major_alleles & ori_mut_nuc)) {
                    //perserve original state if possible
                    altered.set_mut_one_hot(major_alleles.choose_first());
                }
                if (altered.get_mut_one_hot() != ori_mut_nuc) {
                    changed_positions.emplace_back(altered, ori_mut_nuc);
                }
            }
            altered.set_auxillary(major_alleles, boundary1_major_allele >> 4);
            //add if have boundary allele or major allele is different from parent allele
            if ((altered.get_par_one_hot() != altered.get_all_major_allele()) ||
                (altered.get_boundary1_one_hot() &&
                 node->children.size() > 1)) {
                new_major_alleles_out.push_back(altered);
            }
            iter++;
        } else {
            nuc_one_hot par_nuc = allele.base.get_par_one_hot();
            nuc_one_hot boundary1_allele = boundary1_major_allele >> 4;
            changed |= (major_alleles != par_nuc);
            //basically the same as above, just don't need to try to perserve original state
            if (boundary1_allele || major_alleles != par_nuc) {
                new_major_alleles_out.push_back(allele.base);
                MAT::Mutation &altered = new_major_alleles_out.back();
                if (par_nuc & major_alleles) {
                    altered.set_mut_one_hot(par_nuc & major_alleles);
                    altered.set_auxillary(major_alleles, boundary1_allele);
                } else {
                    altered.set_mut_one_hot(major_alleles.choose_first());
                    altered.set_auxillary(major_alleles, boundary1_allele);
                    changed_positions.emplace_back(altered, par_nuc);
                }
            }
        }
    }
    while (iter != end) {
        if (iter->get_all_major_allele() != iter->get_mut_one_hot()) {
            changed = true;
        }
        if (iter->is_valid()) {
            new_major_alleles_out.push_back(*iter);
            new_major_alleles_out.back().set_auxillary(iter->get_mut_one_hot(),
                                                       0);
        }
        iter++;
    }
    return changed;
}
