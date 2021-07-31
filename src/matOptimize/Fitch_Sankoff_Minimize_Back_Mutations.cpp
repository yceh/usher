#include "Fitch_Sankoff.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <climits>
#include <cstdint>
#include <tbb/concurrent_vector.h>
#include <utility>
#include <vector>
#include "min_back.hpp"
namespace MAT = Mutation_Annotated_Tree;
// Update score entry at node and assumed allele (par_allele) refered by 'out'
// from contribution of a particular child referred in 'child_info'
static int update_info_per_allele(const backward_info &child_info,
                                  int par_allele,
                                  backward_info_per_allele &out) {
    // Initialization, worst values possible, better is smaller
    int best_allele = 5;
    unsigned int best_parsimony_score = INT_MAX;
    unsigned int best_bac_mutation_count = INT_MAX;
    // See which is the best allele assignment 'child_node_allele_2_bit' for
    // child_info, given 'par_allele'
    for (int child_node_allele_2_bit = 0; child_node_allele_2_bit < 4;
         child_node_allele_2_bit++) {
        const auto &child_info_this_allele =
            child_info[child_node_allele_2_bit];
        // Increment parsimony score over parsimony score at the subtree if
        // don't match
        auto this_par_score = child_info_this_allele.parsimony +
                              (child_node_allele_2_bit != par_allele);
        // No need to continue if parsimony is suboptimal
        if (this_par_score > best_parsimony_score) {
            continue;
        }
        // back mutations in the subtree rooted at child node
        // back mutation: A(par_allele)->B->A(par_allele), when par_allele is
        // the same as some mut_nuc in some mutations in the subtree, regardless
        // of B (the par_nuc in that mutation) + existing back mutations
        auto this_back_mut_count =
            child_info_this_allele.back_mutations_count +
            child_info_this_allele.get_total_muts_to_particular_allele(
                par_allele);

        // update if better parsimony
        if (this_par_score < best_parsimony_score) {
            best_parsimony_score = this_par_score;
            best_bac_mutation_count = this_back_mut_count;
            best_allele = child_node_allele_2_bit;
        } else if (this_back_mut_count < best_bac_mutation_count) {
            // update if less back mutations
            assert(this_par_score == best_parsimony_score);
            best_bac_mutation_count = this_back_mut_count;
            best_allele = child_node_allele_2_bit;
        }
    }
    assert(best_parsimony_score<=child_info[par_allele].parsimony);
    assert(best_allele != 5);
    assert(best_bac_mutation_count!=INT_MAX);
    out.parsimony += best_parsimony_score;
    out.back_mutations_count += best_bac_mutation_count;
    // Accumulate mutations count
    out.add_other_mutations(par_allele, child_info[best_allele]);
    // Add one mutation for parent node to child node
    if (par_allele != best_allele) {
        out.get_mut_count(best_allele, par_allele)++;
    }
    return best_allele;
}
static void individual_backward_pass(size_t node_idx,
                                     const backward_pass_range &children_range,
                                     std::vector<backward_info> &score_vec,
                                     std::vector<uint8_t> &best_allele) {
    auto child_end_idx =
        children_range.first_child_bfs_idx + children_range.child_size;
    for (size_t child_node_idx = children_range.first_child_bfs_idx;
         child_node_idx < child_end_idx; child_node_idx++) {
        assert(child_node_idx>node_idx);
        int best_allele_sel = 0;
        for (int this_node_allele_2_bit = 0; this_node_allele_2_bit < 4;
             this_node_allele_2_bit++) {
            int this_best_allele = update_info_per_allele(
                score_vec[child_node_idx], this_node_allele_2_bit,
                score_vec[node_idx][this_node_allele_2_bit]);
            best_allele_sel |= (this_best_allele << (2 * this_node_allele_2_bit));
            assert(((best_allele_sel >> (2 * this_node_allele_2_bit)) & 3)==this_best_allele);
        }
        best_allele[child_node_idx] = best_allele_sel;
    }
}

std::pair<size_t, size_t>
backward_pass(const std::vector<backward_pass_range> &child_idx_range,
              std::vector<uint8_t> &best_allele, const mutated_t &mutated,
              uint8_t ref_nuc,std::vector<backward_info>& score_vec) {
    best_allele.resize(child_idx_range.size());
    auto mutated_iter = mutated.begin();
    for (long node_idx = child_idx_range.size() - 1; node_idx >= 0;
         node_idx--) {
        auto child_size = child_idx_range[node_idx].child_size;

        // leaf node
        if (child_size == 0) {
            uint8_t allele;
            if (mutated_iter->first != node_idx) {
                allele = ref_nuc;
                assert(mutated_iter->first < node_idx);
            } else {
                allele = one_hot_to_two_bit(mutated_iter->second);
                mutated_iter++;
            }
            score_vec[node_idx].clear_leaf(allele);
            uint8_t this_allele = 0;
            for (int allele_idx = 0; allele_idx < 4; allele_idx++) {
                this_allele = (this_allele << 2) | allele;
            }
            assert((this_allele&3)==allele);
            assert(((this_allele>>2)&3)==allele);
            assert(((this_allele>>4)&3)==allele);
            assert(((this_allele>>6)&3)==allele);
            best_allele[node_idx] = this_allele;
        } else {
            score_vec[node_idx].clear();
            individual_backward_pass(node_idx, child_idx_range[node_idx],
                                     score_vec, best_allele);
        }
    }
    backward_info_per_allele ignored;
    auto root_best_allele=update_info_per_allele(score_vec[0], ref_nuc,ignored);
    best_allele[0]=root_best_allele<<(2*ref_nuc);
    auto root_allele = score_vec[0][ref_nuc];
    return std::make_pair(root_allele.parsimony,
                          root_allele.back_mutations_count);
}

static void
individual_forward_pass(uint8_t parent_state, uint8_t &this_state,
                        const MAT::Mutation &base,
                        tbb::concurrent_vector<MAT::Mutation> &output) {
    assert(parent_state<=3);
    this_state = (this_state >> (2 * parent_state)) & 3;
    assert(this_state<=3);
    if (this_state != parent_state) {
        MAT::Mutation this_mut(base);
        this_mut.set_par_mut(two_bit_to_one_hot(parent_state), two_bit_to_one_hot(this_state));
        this_mut.set_auxillary(two_bit_to_one_hot(this_state), 0);
        output.push_back(this_mut);
    }
}

static void forward_pass(
    const std::vector<forward_pass_range> &parent_idx,
    const Mutation_Annotated_Tree::Mutation &base,
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>
        &output,
    uint8_t ref_nuc, std::vector<uint8_t> &best_allele) {
    individual_forward_pass(ref_nuc, best_allele[0], base, output[0]);
    for (size_t node_idx = 1; node_idx < parent_idx.size(); node_idx++) {
        individual_forward_pass(
            best_allele[parent_idx[node_idx].parent_bfs_idx],
            best_allele[node_idx], base, output[node_idx]);
    }
}

std::pair<size_t, size_t> Fitch_Sankoff_Minimize_back_mutation(
    const std::vector<backward_pass_range> &child_idx_range,
    const std::vector<forward_pass_range> &parent_idx,
    const Mutation_Annotated_Tree::Mutation &base,
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>
        &output,
    const mutated_t &mutated,std::vector<backward_info>& score_vec,std::vector<uint8_t>& best_allele) {
    uint8_t ref_nuc=one_hot_to_two_bit(base.get_ref_one_hot());
    auto expected =
    backward_pass(child_idx_range, best_allele, mutated, ref_nuc,score_vec);
    forward_pass(parent_idx, base, output, ref_nuc, best_allele);
    return expected;
}