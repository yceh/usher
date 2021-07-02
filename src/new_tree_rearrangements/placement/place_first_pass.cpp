#include "../Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include <limits>
#include <mutex>
#include <tbb/concurrent_vector.h>
#include <tbb/task.h>
#include <vector>
#include "placement.hpp"
extern std::atomic<size_t> estimate_par_increase;
/**
 * @param[in] node
 * @param[out] lower_bound lower_bound of parsimony possible if placed at subtree rooted at node
 * @param[out] parsimony_score_if_split_here parsimony score if placed between node and its parent
 * @param[in] par_mut_iter Mutations if placed as child of parent of node
 * @param mutations_merged Mutation relative to node, for trying children of node
 * @return Whether there are mutations splitted under the new branching node, if not, for next level of recursion
 */
bool merged_mutation_if_continue_down(
    const MAT::Node *node, int &lower_bound, int &parsimony_score_if_split_here,
    range<MAT::Valid_Mutation> par_mut_iter,
    std::vector<MAT::Valid_Mutation> &mutations_merged) {
    bool have_shared=false;
    for (const auto &mut : node->valid_mutations) {
        while (par_mut_iter && (par_mut_iter->position < mut.position)) {
            // mutations needed if placed at parent and no change
            mutations_merged.push_back(*par_mut_iter);
            if(!(par_mut_iter->get_mut_one_hot()&par_mut_iter->get_par_one_hot())){
                parsimony_score_if_split_here++;
                // if none of the descendent mutation can mutate to the desired
                // allele, increment lower bound
                if (!(par_mut_iter->allele_present_among_descendents &
                      par_mut_iter->get_mut_one_hot())) {
                    lower_bound++;
                }
            }
            par_mut_iter++;
        }
        mutations_merged.push_back(mut);
        auto &inserted = mutations_merged.back();
        if (par_mut_iter && (par_mut_iter->position == mut.position)) {
            if (par_mut_iter->get_mut_one_hot() == mut.get_mut_one_hot()) {
                // coincide, and a mutation at node consumed that mutation
                mutations_merged.pop_back();
                have_shared=true;
                par_mut_iter++;
                continue;
            } else {
                // coincide, but split
                inserted.set_par_mut(mut.get_mut_one_hot(),par_mut_iter->get_mut_one_hot());
                par_mut_iter++;
            }
        } else {
            //back mutation
            inserted.set_par_mut(mut.get_mut_one_hot(),mut.get_par_one_hot());
        }
        //cannot test equality directly, as the sample may contain ambiguous mutations
        if(!(inserted.get_mut_one_hot()&inserted.get_par_one_hot())){
            parsimony_score_if_split_here++;
            if (!(inserted.get_mut_one_hot() &
                  inserted.allele_present_among_descendents)) {
                lower_bound++;
            }
        }
    }
    while (par_mut_iter) {
        mutations_merged.push_back(*par_mut_iter);
        if(!(par_mut_iter->get_mut_one_hot()&par_mut_iter->get_par_one_hot())){
            parsimony_score_if_split_here++;
            // if none of the descendent mutation can mutate to the desired
            // allele, increment lower bound
            if (!(par_mut_iter->allele_present_among_descendents &
                  par_mut_iter->get_mut_one_hot())) {
                lower_bound++;
            }
        }
        par_mut_iter++;
    }
    return have_shared;
}

int parsimony_score_if_merge_with_another_children(
    range<MAT::Valid_Mutation> par_mut_iter,
    range<MAT::Valid_Mutation> other_child_mut_iter) {
    int to_return = 0;
    for (; par_mut_iter; par_mut_iter++) {
        while (other_child_mut_iter &&
               (other_child_mut_iter->position < par_mut_iter->position)) {
            other_child_mut_iter++;
        }
        if (other_child_mut_iter &&
            (other_child_mut_iter->position == par_mut_iter->position)) {
            if (!(other_child_mut_iter->get_mut_one_hot() &
                par_mut_iter->get_mut_one_hot()) ){
                to_return++;
            }
            other_child_mut_iter++;
        } else {
            to_return++;
        }
    }
    return to_return;
}
struct placement_result_ifc {
    int best_score;
    // placing the sample between best node and its parent, or
    // between its parent and a placed sample that alreday splitted the edge
    // between best node and its parent gives lowest parsimony score
    const MAT::Node *best_node;
    std::vector<MAT::Valid_Mutation> mutations_relative_to_parent;
    std::mutex mutex;
};

static int place_first_pass_helper(
    placement_result_ifc &output,
    const std::vector<MAT::Valid_Mutation> &mutations, const MAT::Node *node,
    const std::vector<first_pass_per_node_t>
        &already_placed_samples,
    std::vector<MAT::Valid_Mutation> &mutations_merged) {
    int lower_bound =
        0; // lower bound of parsimony score if placed at subtree rooted at node
    int parsimony_score_if_split_here =
        0; // parsimony score if placed between node and its parent
    range<MAT::Valid_Mutation> par_mut_iter(
        mutations); // iterator for would-be mutations if placed at parent
    // mutations needed if placed at node
    bool have_shared=merged_mutation_if_continue_down(node, lower_bound,
                                     parsimony_score_if_split_here,
                                     par_mut_iter, mutations_merged);
    auto min_score = std::min(output.best_score, parsimony_score_if_split_here);
    // searching pending children
    /*for (const auto &placed_child : already_placed_samples[node->dfs_index]) {
        if (!placed_child.ready) {
            continue;
        }
        // Best estimate is all mutations in child consume one mutation from
        // parent
        if ((int)mutations.size() -
                (int)placed_child.mutations_relative_to_this_node.size() <
            min_score) {
            auto score = parsimony_score_if_merge_with_another_children(
                mutations, placed_child.mutations_relative_to_this_node);
            parsimony_score_if_split_here =
                std::min(parsimony_score_if_split_here, score);
            min_score =
                std::min(output.best_score, parsimony_score_if_split_here);
        }
    }*/
    // output
    if (parsimony_score_if_split_here < output.best_score) {
        std::lock_guard<std::mutex> lock(output.mutex);
        if (parsimony_score_if_split_here < output.best_score) {
            output.best_score = parsimony_score_if_split_here;
            output.best_node = node;
            output.mutations_relative_to_parent = std::move(mutations);
        }
    }
    return lower_bound;
}
/**
 * @brief Serial DFS placement
 * @param best_node
 * @param best_score
 * @param mutations mutation needed if placed at parent of "node"
 * @param node start searching from subtree rooted at node
 */
static void place_first_pass_serial_each_level(
    placement_result_ifc &output,
    const std::vector<MAT::Valid_Mutation> &mutations, MAT::Node *node,
    const std::vector<first_pass_per_node_t>
        &already_placed_samples) {
    std::vector<MAT::Valid_Mutation> mutations_merged;
    auto lower_bound = place_first_pass_helper(
        output, mutations, node, already_placed_samples, mutations_merged);
    // go down if can reduce parsimony score
    if (lower_bound < output.best_score) {
        for (auto c : node->children) {
            place_first_pass_serial_each_level(output, mutations_merged, c,
                                               already_placed_samples);
        }
    }
}
struct first_pass_functor_cont : public tbb::task {
    std::vector<MAT::Valid_Mutation> mutations;
    task *execute() override { return nullptr; }
};
#define THRESHOLD 30
struct first_pass_functor : public tbb::task {
    const MAT::Node *node;
    const std::vector<MAT::Valid_Mutation> &mutations;
    placement_result_ifc &output;
    const std::vector<first_pass_per_node_t>
        &already_placed_samples;
    first_pass_functor(
        placement_result_ifc &output,
        const std::vector<MAT::Valid_Mutation> &mutations, MAT::Node *node,
        const std::vector<first_pass_per_node_t>
            &already_placed_samples)
        : node(node), mutations(mutations), output(output),
          already_placed_samples(already_placed_samples) {}
    task *execute() override {
        std::vector<MAT::Valid_Mutation> mutations_merged;
        auto lower_bound = place_first_pass_helper(
            output, mutations, node, already_placed_samples, mutations_merged);
        // go down if can reduce parsimony score
        if (lower_bound < output.best_score) {
            //agglomerative pattern, if there are less work, do not parallelize
            if (node->dfs_end_index - node->dfs_index < THRESHOLD) {
                for (auto c : node->children) {
                    place_first_pass_serial_each_level(
                        output, mutations_merged, c, already_placed_samples);
                }
            } else {
                auto cont =
                    new (allocate_continuation()) first_pass_functor_cont();
                cont->mutations = std::move(mutations_merged);
                cont->set_ref_count(node->children.size());
                for (auto c : node->children) {
                    cont->spawn(*new (cont->allocate_child())
                                    first_pass_functor(output, cont->mutations,
                                                       c,
                                                       already_placed_samples));
                }
            }
        }
        return nullptr;
    }
};
/**
 * @brief First step of placement, cluster samples into their most parsimonious
 * placement under existing tree, may take some of already placed sample into
 * consideration
 * @param sample_name
 * @param mutations Mutations of the sample to place relative to ref
 * @param tree tree to place sample in
 * @param output samples placed between the node of correspinding bfs index and
 * its parent
 * @param above_root samples placed above root
 */
void place_first_pass(
    std::string &sample_name, std::vector<MAT::Valid_Mutation> &mutations,
    MAT::Tree *tree,
    std::vector<first_pass_per_node_t> &output) {
    placement_result_ifc result;
    result.best_node = nullptr;
    result.best_score = mutations.size()+1;
    tbb::task::spawn_root_and_wait(
        *new (tbb::task::allocate_root())
            first_pass_functor(result, mutations, tree->root, output));
    assert (result.best_node);
    estimate_par_increase.fetch_add(result.best_score);
    output[result.best_node->dfs_index].emplace_back(
        sample_name, result.mutations_relative_to_parent);
}
