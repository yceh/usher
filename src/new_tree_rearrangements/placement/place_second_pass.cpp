#include "placement.hpp"
#include "src/new_tree_rearrangements/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <cstdio>
#include <limits>
#include <mutex>
#include <string>
#include <tbb/concurrent_vector.h>
#include <tbb/mutex.h>
#include <tbb/null_mutex.h>
#include <tbb/queuing_mutex.h>
#include <tbb/queuing_rw_mutex.h>
#include <utility>
#include <vector>
extern std::atomic<size_t> estimate_par_increase;
static void get_mutation_set_from_root(MAT::Node *node,
                                       MAT::Mutations_Collection &mutations,
                                       MAT::Node *stop_node) {
    MAT::Mutations_Collection temp(mutations);
    while (node && node != stop_node) {
        temp.merge(node->mutations, MAT::Mutations_Collection::KEEP_SELF);
        node = node->parent;
    }
    mutations.clear();
    mutations.reserve(temp.size());
    for (const auto &mut : temp) {
        if (mut.get_par_one_hot() != mut.get_all_major_allele()) {
            mutations.push_back(mut);
        }
    }
}
template <typename mutex_type> struct place_result {
    mutex_type mutex; // tbb::null_mutex for small bucket, actual mutex for
                      // large bucket
    MAT::Node *sibling_node; // sibling node or parent node of best placement
    int parsimoy_score;
    /*
            sibling_node->parent
                |
                |(shared)
                |
                new (MAT::Node)
            /           \
  (sibling /mut)     (new\ insert mut)
          /               \
        sibling_node     placed sample
    */
    MAT::Mutations_Collection shared;
    MAT::Mutations_Collection sibling_mut;
    std::vector<MAT::Valid_Mutation> new_insert_mut;
};
/**
 * @brief Try to place the new sample as sibling of sibling_node
 * @param sibling_node
 * @param[in] mutations mutations of the new sample to place relative to parent
 * of sibling_node
 * @param[out] mutations_if_placed_under_sibling mutations of the new sample to
 * place relative to sibling_node, for next recursion when trying children of
 * sibling node
 * @param[inout] result current best placement and corresponding branch
 * splitting
 * @return (void)
 */
template <typename mutex_type>
static void try_place_as_sibling(
    MAT::Node *sibling_node,
    range<MAT::Valid_Mutation> insert_mut_relative_to_parent,
    std::vector<MAT::Valid_Mutation> &mutations_if_placed_under_sibling,
    place_result<mutex_type> &result) {
    // same as in output
    MAT::Mutations_Collection shared_out;
    MAT::Mutations_Collection sibling_muts_out;
    std::vector<MAT::Valid_Mutation> new_insert_mut_out;
    int score =
        0; // parsimony score if placed as sibling or child of sibling_node
    for (const auto &sibling_mut : sibling_node->mutations) {
        if (sibling_mut.get_par_one_hot() ==
            sibling_mut.get_all_major_allele()) {
            continue;
        }
        // From parent, no change
        while (insert_mut_relative_to_parent &&
               insert_mut_relative_to_parent->position <
                   sibling_mut.get_position()) {
            if (!(insert_mut_relative_to_parent->get_par_one_hot() &
                  insert_mut_relative_to_parent->get_mut_one_hot())) {
                score++;
            }
            // Would be unique to placed sample
            new_insert_mut_out.push_back(*insert_mut_relative_to_parent);
            mutations_if_placed_under_sibling.push_back(
                *insert_mut_relative_to_parent);
            insert_mut_relative_to_parent++;
        }
        if (insert_mut_relative_to_parent &&
            (insert_mut_relative_to_parent->position ==
             sibling_mut.get_position())) {
            // coincide, split
            if (insert_mut_relative_to_parent->get_mut_one_hot() !=
                    sibling_mut.get_all_major_allele() ||
                insert_mut_relative_to_parent->get_mut_one_hot()
                    .is_ambiguous() ||
                sibling_mut.get_all_major_allele().is_ambiguous()) {
                // increment, as added a valid mutation to newly inserted sample
                // branch
                nuc_one_hot shared_nuc =
                    sibling_mut.get_all_major_allele() &
                    insert_mut_relative_to_parent->get_mut_one_hot();
                if (shared_nuc) {
                    shared_nuc = shared_nuc.choose_first();
                    shared_out.push_back(sibling_mut);
                    shared_out.back().set_mut_one_hot(shared_nuc);
                    if (insert_mut_relative_to_parent->get_mut_one_hot() !=
                        shared_nuc) {
                        new_insert_mut_out.push_back(
                            *insert_mut_relative_to_parent);
                        new_insert_mut_out.back().set_par_mut(
                            shared_nuc,
                            insert_mut_relative_to_parent->get_mut_one_hot());
                    }
                    if (shared_nuc != sibling_mut.get_all_major_allele()) {
                        sibling_muts_out.push_back(sibling_mut);
                        sibling_muts_out.back().set_par_one_hot(shared_nuc);
                    }
                } else {
                    score++;
                    // split
                    sibling_muts_out.push_back(sibling_mut);
                    new_insert_mut_out.push_back(
                        *insert_mut_relative_to_parent);
                }
                if(sibling_mut.get_mut_one_hot()!=insert_mut_relative_to_parent->get_mut_one_hot()){
                    mutations_if_placed_under_sibling.push_back(
                        *insert_mut_relative_to_parent);
                    mutations_if_placed_under_sibling.back().set_par_mut(
                        sibling_mut.get_mut_one_hot(),
                        insert_mut_relative_to_parent->get_mut_one_hot());
                }

            } else {
                // coincide
                shared_out.push_back(sibling_mut);
            }
            insert_mut_relative_to_parent++;
        } else {
            if (sibling_mut.is_valid()) {
                // back mutation
                mutations_if_placed_under_sibling.push_back(sibling_mut);
                mutations_if_placed_under_sibling.back().set_par_mut(
                    sibling_mut.get_mut_one_hot(),
                    sibling_mut.get_par_one_hot());
            }
            // don't increment score, as it is spliting existing mutations
            sibling_muts_out.push_back(sibling_mut);
        }
    }
    while (insert_mut_relative_to_parent) {
        if (!(insert_mut_relative_to_parent->get_par_one_hot() &
              insert_mut_relative_to_parent->get_mut_one_hot())) {
            score++;
        }
        new_insert_mut_out.push_back(*insert_mut_relative_to_parent);
        mutations_if_placed_under_sibling.push_back(
            *insert_mut_relative_to_parent);
        insert_mut_relative_to_parent++;
    }
    if (score < result.parsimoy_score) {
        // std::lock_guard<mutex_type> lk(result.mutex);
        if (score < result.parsimoy_score) {
            result.parsimoy_score = score;
            result.sibling_node = sibling_node;
            result.shared = std::move(shared_out);
            result.sibling_mut = std::move(sibling_muts_out);
            result.new_insert_mut = std::move(new_insert_mut_out);
        }
    }
}
// Trying to place as sibling of already placed samples in the bucket in DFS
// order, until reaching stop node, where it is another bucket
template <typename mutex_type>
static void find_best_sibling_node_truncated_serial(
    MAT::Node *start_node, const MAT::Node *stop_node,
    const std::vector<MAT::Valid_Mutation> &mutations,
    place_result<mutex_type> &out) {
    std::vector<MAT::Valid_Mutation> mutations_if_placed_under_start_node;
    try_place_as_sibling(start_node, mutations,
                         mutations_if_placed_under_start_node, out);
    if (start_node != stop_node) {
        for (const auto c : start_node->children) {
            find_best_sibling_node_truncated_serial(
                c, stop_node, mutations_if_placed_under_start_node, out);
        }
    }
}
// Forward pass on to_fill
void fill_mutation_collection(MAT::Node *to_fill,
                              const std::vector<MAT::Valid_Mutation> &filling) {
    to_fill->mutations.clear();
    to_fill->mutations.reserve(filling.size());
    for (const auto &mut : filling) {
        auto mut_nuc = mut.get_mut_one_hot() & mut.get_par_one_hot();
        if (!mut_nuc) {
            mut_nuc = mut.get_mut_one_hot().choose_first();
        }
        to_fill->mutations.push_back(MAT::Mutation(
            mut.chrom_idx, mut.position, mut.get_par_one_hot(), mut_nuc));
        to_fill->mutations.back().set_auxillary(mut.get_mut_one_hot(),
                                                0xf ^ (mut.get_mut_one_hot()));
    }
}
// Only fill vaid mutations for internal nodes, and ambiguous bases for leaf
// nodes
template <typename mutex_type>
static std::pair<MAT::Node *, MAT::Node *>
place_node_partial(std::string &sample_name, place_result<mutex_type> &to_place,
                   tbb::concurrent_vector<MAT::Node *> &nodes_need_identifier,MAT::Node* end_node,MAT::Node* ori_parent,tbb::queuing_rw_mutex &end_node_mutexes, tbb::queuing_rw_mutex &ori_parent_mutex) {
    auto new_sample_node = new MAT::Node();
    new_sample_node->identifier = std::move(sample_name);
    nodes_need_identifier.push_back(new_sample_node);
    fill_mutation_collection(new_sample_node, to_place.new_insert_mut);

    // Not necessary to split branch, as no sibling mutation, so as direct
    // children
    if (to_place.sibling_mut.empty() && (!to_place.sibling_node->is_leaf())) {
        if(to_place.sibling_node==end_node){
            tbb::queuing_rw_mutex::scoped_lock lk(end_node_mutexes,true);
            to_place.sibling_node->children.push_back(new_sample_node);
        }else{
            to_place.sibling_node->children.push_back(new_sample_node);
            assert(to_place.sibling_node->identifier=="");
        }
        new_sample_node->parent = to_place.sibling_node;
        return std::make_pair(nullptr, new_sample_node);
    } else {
        auto sibling_parent = to_place.sibling_node->parent;
        // Split the branch
        auto new_internal_node = new MAT::Node();
        nodes_need_identifier.push_back(new_internal_node);
        new_internal_node->mutations = std::move(to_place.shared);
        // Replace sibling node by internal node
        if (ori_parent==sibling_parent) {
            tbb::queuing_rw_mutex::scoped_lock lk(ori_parent_mutex,false);
            auto iter =
            std::find(sibling_parent->children.begin(),
                      sibling_parent->children.end(), to_place.sibling_node);
            *iter = new_internal_node;    
        }else{
            assert(sibling_parent->identifier=="");
        auto iter =
            std::find(sibling_parent->children.begin(),
                      sibling_parent->children.end(), to_place.sibling_node);
        *iter = new_internal_node;
        }
        new_internal_node->parent = sibling_parent;
        // Add sibling node
        new_internal_node->children.push_back(to_place.sibling_node);
        to_place.sibling_node->parent = new_internal_node;
        // Add sample node
        new_internal_node->children.push_back(new_sample_node);
        new_sample_node->parent = new_internal_node;
        // Fill mutations
        to_place.sibling_node->mutations = std::move(to_place.sibling_mut);
        assert(!new_internal_node->mutations.empty() ||
               !to_place.sibling_node->mutations.empty() ||
               to_place.sibling_node->is_leaf());
        return std::make_pair(new_internal_node, new_sample_node);
    }
}
// Sequential placement within each bucket
void process_sample_serial(
    MAT::Node *node, first_pass_per_node_t &samples,
    tbb::concurrent_vector<MAT::Node *> &nodes_need_identifier,
    std::vector<tbb::queuing_rw_mutex> &mutexes, tbb::queuing_rw_mutex &root_mutex) {
    if (node->is_root()) {
        for (auto &in : samples) {
            auto new_sample_node = new MAT::Node();
            new_sample_node->identifier = std::move(in.sample_name);
            nodes_need_identifier.push_back(new_sample_node);
            fill_mutation_collection(new_sample_node,
                                     in.mutations_relative_to_this_node);
            {
                tbb::queuing_rw_mutex::scoped_lock lk(mutexes[0],true);
                node->children.push_back(new_sample_node);
            }
            new_sample_node->parent=node;
        }
        return;
    }
    place_result<tbb::null_mutex> result;
    auto start_node = node;
    auto ori_parent = node->parent;
    for (auto &in : samples) {
        std::vector<MAT::Valid_Mutation> ori_sample_mutations(
            in.mutations_relative_to_this_node);

        result.parsimoy_score =
            std::numeric_limits<decltype(result.parsimoy_score)>::max();
        assert(start_node->parent == ori_parent);
        find_best_sibling_node_truncated_serial(
            start_node, node, in.mutations_relative_to_this_node, result);
        assert(result.parsimoy_score !=
               std::numeric_limits<decltype(result.parsimoy_score)>::max());
        estimate_par_increase.fetch_add(result.parsimoy_score);

#ifdef CHECK_PLACEMENT
        MAT::Mutations_Collection sibling_state_old;
        sibling_state_old.reserve((result.sibling_node->mutations).size());
        for (const auto &old_mut : result.sibling_node->mutations) {
            if (old_mut.get_all_major_allele() != old_mut.get_par_one_hot()) {
                sibling_state_old.push_back(old_mut);
            }
        }
#endif

        auto temp =
            place_node_partial(in.sample_name, result, nodes_need_identifier,node,ori_parent,mutexes[node->dfs_index],mutexes[ori_parent->dfs_index]);
        auto new_internal_node = temp.first;

#ifdef CHECK_PLACEMENT
        if (new_internal_node) {
            MAT::Mutations_Collection sibling_state_new(
                result.sibling_node->mutations);
            assert(result.sibling_node->parent == new_internal_node);
            sibling_state_new.merge(new_internal_node->mutations,
                                    MAT::Mutations_Collection::KEEP_SELF);
            auto old_iter = sibling_state_old.begin();
            for (const auto &mut : sibling_state_new) {
                if (mut.get_all_major_allele() == mut.get_par_one_hot()) {
                    continue;
                }
                assert(mut.get_position() == old_iter->get_position() &&
                       mut.get_all_major_allele() ==
                           old_iter->get_all_major_allele());
                old_iter++;
            }
            assert(old_iter == sibling_state_old.end());
        }
        MAT::Mutations_Collection sample_mutations_placed(
            temp.second->mutations);
        auto old_iter = ori_sample_mutations.begin();
        get_mutation_set_from_root(temp.second->parent, sample_mutations_placed,
                                   ori_parent);
        for (const auto &mut : sample_mutations_placed) {
            assert(mut.get_position() == old_iter->get_position() &&
                   mut.get_all_major_allele() == old_iter->get_mut_one_hot() &&
                   old_iter->get_par_one_hot() == mut.get_par_one_hot());
            old_iter++;
        }
        assert(old_iter == ori_sample_mutations.end());
#endif

        // Made a new internal node between start node and its parent (the
        // bucket), so start searching from this new internal node
        if (result.sibling_node == start_node && new_internal_node) {
            start_node = new_internal_node;
        }
    }
}