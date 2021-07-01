#include "src/new_tree_rearrangements/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include <limits>
#include <mutex>
#include <string>
#include <tbb/concurrent_vector.h>
#include <tbb/null_mutex.h>
#include <utility>
#include <vector>
#include "placement.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"

template <typename mutex_type> struct place_result {
    mutex_type mutex; //tbb::null_mutex for small bucket, actual mutex for large bucket
    MAT::Node *sibling_node; //sibling node or parent node of best placement
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
    std::vector<MAT::Valid_Mutation> shared;
    std::vector<MAT::Valid_Mutation> sibling_mut;
    std::vector<MAT::Valid_Mutation> new_insert_mut;
};
/**
 * @brief Try to place the new sample as sibling of sibling_node
 * @param sibling_node
 * @param[in] mutations mutations of the new sample to place relative to parent of sibling_node
 * @param[out] mutations_if_placed_under_sibling mutations of the new sample to place relative to sibling_node, for next recursion when trying children of sibling node
 * @param[inout] result current best placement and corresponding branch splitting
 * @return (void)
 */
template <typename mutex_type>
static void try_place_as_sibling(
    MAT::Node *sibling_node, range<MAT::Valid_Mutation> mutations,
    std::vector<MAT::Valid_Mutation> &mutations_if_placed_under_sibling,
    place_result<mutex_type> &result) {
    //same as in output
    std::vector<MAT::Valid_Mutation> shared;
    std::vector<MAT::Valid_Mutation> sibling_mut;
    std::vector<MAT::Valid_Mutation> new_insert_mut;
    int score=0; //parsimony score if placed as sibling or child of sibling_node
    for (const auto &mut : sibling_node->valid_mutations) {
        //From parent, no change
        while (mutations && mutations->position < mut.position) {
            if(!(mutations->get_par_one_hot()&mutations->get_mut_one_hot())){
                score++;
            }
            //Would be unique to placed sample
            new_insert_mut.push_back(*mutations);
            mutations_if_placed_under_sibling.push_back(*mutations);
            mutations++;
        }
        if (mutations && (mutations->position == mut.position)) {
            //coincide, split
            if (mutations->get_mut_one_hot() != mut.get_mut_one_hot()||mut.get_mut_one_hot().is_ambiguous()) {
                //increment, as added a valid mutation to newly inserted sample branch
                if(!(mutations->get_par_one_hot()&mutations->get_mut_one_hot())){
                    score++;
                }
                mutations_if_placed_under_sibling.push_back(*mutations);
                mutations_if_placed_under_sibling.back().set_par_mut(mut.get_mut_one_hot(),mutations->get_mut_one_hot());
                //split
                sibling_mut.push_back(mut);
                new_insert_mut.push_back(*mutations);
            } else {
                //coincide
                shared.push_back(mut);
            }
            mutations++;
        } else {
            //back mutation
            mutations_if_placed_under_sibling.push_back(mut);
            mutations_if_placed_under_sibling.back().set_par_mut(mut.get_mut_one_hot(),
                mut.get_par_one_hot());
            //don't increment score, as it is spliting existing mutations
            sibling_mut.push_back(mut);
        }
    }
    while (mutations) {
        if(!(mutations->get_par_one_hot()&mutations->get_mut_one_hot())){
            score++;
        }
        new_insert_mut.push_back(*mutations);
        mutations_if_placed_under_sibling.push_back(*mutations);
        mutations++;
    }
    //placing as children. next level
    if (sibling_mut.empty()) {
        return;
    }
    if (score < result.parsimoy_score) {
        std::lock_guard<mutex_type> lk(result.mutex);
        if (score < result.parsimoy_score) {
            result.parsimoy_score=score;
            result.sibling_node = sibling_node;
            result.shared = std::move(shared);
            result.sibling_mut = std::move(sibling_mut);
            result.new_insert_mut = std::move(new_insert_mut);
        }
    }
}
//Trying to place as sibling of already placed samples in the bucket in DFS order, until reaching stop node, where it is another bucket
template <typename mutex_type>
static void find_best_sibling_node_truncated_serial(
    MAT::Node *start_node, const MAT::Node *stop_node,
    const std::vector<MAT::Valid_Mutation> &mutations,place_result<mutex_type>& out){
    std::vector<MAT::Valid_Mutation> mutations_if_placed_under_start_node;
    try_place_as_sibling(start_node, mutations, mutations_if_placed_under_start_node, out);
    if (start_node!=stop_node) {
        for(const auto c:start_node->children){
            find_best_sibling_node_truncated_serial(c, stop_node, mutations_if_placed_under_start_node, out);
        }
    }
}
//Forward pass on to_fill 
void fill_mutation_collection(MAT::Node* to_fill, const std::vector<MAT::Valid_Mutation>& filling){
    to_fill->mutations.clear();
    to_fill->mutations.reserve(filling.size());
    for(const auto & mut:filling){
        auto mut_nuc=mut.get_mut_one_hot()&mut.get_par_one_hot();
        if (!mut_nuc) {
            mut_nuc=mut.get_mut_one_hot().choose_first();
        }
        to_fill->mutations.push_back(MAT::Mutation(mut.chrom_idx,mut.position,mut.get_par_one_hot(),mut_nuc));
        to_fill->mutations.back().set_auxillary(mut.get_mut_one_hot(), 0xf^(mut.get_mut_one_hot()));
    }
}
//Only fill vaid mutations for internal nodes, and ambiguous bases for leaf nodes
template <typename mutex_type>
static MAT::Node* place_node_partial(std::string& sample_name,place_result<mutex_type>& to_place,tbb::concurrent_vector<MAT::Node*>& nodes_need_identifier){
    auto new_sample_node=new MAT::Node();
    new_sample_node->identifier=std::move(sample_name);
    fill_mutation_collection(new_sample_node, to_place.new_insert_mut);
    nodes_need_identifier.push_back(new_sample_node);
    
    auto sibling_parent=to_place.sibling_node->parent;
    //Not necessary to split branch, as no shared mutation, so as direct children
    if (to_place.shared.empty()) {
        sibling_parent->children.push_back(new_sample_node);
        new_sample_node->parent=sibling_parent;
        return nullptr;
    }else {
        //Split the branch
        auto new_internal_node=new MAT::Node();
        nodes_need_identifier.push_back(new_internal_node);
        fill_mutation_collection(new_internal_node, to_place.shared);
        //Replace sibling node by internal node
        auto iter=std::find(sibling_parent->children.begin(),sibling_parent->children.end(),to_place.sibling_node);
        *iter=new_internal_node;
        new_internal_node->parent=sibling_parent;
        //Add sibling node
        new_internal_node->children.push_back(to_place.sibling_node);
        to_place.sibling_node->parent=new_internal_node;
        //Add sample node
        new_internal_node->children.push_back(new_sample_node);
        new_sample_node->parent=new_internal_node;
        //Fill mutations
        fill_mutation_collection(to_place.sibling_node, to_place.sibling_mut);
        return new_internal_node;
    }
}
//Sequential placement within each bucket
void process_sample_serial(MAT::Node* node,first_pass_per_node_t& samples,tbb::concurrent_vector<MAT::Node*>& nodes_need_identifier){
    place_result<tbb::null_mutex> result;
    auto start_node=node;
    for(auto& in:samples){
        result.parsimoy_score=std::numeric_limits<decltype(result.parsimoy_score)>::max();
        find_best_sibling_node_truncated_serial(start_node, node, in.mutations_relative_to_this_node, result);
        assert(result.parsimoy_score!=std::numeric_limits<decltype(result.parsimoy_score)>::max());
        assert(!result.sibling_mut.empty());
        auto new_internal_node=place_node_partial(in.sample_name, result, nodes_need_identifier);
        //Made a new internal node between start node and its parent (the bucket), so start searching from this new internal node
        if (result.sibling_node==start_node) {
            start_node=new_internal_node;
        }
    }
}