#include "../mutation_annotated_tree.hpp"
#include "src/new_tree_rearrangements/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include <mutex>
#include <string>
#include <tbb/concurrent_vector.h>
#include <tbb/null_mutex.h>
#include <utility>
#include <vector>
struct placed_sample_info {
    std::string sample_name;
    std::vector<MAT::Valid_Mutation> mutations_relative_to_this_node;
    placed_sample_info(
        std::string sample_name,
        std::vector<MAT::Valid_Mutation> mutations_relative_to_this_node)
        : sample_name(sample_name),
          mutations_relative_to_this_node(mutations_relative_to_this_node) {}
};
namespace MAT = Mutation_Annotated_Tree;
template <typename mutex_type> struct place_result {
    mutex_type mutex;
    MAT::Node *sibling_node;
    int parsimoy_score;
    std::vector<MAT::Valid_Mutation> shared;
    std::vector<MAT::Valid_Mutation> sibling_mut;
    std::vector<MAT::Valid_Mutation> new_insert_mut;
};
template <typename mutex_type>
static void try_place_as_sibling(
    MAT::Node *sibling_node, range<MAT::Valid_Mutation> mutations,
    std::vector<MAT::Valid_Mutation> &mutations_if_placed_under_sibling,
    place_result<mutex_type> &result) {
    std::vector<MAT::Valid_Mutation> shared;
    std::vector<MAT::Valid_Mutation> sibling_mut;
    std::vector<MAT::Valid_Mutation> new_insert_mut;
    for (const auto &mut : sibling_node->valid_mutations) {
        while (mutations && mutations->position < mut.position) {
            new_insert_mut.push_back(*mutations);
            mutations_if_placed_under_sibling.push_back(*mutations);
            mutations++;
        }
        if (mutations && (mutations->position == mut.position)) {
            if (mutations->get_mut_one_hot() != mut.get_mut_one_hot()) {
                mutations_if_placed_under_sibling.push_back(*mutations);
                sibling_mut.push_back(mut);
                new_insert_mut.push_back(*mutations);
            } else {
                shared.push_back(mut);
            }
            mutations++;
        } else {
            mutations_if_placed_under_sibling.push_back(mut);
            mutations_if_placed_under_sibling.back().set_mut_one_hot(
                mut.get_ref_one_hot());
            sibling_mut.push_back(mut);
        }
    }
    while (mutations) {
        new_insert_mut.push_back(*mutations);
        mutations_if_placed_under_sibling.push_back(*mutations);
        mutations++;
    }
    if (new_insert_mut.size() < result.parsimoy_score) {
        std::lock_guard<mutex_type> lk(result.mutex);
        if (new_insert_mut.size() < result.parsimoy_score) {
            result.sibling_node = sibling_node;
            result.shared = std::move(shared);
            result.sibling_mut = std::move(sibling_mut);
            result.new_insert_mut = std::move(new_insert_mut);
        }
    }
}
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
template <typename mutex_type>
static MAT::Node* place_node_partial(std::string& sample_name,place_result<mutex_type>& to_place,tbb::concurrent_vector<MAT::Node*>& nodes_need_identifier){
    auto new_sample_node=new MAT::Node();
    new_sample_node->identifier=std::move(sample_name);
    new_sample_node->valid_mutations=std::move(to_place.new_insert_mut);
    nodes_need_identifier.push_back(new_sample_node);
    auto sibling_parent=to_place.sibling_node->parent;
    if (to_place.shared.empty()) {
        sibling_parent->children.push_back(new_sample_node);
        new_sample_node->parent=sibling_parent;
        return nullptr;
    }else {
        auto new_internal_node=new MAT::Node();
        nodes_need_identifier.push_back(new_internal_node);
        new_internal_node->valid_mutations=std::move(to_place.shared);
        auto iter=std::find(sibling_parent->children.begin(),sibling_parent->children.end(),to_place.sibling_node);
        *iter=new_internal_node;
        new_internal_node->parent=sibling_parent;
        new_internal_node->children.push_back(to_place.sibling_node);
        to_place.sibling_node->parent=new_internal_node;
        to_place.sibling_node->valid_mutations=std::move(to_place.sibling_mut);
        new_internal_node->children.push_back(new_sample_node);
        return new_internal_node;
    }
}
static void process_sample_serial(MAT::Node* node,placed_sample_info& in,tbb::concurrent_vector<MAT::Node*>& nodes_need_identifier){
    place_result<tbb::null_mutex> result;
    
    find_best_sibling_node_truncated_serial(node, node, in.mutations_relative_to_this_node, result);
    place_node_partial(in.sample_name, result, nodes_need_identifier);

}

    for(const MAT::Valid_Mutation& mut:to_place.new_insert_mut){
        new_sample_node->mutations.push_back(MAT::Mutation(mut.chrom_idx,mut.position,0,mut.get_mut_one_hot()));
        new_sample_node->mutations.back().set_auxillary(mut.get_mut_one_hot(), 0xf^(mut.get_mut_one_hot()));
    }