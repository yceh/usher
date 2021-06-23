#include "src/new_tree_rearrangements/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#include <mutex>
#include <string>
#include <tbb/concurrent_vector.h>
#include <tbb/null_mutex.h>
#include <utility>
#include <vector>
#include "placement.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"

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
    int score=0;
    for (const auto &mut : sibling_node->valid_mutations) {
        while (mutations && mutations->position < mut.position) {
            if(!(mutations->get_par_one_hot()&mutations->get_mut_one_hot())){
                score++;
            }
            new_insert_mut.push_back(*mutations);
            mutations_if_placed_under_sibling.push_back(*mutations);
            mutations++;
        }
        if (mutations && (mutations->position == mut.position)) {
            if (mutations->get_mut_one_hot() != mut.get_mut_one_hot()||mut.get_mut_one_hot().is_ambiguous()) {
                if(!(mutations->get_par_one_hot()&mutations->get_mut_one_hot())){
                    score++;
                }
                mutations_if_placed_under_sibling.push_back(*mutations);
                mutations_if_placed_under_sibling.back().set_par_mut(mut.get_mut_one_hot(),mutations->get_mut_one_hot());
                sibling_mut.push_back(mut);
                new_insert_mut.push_back(*mutations);
            } else {
                shared.push_back(mut);
            }
            mutations++;
        } else {
            score++;
            mutations_if_placed_under_sibling.push_back(mut);
            mutations_if_placed_under_sibling.back().set_par_mut(mut.get_mut_one_hot(),
                mut.get_par_one_hot());
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
template <typename mutex_type>
static MAT::Node* place_node_partial(std::string& sample_name,place_result<mutex_type>& to_place,tbb::concurrent_vector<MAT::Node*>& nodes_need_identifier){
    auto new_sample_node=new MAT::Node();
    new_sample_node->identifier=std::move(sample_name);
    fill_mutation_collection(new_sample_node, to_place.new_insert_mut);
    nodes_need_identifier.push_back(new_sample_node);
    auto sibling_parent=to_place.sibling_node->parent;
    if (to_place.shared.empty()) {
        sibling_parent->children.push_back(new_sample_node);
        new_sample_node->parent=sibling_parent;
        return nullptr;
    }else {
        auto new_internal_node=new MAT::Node();
        nodes_need_identifier.push_back(new_internal_node);
        fill_mutation_collection(new_internal_node, to_place.shared);
        auto iter=std::find(sibling_parent->children.begin(),sibling_parent->children.end(),to_place.sibling_node);
        *iter=new_internal_node;
        new_internal_node->parent=sibling_parent;
        new_internal_node->children.push_back(to_place.sibling_node);
        to_place.sibling_node->parent=new_internal_node;
        fill_mutation_collection(to_place.sibling_node, to_place.sibling_mut);
        new_internal_node->children.push_back(new_sample_node);
        return new_internal_node;
    }
}
static void process_sample_serial(MAT::Node* node,first_pass_per_node_t& samples,tbb::concurrent_vector<MAT::Node*>& nodes_need_identifier){
    place_result<tbb::null_mutex> result;
    auto start_node=node;
    for(auto& in:samples){
        find_best_sibling_node_truncated_serial(start_node, node, in.mutations_relative_to_this_node, result);
        auto new_internal_node=place_node_partial(in.sample_name, result, nodes_need_identifier);
        if (result.sibling_node==start_node) {
            start_node=new_internal_node;
        }
    }
}