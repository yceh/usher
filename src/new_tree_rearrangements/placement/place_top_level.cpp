#include "placement.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include "../apply_move/apply_move.hpp"
struct New_Sample_t{
    std::string name;
    std::vector<Mutation_Annotated_Tree::Valid_Mutation> mutations;
};
void place_samples(std::vector<New_Sample_t>& new_samples,MAT::Tree *tree){
    placement_prep(tree);
    auto dfs_ordered_nodes=tree->depth_first_expansion();
    std::vector<first_pass_per_node_t> output(dfs_ordered_nodes.size());
    first_pass_per_node_t above_root;
    tbb::parallel_for(tbb::blocked_range<size_t>(0,new_samples.size()),[&new_samples,tree,&output,&above_root](const tbb::blocked_range<size_t>& range){
        for (size_t idx=range.begin(); idx<range.end(); idx++) {
            place_first_pass(new_samples[idx].name, new_samples[idx].mutations,tree,output,above_root);
        }
    });
    std::vector<std::pair<MAT::Node*,first_pass_per_node_t>> non_empty_buckets;
    non_empty_buckets.reserve(dfs_ordered_nodes.size());
    for (size_t idx=0; idx<dfs_ordered_nodes.size(); idx++) {
        if (output[idx].size()) {
            non_empty_buckets.emplace_back(dfs_ordered_nodes[idx],output[idx]);
        }
    }
    if (!above_root.empty()) {
        auto new_root=new MAT::Node();
        new_root->identifier=std::to_string(++tree->curr_internal_node);
        new_root->children.push_back(tree->root);
        tree->root->parent=new_root;
        tree->root=new_root;
        non_empty_buckets.emplace_back(new_root,above_root);
    }
    tbb::concurrent_vector<MAT::Node*> nodes_need_identifier;
    tbb::parallel_for(tbb::blocked_range<size_t>(0,non_empty_buckets.size()+1),[&non_empty_buckets,nodes_need_identifier](const tbb::blocked_range<size_t>& range){
        for (size_t idx=range.begin(); idx<range.end()-1; idx++) {
            process_sample_serial(non_empty_buckets[idx].first,non_empty_buckets[idx].second,nodes_need_identifier);
        }
    });
    std::vector<MAT::Node *> updated_nodes;
    updated_nodes.reserve(nodes_need_identifier.size());
    for(auto node:nodes_need_identifier){
        if (node->identifier=="") {
            node->identifier=std::to_string(++tree->curr_internal_node);
        }
        tree->all_nodes.emplace(node->identifier,node);
    }
    std::vector<Altered_Node_t> nodes_with_changed_states;
    reassign_backward_pass(updated_nodes, nodes_with_changed_states);
    forward_pass(nodes_with_changed_states);
}