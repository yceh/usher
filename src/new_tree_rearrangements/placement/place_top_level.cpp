#include "placement.hpp"
#include "src/new_tree_rearrangements/check_samples.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <chrono>
#include <cstdio>
#include <mutex>
#include <string>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/queuing_rw_mutex.h>
#include <unordered_map>
#include "../apply_move/apply_move.hpp"
#include "../tree_rearrangement_internal.hpp"
std::atomic<size_t> estimate_par_increase(0);
#ifdef CHECK_PLACEMENT
static void get_mutation_set_from_root(MAT::Node *node,
                                       MAT::Mutations_Collection &mutations) {
    MAT::Mutations_Collection temp;
    while (node) {
        temp.merge(node->mutations, MAT::Mutations_Collection::KEEP_SELF);
        node = node->parent;
    }
    for (auto &m : temp) {
        if (m.get_ref_one_hot() != m.get_mut_one_hot()) {
            mutations.push_back(m);
        }
    }
}
static void get_mutation_set_from_root_ambiguious(MAT::Node *node,
                                       MAT::Mutations_Collection &mutations) {
    MAT::Mutations_Collection temp(node->mutations);
    auto last_node=node;
    node=node->parent;
    while (node) {
        temp.merge(node->mutations, MAT::Mutations_Collection::KEEP_SELF);
        assert(std::find(node->children.begin(),node->children.end(),last_node)!=node->children.end());
        last_node=node;
        node = node->parent;
    }
    for (auto &m : temp) {
        if (m.get_ref_one_hot() != m.get_all_major_allele()) {
            mutations.push_back(m);
        }
    }
}
#endif
void place_samples(std::vector<New_Sample_t>& new_samples,MAT::Tree *tree,Original_State_t& ori_state){
    auto sample_count=new_samples.size();
    auto start_time=std::chrono::steady_clock::now();
    tbb::parallel_for(tbb::blocked_range<size_t>(0,new_samples.size()),[&new_samples,&ori_state](const tbb::blocked_range<size_t>& range){
        for (size_t idx=range.begin(); idx<range.end(); idx++) {
            Mutation_Set mut_set;
            const auto samp=new_samples[idx];
            for (const auto& temp : samp.mutations) {
                MAT::Mutation mut(temp.chrom_idx,temp.position,temp.get_par_one_hot(),temp.get_mut_one_hot());
                mut_set.insert(mut);
            }
            ori_state.emplace(samp.name,std::move(mut_set));
        }
    });
    placement_prep(tree);
#ifdef CHECK_PLACEMENT
    for (const auto& sample : new_samples) {
        for (const auto & mut : sample.mutations) {
            assert(mut.get_par_one_hot()==MAT::Mutation::refs[mut.position]);
        }
    }
    std::vector<New_Sample_t> new_samples_ori(new_samples);
#endif
    auto dfs_ordered_nodes=tree->depth_first_expansion();
    std::vector<first_pass_per_node_t> output(dfs_ordered_nodes.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0,new_samples.size()),[&new_samples,tree,&output](const tbb::blocked_range<size_t>& range){
        for (size_t idx=range.begin(); idx<range.end(); idx++) {
            place_first_pass(new_samples[idx].name, new_samples[idx].mutations,tree,output);
        }
    });
    fprintf(stderr, "estimate increase: %zu\n",estimate_par_increase.load());
    std::vector<std::pair<MAT::Node*,first_pass_per_node_t>> non_empty_buckets;
    non_empty_buckets.reserve(dfs_ordered_nodes.size());
    for (size_t idx=0; idx<dfs_ordered_nodes.size(); idx++) {
        if (output[idx].size()) {
            non_empty_buckets.emplace_back(dfs_ordered_nodes[idx],std::move(output[idx]));
        }
    }
#ifdef CHECK_PLACEMENT
    std::unordered_map<std::string, const std::vector<MAT::Valid_Mutation>*> sample_mut_map;
    for (const auto & temp : new_samples_ori) {
        sample_mut_map.insert(std::make_pair(temp.name,&temp.mutations));
    }
    for (const auto & bucket : non_empty_buckets) {
        MAT::Mutations_Collection sibling_to_root;
        if (bucket.first->parent) {
            get_mutation_set_from_root(bucket.first->parent, sibling_to_root);
        }
        for(const auto& sample:bucket.second){
            MAT::Mutations_Collection last_step;
            last_step.reserve(sample.mutations_relative_to_this_node.size());
            for(const auto& mut: sample.mutations_relative_to_this_node){
                last_step.push_back(MAT::Mutation(mut.chrom_idx,mut.position,mut.get_par_one_hot(),mut.get_mut_one_hot()));
            }
            last_step.merge(sibling_to_root, MAT::Mutations_Collection::KEEP_SELF);
            auto ori_mut_vect=sample_mut_map[sample.sample_name];
            auto ori_iter=ori_mut_vect->begin();
            for(const auto& placed_mut:last_step){
                assert(placed_mut.get_position()==ori_iter->get_position());
                assert(placed_mut.get_mut_one_hot()==ori_iter->get_mut_one_hot());
                ori_iter++;
            }
            assert(ori_iter==ori_mut_vect->end());
        }
    }
#endif
    tbb::concurrent_vector<MAT::Node*> nodes_need_identifier;
    std::vector<tbb::queuing_rw_mutex> mutexes(dfs_ordered_nodes.size());
    tbb::queuing_rw_mutex root_mutex;
    estimate_par_increase=0;
    tbb::parallel_for(tbb::blocked_range<size_t>(0,non_empty_buckets.size()),[&non_empty_buckets,&nodes_need_identifier,&mutexes,&root_mutex](const tbb::blocked_range<size_t>& range){
        for (size_t idx=range.begin(); idx<range.end(); idx++) {
            process_sample_serial(non_empty_buckets[idx].first,non_empty_buckets[idx].second,nodes_need_identifier,mutexes,root_mutex);
        }
    });
    fprintf(stderr, "place increase: %zu\n",estimate_par_increase.load());
    std::vector<MAT::Node *> updated_nodes;
    updated_nodes.reserve(nodes_need_identifier.size());
    for(auto node:nodes_need_identifier){
        if (node->identifier=="") {
            node->identifier=std::to_string(++tree->curr_internal_node);
        }
        tree->all_nodes.emplace(node->identifier,node);
        if (!node->children.empty()) {
            updated_nodes.push_back(node);
        }
    }
#ifdef CHECK_PLACEMENT
    for (const auto & temp : new_samples_ori) {
        MAT::Mutations_Collection placed;
        auto node=tree->get_node(temp.name);
        assert(node);
        get_mutation_set_from_root_ambiguious(node, placed);
        auto ori_iter=temp.mutations.begin();
        for(const auto& placed_mut:placed){
            if(!(placed_mut.get_position()==ori_iter->get_position()&&placed_mut.get_all_major_allele()==ori_iter->get_mut_one_hot())){
                fprintf(stderr, "%s \n",temp.name.c_str());
                assert(false);
            }
            ori_iter++;
        }
        if(ori_iter!=temp.mutations.end()){
            fprintf(stderr, "%s \n",temp.name.c_str());
            assert(false);
        }
    }
#endif
    std::vector<Altered_Node_t> nodes_with_changed_states;
    check_samples(tree->root, ori_state, tree);
    tree->depth_first_expansion();
    reassign_backward_pass(updated_nodes, nodes_with_changed_states);
    forward_pass(nodes_with_changed_states);
    clean_tree(*tree);
    check_samples(tree->root, ori_state, tree);
    fprintf(stderr, "Take %ld second to place %zu nodes",std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now()-start_time).count(),sample_count);
}