#include "placement.hpp"
#include "src/new_tree_rearrangements/check_samples.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <cstddef>
#include <cstdio>
#include <string>
#include <unordered_set>
#include <utility>
void merge_back(const std::vector<MAT::Mutation>& lower,const std::vector<MAT::Mutation>& upper,  std::vector<MAT::Mutation>& out){
    out.reserve(upper.size()+lower.size());
    auto lower_iter=lower.begin();
    for (const auto& upp : upper) {
        while(lower_iter!=lower.end()&&lower_iter->get_position()<upp.get_position()){
            out.push_back(*lower_iter);
            lower_iter++;
        }
        if (lower_iter!=lower.end()&&lower_iter->get_position()==upp.get_position()) {
            if (upp.get_par_one_hot()!=lower_iter->get_mut_one_hot()) {
                out.push_back(*lower_iter);
                out.back().set_par_one_hot(upp.get_par_one_hot());
            }
            lower_iter++;
        }else {
            out.push_back(upp);
        }
    }
    while(lower_iter!=lower.end()){
        out.push_back(*lower_iter);
        lower_iter++;
    }
}
#ifdef LITE_DETAIL
void get_mutations_relative_to_root(MAT::Node* node,std::vector<MAT::Mutation>& muts){
    while (node) {
        std::vector<MAT::Mutation> output;
        merge_back(muts, node->mutations, output);
        muts=std::move(output);
        node=node->parent;
    }
}
#endif
void merge_children(MAT::Node* node,MAT::Tree& tree,std::unordered_set<size_t>& removed_nodes);
void clean_up_src_par(MAT::Node *par,MAT::Tree* tree,std::unordered_set<size_t>& removed_nodes){
    if (par->children.size()==1&&(!par->is_root())) {
        std::vector<MAT::Mutation> out;
        auto only_child = par->children[0];
        merge_back(only_child->mutations, par->mutations,out );
        only_child->mutations = std::move(out);
        only_child->parent=par->parent;
        auto& par_children=par->parent->children;
        *(std::find(par_children.begin(),par_children.end(),par))=only_child;
        tree->all_nodes.erase(par->identifier);
        removed_nodes.insert((size_t) par);
        merge_children(only_child, *tree, removed_nodes);
    }
}
void merge_children(MAT::Node* node,MAT::Tree& tree,std::unordered_set<size_t>& removed_nodes){
    if (node->mutations.empty()&&(!node->is_leaf())&&(!node->is_root())) {
        auto& par_mut=node->parent->children;
        auto this_iter=std::find(par_mut.begin(),par_mut.end(),node);
        *this_iter=node->children[0];
        node->children[0]->parent=node->parent;
        for(int idx=1;idx<node->children.size();idx++){
            node->children[idx]->parent=node->parent;
            par_mut.push_back(node->children[idx]);
        }
        tree.all_nodes.erase(node->identifier);
        removed_nodes.insert((size_t)node);
        assert(par_mut.end()==std::find(par_mut.begin(),par_mut.end(),node));
    }
}
void compare_mutations(
    std::vector<MAT::Mutation> &original_mutations,
    std::vector<MAT::Mutation> &new_mutations) {
    auto new_mut_iter=new_mutations.begin();
    for (const auto &mut : original_mutations) {
        assert(mut.get_position() == new_mut_iter->get_position() &&
               mut.get_mut_one_hot() == new_mut_iter->get_mut_one_hot() &&
               mut.get_par_one_hot() == new_mut_iter->get_par_one_hot());
        new_mut_iter++;
    }
    assert(new_mut_iter==new_mutations.end());
}
void move_node(MAT::Node *src, MAT::Node *dst,
               const std::vector<MAT::Mutation> &mutations_relative_to_parent,
               size_t parsimony_score, MAT::Tree &tree,std::vector<MAT::Node*>& nodes_may_changed) {
#ifdef LITE_DETAIL
    std::vector<MAT::Mutation> original_mutations;
    get_mutations_relative_to_root(src, original_mutations);
    std::vector<MAT::Mutation> input_mutations(mutations_relative_to_parent);
    get_mutations_relative_to_root(dst->parent, input_mutations);
    compare_mutations(original_mutations, input_mutations);
    std::vector<MAT::Mutation> dst_mutations;
    get_mutations_relative_to_root(dst, dst_mutations);
#endif
    std::vector<MAT::Mutation> shared_mutations;
    std::vector<MAT::Mutation> new_insert_mut;
    std::vector<MAT::Mutation> sibling_mut;
    parsimony_score+=src->mutations.size();
    auto par_mut_iter=mutations_relative_to_parent.begin();
    for (const auto& sib : dst->mutations) {
        while (par_mut_iter!=mutations_relative_to_parent.end()&&par_mut_iter->get_position()<sib.get_position()) {
            new_insert_mut.push_back(*par_mut_iter);
            par_mut_iter++;
        }
        if (par_mut_iter!=mutations_relative_to_parent.end()&&par_mut_iter->get_position()==sib.get_position()) {
            if (par_mut_iter->get_mut_one_hot()==sib.get_mut_one_hot()) {
                shared_mutations.push_back(*par_mut_iter);
            }else {
                //split
                new_insert_mut.push_back(*par_mut_iter);
                sibling_mut.push_back(sib);
            }
            par_mut_iter++;
        }else {
            sibling_mut.push_back(sib);
        }
    }
    while (par_mut_iter!=mutations_relative_to_parent.end()) {
            new_insert_mut.push_back(*par_mut_iter);
            par_mut_iter++;
    }
    assert(new_insert_mut.size()==parsimony_score);
    auto& src_par_children=src->parent->children;
    src_par_children.erase(std::find(src_par_children.begin(),src_par_children.end(),src));
    nodes_may_changed.push_back(src->parent);
    src->mutations=std::move(new_insert_mut);
    if (shared_mutations.empty()) {
        src->parent=dst->parent;
        dst->parent->children.push_back(src);
    }
    else if ((!sibling_mut.empty())||dst->is_leaf()) {
        auto& dst_par_children=dst->parent->children;
        auto dst_iter=std::find(dst_par_children.begin(),dst_par_children.end(),dst);
        auto new_node=new MAT::Node;
        new_node->identifier=std::to_string(++tree.curr_internal_node);
        tree.all_nodes.emplace(new_node->identifier,new_node);
        *dst_iter=new_node;
        new_node->parent=dst->parent;
        new_node->mutations=std::move(shared_mutations);
        dst->parent=new_node;
        dst->mutations=std::move(sibling_mut);
        new_node->children.push_back(dst);
        src->parent=new_node;
        new_node->children.push_back(src);
    }else {
        dst->children.push_back(src);
        src->parent=dst;
    }
    if(src->mutations.size()==0&&(!src->is_leaf())){
        nodes_may_changed.push_back(src);
    }
    assert(dst->mutations.size()||dst->is_leaf());
    src->changed=true;
#ifdef LITE_DETAIL
    std::vector<MAT::Mutation> new_mutations;
    get_mutations_relative_to_root(src, new_mutations);
    compare_mutations(original_mutations, new_mutations);
    std::vector<MAT::Mutation> new_dst_mutations;
    get_mutations_relative_to_root(dst, new_dst_mutations);
    compare_mutations(dst_mutations, new_dst_mutations);

#endif
}
size_t single_move_last_stage(
    const MAT::Node *node,const std::vector<MAT::Mutation>& par_mut) {
    auto par_mut_iter=par_mut.begin();
    size_t parsimony_score_if_split_here=0;
    for (const auto &mut : node->mutations) {
        while (par_mut_iter!=par_mut.end() && (par_mut_iter->get_position() < mut.get_position())) {
            // mutations needed if placed at parent and no change
            assert(!(par_mut_iter->get_mut_one_hot()&par_mut_iter->get_par_one_hot()));
            parsimony_score_if_split_here++;
            par_mut_iter++;
        }
        if (par_mut_iter!=par_mut.end() && (par_mut_iter->get_position() == mut.get_position())) {
            assert(par_mut_iter->get_par_one_hot()==mut.get_par_one_hot());
            if (par_mut_iter->get_mut_one_hot() != mut.get_mut_one_hot()) {
                parsimony_score_if_split_here++;
            }
            par_mut_iter++;
        }
    }
    while (par_mut_iter!=par_mut.end()) {
        assert(!(par_mut_iter->get_mut_one_hot()&par_mut_iter->get_par_one_hot()));
        parsimony_score_if_split_here++;
        par_mut_iter++;
    }
    return parsimony_score_if_split_here;
}

void apply_moves(std::vector<Profitable_Moves_ptr_t> &all_moves, MAT::Tree &t,
                 std::vector<MAT::Node *> &bfs_ordered_nodes,
                 tbb::concurrent_vector<MAT::Node *> &to_filter
#ifdef CHECK_STATE_REASSIGN
                 ,
                Original_State_t& original_state
#endif
){
    std::unordered_set<size_t> deleted_nodes;
    std::vector<MAT::Node*> changed_nodes;
    for (auto move : all_moves) {
        move_node(move->src, move->get_dst(),move->mutations_relative_to_parent,move->score_change,t,changed_nodes);
    }
    for(auto node:changed_nodes){
        if (deleted_nodes.count((size_t) node)) {
            continue;
        }
        clean_up_src_par(node, &t, deleted_nodes);
        merge_children(node, t, deleted_nodes);
    }
    tbb::concurrent_vector<MAT::Node *> to_filter_new;
    to_filter_new.reserve(to_filter.size());
    for (const auto node : to_filter) {
        if (!deleted_nodes.count((size_t)node) ){
            to_filter_new.push_back(node);
        }
    }
    to_filter=std::move(to_filter_new);
    for (auto node : deleted_nodes) {
        delete (MAT::Node*)node;
    }
#ifdef CHECK_STATE_REASSIGN
    check_samples(t.root, original_state, &t);
    check_clean(t.breadth_first_expansion());
#endif
}
void remove_nodes(MAT::Node* const root,size_t idx,MAT::Tree* tree){
    if (!(root->is_leaf()||root->is_root())) {
    if (root->mutations.empty()) {
        root->parent->children[idx]=root->children[0];
        root->children[0]->parent=root->parent;
        for(size_t idx=1;idx<root->children.size();idx++){
            root->children[idx]->parent=root->parent;
            root->parent->children.push_back(root->children[idx]);
        }
        tree->all_nodes.erase(root->identifier);
        remove_nodes(root->children[0], idx,tree);
        delete root;
        return;
    }
    if (root->children.size()==1) {
        std::vector<MAT::Mutation> merged;
        merge_back(root->children[0]->mutations, root->mutations, merged);
        root->children[0]->mutations=std::move(merged);
        root->children[0]->parent=root->parent;
        root->parent->children[idx]=root->children[0];
        tree->all_nodes.erase(root->identifier);
        delete root;
        return;
    }
    }
    for(size_t idx=0;idx<root->children.size();idx++){
        remove_nodes(root->children[idx],idx,tree);
    }
}
void check_clean(const std::vector<MAT::Node *> &bfs_ordered_nodes){
    for (const auto node : bfs_ordered_nodes) {
        if (node->is_leaf()||node->is_root()) {
            continue;
        }
        if (node->mutations.empty()) {
            fprintf(stderr, "node %zu no mutation\n",node);
            assert(false);
        }
        if (node->children.size()==1) {
            fprintf(stderr, "node %zu no one children\n",node);
            assert(false);
        }
    }
}
void merge_down(const std::vector<MAT::Mutation>& upper,const std::vector<MAT::Mutation>& lower,std::vector<MAT::Mutation>& out){
    auto lower_iter=lower.begin();
    auto lower_end=lower.end();
    for (const auto&  el:upper) {
        while (lower_iter!=lower.end()&&lower_iter->get_position()<el.get_position()) {
            out.push_back(lower_iter->invert());
            lower_iter++;
        }
        if (lower_iter!=lower_end&&lower_iter->get_position()==el.get_position()) {
            if (lower_iter->get_mut_one_hot()!=el.get_mut_one_hot()) {
                out.push_back(el);
                out.back().set_par_mut(lower_iter->get_mut_one_hot(), el.get_mut_one_hot());
            }
            lower_iter++;
        }else {
            out.push_back(el);
        }
    }
    while (lower_iter!=lower.end()) {
        out.push_back(lower_iter->invert());
        lower_iter++;
    }
}
void individual_move(MAT::Node* src,MAT::Node* dst,MAT::Node* LCA,output_t& out){
    if (dst==LCA) {
        LCA=LCA->parent;
    }
    std::vector<MAT::Mutation> mutations(src->mutations);
    MAT::Node* node=src->parent;
    std::vector<MAT::Node*> src_to_LCA;
    while (node!=LCA) {
        std::vector<MAT::Mutation> mutations_new;
        merge_back(mutations,node->mutations ,mutations_new);
        mutations=std::move(mutations_new);
        src_to_LCA.push_back(node);
        node=node->parent;
    }
    assert(src_to_LCA.empty()||src_to_LCA.back()!=LCA);
    node=dst;
    std::vector<MAT::Node*> dst_to_LCA;
    while (node!=LCA) {
        dst_to_LCA.push_back(node);
        node=node->parent;
    }
    for(auto iter=dst_to_LCA.rbegin();iter<dst_to_LCA.rend()-1;iter++){
        std::vector<MAT::Mutation> mutations_new;
        merge_down(mutations, (*iter)->mutations, mutations_new);
        mutations=std::move(mutations_new);
    }
    size_t new_par_score=single_move_last_stage(dst,mutations);
    if (new_par_score<src->mutations.size()) {
        auto new_move=new Profitable_Moves(new_par_score-src->mutations.size(),LCA,dst,mutations);
        new_move->src_to_LCA=std::move(src_to_LCA);
        new_move->src=src;
        new_move->dst_to_LCA=std::move(dst_to_LCA);
        if ((!new_move->src_to_LCA.empty())&&new_move->src_to_LCA.front()==new_move->dst_to_LCA.back()) {
            new_move->src_to_LCA.erase(new_move->src_to_LCA.begin());
        }
        if ((!new_move->src_to_LCA.empty())&&new_move->src_to_LCA.back()==new_move->dst_to_LCA.back()) {
            new_move->src_to_LCA.pop_back();
        }
        std::unordered_set<MAT::Node*> nodes_in_path;
        new_move->apply_nodes([&nodes_in_path](MAT::Node* node){
            assert(nodes_in_path.insert(node).second);
        });
        out.push_back(new_move);
    }
}
