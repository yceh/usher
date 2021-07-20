#include <atomic>
#include <cstddef>
#include <tuple>
#include <utility>
#include <vector>
#include "placement.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#ifdef PROFILE
extern std::vector<std::atomic<size_t>> stoped_radius;
#endif
typedef std::vector<Profitable_Moves_ptr_t> output_t;
/**
 * @param[in] node
 * @param[out] lower_bound lower_bound of parsimony possible if placed at subtree rooted at node
 * @param[out] parsimony_score_if_split_here parsimony score if placed between node and its parent
 * @param[in] par_mut_iter Mutations if placed as child of parent of node
 * @param mutations_merged Mutation relative to node, for trying children of node
 * @return Whether there are mutations splitted under the new branching node, if not, for next level of recursion
 */
std::pair<unsigned int,unsigned int> merged_mutation_if_continue_down(
    const MAT::Node *node,
    const std::vector<MAT::Mutation>& par_mut,
    std::vector<MAT::Mutation>& mutations_merged) {
    auto par_mut_iter=par_mut.begin();
    mutations_merged.reserve(node->mutations.size()+par_mut.size());
    unsigned int parsimony_score_if_split_here=0;
    unsigned int lower_bound=0;
    for (const auto &mut : node->mutations) {
        while (par_mut_iter!=par_mut.end() && (par_mut_iter->get_position() < mut.get_position())) {
            // mutations needed if placed at parent and no change
            mutations_merged.push_back(*par_mut_iter);
            assert(!(par_mut_iter->get_mut_one_hot()&par_mut_iter->get_par_one_hot()));
            parsimony_score_if_split_here++;
            // if none of the descendent mutation can mutate to the desired
            // allele, increment lower bound
            if (!(par_mut_iter->get_descendent_mut() &
                  par_mut_iter->get_mut_one_hot())) {
                lower_bound++;
            }
            par_mut_iter++;
        }
        mutations_merged.push_back(mut);
        auto &inserted = mutations_merged.back();
        if (par_mut_iter!=par_mut.end() && (par_mut_iter->get_position() == mut.get_position())) {
            assert(par_mut_iter->get_par_one_hot()==mut.get_par_one_hot());
            //assert((par_mut_iter->get_descendent_mut()&par_mut_iter->get_mut_one_hot())==(par_mut_iter->get_descendent_mut()&(mut.get_descendent_mut()|mut.get_mut_one_hot())&par_mut_iter->get_mut_one_hot()));
            if (par_mut_iter->get_mut_one_hot() == mut.get_mut_one_hot()) {
                // coincide, and a mutation at node consumed that mutation
                assert(par_mut_iter->get_descendent_mut()&mut.get_mut_one_hot());
                mutations_merged.pop_back();
                par_mut_iter++;
                continue;
            } else {
                // coincide, but split
                parsimony_score_if_split_here++;
                inserted.set_par_mut(mut.get_mut_one_hot(),par_mut_iter->get_mut_one_hot());
                par_mut_iter++;
            }
        } else {
            //back mutation
            inserted.set_par_mut(mut.get_mut_one_hot(),mut.get_par_one_hot());
        }
        assert(inserted.get_mut_one_hot()!=inserted.get_par_one_hot());
        if (!(inserted.get_mut_one_hot() &
              inserted.get_descendent_mut())) {
            lower_bound++;
        }
    }
    while (par_mut_iter!=par_mut.end()) {
        mutations_merged.push_back(*par_mut_iter);
        assert(!(par_mut_iter->get_mut_one_hot()&par_mut_iter->get_par_one_hot()));
        parsimony_score_if_split_here++;
        // if none of the descendent mutation can mutate to the desired
        // allele, increment lower bound
        if (!(par_mut_iter->get_descendent_mut() &
              par_mut_iter->get_mut_one_hot())) {
            lower_bound++;
        }
        par_mut_iter++;
    }
    return std::make_pair(parsimony_score_if_split_here, lower_bound);
}



static std::tuple<size_t,size_t,size_t> upward( std::vector<MAT::Mutation>& par_mut, MAT::Node *node,std::vector<MAT::Mutation> & mutations_merged){
    size_t lower_bound=0;
    auto par_mut_iter=par_mut.begin();
    mutations_merged.reserve(par_mut.size()+node->mutations.size());
    size_t parsimony_score_if_split_here=0;
    size_t side_way_lower_bound=0;
    for (const auto &mut : node->mutations) {
        while (par_mut_iter!=par_mut.end() && (par_mut_iter->get_position() < mut.get_position())) {
            // mutations needed if placed at parent and no change
            mutations_merged.push_back(*par_mut_iter);
            if (mutations_merged.back().decrement_tip_distance()) {
                lower_bound++;
                side_way_lower_bound++;
            }
            assert(!(par_mut_iter->get_mut_one_hot()&par_mut_iter->get_par_one_hot()));
            parsimony_score_if_split_here++;
            // Already covered by lower bound of going up, NOP
            par_mut_iter++;
        }
        if (par_mut_iter!=par_mut.end() && (par_mut_iter->get_position() == mut.get_position())) {
            assert(par_mut_iter->get_par_one_hot()==mut.get_mut_one_hot());
            if (par_mut_iter->get_mut_one_hot() == mut.get_par_one_hot()) {
                assert(!par_mut_iter->is_last_mut());
                // coincide, and a mutation at node consumed that mutation
                par_mut_iter++;
                continue;
            } else {
                // coincide, but split
                parsimony_score_if_split_here++;
                //mut nuc not the same, so not change whether it is tip or descendent alleles
                mutations_merged.push_back(*par_mut_iter);
                auto& inserted=mutations_merged.back();
                inserted.set_par_mut(mut.get_par_one_hot(),par_mut_iter->get_mut_one_hot());
                //If the last pushed added mut is already the tip carraying the allele, it have been accounted in the lower bound of going up, if not check the sibling state set to see whether it worth to go side way.
                if(!(mut.get_sibling_mut()&par_mut_iter->get_mut_one_hot())){
                    side_way_lower_bound++;
                }
                par_mut_iter++;
            }
        } else {
            mutations_merged.push_back(mut);
            auto& inserted=mutations_merged.back();
            //It may be possible that it is a tip, yet some sibling have that state, so not increment sideway lower bound so fast
            if(!(mut.get_mut_one_hot()&mut.get_sibling_mut())){
                side_way_lower_bound++;
            }
            inserted.set_descendent(0xf);
        }
        assert(mutations_merged.back().get_mut_one_hot()!=mutations_merged.back().get_par_one_hot());
    }
    while (par_mut_iter!=par_mut.end()) {
        mutations_merged.push_back(*par_mut_iter);
        assert(!(par_mut_iter->get_mut_one_hot()&par_mut_iter->get_par_one_hot()));
        if (mutations_merged.back().decrement_tip_distance()) {
            lower_bound++;
            side_way_lower_bound++;
        }
        parsimony_score_if_split_here++;
        par_mut_iter++;
    }

    assert(side_way_lower_bound>=lower_bound);
    return std::make_tuple(lower_bound,side_way_lower_bound,parsimony_score_if_split_here);
}
void add_output(int parsimony_score,output_t &output,MAT::Node* dst,MAT::Node* LCA,const std::vector<MAT::Mutation> & mut_relative_to_parent){
    assert(LCA!=dst);
    if (parsimony_score<output[0]->score_change) {
        for(const auto change:output){
            delete change;
        }
        output.clear();
    }
    if (output.empty()||output[0]->score_change==parsimony_score) {
        auto this_move=new Profitable_Moves(parsimony_score,LCA,dst,mut_relative_to_parent);
        output.push_back(this_move);
    }
}
static int place_first_pass_helper(
    output_t &output,
    const std::vector<MAT::Mutation> &mutations, MAT::Node *node,
    std::vector<MAT::Mutation> &mutations_merged, MAT::Node* LCA
    #ifndef NDEBUG
    ,size_t old_lower_bound
    #endif
    ) {
    // mutations needed if placed at node
    auto res=merged_mutation_if_continue_down(node,
                                     mutations, mutations_merged);
    auto lower_bound =
        res.second; // lower bound of parsimony score if placed at subtree rooted at node
    auto parsimony_score_if_split_here =
        res.first; // parsimony score if placed between node and its parent
    assert(parsimony_score_if_split_here>=old_lower_bound);
    // output
    add_output(parsimony_score_if_split_here, output, node, LCA, mutations);
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
    output_t &output,
    const std::vector<MAT::Mutation> &mutations, MAT::Node *node,MAT::Node* LCA,unsigned int radius_left
    #ifndef NDEBUG
    ,size_t lower_bound_prev
    #endif
    ) {
    std::vector<MAT::Mutation> mutations_merged;
    auto lower_bound = place_first_pass_helper(
        output, mutations, node, mutations_merged,LCA
        #ifndef NDEBUG
        ,lower_bound_prev
        #endif
        );
    // go down if can reduce parsimony score
    if (!radius_left) {
        #ifdef PROFILE
            stoped_radius[radius_left]++;
        #endif
        return;
    }
    #ifndef CHECK_BOUND
    if (lower_bound >= output[0]->score_change) {
        #ifdef PROFILE
            stoped_radius[radius_left]++;
        #endif
        return;
    }
    #endif
    for (auto c : node->children) {
        place_first_pass_serial_each_level(output, mutations_merged, c,LCA,radius_left-1
        #ifndef NDEBUG
        ,lower_bound
        #endif
        );
    }
}
std::tuple<size_t,size_t,size_t> init_mut_vector(MAT::Node* src,std::vector<MAT::Mutation>& out){
    out.reserve(src->mutations.size());
    size_t lower_bound=0;
    size_t sibling_bound=0;
    for (const auto& mut : src->mutations) {
        out.push_back(mut);
        if (!(mut.get_mut_one_hot()&mut.get_sibling_mut())) {
            sibling_bound++;
        }
    }
    return std::make_tuple(lower_bound,sibling_bound,0);
}
void find_place(MAT::Node* src,output_t &output,unsigned int radius_left){
    if (src->is_root()) {
        return;
    }
    #ifdef LITE_DETAIL_SEARCH
    std::vector<MAT::Mutation> ori_src;
    get_mutations_relative_to_root(src,ori_src);
    #endif
    std::vector<MAT::Mutation> mutations;
    auto bounds=init_mut_vector(src, mutations);
    MAT::Node* this_node=src->parent;
    MAT::Node* last_node=src;
    std::vector<MAT::Mutation> mutations_next;
    auto dummy=new Profitable_Moves;
    dummy->score_change=src->mutations.size();
    output.push_back(dummy);
    while (this_node&&radius_left>0) {
        #ifndef CHECK_BOUND
        if (std::get<1>(bounds)<output[0]->score_change) {
        #endif
            for (auto c : this_node->children) {
                if (c!=last_node) {
                    place_first_pass_serial_each_level(output, mutations, c,this_node,radius_left
        #ifndef NDEBUG
                    ,std::get<1>(bounds)
        #endif
                    );
                }
            }
        #ifndef CHECK_BOUND
        #ifdef PROFILE
        }else {
            stoped_radius[radius_left]++;
        #endif
        }
        #endif
        mutations_next.clear();
        auto old_lower_bound=std::get<0>(bounds);
        bounds=upward(mutations, this_node,mutations_next);
        assert(std::get<2>(bounds)>=old_lower_bound);
        add_output(std::get<2>(bounds), output, this_node, this_node->parent, mutations_next);
    #ifndef CHECK_BOUND
        if(std::get<0>(bounds)>output[0]->score_change){
            #ifdef PROFILE
                stoped_radius[radius_left]++;
            #endif
            break;
        }
    #endif
        mutations=std::move(mutations_next);
        last_node=this_node;
        this_node=this_node->parent;
        #ifdef LITE_DETAIL_SEARCH
        std::vector<MAT::Mutation> new_muts(mutations);
        get_mutations_relative_to_root(this_node, new_muts);
        compare_mutations(ori_src,new_muts);
        #endif
        radius_left--;
    }
    #ifdef PROFILE
    if(!radius_left){
        stoped_radius[radius_left]++;
    }
    #endif
    if (output[0]->score_change==src->mutations.size()) {
        for (auto move : output) {
            delete move;
        }
        output.clear();
    }else {
        for (auto & out : output) {
            out->src=src;
            out->score_change=out->score_change-src->mutations.size();
            auto node=src->parent;
            while (node!=out->LCA) {
                out->src_to_LCA.push_back(node);
                node=node->parent;
            }
            node=out->get_dst()->parent;
            while (node!=out->LCA) {
                out->dst_to_LCA.push_back(node);
                node=node->parent;
            }
            if ((!out->src_to_LCA.empty())&&out->dst_to_LCA.back()==out->src_to_LCA.front()) {
                out->src_to_LCA.erase(out->src_to_LCA.begin());
            }
            if ((!out->src_to_LCA.empty())&&out->src_to_LCA.back()==out->dst_to_LCA.back()) {
            out->src_to_LCA.pop_back();
        }
        }
    }
}