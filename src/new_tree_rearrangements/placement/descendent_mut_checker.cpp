#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <cstdio>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <utility>
#include <vector>
namespace MAT=Mutation_Annotated_Tree;
void check_descendent(MAT::Tree& tree){
    struct descendant_allele_checker{
        const MAT::Tree& tree;
        std::pair<const MAT::Mutation*,unsigned char> check_pos(int pos,const MAT::Node* start_node, char par_state)const{
            const MAT::Mutation* pos_found=nullptr;
            auto this_state=par_state;
            for (const auto& mut : start_node->mutations) {
                if (mut.get_position()==pos) {
                    this_state=mut.get_mut_one_hot();
                    pos_found=&mut;
                }
            }
            std::unordered_map<char,std::vector<const MAT::Mutation*>> children_states;
            char tip_set=0;
            for (const auto child : start_node->children) {
                auto res=check_pos(pos, child, this_state);
                auto insert_res=children_states.emplace(res.second&0xf,std::vector<const MAT::Mutation*>());
                insert_res.first->second.push_back(res.first);
                tip_set|=(res.second>>4);
            }
            if (pos_found&&pos_found->is_mut_tip()) {
                assert(!(pos_found->get_mut_one_hot()&tip_set));
                tip_set|=pos_found->get_mut_one_hot();
            }
            char descendant_states_children=0;
            for (auto& p : children_states) {
                descendant_states_children|=p.first;
            }
            auto descendant_states=descendant_states_children|this_state;
            assert((descendant_states|tip_set)==descendant_states);
            if (pos_found) {
                assert(pos_found->get_descendent_mut()==(descendant_states));
            }
            for (const auto& p : children_states) {
                if (p.second.size()>1) {
                    for (const auto& mut_p : p.second) {
                        if (mut_p) {
                            assert(mut_p->get_sibling_mut()==descendant_states_children);
                        }
                    }
                }else {
                    assert(p.second.size()==1);
                    char other_state=0;
                    for (const auto& p2 : children_states) {
                         if(p2.first!=p.first){
                             other_state|=p2.first;
                         }
                    }
                    if(p.second[0]){
                        if (p.second[0]->get_sibling_mut()!=other_state&&other_state) {
                            fprintf(stderr, "At position %d of %zu, sibling muts expect %d but get %d, children table at %zu\n",pos,start_node,other_state,p.second[0]->get_sibling_mut().get_nuc_no_check(),&children_states);
                            assert(false);
                        }
                    }
                }
            }
            return std::make_pair(pos_found, (tip_set<<4)|descendant_states);
        }
        void operator()(tbb::blocked_range<size_t>& range) const{
            for (auto pos=range.begin(); pos<range.end(); pos++) {
                if (MAT::Mutation::refs[pos]) {
                    auto ret=check_pos(pos,tree.root,MAT::Mutation::refs[pos]);
                    //assert((ret.second>>4|MAT::Mutation::refs[pos])==(ret.second&0xf));
                }
            }
        }
    };
    tbb::parallel_for(tbb::blocked_range<size_t>(0,MAT::Mutation::refs.size()),descendant_allele_checker{tree});
}