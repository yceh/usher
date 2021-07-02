#include "mutation_annotated_tree.hpp"
#include "tbb/task.h"
#include <atomic>
#include <vector>
//find mutated alleles shared by both of the input
static void intersect_allele(const Mutation_Annotated_Tree::Mutations_Collection& in1, const Mutation_Annotated_Tree::Mutations_Collection& in2, Mutation_Annotated_Tree::Mutations_Collection& out){
    if (in1.empty()||in2.empty()) {
        return;
    }
    auto in1_iter=in1.begin();
    auto in1_end=in1.end();
    for(const auto& mut:in2){
        while (in1_iter!=in1_end&&in1_iter->get_position()<mut.get_position()) {
            in1_iter++;
            //skip all muts present in only one of the inout
        }
        if (in1_iter!=in1_end&&in1_iter->get_position()==mut.get_position()){
            //output common major allele
            nuc_one_hot common=in1_iter->get_all_major_allele()&mut.get_all_major_allele();
            if (common!=mut.get_par_one_hot()) {
                out.push_back(mut);
                out.back().set_auxillary(common, 0xf&(~common));
            }
            in1_iter++;
        }
    }
}
//merge sortish way to find major alleles shared by all nodes from begin to end
static void get_common_mutations(std::vector<Mutation_Annotated_Tree::Node*>::const_iterator begin,std::vector<Mutation_Annotated_Tree::Node*>::const_iterator end,Mutation_Annotated_Tree::Mutations_Collection& out){
    if (end-begin==3) {
        Mutation_Annotated_Tree::Mutations_Collection temp;
        intersect_allele((*begin)->mutations, (*(begin+1))->mutations, temp);
        intersect_allele(temp, (*(begin+2))->mutations, out);
    } else if (end-begin==2) {
        intersect_allele((*begin)->mutations, (*(begin+1))->mutations, out);
    } else {
        Mutation_Annotated_Tree::Mutations_Collection temp1;
        Mutation_Annotated_Tree::Mutations_Collection temp2;
        auto middle=begin+(end-begin)/2;
        get_common_mutations(begin,middle , temp1);
        get_common_mutations(middle,end,temp2);
        intersect_allele(temp1,temp2, out);
    }
}
//Condense nodes, this parallelization is more 
//for getting around stack overflow than performance, 
//as condensing nodes doesn't take too long anyway
struct Node_Condenser:public tbb::task{
    Mutation_Annotated_Tree::Tree::condensed_node_t & condensed_nodes;
    Mutation_Annotated_Tree::Node* root;
    std::atomic<size_t>& condensed_nodes_count;
    Node_Condenser(Mutation_Annotated_Tree::Tree::condensed_node_t & removed_nodes,Mutation_Annotated_Tree::Node* root,std::atomic<size_t>& condensed_nodes_count):condensed_nodes(removed_nodes),root(root),condensed_nodes_count(condensed_nodes_count){}
    tbb::task* execute() override{
        //condensed children of root
        std::vector<Mutation_Annotated_Tree::Node*> polytomy_nodes;
        //Task for each non-leaf node 
        std::vector<tbb::task*> to_spawn;
        //continuation
        tbb::empty_task* empty=new(allocate_continuation()) tbb::empty_task();
        //instead of removing,just copy children over 
        std::vector<Mutation_Annotated_Tree::Node*> new_children;
        new_children.reserve(root->children.size());
        for(auto node:root->children){
            if (node->is_leaf()) {
                if (node->no_valid_mutation()) {
                    //have no valid mutation, so condense it
                    polytomy_nodes.push_back(node);
                }else {
                    new_children.push_back(node);
                }
            }else {
                // a new task for non-leaf child
                to_spawn.push_back(new (
                    empty->allocate_child()) Node_Condenser(condensed_nodes,node,condensed_nodes_count));
                new_children.push_back(node);
            }
        }
        empty->set_ref_count(to_spawn.size());
        for(auto task:to_spawn){
            empty->spawn(*task);
        }
        if (polytomy_nodes.size()>1) {
            Mutation_Annotated_Tree::Mutations_Collection new_shared_muts;
            //get all ambiguous mutation, to place on the condensed node
            get_common_mutations(polytomy_nodes.begin(), polytomy_nodes.end(), new_shared_muts);
#ifdef CHECK_CONDENSE
            for(const auto node:polytomy_nodes){
                for(const auto& mut:new_shared_muts){
                    auto iter=node->mutations.find(mut);
                    assert(iter!=node->mutations.end());
                    assert(iter->get_all_major_allele()&mut.get_all_major_allele());
                }
            }
#endif
            //replace the first node to condense with the condensed node in place
            std::vector<std::string> children_id{std::move(polytomy_nodes[0]->identifier)};
            children_id.reserve(polytomy_nodes.size());
            new_children.push_back(polytomy_nodes[0]);
            //get a new name for condensed node
            polytomy_nodes[0]->identifier="node_" + std::to_string(++condensed_nodes_count) + "_condensed_" + std::to_string(polytomy_nodes.size()) + "_leaves";
            polytomy_nodes[0]->mutations.swap(new_shared_muts);
            //record original sample name
            for (size_t replaced_child_idx=1; replaced_child_idx<polytomy_nodes.size(); replaced_child_idx++) {
                children_id.push_back(std::move(polytomy_nodes[replaced_child_idx]->identifier));
                delete polytomy_nodes[replaced_child_idx];
            }
            //assert(res.second);
            root->children.swap(new_children);
            std::vector<std::string> temp;
            auto condensed_node_insert_result=condensed_nodes.emplace(polytomy_nodes[0]->identifier,temp);
            //retry if failed to insert (shouldn't happen)
            while (!condensed_node_insert_result.second) {
                polytomy_nodes[0]->identifier="node_" + std::to_string(++condensed_nodes_count) + "_condensed_" + std::to_string(polytomy_nodes.size()) + "_leaves";
                condensed_node_insert_result=condensed_nodes.emplace(polytomy_nodes[0]->identifier,temp);
            }
            condensed_node_insert_result.first->second=std::move(children_id);
        }
        //scheduler bypass to run continuation task if no task spawned
        return to_spawn.empty()?empty:nullptr;
    }
};
static Mutation_Annotated_Tree::Node* update_node_map(const std::pair<std::string,std::vector<std::string>>& condensed,Mutation_Annotated_Tree::Tree& tree){
        Mutation_Annotated_Tree::Node* ori_node=tree.get_node(condensed.second[0]);
        //update all node, as the first condensed children is replaced with resulting condensed node in place, 
        //can get the pointer to condensed node with the name of first condensed children
        tree.all_nodes.emplace(ori_node->identifier,ori_node);
        for(auto node_id:condensed.second){
            tree.all_nodes.erase(node_id);
        }
        return ori_node;
}
void Mutation_Annotated_Tree::Tree::condense_leaves() {
    bool contain_condensed_nodes=condensed_nodes.size() > 0;
    std::atomic<size_t> condensed_nodes_count(condensed_nodes.size());
    if (contain_condensed_nodes) {
        Mutation_Annotated_Tree::Tree::condensed_node_t temp;
        tbb::task::spawn_root_and_wait(*new(tbb::task::allocate_root()) Node_Condenser(temp,root,condensed_nodes_count));
        for(auto& condensed:temp){
            auto ori_node=update_node_map(condensed,*this);
            for (const auto& removed_node : condensed.second) {
                auto iter=condensed_nodes.find(removed_node);
                if (iter!=condensed_nodes.end()) {
                    condensed.second.insert(condensed.second.end(),iter->second.begin(),iter->second.end());
                    iter->second.clear();
                }
                condensed_nodes.unsafe_erase(iter);
            }
        }
    }else{
    tbb::task::spawn_root_and_wait(*new(tbb::task::allocate_root()) Node_Condenser(condensed_nodes,root,condensed_nodes_count));
    //assert(condensed_nodes_count.load()==condensed_nodes.size());
    for(const auto& condensed:condensed_nodes){
        update_node_map(condensed,*this);
    }
    }
}