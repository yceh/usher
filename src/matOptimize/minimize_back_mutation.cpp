#include "mutation_annotated_tree.hpp"
#include "src/matOptimize/check_samples.hpp"
#include <atomic>
#include <cstddef>
#include <cstdio>
#include <tbb/blocked_range.h>
#include <tbb/flow_graph.h>
#include <tbb/parallel_for_each.h>
#include "Fitch_Sankoff.hpp"
#include <tbb/task.h>
#include <vector>
#include "min_back.hpp"
namespace MAT = Mutation_Annotated_Tree;
std::atomic<size_t> back_mutation_count;
static MAT::Mutations_Collection merge_mutation(const MAT::Mutations_Collection& parent,const MAT::Mutations_Collection& this_node){
    MAT::Mutations_Collection merged;
    merged.reserve(parent.size()+this_node.size());
    auto iter=parent.begin();
    auto end=parent.end();
    for (const auto& this_mut : this_node) {
        while (iter->get_position()<this_mut.get_position()) {
            merged.push_back(*iter);
            iter++;
        }
        if (this_mut.get_position()==iter->get_position()) {
            if (this_mut.get_mut_one_hot()==iter->get_par_one_hot()) {
                back_mutation_count++;
            }else {
                merged.push_back(this_mut);
            }
            iter++;
        }else {
            merged.push_back(this_mut);    
        }
    }
    merged.mutations.insert(merged.end(),iter,end);
    return merged;
}
struct Get_Leaf_State:public tbb::task{
    MAT::Mutations_Collection parent_mutations;
    MAT::Node* node;
    Get_Leaf_State(MAT::Mutations_Collection&& parent_mutations,MAT::Node* node):parent_mutations(parent_mutations),node(node){}
    task * execute() override{
        auto empty=new (allocate_continuation()) tbb::empty_task;
        empty->set_ref_count(node->children.size());
        std::vector<tbb::task*> tasks;
        for (const auto child : node->children) {
            auto merged_mutations=merge_mutation(parent_mutations, child->mutations);
            if (!child->is_leaf()) {
                tasks.push_back(new(empty->allocate_child()) Get_Leaf_State(std::move(merged_mutations),child));
            }
        }
        empty->set_ref_count(tasks.size());
        for(auto task:tasks){
            empty->spawn(*task);
        }
        return tasks.empty()?empty:nullptr;
    }
};
size_t get_back_mutations_count(const MAT::Tree& tree){
    back_mutation_count.store(0);
    MAT::Mutations_Collection root_muts(tree.root->mutations);
    root_muts.push_back(MAT::Mutation(INT_MAX));
    tbb::task::spawn_root_and_wait(*new(tbb::task::allocate_root()) Get_Leaf_State(std::move(root_muts),tree.root));
    return back_mutation_count.load();
}
std::pair<size_t, size_t> Fitch_Sankoff_Minimize_back_mutation(
    const std::vector<backward_pass_range> &child_idx_range,
    const std::vector<forward_pass_range> &parent_idx,
    const Mutation_Annotated_Tree::Mutation &base,
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>
        &output,
    const mutated_t &mutated,std::vector<backward_info>& score_vec,std::vector<uint8_t>& best_allele) ;
thread_local struct {
    std::vector<backward_info> score_vec;
    std::vector<uint8_t> best_allele;
} FS_containers;
int main (int argc, char** argv){
    MAT::Tree tree=MAT::load_mutation_annotated_tree(argv[1]);
    size_t old_back_mutation_count=get_back_mutations_count(tree);
    Original_State_t origin_states;
    check_samples(tree.root, origin_states, &tree);
    size_t old_par=tree.get_parsimony_score();
    fprintf(stderr, "Old parsimony score: %zu, old mutation count %zu\n ",old_par,old_back_mutation_count);
    auto bfs_ordered_nodes=tree.breadth_first_expansion();
    std::vector<backward_pass_range> child_idx_range;
    std::vector<forward_pass_range> parent_idx;
    std::vector<std::pair<MAT::Mutation,tbb::concurrent_vector<std::pair<size_t,char>>>> pos_mutated(MAT::Mutation::refs.size());
    tbb::parallel_for_each(origin_states.begin(),origin_states.end(),[&pos_mutated,tree](const std::pair<std::string, Mutation_Set>& sample_mutations){
        for (const MAT::Mutation &m : sample_mutations.second) {
            pos_mutated[m.get_position()].first=m;
            pos_mutated[m.get_position()].second.emplace_back(tree.get_node(sample_mutations.first)->bfs_index,m.get_all_major_allele());
        }
    });
    Fitch_Sankoff_prep(bfs_ordered_nodes,child_idx_range, parent_idx);
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>
        output(bfs_ordered_nodes.size());
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0,pos_mutated.size()),
        [&output,&child_idx_range,&parent_idx,&pos_mutated](const tbb::blocked_range<size_t>& in) {
            FS_containers.score_vec.resize(child_idx_range.size());
            FS_containers.best_allele.resize(child_idx_range.size());
            //std::vector<backward_info> score_vec(child_idx_range.size());
            //std::vector<uint8_t> best_allele(child_idx_range.size());
            
            for (size_t idx=in.begin(); idx<in.end(); idx++) {
                if (pos_mutated[idx].second.empty()) {
                    continue;
                }
                std::vector<std::pair<long,nuc_one_hot>> mutated_nodes_idx(pos_mutated[idx].second.begin(),pos_mutated[idx].second.end());
                std::sort(mutated_nodes_idx.begin(),mutated_nodes_idx.end(),mutated_t_comparator());
                mutated_nodes_idx.emplace_back(0,0xf);
                Fitch_Sankoff_Minimize_back_mutation(child_idx_range,parent_idx, pos_mutated[idx].first,output, mutated_nodes_idx
                ,FS_containers.score_vec,FS_containers.best_allele);
                //,score_vec,best_allele);
    
            }
        });
    tbb::affinity_partitioner ap;
    //sort and fill
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, bfs_ordered_nodes.size()),
        [&bfs_ordered_nodes, &output](tbb::blocked_range<size_t> r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                const auto &to_refill = output[i];
                bfs_ordered_nodes[i]->mutations.mutations.clear();
                bfs_ordered_nodes[i]->refill(to_refill.begin(), to_refill.end(),
                                             to_refill.size());
            }
        },
        ap);
    check_samples(tree.root, origin_states, &tree);
    auto new_mut_count=tree.get_parsimony_score();
    auto new_back_mutation_count=get_back_mutations_count(tree);
    fprintf(stderr, "New parsimony score: %zu, new mutation count %zu\n ",new_mut_count,new_back_mutation_count);
    MAT::save_mutation_annotated_tree(tree, argv[2]);

}