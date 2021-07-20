#include "../check_samples.hpp"
#include "../priority_conflict_resolver.hpp"
#include "placement.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <cstdio>
#include <cstdlib>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#include <utility>
#include <vector>
unsigned int search_radius;
MAT::Node* get_LCA(MAT::Node* src,MAT::Node* dst);
int main(int argc, char **argv) {
    MAT::Tree tree = MAT::load_mutation_annotated_tree(argv[1]);
    search_radius=atoi(argv[2]);
    Original_State_t ori_state;
    check_samples(tree.root, ori_state, &tree);
    remove_nodes(tree.root,0,&tree);
    check_samples(tree.root, ori_state, &tree);
    FILE *log = fopen("/dev/null", "w");
    auto bfs_ordered_nodes = tree.breadth_first_expansion();
    /*auto src=tree.get_node("Hangzhou_ZJU-01_2020");
    auto dst=tree.get_node("49");
     output_t out;
    individual_move(src, dst, get_LCA(src, dst), out);*/
    tbb::concurrent_vector<MAT::Node *> nodes_to_search(
        bfs_ordered_nodes.begin(), bfs_ordered_nodes.end());
    size_t last_parsimony_score=tree.get_parsimony_score();
    while(true){
    size_t par_score;
    while (!nodes_to_search.empty()) {
        par_score=optimize_tree(bfs_ordered_nodes, nodes_to_search, tree, search_radius, log
#ifndef NDEBUG
                      ,
                      ori_state
#endif
        );
        check_samples(tree.root, ori_state, &tree);
        fprintf(stderr, "%zu\n",par_score);
        bfs_ordered_nodes=tree.breadth_first_expansion();
    }
    if (par_score<last_parsimony_score) {
        last_parsimony_score=par_score;
    }else{
        break;
    }
    nodes_to_search=tbb::concurrent_vector<MAT::Node *> (bfs_ordered_nodes.begin(), bfs_ordered_nodes.end());
    }
    save_final_tree(tree, ori_state, "out.pb");
    tree.delete_nodes();
}