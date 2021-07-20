#include "../check_samples.hpp"
#include "../priority_conflict_resolver.hpp"
#include "placement.hpp"
#include "src/new_tree_rearrangements/mutation_annotated_tree.hpp"
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#include <utility>
#include <vector>
#ifdef PROFILE
std::vector<std::atomic<size_t>> stoped_radius;
#endif
//see whether within radius of root have nodes that have same mutation as root and have changed in previous iteration (or this is the first iteration)
static void check_changed_neighbor(int radius, MAT::Node* root, MAT::Node* exclude,bool& found,bool is_first,const std::unordered_set<int>& pos){
    if (root->changed) {
        found=true;
        return;
    }
    if (!radius||found) {
        return;
    }
    if (root->parent&&root->parent!=exclude) {
        check_changed_neighbor(radius-1, root->parent, root, found,is_first,pos);
    }
    for(auto node:root->children){
        if (found) {
            return;
        }
        if(node!=exclude){
            check_changed_neighbor(radius-1, node, root, found,is_first,pos);
        }
    }
}
//find src node to search
void find_nodes_to_move(const std::vector<MAT::Node *> &bfs_ordered_nodes,
                        tbb::concurrent_vector<MAT::Node *> &output,
                        bool is_first, int radius) {
    auto start=std::chrono::steady_clock::now();
    if (is_first) {
        output=tbb::concurrent_vector<MAT::Node*>(bfs_ordered_nodes.begin(),bfs_ordered_nodes.end());
    }else{
    output.reserve(bfs_ordered_nodes.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0,bfs_ordered_nodes.size()),[&bfs_ordered_nodes,is_first,radius,&output](const tbb::blocked_range<size_t>& r){
        for (size_t idx=r.begin(); idx<r.end(); idx++) {
            auto node=bfs_ordered_nodes[idx];
            std::unordered_set<int> pos;
            pos.reserve(node->mutations.size()*4);
            for(const auto& mut:node->mutations){
                pos.insert(mut.get_position());
            }
            bool found=false;
            check_changed_neighbor(radius, node, nullptr, found, is_first, pos);
            if (found) {
                output.push_back(node);
            }
        }
    });
    }
    std::minstd_rand g(std::chrono::steady_clock::now().time_since_epoch().count());
    std::shuffle(output.begin(), output.end(), g);
    fprintf(stderr, "Will search %f of nodes\n",(double)output.size()/(double)bfs_ordered_nodes.size());
    std::chrono::duration<double> elapsed_seconds = std::chrono::steady_clock::now()-start;
    fprintf(stderr, "Took %f s to find nodes to move\n",elapsed_seconds.count());
}
unsigned int search_radius;
MAT::Node* get_LCA(MAT::Node* src,MAT::Node* dst);
int main(int argc, char **argv) {
    MAT::Tree tree = MAT::load_mutation_annotated_tree(argv[1]);
    search_radius=atoi(argv[2]);
    stoped_radius=std::vector<std::atomic<size_t>>(search_radius+1);
    Original_State_t ori_state;
    check_samples(tree.root, ori_state, &tree);
    remove_nodes(tree.root,0,&tree);
    check_samples(tree.root, ori_state, &tree);
    FILE *log = fopen("/dev/null", "w");
    auto bfs_ordered_nodes = tree.breadth_first_expansion();
    /*auto src=tree.get_node("cov-100001875387|OD986787.1|2021-03-12");
     output_t out;
    find_place(src, out, 5);*/
    tbb::concurrent_vector<MAT::Node *> nodes_to_search;
    find_nodes_to_move(bfs_ordered_nodes, nodes_to_search, true, search_radius);
    size_t last_parsimony_score=tree.get_parsimony_score();
    while(true){
    for (auto node : bfs_ordered_nodes) {
        node->changed=false;
    }
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
    find_nodes_to_move(bfs_ordered_nodes, nodes_to_search, false, search_radius);
    }
    save_final_tree(tree, ori_state, "out.pb");
    tree.delete_nodes();
#ifdef PROFILE
    for (int i=0; i<=search_radius; i++) {
        fprintf(stderr, "%zu\n",stoped_radius[i].load());
    }
#endif
}