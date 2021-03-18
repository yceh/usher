#include "src/mutation_annotated_tree.hpp"
#include "tree_rearrangement_internal.hpp"
#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdio>
#include <deque>
#include <mutex>
#include <string>
#include <tbb/blocked_range.h>
#include "check_samples.hpp"
#include <tbb/flow_graph.h>
#include <thread>
#include <utility>
#include <vector>
#include "src/tree_rearrangement.hpp"
#include "priority_conflict_resolver.hpp"
#include <chrono>
extern uint32_t num_cores;
std::atomic_bool interrupted(false);
typedef std::vector<std::pair<MAT::Node *, MAT::Node *>> Deleted_Map_t;

void find_nodes_with_recurrent_mutations(std::vector<MAT::Node *>& all_nodes, std::vector<MAT::Node *>& output);
static void fix_condensed_nodes(MAT::Tree* tree){
    std::vector<MAT::Node*> nodes_to_fix;
    for(auto iter:tree->all_nodes){
        if (tree->condensed_nodes.count(iter.first)&&(!iter.second->mutations.empty())) {
            nodes_to_fix.push_back(iter.second);
        }
    }
    for(auto node:nodes_to_fix){
        std::string ori_identifier(node->identifier);
        tree->rename_node(ori_identifier, std::to_string(++tree->curr_internal_node));
        tree->create_node(ori_identifier,node);
    }
}
/*
static void push_all_nodes(MAT::Tree* tree,std::vector<MAT::Node*>& nodes){
    nodes.clear();
    for(auto a:tree->all_nodes){
        nodes.push_back(a.second);
    }
}*/
class save_output_thread{
    int output_count;
    std::atomic_bool done;
    std::thread* saving_thread;
    const std::string outname;
    MAT::Tree& to_save;
    void run(){
        output_count++;
        std::string this_out_name(outname);
        this_out_name.append(std::to_string(output_count)).append(".pb");
        MAT::save_mutation_annotated_tree(to_save, this_out_name);
        done.store(true);
    }
    public:
    save_output_thread(MAT::Tree& to_save, const std::string& outname):outname(outname),to_save(to_save){
        output_count=0;
        done.store(true);
    }
    void wait(){
        /*if (!done.load()) {
            saving_thread->join();
            delete saving_thread;
        }*/
    }
    void start(){
        run();
        //done.store(false);
        //saving_thread=new std::thread(&save_output_thread::run,this);
    }

};

void feed_nodes(int radius,std::vector<MAT::Node *> &to_feed,
                         std::deque<MAT::Node*>& out_queue,
                         std::mutex& out_mutex,
                         std::condition_variable& out_pushed,
                         const std::vector<MAT::Node *> &dfs_ordered_nodes) ;
typedef tbb::flow::multifunction_node<char,std::tuple<MAT::Node*>> Seach_Source_Node_t;
struct Search_Source{
    std::deque<MAT::Node*>& buffer;
    std::mutex& queue_mutex;
    std::condition_variable& ready;
    unsigned int conflict_pct_limit;
    unsigned int min_batch_size;
    mutable unsigned int conflicting_count;
    mutable int inflight_left;
    unsigned int inflight_limit;
    mutable size_t all_count;
    std::atomic_bool& all_nodes_done;
    std::atomic_long& states_in_flight;
    std::atomic_long& states_to_calculate;
    Search_Source*& handle;
    mutable size_t processed;
    std::chrono::time_point<std::chrono::steady_clock> start_time;
    Search_Source(std::deque<MAT::Node*>& buff, std::mutex& queue_mutex, std::condition_variable& ready,unsigned int inflight_limit,unsigned int conflict_pct_limit,unsigned int min_batch_size,std::atomic_bool& all_nodes_done,std::atomic_long& states_in_flight,std::atomic_long& states_to_calculate, Search_Source*& handle):buffer(buff),queue_mutex(queue_mutex),ready(ready),conflict_pct_limit(conflict_pct_limit),min_batch_size(min_batch_size),conflicting_count(0),inflight_left(inflight_limit),inflight_limit(inflight_limit),all_count(0),all_nodes_done(all_nodes_done),states_in_flight(states_in_flight),states_to_calculate(states_to_calculate),handle(handle){
        handle=this;
        processed=0;
        start_time=std::chrono::steady_clock::now();
    }
    void reset(){
        fprintf(stderr,"%d conflicting moves over %zu moves before\n",conflicting_count,all_count);
        all_count=0;
        start_time=std::chrono::steady_clock::now();
        processed=0;
        inflight_left=inflight_limit;
        conflicting_count=0;
        states_to_calculate.store(0);
    }
    void operator()(char in,Seach_Source_Node_t::output_ports_type& out) const{
        handle=const_cast<Search_Source*>(this);
        inflight_left++;
        processed++;
        //fprintf(stderr,"%d conflicting moves over %zu moves so far\n",conflicting_count,all_count);
        if(interrupted.load()){
            return;
        }
        if (in & MOVE_FOUND_MASK) {
            all_count++;
            if (!(in & NONE_CONFLICT_MOVE_MASK)) {
                conflicting_count++;
            }
        }
        fprintf(stderr, "state in flight: %ld\n",states_in_flight.load());
        if (states_in_flight.load()>0x8000000||((conflicting_count > min_batch_size) &&
            (all_count * conflict_pct_limit < conflicting_count * 100))) {
            return;
        }
        if(all_nodes_done){
           //fprintf(stderr,"all_nodes_done");
            return;
        }
        while(states_to_calculate.load()<3){
            if(inflight_left<=0){
                return;
            }
            std::unique_lock<std::mutex> lk(queue_mutex);
            while (buffer.empty()) {
                ready.wait(lk);
            }
            if (buffer.front()==nullptr) {
                buffer.pop_front();
                all_nodes_done.store(true);
           	fprintf(stderr,"all_nodes_done\n");
                return;
            }else {
                std::get<0>(out).try_put(buffer.front());
                states_to_calculate.fetch_add(1);
                buffer.pop_front();
                inflight_left--;
                auto time_now=std::chrono::steady_clock::now();
                auto time_elapsed=std::chrono::duration_cast<std::chrono::seconds>(time_now-start_time);
                float curr_rate=(60*(float)processed)/time_elapsed.count();
                auto to_process=buffer.size();
                float time_left=(to_process/curr_rate)/60;
                fprintf(stderr, "Processing %f nodes/minute, %zu nodes left, estimated %f hours left\n",curr_rate,to_process,time_left);
            }
        }
        //fprintf(stderr,"in flight quota exceeded\n");
        
    }
};
void clean_internal_nodes_with_no_mutation(MAT::Node *this_node,
                                           MAT::Tree *this_tree,Deleted_Map_t& deleted_map) {
    if (this_node->is_leaf()||this_node->is_root()) {
        return;
    }
    std::vector<MAT::Node *> &parent_children = this_node->parent->children;
    std::vector<MAT::Node *> this_node_ori_children = this_node->children;

    if (this_node->mutations.empty()) {
        auto iter = std::find(parent_children.begin(), parent_children.end(),
                              this_node);
        assert(iter != parent_children.end());
        parent_children.erase(iter);
        this_tree->all_nodes.erase(this_node->identifier);
        for (MAT::Node *child : this_node_ori_children) {
            child->parent = this_node->parent;
            parent_children.push_back(child);
        }
        deleted_map.emplace_back(this_node,this_node->parent);
        delete this_node;
    }

    for (MAT::Node *child : this_node_ori_children) {
        clean_internal_nodes_with_no_mutation(child, this_tree,deleted_map);
    }
}
void Tree_Rearrangement::refine_trees(MAT::Tree& this_tree,int radius,int min_batch,int conflict_pct,const std::string& out_name) {

        fprintf(stderr, "Before refinement: %zu \n",
                this_tree.get_parsimony_score());
        Original_State_t ori;
        std::vector<MAT::Node*>& to_optimize=this_tree.new_nodes;
        check_samples(this_tree.root, ori,&this_tree);
        Pending_Moves_t pending_moves;
        tbb::concurrent_vector<Profitable_Move*> profitable_moves;
        {
        Deleted_Map_t ignored;
        clean_internal_nodes_with_no_mutation(this_tree.root, &this_tree,ignored);
        }
        Original_State_t copy(ori);
        check_samples(this_tree.root, copy, &this_tree);

        auto dfs_ordered_nodes = this_tree.depth_first_expansion();
        for(auto node:dfs_ordered_nodes){
            node->tree=&this_tree;
        }
        std::atomic_long states_in_flight(0);
        std::atomic_long states_to_calculate(0);
        //building pipeline
        tbb::queuing_rw_mutex mutex;
        tbb::flow::graph search_graph;
        std::deque<MAT::Node*> search_queue;
        std::mutex queue_mutex;
        std::condition_variable queue_ready;
        std::atomic_bool all_nodes_done;
        Search_Source* handle;
        
        Search_Source input_functor(search_queue,queue_mutex,queue_ready,num_cores,conflict_pct,min_batch,all_nodes_done,states_in_flight,states_to_calculate,handle);
        Seach_Source_Node_t input(search_graph,1,input_functor);

        tbb::flow::function_node<MAT::Node*,Possible_Moves*> neighbors_finder(search_graph,tbb::flow::unlimited,Neighbors_Finder(radius));
        tbb::flow::make_edge(input,neighbors_finder);

        Parsimony_Score_Calculator_t parsimony_score_calculator(search_graph,tbb::flow::unlimited,Parsimony_Score_Calculator{ori,dfs_ordered_nodes,states_in_flight,states_to_calculate});
        tbb::flow::make_edge(neighbors_finder,parsimony_score_calculator);
        
        tbb::flow::function_node<Candidate_Moves*,Profitable_Moves_From_One_Source*> profitable_move_enumerator(search_graph,tbb::flow::unlimited,Profitable_Moves_Enumerator{dfs_ordered_nodes,ori,states_in_flight});
        tbb::flow::make_edge(parsimony_score_calculator,profitable_move_enumerator);

        Cross_t potential_crosses;
        int nodes_inside=0;
        tbb::flow::function_node<Profitable_Moves_From_One_Source*,char> conflict_resolver(search_graph,1,Conflict_Resolver{potential_crosses,to_optimize,nodes_inside,states_in_flight});
        tbb::flow::make_edge(profitable_move_enumerator,conflict_resolver);
        tbb::flow::make_edge(conflict_resolver,input);

        bool have_improvement=true;
        int old_pasimony_score=this_tree.get_parsimony_score();
        int moves_found;
        save_output_thread output_saver(this_tree,out_name);
        while (have_improvement) {
            have_improvement = false;
            find_nodes_with_recurrent_mutations(dfs_ordered_nodes, to_optimize);
            fprintf(stderr, "next_round\n");
            // push_all_nodes(&this_tree, to_optimize);
            while (!to_optimize.empty()) {
                input.try_put(0);
                feed_nodes(radius, to_optimize, search_queue, queue_mutex,
                           queue_ready, dfs_ordered_nodes);
                all_nodes_done.store(false);
                to_optimize.clear();
                while (!all_nodes_done.load()) {
                    std::vector<Profitable_Moves_ptr_t> non_conflicting_moves;
                    handle->reset();
                    potential_crosses.clear();

                    input.try_put(0);
                    search_graph.wait_for_all();
                    states_in_flight.store(0);
                    schedule_moves(potential_crosses, non_conflicting_moves,nodes_inside);
                    moves_found = non_conflicting_moves.size();
                    fprintf(stderr, "%zu moves profitable\n",
                            non_conflicting_moves.size());
                    if (!non_conflicting_moves.empty()) {
                        Pending_Moves_t tree_edits;
                        // tbb::parallel_for(tbb::blocked_range<size_t>(0,
                        // non_conflicting_moves.size()),
                        // Move_Executor{dfs_ordered_nodes,this_tree,non_conflicting_moves,tree_edits});
                        output_saver.wait();
                        Move_Executor temp{dfs_ordered_nodes, this_tree,
                                           non_conflicting_moves, tree_edits,
                                           ori};
                        tbb::blocked_range<size_t> range_temp(
                            0, non_conflicting_moves.size());
                        temp(range_temp);
#ifndef NDEBUG
                        for (auto m : non_conflicting_moves) {
                            #ifdef CHECK_LEAK
                            m->destructed=true;
                            #else
                            delete m;
                            #endif
                        }
#endif
                        non_conflicting_moves.clear();
                        // this_tree.finalize();
                        /*
                        tbb::parallel_for_each(
                            tree_edits.begin(), tree_edits.end(),
                            [&this_tree,
                             &ori](const std::pair<MAT::Node *, ConfirmedMove>
                        &in) {
                           });
                        */
                        
                        Deleted_Map_t deleted_map;
                        for (Pending_Moves_t::const_iterator iter =
                                 tree_edits.begin();
                             iter != tree_edits.end(); iter++) {
                            finalize_children(
                                reinterpret_cast<MAT::Node *>(iter->first),
                                const_cast<ConfirmedMove &>(iter->second),
                                &this_tree, ori, deleted_map);
                        }
                        clean_internal_nodes_with_no_mutation(this_tree.root, &this_tree,deleted_map);
                        Original_State_t copy(ori);
                        check_samples(this_tree.root, copy, &this_tree);
                        fix_condensed_nodes(&this_tree);
                        this_tree.reassign_level();
                        //#ifndef NDEBUG
                        output_saver.start();
                        for (MAT::Node *&next_node : to_optimize) {
                            for (const auto &deleted : deleted_map) {
                                if (next_node == deleted.first) {
                                    next_node = deleted.second;
                                }
                            }
                        }
                        for (MAT::Node *&next_node : search_queue) {
                            for (const auto &deleted : deleted_map) {
                                if (next_node == deleted.first) {
                                    next_node = deleted.second;
                                }
                            }
                        }

                        dfs_ordered_nodes = this_tree.depth_first_expansion();

                        int new_score = this_tree.get_parsimony_score();
                        //assert(old_pasimony_score - new_score >= moves_found);
                        fprintf(stderr, "Between refinement: %d \n",new_score);
                        //#endif
                        if(new_score<old_pasimony_score){
                            have_improvement = true;
                        }
                        old_pasimony_score = new_score;
                    }else{
                        fprintf(stderr, "no moves profitable\n");
                    }
                    if(interrupted.load()){
                        goto END;
                    }
                }
                fprintf(stderr, "%zu nodes deferred next round",
                        to_optimize.size());
            }
        }
END:
        check_samples(this_tree.root, ori,&this_tree);
        fprintf(stderr, "After refinement: %zu \n",
                this_tree.get_parsimony_score());
        fix_condensed_nodes(&this_tree);
        this_tree.reassign_level();
}
