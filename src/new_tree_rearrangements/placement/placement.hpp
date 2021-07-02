#include "../mutation_annotated_tree.hpp"
#include <atomic>
#include <tbb/queuing_rw_mutex.h>
#include <tbb/tbb_allocator.h>
namespace MAT = Mutation_Annotated_Tree;
struct placed_sample_info {
    std::vector<MAT::Valid_Mutation> mutations_relative_to_this_node;
    std::atomic<bool> ready;
    std::string sample_name;
    placed_sample_info(placed_sample_info&& other):mutations_relative_to_this_node(std::move(other.mutations_relative_to_this_node)),ready(true),sample_name(std::move(other.sample_name)){}
    placed_sample_info(const placed_sample_info& other):mutations_relative_to_this_node(std::move(other.mutations_relative_to_this_node)),ready(true),sample_name(std::move(other.sample_name)){}
    placed_sample_info(
        std::string& sample_name,
        std::vector<MAT::Valid_Mutation>& mutations_relative_to_this_node)
        :mutations_relative_to_this_node(std::move(mutations_relative_to_this_node)),ready(true),sample_name(std::move(sample_name)) {
        }
};
//,tbb::zero_allocator<placed_sample_info>
typedef tbb::concurrent_vector<placed_sample_info> first_pass_per_node_t;
void place_first_pass(
    std::string &sample_name, std::vector<MAT::Valid_Mutation> &mutations,
    MAT::Tree *tree,
    std::vector<first_pass_per_node_t> &output); 
void process_sample_serial(MAT::Node* node,
  first_pass_per_node_t& samples,
  tbb::concurrent_vector<MAT::Node*>& nodes_need_identifier,std::vector<tbb::queuing_rw_mutex>& mutexes,tbb::queuing_rw_mutex& root_mutex);
void placement_prep(MAT::Tree *t);
