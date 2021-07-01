#include "../mutation_annotated_tree.hpp"
namespace MAT = Mutation_Annotated_Tree;
struct placed_sample_info {
    std::string sample_name;
    std::vector<MAT::Valid_Mutation> mutations_relative_to_this_node;
    placed_sample_info(
        std::string sample_name,
        std::vector<MAT::Valid_Mutation> mutations_relative_to_this_node)
        : sample_name(sample_name),
          mutations_relative_to_this_node(mutations_relative_to_this_node) {}
};
typedef tbb::concurrent_vector<placed_sample_info> first_pass_per_node_t;
void place_first_pass(
    std::string &sample_name, std::vector<MAT::Valid_Mutation> &mutations,
    MAT::Tree *tree,
    std::vector<first_pass_per_node_t> &output,
    first_pass_per_node_t &above_root); 
void process_sample_serial(MAT::Node* node,
  first_pass_per_node_t& samples,
  tbb::concurrent_vector<MAT::Node*>& nodes_need_identifier);
void placement_prep(MAT::Tree *t);
