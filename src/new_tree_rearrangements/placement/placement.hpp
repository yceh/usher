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