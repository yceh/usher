#ifndef AMPLICON_H
#define AMPLICON_H

#include "src/mutation_annotated_tree.hpp"
#include "src/ripplesUtils/text_parser.hpp"
#include "src/usher_graph.hpp"

namespace MAT = Mutation_Annotated_Tree;

void get_starting_coordinates(std::string sam_file_path,
                              std::vector<int> &start_coordinates);
inline uint64_t str_view_to_uint64(std::string_view str) noexcept;

#endif

