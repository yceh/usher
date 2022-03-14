#ifndef AMPLICON_H
#define AMPLICON_H

#include "src/mutation_annotated_tree.hpp"
#include "src/ripplesUtils/text_parser.hpp"
#include "src/usher_graph.hpp"

namespace MAT = Mutation_Annotated_Tree;

void get_coordinates(
    std::string sam_file_path,
    std::vector<std::tuple<int, int>> &samples_start_end_coordinates);

inline uint64_t str_view_to_uint64(std::string_view str) noexcept;

void placement(MAT::Tree *T, std::string &vcf_filename,
               std::vector<Missing_Sample> &missing_samples,
							 size_t total_nodes,
               std::vector<std::tuple<int, int>> samples_start_end_coordinates,
               bool create_new_mat);

#endif

