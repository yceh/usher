#ifndef EXTRACT_FORMATS_H
#define EXTRACT_FORMATS_H

#include "../mutation_annotated_tree.hpp"

namespace MAT = Mutation_Annotated_Tree;

void get_trios(MAT::Tree T, std::string filepath);

void get_parents(Mutation_Annotated_Tree::Tree *T,
                 std::unordered_set<std::string> &need_parents,
                 std::unordered_set<std::string> &all_nodes);

#endif

