#include "extract_formats.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>

void get_parents(Mutation_Annotated_Tree::Tree *T,
                 std::unordered_set<std::string> &need_parents,
                 std::unordered_set<std::string> &all_nodes) {

    FILE *node_to_parent_fp =
        fopen("recombination/filtering/data/nodeToParent.txt", "w");

    if (node_to_parent_fp == NULL) {
        fprintf(stderr, "Error.  Please call like this: build/ripplesUtils "
                        "<MAT_name.pb>\n");
    }
    fprintf(node_to_parent_fp, "node\tparent\n");

    for (const auto &id : need_parents) {
        auto node = T->get_node(id);
        if (node == NULL) {
            continue;
        }
        auto parent = node->parent;
        if (parent == NULL) {
            continue;
        }
        std::string node_id = node->identifier;
        std::string parent_id = parent->identifier;

        // Insert parent node_id into collection of all relevant nodes
        all_nodes.insert(parent_id);

        // Format nodeToParent w/out "node_" TODO: temporary fix for sims
        bool format_flag = true;
        if (format_flag == true) {
            if (node_id[0] == 'n') {
                node_id = node_id.substr(5, node_id.size());
            }
            if (parent_id[0] == 'n') {
                parent_id = parent_id.substr(5, parent_id.size());
            }
        }

        // Write to nodeToParent.txt file
        fprintf(node_to_parent_fp, "%s\t%s\n", node->identifier.c_str(),
                parent_id.c_str());
    }
    fclose(node_to_parent_fp);
}
