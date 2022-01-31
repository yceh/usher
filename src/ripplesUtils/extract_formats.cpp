#include "extract_formats.h"
#include <fstream>
#include <iostream>
#include <sstream>

void get_parents(Mutation_Annotated_Tree::Tree *T,
                 std::unordered_set<std::string> &need_parents,
                 std::unordered_set<std::string> &all_nodes) {

    FILE *node_to_parent_fp =
        fopen("recombination/filtering/data/nodeToParent.txt", "w");
    fprintf(node_to_parent_fp, "node\tparent\n");

    for (const auto &node_id : need_parents) {
        auto node = T->get_node(node_id);
        if (node == NULL) {
            continue;
        }
        auto parent = node->parent;
        if (parent == NULL) {
            continue;
        }
        std::string parent_id = parent->identifier;
        // Insert parent node_id into collection of all relevant nodes
        all_nodes.insert(parent_id);

        // Output just node ids without "node_" prefix
        // fprintf(node_to_parent_fp, "%s\t", node_id.c_str());
        // fprintf(node_to_parent_fp, "%s\n", parent_id.c_str());

        // Temporary formatting for getABABA.py in filtration pipeline (remove
        // "node_")
        fprintf(node_to_parent_fp, "%s\t",
                node_id.substr(5, node_id.size()).c_str());
        fprintf(node_to_parent_fp, "%s\n",
                parent_id.substr(5, parent_id.size()).c_str());
    }

    fclose(node_to_parent_fp);
}
