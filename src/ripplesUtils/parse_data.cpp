#include "extract_formats.h"
#include "text_parser.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

void get_trios(MAT::Tree T, std::string filepath) {

    // Read entire dataset into memory
    text_parser rec(filepath);

    std::unordered_set<std::string> all_nodes;
    std::unordered_set<std::string> needParent;

    // Skip over header line (starting with #)
    rec.next_line();
    std::string_view yes{"y"};

    for (; !rec.done(); rec.next_line()) {
        auto recomb = "node_" + std::string{rec.get_value(0)};
        auto donor = "node_" + std::string{rec.get_value(3)};
        auto acceptor = "node_" + std::string{rec.get_value(6)};

        // Get the recomb, donor, acceptor trios from each line
        all_nodes.insert(recomb);
        all_nodes.insert(donor);
        all_nodes.insert(acceptor);

        // Check if donor is placed as sibling
        if (rec.get_value(4) == yes) {
            needParent.insert(donor);
        }
        // Check if acceptor is placed as sibling
        if (rec.get_value(7) == yes) {
            needParent.insert(acceptor);
        }
    }
    // Call get_parents to retrieve the parents for all the donors/acceptors
    // placed as siblings
    get_parents(&T, needParent, all_nodes);

    // Create allRelevantNodeNames.txt file
    FILE *allRelevantNodeNames_fp =
        fopen("recombination/filtering/data/allRelevantNodeNames.txt", "w");
    for (const auto &n : all_nodes) {
        fprintf(allRelevantNodeNames_fp, "%s\n", n.c_str());
    }
    fclose(allRelevantNodeNames_fp);
}

