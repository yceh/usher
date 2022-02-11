#include "extract_formats.hpp"
#include <iostream>
#include <string>
#include <string_view>

int main(int argc, char **argv) {

    std::string input_mat_filename;
    const std::string &combined_pvals_filename =
        "recombination/filtering/data/combinedCatOnlyBestWithPVals.txt";

    if (argc == 2) {
        input_mat_filename = std::string{argv[1]};
    } else {
        fprintf(stderr, "Error.  Please call like this: build/ripplesUtils "
                        "<MAT_name.pb>\n");
        exit(1);
    }

    fprintf(stdout, "Loading input MAT file %s.\n", input_mat_filename.c_str());
    // Load input MAT and uncondense tree
    MAT::Tree T;
    if (input_mat_filename.find(".pb\0") != std::string::npos) {
        T = MAT::load_mutation_annotated_tree(input_mat_filename);
        T.uncondense_leaves();
    } else {
        fprintf(
            stderr,
            "ERROR: Input file ending not recognized. Must be .json or .pb\n");
        exit(1);
    }
    fprintf(stdout, "Completed loading MAT, getting trios now.\n");

    // Ouputs two files:  allRelevantNodeNames.txt and nodeToParent.txt
    get_trios(T, combined_pvals_filename);

    return 0;
}
