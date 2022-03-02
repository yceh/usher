#include "fmt/core.h"
#include "fmt/format.h"
#include <memory>
#include <random>
#include <vector>
//#include "tbb/concurrent_unordered_set.h"
#include "amplicon.hpp"
#include "src/ripplesUtils/text_parser.hpp"

void get_starting_coordinates(std::string sam_file_path,
                              std::vector<int> &start_coordinates) {

    // Read entire SAM file into memory
    text_parser sam(sam_file_path);

    char at = '@';
    while (1) {
        // Get first character in current line, and check if '@'
        std::string first_char = std::string{sam.get_value(0).at(0)};
        if (sam.get_value(0)[0] != at) {
            //  Once moved past metadata headers
            break;
        }
        // If current line still metadata, move to the next line
        sam.next_line();
    }

    //  Current line now first line of actual alignment data
    for (; !sam.done(); sam.next_line()) {
        // Get the start position of the aligment read, which is column 4
        start_coordinates.push_back(str_view_to_uint64(sam.get_value(3)));
    }
}

inline uint64_t str_view_to_uint64(std::string_view str) noexcept {
    uint64_t num = 0;
    for (auto &c : str)
        (num *= 10) += c - '0';
    return num;
}

void placement(MAT::Tree *T, std::string &vcf_filename,
               std::vector<Missing_Sample> &missing_samples,
               bool create_new_mat) {}
