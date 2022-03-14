#include "amplicon.hpp"
#include "fmt/core.h"
#include "fmt/format.h"
#include "src/ripplesUtils/text_parser.hpp"
#include "tbb/concurrent_unordered_set.h"
#include <memory>
#include <numeric>
#include <random>
#include <tuple>
#include <vector>

// Get starting and ending genomic coordinates of where the amplicon aligned to
// reference
void get_coordinates(
    std::string fasta_file_path,
    std::vector<std::tuple<int, int>> &samples_start_end_coordinates) {

    //  Read entire fasta file into memory
    text_parser fasta(fasta_file_path);

    // Skip over reference sequence
    fasta.next_line();

    char at = '>';
    for (; !fasta.done(); fasta.next_line()) {
        if (fasta.get_value(0)[0] != at) {
            continue;
        }
        // Get coordinates debugging
        // fmt::print("Name of amplicon sample: {}\t", fasta.get_value(0));
        // fmt::print("Size of amplicon: {}\t", fasta.get_value(1));
        // fmt::print("Starting position: {}\t", fasta.get_value(2));
        // fmt::print("Ending position: {}\n", fasta.get_value(3));

        int start = str_view_to_uint64(fasta.get_value(2));
        int end = str_view_to_uint64(fasta.get_value(3));
        std::tuple<int, int> tup = std::make_tuple(start, end);
				// Add each start and end coordinate tuple for each amplicon
        samples_start_end_coordinates.push_back(tup);
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
               std::vector<std::tuple<int, int>> samples_start_end_coordinates,
               bool create_new_mat) {}
