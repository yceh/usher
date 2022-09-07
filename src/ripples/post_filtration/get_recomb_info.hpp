#ifndef GET_RECOMB_INFO_H
#define GET_RECOMB_INFO_H

#include "src/mutation_annotated_tree.hpp"
#include "src/ripples/util/text_parser.hpp"

namespace MAT = Mutation_Annotated_Tree;

struct Recombinant {
    float recomb_rank;
    std::string recomb_node_id;
    std::string donor_node_id;
    std::string acceptor_node_id;
    std::tuple<std::string, std::string>
        breakpoint_intervals; // Breakpoint interval <1,2>

    Recombinant(std::string recomb_node_id) { recomb_node_id = recomb_node_id; }
};

struct Ranked_Recombinant {
    float recomb_rank;
    std::string recomb_node_id;

    Ranked_Recombinant(std::string id) { recomb_node_id = id; }
};

void write_recombination_list(
    MAT::Tree &T, std::unordered_map<std::string, Recombinant> &recombinants,
    std::vector<Ranked_Recombinant> &ranked_recombs, std::ofstream &outfile,
    std::vector<std::string> &header_list);

void get_recombination_info(
    MAT::Tree &T, std::string tree_date,
    std::unordered_map<std::string_view, std::string_view>
        &node_to_inferred_date,
    std::string filtered_recomb_file, std::ofstream &outfile,
    std::vector<std::string> header_list);

void get_recombination_info_using_descendants(
    MAT::Tree &T, std::string tree_date, std::string filtered_recomb_file,
    std::unordered_map<std::string_view, std::string_view> &descendant_to_date,
    std::ofstream &outfile, std::vector<std::string> header_list);

void chron_id_mapping(MAT::Tree &T,
                      std::unordered_map<std::string, std::string> &id_map);

void tsv_to_dict(std::string tsv_file,
                 std::unordered_map<std::string_view, std::string_view> &map,
                 int key_col, int val_col, bool header);

inline float recombinant_rank(int days, int num_descendants) noexcept {
    return days / static_cast<float>(num_descendants);
}

int elapsed_days(std::string tree_date,
                 std::string inferred_recomb_date) noexcept;
std::vector<std::string> format_date(std::string date);

#endif
