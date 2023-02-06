#ifndef GET_RECOMB_INFO_H
#define GET_RECOMB_INFO_H

#include "src/mutation_annotated_tree.hpp"
#include "src/ripples/util/text_parser.hpp"
#include <cstddef>
#include <math.h>

namespace MAT = Mutation_Annotated_Tree;

struct Recombinant {
    int original_parsimony;
    int parsimony_improvement;
    float recomb_rank;
    std::string recomb_node_id;
    std::string donor_node_id;
    std::string acceptor_node_id;
    std::string informative_seq;      // eg) "AAAAAAAAAAABBB"
    std::string informative_position; // eg) "1,2,..."
    std::string mnk_3seq_values;      // eg) "(M,N,K)"
    std::string p_value_3seq;
		std::string inferred_date;
    std::string descendants;
    std::tuple<std::string, std::string>
        breakpoint_intervals; // Breakpoint interval <1,2>
    std::string filter;
    Recombinant(std::string recomb_node_id) :recomb_node_id(recomb_node_id){} 
};

struct Ranked_Recombinant {
    float recomb_rank;
    size_t id;

    Ranked_Recombinant(size_t id):id(id){}
};

struct RepresentativeSample {
    int num_novel_mutations;
    std::string node_id;

    RepresentativeSample(std::string id) {
        node_id = id;
        num_novel_mutations = 0;
    }
};

void write_recombination_list(
    MAT::Tree &T, std::unordered_map<std::string, Recombinant> &recombinants,
    std::vector<Ranked_Recombinant> &ranked_recombs, std::ofstream &outfile,
    std::vector<std::string> &header_list,std::vector<Recombinant>& filtered_out_recombs);

std::vector<std::string> get_recombination_info(
    MAT::Tree &T, std::string tree_date,
    std::unordered_map<std::string, std::string> &node_to_inferred_date,
    std::string filtered_recomb_file, std::ofstream &outfile,
    std::vector<std::string> &header_list);

void get_recombination_info_using_descendants(
    MAT::Tree &T, std::string tree_date, std::string filtered_recomb_file,
    std::unordered_map<std::string, std::string> &descendant_to_date,
    std::ofstream &outfile, std::vector<std::string> &header_list);

void chron_id_mapping(MAT::Tree &T,
                      std::unordered_map<std::string, std::string> &id_map);

void tsv_to_dict(std::string tsv_file,
                 std::unordered_map<std::string, std::string> &map, int key_col,
                 int val_col, bool header);

inline float recombinant_rank(int days, int num_descendants, int weight) noexcept {
    return (static_cast<float>(num_descendants) / days) * pow(2, -weight);
}

int elapsed_days(std::string tree_date,
                 std::string inferred_recomb_date) noexcept;

std::vector<std::string> format_date(std::string date);

std::string parse_string(const std::string s, int n, std::string delimiter);

void internal_node_msa(MAT::Tree T, std::string sample_paths,
                       std::vector<std::string> &trio_node_ids,
                       std::string reference_file, std::ofstream &trio_msa);

std::string find_representative_sample(MAT::Tree &T,
                                       std::string internal_node_id);

std::vector<std::string>
get_node_descendants(MAT::Tree &T, std::string internal_node_id, int max_desc);

#endif
