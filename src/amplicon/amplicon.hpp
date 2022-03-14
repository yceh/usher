#ifndef AMPLICON_H
#define AMPLICON_H

#include "src/mutation_annotated_tree.hpp"
#include "src/ripplesUtils/text_parser.hpp"
#include "src/usher_graph.hpp"

namespace MAT = Mutation_Annotated_Tree;

struct Pruned_Sample {
    std::string sample_name;
    std::vector<MAT::Mutation> sample_mutations;
    std::unordered_set<uint32_t> positions;

    // Assumes mutations are added in reverse chrono order
    void add_mutation(MAT::Mutation mut) {
        // If not reversal to reference allele
        if ((mut.ref_nuc != mut.mut_nuc) &&
            (positions.find(mut.position) == positions.end())) {
            auto iter = std::lower_bound(sample_mutations.begin(),
                                         sample_mutations.end(), mut);
            auto m = mut.copy();
            m.par_nuc = m.ref_nuc;
            sample_mutations.insert(iter, m);
        }
        positions.insert(mut.position);
    }

    Pruned_Sample(std::string name) {
        sample_name = name;
        sample_mutations.clear();
        positions.clear();
    }
};

void get_coordinates(
    std::string sam_file_path,
    std::vector<std::tuple<int, int>> &samples_start_end_coordinates);

inline uint64_t str_view_to_uint64(std::string_view str) noexcept;

void placement(MAT::Tree *T, std::string &vcf_filename,
               std::vector<Missing_Sample> &missing_samples,
               std::vector<std::tuple<int, int>> samples_start_end_coordinates,
               bool create_new_mat);

#endif

