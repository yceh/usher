#include "amplicon.hpp"
#include "tbb/concurrent_unordered_set.h"
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <random>
#include <stdio.h>
#include <time.h>
#include <tuple>
#include <unordered_map>

namespace po = boost::program_options;

int main(int argc, char **argv) {
    Timer timer;

    std::string outdir;

    po::options_description desc("optimize options");
    desc.add_options()(
        "input-mat,i", po::value<std::string>()->required(),
        "Input mutation-annotated tree file to optimize [REQUIRED].")(
        "outdir,d", po::value<std::string>(&outdir)->default_value("."),
        "Output directory to dump output and log files [DEFAULT uses current "
        "directory]")(
        "vcf,v", po::value<std::string>()->required(),
        "Input VCF file (in uncompressed or gzip-compressed .gz format) "
        "[REQUIRED]")("fasta,f", po::value<std::string>()->required(),
                      "Input aligned fasta alignment file) "
                      "[REQUIRED]");

    po::options_description all_options;
    all_options.add(desc);
    po::positional_options_description p;
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv)
                      .options(all_options)
                      .positional(p)
                      .run(),
                  vm);
        po::notify(vm);
    } catch (std::exception &e) {
        std::cerr << desc << std::endl;
        // Return with error code 1 unless
        // the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }

    MAT::Tree T;
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string input_vcf_filename = vm["vcf"].as<std::string>();
    std::string input_fasta_filename = vm["fasta"].as<std::string>();

    timer.Start();
    fprintf(stdout, "Loading input MAT file = %s\n",
            input_mat_filename.c_str());
    fprintf(stdout, "Loading input FASTA file = %s\n",
            input_fasta_filename.c_str());

    // Load input MAT to place amplicons onto and uncondense tree
    T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld msec \n", timer.Stop());

    auto bfs = T.breadth_first_expansion();
    size_t total_nodes = bfs.size();
    // tbb::concurrent_unordered_set<std::string> nodes_to_consider;

    std::unordered_map<MAT::Node *, size_t> tree_num_leaves;
    for (int i = int(bfs.size()) - 1; i >= 0; i--) {
        auto n = bfs[i];
        size_t desc = 1;
        for (auto child : n->children) {
            desc += tree_num_leaves[child];
        }
        tree_num_leaves[n] = desc;
    }

    // All input amplicon sequences to place, not found already on input tree
    std::vector<Missing_Sample> missing_samples;

    // Read VCF and extract new amplicon sequences to place
    MAT::read_vcf(&T, input_vcf_filename, missing_samples, false);

    auto num_samples = missing_samples.size();
    // Check that there are actually missing samples from VCF to add to tree
    if (num_samples <= 0) {
        printf("No new samples found in input VCF to add to input MAT. "
               "Exiting program.\n");
        return 0;
    }
    // Amplicon samples to place
    for (auto &sample : missing_samples) {
        printf("%s\n", sample.name.c_str());
    }
    printf("Found %lu missing amplicon samples.\n", num_samples);

    boost::filesystem::path path(outdir);
    if (!boost::filesystem::exists(path)) {
        boost::filesystem::create_directory(path);
    }
    path = boost::filesystem::canonical(outdir);
    outdir = path.generic_string();

    static tbb::affinity_partitioner ap;

    std::vector<std::tuple<int, int>> samples_start_end_coordinates;
    samples_start_end_coordinates.reserve(num_samples);

    // Get start and ending coordinates of each amplicon aligned to reference
    // from aligned multifasta file
    get_coordinates(input_fasta_filename, samples_start_end_coordinates);
    /*
    for (auto &tup : sample_start_end_coordinates) {
      printf("START: %d\t", std::get<0>(tup));
      printf("END: %d\n", std::get<1>(tup));
    }
    */
    printf("All amplicon sample alignment coordinates retrieved.\n");

    std::vector<std::vector<MAT::Mutation>> node_excess_mutations(total_nodes);
    std::vector<std::vector<MAT::Mutation>> imputed_mutations(total_nodes);
    std::vector<int> node_set_difference(total_nodes);

    size_t best_node_num_leaves = 0;
    int best_set_difference = 1e9;

    std::vector<bool> node_has_unique(total_nodes);
    size_t best_j = 0;
    bool best_node_has_unique = false;

    size_t best_distance = 1e9;

    std::vector<size_t> best_j_vec;

    size_t num_best = 1;
    MAT::Node *best_node = T.root;
    best_j_vec.emplace_back(0);

    std::vector<size_t> node_distance(total_nodes, 0);

    // Iterate over all the missing amplicon samples, and place them one at a
    // time
    for (size_t i = 0; i < missing_samples.size(); ++i) {

        auto missing_sample = missing_samples[i];
        size_t num_mutations = missing_sample.mutations.size();

        fprintf(stderr, "Currently placing amplicon: %s\n",
                missing_sample.name.c_str());

        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, total_nodes),
            [&](tbb::blocked_range<size_t> r) {
                for (size_t k = r.begin(); k < r.end(); ++k) {
                    // TODO:
                    /*
                    if (tree_num_leaves[bfs[k]] < num_descendants) {
                       continue;
                     }
                    */

                    node_has_unique[k] = false;
                    mapper2_input inp;
                    inp.T = &T;
                    inp.node = bfs[k];
                    // TODO
                    // inp.missing_sample_mutations =
                    //&pruned_sample.sample_mutations;
                    inp.missing_sample_mutations = &missing_sample.mutations;
                    inp.excess_mutations = &node_excess_mutations[k];
                    inp.imputed_mutations = &imputed_mutations[k];
                    inp.best_node_num_leaves = &best_node_num_leaves;
                    inp.best_set_difference = &best_set_difference;
                    inp.best_node = &best_node;
                    inp.best_j = &best_j;
                    inp.num_best = &num_best;
                    inp.j = k;
                    inp.has_unique = &best_node_has_unique;

                    inp.set_difference = &node_set_difference[k];

                    inp.distance = node_distance[k];
                    inp.best_distance = &best_distance;

                    inp.best_j_vec = &best_j_vec;
                    inp.node_has_unique = &(node_has_unique);

                    // Placement
                    mapper2_body(inp, true);
                }
            },
            ap);

        /*
        for(auto &mat_mutation: missing_sample.mutations){
          fprintf(stderr, "Currently trying to place amplicon: %d\n",
        mat_mutation.mut_nuc);
        }
        break;
        */
    }
}

