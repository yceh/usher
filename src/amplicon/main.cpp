#include "amplicon.hpp"
#include "fmt/core.h"
#include "fmt/format.h"
#include "tbb/concurrent_unordered_set.h"
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <random>
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
        "directory]")
        /*
        ("start-coordinate-range,r", po::value<int>()->default_value(1e3),
        "Starting genomic coordinate position of the amplicon sequence "
        "in reference sequence [REQUIRED].")
        ("end-coordinate-range,R", po::value<int>()->default_value(1e7),
        "Ending genomic coordinate position of the amplicon sequence "
        "in the reference sequence.")
        */
        ("vcf,v", po::value<std::string>()->required(),
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
    fmt::print(stderr, "Loading input MAT file = {}\n", input_mat_filename);
    fmt::print(stderr, "Loading input FASTA file = {}\n", input_fasta_filename);

    // Load input MAT to place amplicons onto and uncondense tree
    T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fmt::print(stderr, "Completed in {} msec \n", timer.Stop());

    auto bfs = T.breadth_first_expansion();
    tbb::concurrent_unordered_set<std::string> nodes_to_consider;
		//fmt::print("{}\n", bfs.size());

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
        fmt::print("No new samples found in input VCF to add to input MAT. "
                   "Exiting program.\n");
        return 0;
    }
    // Amplicon samples to place
    for (auto &sample : missing_samples) {
        fmt::print("{}\n", sample.name);
    }
    fmt::print("Found {} missing amplicon samples.\n", num_samples);

    boost::filesystem::path path(outdir);
    if (!boost::filesystem::exists(path)) {
        fmt::print(stderr, "Creating output directory.\n", timer.Stop());
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
        fmt::print("START: {}\t", std::get<0>(tup));
        fmt::print("END: {}\n", std::get<1>(tup));
    }
		*/
    fmt::print("All amplicon sample alignment coordinates retrieved.\n");
}

