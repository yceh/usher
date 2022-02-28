#include "amplicon.hpp"
#include "fmt/core.h"
#include "fmt/format.h"
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <time.h>

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

    timer.Start();
    fmt::print(stderr, "Loading input MAT file = {}\n", input_mat_filename);

    // Load input MAT to place amplicons onto and uncondense tree
    T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fmt::print(stderr, "Completed in {} msec \n", timer.Stop());

    // All input amplicon sequences to place, not found already on input tree
    std::vector<Missing_Sample> missing_samples;

    // Read VCF and extract new amplicon sequences to place
    MAT::read_vcf(&T, input_vcf_filename, missing_samples, false);

    // Print out names of new amplicon samples to place
    for (auto sample : missing_samples) {
        fmt::print("{}\n", sample.name);

        for (auto mutation : sample.mutations) {
            fmt::print("{}\n", mutation.position);
        }
    }
    fmt::print("Finished printing new amplicon samples from input VCF\n");
}

