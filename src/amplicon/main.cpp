#include "amplicon.hpp"
#include "fmt/core.h"
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <time.h>

namespace po = boost::program_options;

int main(int argc, char **argv) {
    Timer timer;

    po::options_description desc("optimize options");
    desc.add_options()(
        "input-mat,i", po::value<std::string>()->required(),
        "Input mutation-annotated tree file to optimize [REQUIRED].")(
        "min-coordinate-range,r", po::value<int>()->default_value(1e3),
        "Starting position of the genomic coordinates of the amplicon sequence "
        "in reference sequence [REQUIRED].")(
        "max-coordinate-range,R", po::value<int>()->default_value(1e7),
        "Ending position of the genomic coordinates of the amplicon sequence "
        "in the reference sequence.");

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
    std::string input_mat_filename = vm["input-mat"].as<std::string>();

    timer.Start();
    fmt::print(stderr, "Loading input MAT file = {}\n", input_mat_filename);

    MAT::Tree T;
    // Load input MAT to place amplicons onto and uncondense tree
    T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fmt::print(stderr, "Completed in {} msec \n", timer.Stop());
}

