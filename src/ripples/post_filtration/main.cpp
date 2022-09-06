#include "get_recomb_info.hpp"
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>

namespace po = boost::program_options;

int main(int argc, char **argv) {

    po::options_description desc("optimize options");
    desc.add_options()("input-mat,i", po::value<std::string>()->required(),
                       "Input mutation-annotated tree file [REQUIRED].")(
        "filtered-recombinants,f", po::value<std::string>()->required(),
        "Input file containing filtered recombinants from filtration "
        "pipeline [REQUIRED].")(
        "date,d", po::value<std::string>()->required(),
        "MAT tree date (format: year-month-day, eg. 2022-08-14) [REQUIRED].")(
        "chronumental-dates,c", po::value<std::string>()->required(),
        "Output inferred dates file from running Chronumental [REQUIRED].")(
        "final-recombinants,r", po::value<std::string>()->required(),
        "Output file containing filtered recombinants with all "
        "information needed for Recombination Tracker UI [REQUIRED].");

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
    std::string filtered_recomb_file =
        vm["filtered-recombinants"].as<std::string>();
    std::string final_recomb_file = vm["final-recombinants"].as<std::string>();
    std::string tree_date = vm["date"].as<std::string>();
    std::string chron_dates_file = vm["chronumental-dates"].as<std::string>();

    // Load input MAT and uncondense tree
    printf("Loading input MAT file\n");
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();

    // Create new recombinant output file
    std::ofstream outfile{final_recomb_file};
    if (!outfile) {
        throw std::runtime_error(
            "ERROR: Cannot find filtered recombination results file");
    }
    std::cout << "Retrieving recombinant node parent clade assignments"
              << "\n";
    std::cout << "Outfile given: " << final_recomb_file << "\n";

    // Load inferred dates for internal nodes from Chronumental output
    std::unordered_map<std::string, std::string> node_to_inferred_date;
    tsv_to_dict(chron_dates_file, node_to_inferred_date, 0, 1, true);

    // Output file columns
    std::vector<std::string> header_list = {
        "Recombinant Node ID",   "Breakpoint Interval 1",
        "Breakpoint Interval 2", "Donor Node ID",
        "Donor Clade",           "Donor Lineage",
        "Acceptor Node ID",      "Acceptor Clade",
        "Acceptor Lineage",      "Recombinant Ranking Score"};

    // NOTE: Chronumental will preserve internal node id naming using Newick
    // generated from  matUtils extract Get information for each column for all
    // filtered recombinants, including rank score, and output to outfile
    get_recombination_info(T, tree_date, node_to_inferred_date,
                           filtered_recomb_file, outfile, header_list);

    std::cout << "Final recombination results written to:  " << final_recomb_file
              << "\n";

    outfile.close();

    return 0;
}

