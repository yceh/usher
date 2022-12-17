#include "get_recomb_info.hpp"
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <string>

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
        "chronumental-dates,c", po::value<std::string>()->default_value(""),
        "If using Chronumental, give output inferred dates file from running "
        "Chronumental. Otherwise earliest date from recombinant node "
        "descendants will be used. ")(
        "metadata,m", po::value<std::string>()->default_value(""),
        "If not using Chronumental, give MAT metadata file which contains "
        "dates for all descendants.")(
        "final-recombinants,r", po::value<std::string>()->required(),
        "Output file containing filtered recombinants with all "
        "information needed for RIVET[REQUIRED].");

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
    std::string metadata_file = vm["metadata"].as<std::string>();

    // Load input MAT and uncondense tree
    printf("Loading input MAT file\n");
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();

    // Get number of leaves in the tree
    auto num_leaves = T.get_num_leaves();
    std::cout << "Tree contains: " << num_leaves << " leaves."
              << "\n";

    // Create new recombinant output file
    std::ofstream outfile{final_recomb_file};
    if (!outfile) {
        throw std::runtime_error(
            "ERROR: Cannot create final recombination output file.");
    }
    std::cout << "Retrieving recombinant node parent clade assignments"
              << "\n";
    std::cout << "Outfile given: " << final_recomb_file << "\n";

    // Output file columns
    std::vector<std::string> header_list = {"Recombinant Node ID",
                                            "Donor Node ID",
                                            "Acceptor Node ID",
                                            "Breakpoint 1 Interval",
                                            "Breakpoint 2 Interval",
                                            "Recombinant Clade",
                                            "Recombinant Lineage",
                                            "Donor Clade",
                                            "Donor Lineage",
                                            "Acceptor Clade",
                                            "Acceptor Lineage",
                                            "Descendant",
                                            "Informative Site Sequence",
                                            "3SEQ (M, N, K)",
                                            "3SEQ P-Value",
                                            "Recombinant Ranking Score",
                                            "Original Parsimony Score",
                                            "Parsimony Score Improvement"};

    std::vector<std::string> trio_node_ids;

    // If Chronumental inferred internal dates file provided, use this method
    if (chron_dates_file != "") {
        std::cout << "Chronumental inferred dates file given: "
                  << chron_dates_file << "\n";
        std::cout << "Using Chronumental for recombinant node ranking."
                  << "\n";

        // Load inferred dates for internal nodes from Chronumental output
        std::unordered_map<std::string_view, std::string_view>
            node_to_inferred_date;
        node_to_inferred_date.reserve(num_leaves);

        tsv_to_dict(chron_dates_file, node_to_inferred_date, 0, 1, true);

        // NOTE: Chronumental will preserve internal node id naming using Newick
        // generated from  matUtils extract. Get information for each column for
        // all filtered recombinants, including rank score, and output to
        // outfile.  Return a vector of string node ids for all trio nodes
        // (recomb, donor, acceptor)
        trio_node_ids =
            get_recombination_info(T, tree_date, node_to_inferred_date,
                                   filtered_recomb_file, outfile, header_list);
    }
    // If no Chronumental inferred dates file given, use alternate method
    // of chosing recombinant node descendant with earliest date
    else {
        // If Chronumental not used, check to make sure metadata input file
        // given
        if (metadata_file == "") {
            throw std::runtime_error(
                "ERROR: If not using Chronumental (-c flag), then metadata "
                "file must be provided through --metadata (-m) flag");
        }
        if (metadata_file.substr(metadata_file.find_last_of(".") + 1) == "gz") {
            throw std::runtime_error("Input metadata file must be unzipped and "
                                     "provided through --metadata (-m) flag");
        }
        std::cout << "Using alternate earliest descendant date method for "
                     "recombinant node ranking."
                  << "\n";

        // Load MAT metadata into dictionary to get dates for descendants
        std::unordered_map<std::string_view, std::string_view>
            descendant_to_date;
        descendant_to_date.reserve(num_leaves);
        tsv_to_dict(metadata_file, descendant_to_date, 0, 2, true);

        // Get all recombinant node information, rank recombinants and write to
        // given output file, return vector of string node_ids for all trio
        // nodes (recomb, donor, acceptor)
        get_recombination_info_using_descendants(
            T, tree_date, filtered_recomb_file, descendant_to_date, outfile,
            header_list);
    }

    std::cout << "Final recombination results written to:  "
              << final_recomb_file << "\n";

    outfile.close();
    return 0;
}

