#include "get_recomb_info.hpp"
#include <boost/algorithm/string/join.hpp>
#include <ctime>
#include <iostream>
#include <string_view>
#include <time.h>
#include <tuple>
#include <unordered_map>

// Expected format for string date: "2022-08-14"
// Returns int vector of size 3 [year, month, day]
std::vector<std::string> format_date(std::string date) {
    int date_length = 10; // Only accepting dates in format, eg) 2022-11-09
    if (date.length() != date_length) {
        throw std::runtime_error("ERROR: format_date() function. Date not in "
                                 "correct format (date.length() = 10)");
    }
    std::string delimiter = "-";
    std::vector<std::string> year_month_day;
    int i = 0;
    std::string token;
    while ((i = date.find(delimiter)) != std::string::npos) {
        token = date.substr(0, i);
        year_month_day.push_back(token);
        date = date.substr(i + 1, date.length());
    }
    year_month_day.push_back(date);

    return year_month_day;
}

//   Calculate elapsed days from given date to present date.
int elapsed_days(std::string tree_date,
                 std::string inferred_recomb_date) noexcept {

    // Parse format of MAT tree date
    auto year_month_day = format_date(tree_date);
    int year = stoi(year_month_day[0]) - 1900; // Years since 1900
    int month = stoi(year_month_day[1]) - 1; // Months since January – [0, 11]
    int day = stoi(year_month_day[2]);       // Day of the month – [1, 31]

    struct std::tm x = {0, 0, 0, day, month, year};
    std::time_t date_1 = std::mktime(&x);

    // Parse format of recombinant node inferred date
    auto _year_month_day = format_date(inferred_recomb_date);
    int _year = stoi(_year_month_day[0]) - 1900; // Years since 1900
    int _month = stoi(_year_month_day[1]) - 1; // Months since January – [0, 11]
    int _day = stoi(_year_month_day[2]);       // Day of the month – [1, 31]

    struct std::tm y = {0, 0, 0, _day, _month, _year};
    std::time_t date_2 = std::mktime(&y);

    // Implicit conversion to int to get elapsed days only
    return std::difftime(date_1, date_2) / (60 * 60 * 24);
}

void write_recombination_list(
    MAT::Tree &T, std::unordered_map<std::string, Recombinant> &recombinants,
    std::vector<Ranked_Recombinant> &ranked_recombs, std::ofstream &outfile,
    std::vector<std::string> &header_list) {

    // Add header for outfile
    for (std::vector<std::string>::iterator it = header_list.begin();
         it != header_list.end(); ++it) {
        if (it == std::prev(header_list.end())) {
            outfile << *it << "\n";
            break;
        }
        outfile << *it << "\t";
    }
    // Print number of internal nodes we have
    std::cout << "Ranked recomb size *3: " << ranked_recombs.size() * 3 << "\n";

    // Create samples output file with the name of all recombinant nodes
    std::string samples_file = "all_trio_nodes.txt";
    std::ofstream samples_outfile{samples_file};
    if (!samples_outfile) {
        throw std::runtime_error(
            "ERROR: Cannot create sample nodes output file.");
    }
    // Create samples output file with the name and all descendants
    std::string samples_descendants = "samples_descendants.tsv";
    std::ofstream descendants_outfile{samples_descendants};
    if (!descendants_outfile) {
        throw std::runtime_error(
            "ERROR: Cannot create recombinant nodes descendants output file.");
    }
    // Add header
    descendants_outfile << "recombinant node id"
                        << "\t"
                        << "descendants"
                        << "\t"
                        << "count"
                        << "\n";

    std::unordered_set<std::string> internal_nodes;
    internal_nodes.reserve(ranked_recombs.size() * 3);

    for (const auto &rr : ranked_recombs) {
        // Get the Recombinant node and write
        outfile << rr.recomb_node_id << "\t";

        Recombinant r = recombinants.at(rr.recomb_node_id);

        // Write donor node id
        outfile << r.donor_node_id << "\t";

        // Write acceptor node id
        outfile << r.acceptor_node_id << "\t";

        //  Write recombinant node breakpoint intervals
        outfile << std::get<0>(r.breakpoint_intervals) << "\t";
        outfile << std::get<1>(r.breakpoint_intervals) << "\t";

        // Write recomb node clade and lineage
        // Lookup recomb node and check it exists in tree
        auto recomb = T.get_node(rr.recomb_node_id);
        if (recomb == NULL) {
            std::cout << "Recomb node is NULL, not finding recomb node id"
                      << "\n";
            exit(1);
        }

        // Write descendants to sample_descendants.tsv out file
        std::vector<std::string> recomb_capped_descendants;
        recomb_capped_descendants.reserve(10000);
        get_node_descendants(T, rr.recomb_node_id, recomb_capped_descendants,
                             10000);

        descendants_outfile << rr.recomb_node_id << "\t";

        if (recomb_capped_descendants.size() != 0) {
            for (const auto &d : recomb_capped_descendants) {
                if (d != recomb_capped_descendants.back()) {
                    descendants_outfile << d << ", ";
                } else {
                    descendants_outfile << d << "\t";
                }
            }
        }
        descendants_outfile << recomb_capped_descendants.size() << "\n";

        // Get recomb clade (nextstrain) and lineage (pangolin designation)
        auto recomb_clade = T.get_clade_assignment(recomb, 0);
        auto recomb_lineage = T.get_clade_assignment(recomb, 1);
        outfile << recomb_clade << "\t";
        outfile << recomb_lineage << "\t";

        // Lookup donor node and check it exists in tree
        auto donor = T.get_node(r.donor_node_id);
        if (donor == NULL) {
            std::cout << "Donor node is NULL, not finding donor node id"
                      << "\n";
            exit(1);
        }
        // Get donor clade (nextstrain) and lineage (pangolin designation)
        // get_clade_assignment(node, 0) returns nextstrain
        // get_clade_assignment(node, 1) returns pangolin
        auto donor_clade = T.get_clade_assignment(donor, 0);
        auto donor_lineage = T.get_clade_assignment(donor, 1);
        outfile << donor_clade << "\t";
        outfile << donor_lineage << "\t";

        // Lookup acceptor node and check it exists in tree
        auto acceptor = T.get_node(r.acceptor_node_id);
        if (acceptor == NULL) {
            std::cout << "Acceptor node is NULL, not finding acceptor node id"
                      << "\n";
            exit(1);
        }

        // Write all trio nodes to sample nodes file   TODO: Clean up
        if (internal_nodes.find(rr.recomb_node_id) == internal_nodes.end()) {
            samples_outfile << rr.recomb_node_id << "\n";
            internal_nodes.insert({rr.recomb_node_id});
        }
        if (internal_nodes.find(r.donor_node_id) == internal_nodes.end()) {
            samples_outfile << r.donor_node_id << "\n";
            internal_nodes.insert({r.donor_node_id});
        }
        if (internal_nodes.find(r.acceptor_node_id) == internal_nodes.end()) {
            samples_outfile << r.acceptor_node_id << "\n";
            internal_nodes.insert({r.acceptor_node_id});
        }

        // Get acceptor clade (nextstrain) and lineage (pangolin
        // designation)
        auto acceptor_clade = T.get_clade_assignment(acceptor, 0);
        auto acceptor_lineage = T.get_clade_assignment(acceptor, 1);
        outfile << acceptor_clade << "\t";
        outfile << acceptor_lineage << "\t";

        // Find representative sample with no (or minimal) number of
        // additional mutations compared to internal trio node
        auto recomb_rep_sample =
            find_representative_sample(T, rr.recomb_node_id);

        // Only show representative descendant in table
        outfile << recomb_rep_sample << "\t";

        // Write informative sequence
        outfile << r.informative_seq << "\t";

        // Write 3SEQ (M, N, K) values
        outfile << r.mnk_3seq_values << "\t";

        // Write 3SEQ P-value
        outfile << r.p_value_3seq << "\t";

        // Write recombinant node ranking score (increasing order)
        outfile << rr.recomb_rank << "\t";

        outfile << r.original_parsimony << "\t";
        outfile << r.parsimony_improvement << "\n";
    }
    outfile.close();
    samples_outfile.close();
    descendants_outfile.close();
}

std::vector<std::string>
get_recombination_info(MAT::Tree &T, std::string tree_date,
                       std::unordered_map<std::string_view, std::string_view>
                           &node_to_inferred_date,
                       std::string filtered_recomb_file, std::ofstream &outfile,
                       std::vector<std::string> header_list) {

    std::cout << "Opening results file: " << filtered_recomb_file << "\n";
    text_parser results(filtered_recomb_file);

    // Keep track of all trio node_ids, used to build trio nodes MSA later
    std::vector<std::string> trio_node_ids;

    // Keep track of all recomb_node_ids and their associated rank
    std::unordered_map<std::string, Recombinant> recombinants;
    std::vector<Ranked_Recombinant> ranked_recombs;

    // Get each detected recombinant node from filtration pipeline output
    for (; !results.done(); results.next_line()) {

        auto recomb_id = std::string{results.get_value(0)};
        if (std::isdigit(recomb_id.at(0)) == 1) {
            recomb_id = "node_" + recomb_id;
        }
        // Record recombinant node id
        trio_node_ids.push_back(recomb_id);

        // Create new recombinant node
        Recombinant r = Recombinant(recomb_id);

        // Get original parsimony and parsimony score improvement
        int org_parsimony = stoi(std::string{results.get_value(9)});
        int result_parsimony = stoi(std::string{results.get_value(11)});
        int parsimony_improvement = org_parsimony - result_parsimony;
        r.original_parsimony = org_parsimony;
        r.parsimony_improvement = parsimony_improvement;

        // Get breakpoint intervals for the recombinant node id
        std::tuple<std::string_view, std::string_view> bp(results.get_value(1),
                                                          results.get_value(2));
        r.breakpoint_intervals = bp;

        // Get donor/acceptor node ids
        auto donor_id = std::string{results.get_value(3)};
        if (std::isdigit(donor_id.at(0)) == 1) {
            donor_id = "node_" + donor_id;
        }
        // Record acceptor node id
        trio_node_ids.push_back(donor_id);

        auto acceptor_id = std::string{results.get_value(6)};
        if (std::isdigit(acceptor_id.at(0)) == 1) {
            acceptor_id = "node_" + acceptor_id;
        }
        // Record acceptor node id
        trio_node_ids.push_back(acceptor_id);

        // Record donor/acceptor node ids for ranking
        r.donor_node_id = donor_id;
        r.acceptor_node_id = acceptor_id;

        // Get informative site from filtration results file,
        // which is at column 16
        r.informative_seq = std::string{results.get_value(16)};

        // Get 3SEQ M, N, K values from filtration results file
        std::string mnk_values = "(";
        mnk_values += std::string{results.get_value(17)} + ", ";
        mnk_values += std::string{results.get_value(18)} + ", ";
        mnk_values += std::string{results.get_value(19)} + ")";
        r.mnk_3seq_values = mnk_values;

        // Get 3SEQ P-value from filtration results file
        r.p_value_3seq = std::string{results.get_value(20)};

        // Get descendants from filtration results file
        r.descendants = std::string{results.get_value(14)};

        // Get the recombinant node, and make sure id exists in tree
        auto recomb = T.get_node(recomb_id);
        if (recomb == NULL) {
            std::cout << "Recomb node is NULL, not finding recomb node id"
                      << "\n";
            exit(1);
        }
        // Get number of descendants for recombinant node
        size_t recomb_num_descendants = T.get_num_leaves(recomb);

        // Parse dates only from Chronumental inferred dates dictionary
        std::string_view recomb_id_view{recomb_id};
        std::string inferred_recomb_date =
            std::string{node_to_inferred_date.at(recomb_id_view)};
        int space_index = inferred_recomb_date.find(" ", 0);
        inferred_recomb_date = inferred_recomb_date.substr(0, space_index);

        // Calculate number of elapsed days since input tree date
        int days = elapsed_days(tree_date, inferred_recomb_date);

        Ranked_Recombinant rr = Ranked_Recombinant(recomb_id);
        // Generate recombinant ranking score
        auto recomb_rank = recombinant_rank(days, recomb_num_descendants);
        r.recomb_rank = recomb_rank;
        rr.recomb_rank = recomb_rank;
        // Add recombination information to collection of detected
        // recombinants
        recombinants.insert({recomb_id, r});

        // Keep track of rank score for each detected recombinant
        ranked_recombs.push_back(rr);
    }
    // Sort the recombinants by max score
    std::sort(ranked_recombs.begin(), ranked_recombs.end(),
              [](const Ranked_Recombinant &a, const Ranked_Recombinant &b) {
                  return a.recomb_rank > b.recomb_rank;
              });

    // Write all final recombinants to output file, in ranked order
    write_recombination_list(T, recombinants, ranked_recombs, outfile,
                             header_list);

    return trio_node_ids;
}

void get_recombination_info_using_descendants(
    MAT::Tree &T, std::string tree_date, std::string filtered_recomb_file,
    std::unordered_map<std::string_view, std::string_view> &descendant_to_date,
    std::ofstream &outfile, std::vector<std::string> header_list) {

    std::cout << "Opening results file: " << filtered_recomb_file << "\n";
    text_parser results(filtered_recomb_file);

    // Keep track of all recomb_node_ids and their associated rank
    std::unordered_map<std::string, Recombinant> recombinants;
    std::vector<Ranked_Recombinant> ranked_recombs;

    // Get each detected recombinant node from filtration pipeline output
    for (; !results.done(); results.next_line()) {

        auto recomb_id = std::string{results.get_value(0)};
        if (std::isdigit(recomb_id.at(0)) == 1) {
            recomb_id = "node_" + recomb_id;
        }
        // Create new recombinant node
        Recombinant r = Recombinant(recomb_id);

        // Get breakpoint intervals for the recombinant node id
        std::tuple<std::string_view, std::string_view> bp(results.get_value(1),
                                                          results.get_value(2));
        r.breakpoint_intervals = bp;

        // Get donor/acceptor node ids
        auto donor_id = std::string{results.get_value(3)};
        if (std::isdigit(donor_id.at(0)) == 1) {
            donor_id = "node_" + donor_id;
        }

        auto acceptor_id = std::string{results.get_value(6)};
        if (std::isdigit(acceptor_id.at(0)) == 1) {
            acceptor_id = "node_" + acceptor_id;
        }
        r.donor_node_id = donor_id;
        r.acceptor_node_id = acceptor_id;

        // Get the recombinant node, and make sure id exists in tree
        auto recomb = T.get_node(recomb_id);
        if (recomb == NULL) {
            std::cout << "Recomb node is NULL, not finding recomb node id"
                      << "\n";
            exit(1);
        }
        // Get number of descendants for recombinant node
        size_t recomb_num_descendants = T.get_num_leaves(recomb);

        // Get all descendants for recombinant node
        auto descendants_vec = T.get_leaves(recomb_id);
        if (descendants_vec.size() == 0) {
            std::cout << "RECOMB NODE ID with no descendants" << recomb_id
                      << "\n";
            throw std::runtime_error(
                "ERROR: Recombinant node doesn't have any descendants.");
        }
        // Use earliest date of recomb node descendants (earliest_days) as
        // proxy for inferred recomb note date
        int earliest_days = 0;
        std::string earliest_descendant = "";
        for (auto node : descendants_vec) {
            std::string_view n{node->identifier};
            if (descendant_to_date[n].size() != 10) {
                continue;
            }
            int desc_days =
                elapsed_days(tree_date, std::string{descendant_to_date[n]});
            //  Ties don't matter, we just want earliest date
            if (desc_days > earliest_days) {
                earliest_descendant = node->identifier;
                earliest_days = desc_days;
            }
        }
        Ranked_Recombinant rr = Ranked_Recombinant(recomb_id);
        // Generate recombinant ranking score, using earliest date from set
        // of recomb node descendants
        auto recomb_rank =
            recombinant_rank(earliest_days, recomb_num_descendants);

				//TODO: Might still want to report this recombinant, list N/A in rank field?
        if (recomb_rank == 0.0) {
            // Move to next recombinant, incomplete date information for
            // this recombinant node, for all node descendants in metadata
            continue;
        }
        r.recomb_rank = recomb_rank;
        rr.recomb_rank = recomb_rank;
        // Add recombination information to collection of detected
        // recombinants
        recombinants.insert({recomb_id, r});

        // Keep track of rank score for each detected recombinant
        ranked_recombs.push_back(rr);
    }
    // Sort the recombinants by max score
    std::sort(ranked_recombs.begin(), ranked_recombs.end(),
              [](const Ranked_Recombinant &a, const Ranked_Recombinant &b) {
                  return a.recomb_rank > b.recomb_rank;
              });

    // Write all final recombinants to output file, in ranked order
    write_recombination_list(T, recombinants, ranked_recombs, outfile,
                             header_list);
}

// Same preorder traversal as Chronumental performs to map
// from ripples node ids -> Chronumental node ids
void chron_id_mapping(MAT::Tree &T,
                      std::unordered_map<std::string, std::string> &id_map) {
    MAT::Node *root = T.root;
    std::stack<MAT::Node *> s;
    std::vector<MAT::Node *> preorder;
    if (root == NULL) {
        std::cout << "ERROR: Empty tree!"
                  << "\n";
        exit(1);
    }
    s.push(root);
    while (!s.empty()) {
        auto node = s.top();
        s.pop();
        preorder.push_back(node);
        for (auto &child : node->children) {
            s.push(child);
        }
    }
    auto dfs = T.depth_first_expansion();
    if (dfs.size() != preorder.size()) {
        std::cout << "ERROR: Traversal sizes not matching."
                  << "\n";
        exit(1);
    }
    id_map.reserve(dfs.size());
    for (size_t i = 0; i < dfs.size(); ++i) {
        id_map.insert({dfs[i]->identifier, preorder[i]->identifier});
    }
}

//  Extract two columns from a TSV file to act as dictionary,
//  one as key, the other as value
void tsv_to_dict(std::string tsv_file,
                 std::unordered_map<std::string_view, std::string_view> &map,
                 int key_col, int val_col, bool header) {
    text_parser file_handle(tsv_file);

    // If file has header, skip over first header line
    // Assuming header is just a single first line in TSV file
    if (header == true) {
        file_handle.next_line();
    }

    for (; !file_handle.done(); file_handle.next_line()) {
        std::string_view key = file_handle.get_value(key_col);
        std::string_view value = file_handle.get_value(val_col);
        map.insert({key, value});
    }
}

// Find a sample that is representative of the given internal node,
// meaning that it is a descendant of this internal node and shares no, or
// minimal additional mutations with the internal node
std::string find_representative_sample(MAT::Tree &T,
                                       std::string internal_node_id) {

    // Get internal node
    MAT::Node *internal_node = T.get_node(internal_node_id);

    // Get all the descendant nodes for internal node
    std::vector<MAT::Node *> desc_vec;
    desc_vec = T.get_leaves(internal_node_id);

    std::vector<RepresentativeSample> rep_samples;
    // Go through all of the internal_node descendants and find a sample
    // with no or fewest additional mutations
    for (auto node : desc_vec) {
        RepresentativeSample r = RepresentativeSample(node->identifier);
        r.num_novel_mutations = node->mutations.size();
        rep_samples.push_back(r);
    }
    //  Sort representative samples to see which contains the fewest
    //  additional mutations wrt given internal node
    std::sort(rep_samples.begin(), rep_samples.end(),
              [](const RepresentativeSample &a, const RepresentativeSample &b) {
                  return a.num_novel_mutations < b.num_novel_mutations;
              });

    // Return the sample with fewest additional mutations
    return rep_samples[0].node_id;
}


void get_node_descendants(MAT::Tree &T, std::string internal_node_id,
                          std::vector<std::string> &descendants, int max_desc) {

    // Get internal node
    MAT::Node *internal_node = T.get_node(internal_node_id);

    // Get all the descendant nodes for internal node
    std::vector<MAT::Node *> desc_vec;
    desc_vec = T.get_leaves(internal_node_id);
    int i = 0;
    for (const auto &d : desc_vec) {
        if (i == max_desc) {
            break;
        }
        descendants.push_back(d->identifier);
        ++i;
    }
    return;
}
