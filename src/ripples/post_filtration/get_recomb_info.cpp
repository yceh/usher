#include "get_recomb_info.hpp"
#include <ctime>
#include <iostream>
#include <string_view>
#include <time.h>
#include <tuple>
#include <unordered_map>

// Format for string date: "2022-08-14"
// Returns int vector of size 3 [year, month, day]
std::vector<std::string> format_date(std::string date) {
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

//  Calculate elapsed days from given date to present date.
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

    for (const auto &rr : ranked_recombs) {
        // Get the Recombinant node and write
        outfile << rr.recomb_node_id << "\t";

        Recombinant r = recombinants.at(rr.recomb_node_id);

        //  Write recombinant node breakpoint intervals
        outfile << std::get<0>(r.breakpoint_intervals) << "\t";
        outfile << std::get<1>(r.breakpoint_intervals) << "\t";

        // Write donor node id
        outfile << r.donor_node_id << "\t";

        // Lookup donor node and check it exists in tree
        auto donor = T.get_node(r.donor_node_id);
        if (donor == NULL) {
            std::cout << "Donor node is NULL, not finding donor node id"
                      << "\n";
            exit(1);
        }
        // Get donor clade (nextstrain) and lineage (pangolin designation)
        // get_clade_assignment(node, 0) => nextstrain
        // get_clade_assignment(node, 1) => pangolin
        auto donor_clade = T.get_clade_assignment(donor, 0);
        auto donor_lineage = T.get_clade_assignment(donor, 1);
        outfile << donor_clade << "\t";
        outfile << donor_lineage << "\t";

        // Write acceptor node id
        outfile << r.acceptor_node_id << "\t";

        // Lookup acceptor node and check it exists in tree
        auto acceptor = T.get_node(r.acceptor_node_id);
        if (acceptor == NULL) {
            std::cout << "Acceptor node is NULL, not finding acceptor node id"
                      << "\n";
            exit(1);
        }

        // Get acceptor clade (nextstrain) and lineage (pangolin designation)
        auto acceptor_clade = T.get_clade_assignment(acceptor, 0);
        auto acceptor_lineage = T.get_clade_assignment(acceptor, 1);
        outfile << acceptor_clade << "\t";
        outfile << acceptor_lineage << "\t";

        // Write recombinant node ranking score (increasing order)
        outfile << rr.recomb_rank << "\n";
    }
    outfile.close();
}

void get_recombination_info(
    MAT::Tree &T, std::string tree_date,
    std::unordered_map<std::string, std::string> &node_to_inferred_date,
    std::string filtered_recomb_file, std::ofstream &outfile,
    std::vector<std::string> header_list) {

    text_parser results(filtered_recomb_file);

    std::cout << "Opening results file: " << filtered_recomb_file << "\n";

    // Calculate number of elapsed days since input tree date
    // int days = elapsed_days(tree_date);

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

        // Parse only date from Chronumental inferred dates dictionary
        std::string inferred_recomb_date = node_to_inferred_date.at(recomb_id);
        int space_index = inferred_recomb_date.find(" ", 0);
        inferred_recomb_date = inferred_recomb_date.substr(0, space_index);

        // Calculate number of elapsed days since input tree date
        int days = elapsed_days(tree_date, inferred_recomb_date);

        Ranked_Recombinant rr = Ranked_Recombinant(recomb_id);
        // Generate recombinant ranking score
        auto recomb_rank = recombinant_rank(days, recomb_num_descendants);
        r.recomb_rank = recomb_rank;
        rr.recomb_rank = recomb_rank;
        // Add recombination information to collection of detected recombinants
        recombinants.insert({recomb_id, r});

        // Keep track of rank score for each detected recombinant
        ranked_recombs.push_back(rr);
    }
    // Sort the recombinants by min score
    std::sort(ranked_recombs.begin(), ranked_recombs.end(),
              [](const Ranked_Recombinant &a, const Ranked_Recombinant &b) {
                  return a.recomb_rank < b.recomb_rank;
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

void tsv_to_dict(std::string tsv_file,
                 std::unordered_map<std::string, std::string> &map, int key_col,
                 int val_col, bool header) {
    text_parser file_handle(tsv_file);

    // If file has header, skip over first header line
    // Assuming header is just a single first line in TSV file
    if (header == true) {
        file_handle.next_line();
    }

    for (; !file_handle.done(); file_handle.next_line()) {
        auto key = std::string{file_handle.get_value(key_col)};
        auto value = std::string{file_handle.get_value(val_col)};
        map.insert({key, value});
    }
}

