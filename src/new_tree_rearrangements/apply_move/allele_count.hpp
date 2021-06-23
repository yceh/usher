#include "../mutation_annotated_tree.hpp"
#include <array>
namespace MAT = Mutation_Annotated_Tree;
//Clas tracking how many non-parent major allele there are among direct children of a node
struct Allele_Count_t {
    //one example mutation at this loci, for position and par_nuc
    MAT::Mutation base;
    //count of each allele
    std::array<int, 4> count;
    //number of children having any non-parent major allele
    int node_cnt;
    //compatability function for getting position of allele count
    int get_position() const { return base.get_position(); }

    //Merge 2 allele count toghether
    void operator+=(const Allele_Count_t &to_add) {
        assert(base.get_position() == to_add.base.get_position());
        for (int i = 0; i < 4; i++) {
            count[i] += to_add.count[i];
        }
        node_cnt += to_add.node_cnt;
    }
    //add from a raw mutation
    void operator+=(const MAT::Mutation &to_add) {
        assert(base.get_position() == to_add.get_position());
        for (int i = 0; i < 4; i++) {
            if ((1 << i) & to_add.get_all_major_allele()) {
                count[i] += 1;
            }
        }
        node_cnt++;
    }

    Allele_Count_t(){};
    //construct from a mutation, with base set and itsmajor allele counted
    Allele_Count_t(const MAT::Mutation &in)
        : base(in), count({0, 0, 0, 0}), node_cnt(0) {
        (*this) += in;
    };
};