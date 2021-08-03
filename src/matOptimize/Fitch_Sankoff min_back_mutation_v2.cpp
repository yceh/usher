#include "Fitch_Sankoff.hpp"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <immintrin.h>
#include <popcntintrin.h>
#include <tbb/concurrent_vector.h>
#include <utility>
#include <vector>
typedef int v4si __attribute__((vector_size(16)));
namespace MAT = Mutation_Annotated_Tree;
// https://stackoverflow.com/questions/776508/best-practices-for-circular-shift-rotate-operations-in-c
uint32_t rotl(uint32_t x, unsigned int n) {
    const decltype(n) mask = (CHAR_BIT * sizeof(x) - 1);
    // just use unsigned int for the count mask unless you're experimenting with
    // compiler behaviour

    assert(n <= mask && "rotate by more than type width");
    n &= mask; // avoid undef behaviour with NDEBUG.
    return (x << n) | (x >> ((-n) & mask));
}

void FS_backward_pass(const std::vector<backward_pass_range> &child_idx_range,
                      std::vector<uint8_t> &boundary1_major_allele,
                      const mutated_t &mutated, nuc_one_hot ref_nuc);
//return (parent allele possibilities not covered by major allele; possible major alleles at this node)
static std::pair<uint8_t, uint8_t>
get_child_allele_possibility(uint8_t parent_major_allele_possibility,
                             uint8_t this_boundary1_major_allele) {
    uint8_t this_major_allele = this_boundary1_major_allele & 0xf;
    uint8_t parent_allele_not_covered =
        parent_major_allele_possibility & (~this_boundary1_major_allele);
    if (parent_allele_not_covered) {
        // May need mutation to parent allele, so all current major allele and
        // matching boundary alleles are possible
        return std::make_pair(parent_allele_not_covered,
                              this_major_allele |
                                  (parent_allele_not_covered &
                                   (this_boundary1_major_allele >> 4)));
    } else {
        // All alleles can follow parent, so let it be so
        return std::make_pair(parent_allele_not_covered,
                              this_major_allele &
                                  parent_major_allele_possibility);
    }
}

struct back_mut_info {
    std::array<v4si, 4> mutation_to;
    std::array<v4si, 4> this_level_mutation_to;
    v4si back_mutations;
    void clear() {
        for (int nuc_idx = 0; nuc_idx < 4; nuc_idx++) {
            mutation_to = {0, 0, 0, 0};
            this_level_mutation_to={0,0,0,0};
        }
    }
    //add count from another
    void operator+=(const back_mut_info &other) {
        back_mutations += other.back_mutations;
        for (int nuc_idx = 0; nuc_idx < 4; nuc_idx++) {
            mutation_to[nuc_idx] += other.mutation_to[nuc_idx];
            this_level_mutation_to[nuc_idx] += other.this_level_mutation_to[nuc_idx];
        }
    }
};

struct undecided_node_info {
    size_t node_idx;
    size_t dependent_par_node_idx;
    uint8_t allele_choice;
    uint8_t resolved_allele;
    undecided_node_info() {}
    undecided_node_info(size_t node_idx, size_t dependent_par_node_idx)
        : node_idx(node_idx), dependent_par_node_idx(dependent_par_node_idx) {}
};

static void set_mutation(uint8_t parent_state, uint8_t this_state,
                         const MAT::Mutation &base,
                         tbb::concurrent_vector<MAT::Mutation> &output) {
    if (this_state != parent_state) {
        MAT::Mutation this_mut(base);
        this_mut.set_par_mut(parent_state, this_state);
        this_mut.set_auxillary(this_state, 0);
        output.push_back(this_mut);
    }
}

std::pair<uint8_t, uint32_t> get_best_major_allele(int this_node_major_allele,
                                                   v4si back_mutation_count) {
    uint32_t best_back_mutation_count_among_major_allele;
    uint8_t best_major_allele;
    for (int this_allele_idx = 0; this_allele_idx < 4; this_allele_idx++) {
        //only interested in major alleles
        if (!(this_node_major_allele & (two_bit_to_one_hot(this_allele_idx)))) {
            continue;
        }
        int back_mut_count = back_mutation_count[this_allele_idx];
        if (back_mut_count < best_back_mutation_count_among_major_allele) {
            best_back_mutation_count_among_major_allele = back_mut_count;
            best_major_allele = this_allele_idx;
        }
    }
    return std::make_pair(best_major_allele,
                          best_back_mutation_count_among_major_allele);
}

uint8_t update_back_mutation(uint8_t this_node_boundary1_major_allele,
                             uint8_t parent_allele_not_matched,
                             back_mut_info &out) {
    auto this_node_boundary_allele = this_node_boundary1_major_allele >> 4;
    auto this_node_major_allele = this_node_boundary1_major_allele & 0xf;
    //Assume all alleles can follow parent by default
    uint8_t allele_choices = 0xE4;
    //mask for windowing the choice w.r.t. a particular parent allele
    int rot_mask = 0xfc;
    //If matched, following parent optimize parsimony score
    if (parent_allele_not_matched) {
        //Alleles not following parent may be chosen from any major allele or itself if it is in boundary allele (inc parsimony score by 1)
    std::pair<uint8_t, uint32_t> major_allele_least_back_mutation;
        major_allele_least_back_mutation =
            get_best_major_allele(this_node_major_allele, out.back_mutations);
        for (int par_nuc_try_idx = 0; par_nuc_try_idx < 4; par_nuc_try_idx++) {
            uint8_t par_nuc_try = two_bit_to_one_hot(par_nuc_try_idx);
            if (!(parent_allele_not_matched & par_nuc_try)) {
                continue;
            }
            // different from parent only when it is not a boundary/major allele
            // or is a boundary allele, but boundary allele have more back
            // mutation
            if ((!(par_nuc_try & this_node_boundary_allele)) ||
                (out.back_mutations[par_nuc_try_idx] >
                 major_allele_least_back_mutation.second)) {
                //copy over the back mutations field of the best major allele
                out.back_mutations[par_nuc_try_idx] =
                    major_allele_least_back_mutation.second;
                out.mutation_to[par_nuc_try_idx] =
                    out.mutation_to[major_allele_least_back_mutation.first];
                assert(par_nuc_try_idx!=major_allele_least_back_mutation.first);
                //record mutation increment to the best major allele
                out.mutation_to[par_nuc_try_idx][major_allele_least_back_mutation.first]=1;
                allele_choices =
                    (allele_choices & rotl(rot_mask, (2 * par_nuc_try_idx))) |
                    (major_allele_least_back_mutation.first
                     << (2 * par_nuc_try_idx));
            }
        }
    }
    return allele_choices;
}
//shared information down invocations of get_allele_choice
struct pos_info {
    const std::vector<backward_pass_range> &child_idx_range;
    std::vector<undecided_node_info> &undecided_nodes;
    std::vector<tbb::concurrent_vector<MAT::Mutation>> &mut_output;
    const MAT::Mutation &base;
    const std::vector<uint8_t> &boundary1_major_allele;
};

uint8_t resolve_allele(uint8_t this_node_allele_choice,uint8_t possible_parent_state){
    return 0x3 & (this_node_allele_choice >> (2 * one_hot_to_two_bit(possible_parent_state)));
}
//update back mutation count and best choice for different parent allele
std::pair<uint8_t, uint8_t> get_allele_choice(size_t node_idx,
                                              uint8_t possible_parent_state,
                                              pos_info &shared,
                                              uint8_t parent_allele_not_matched,
                                              back_mut_info &back_mut) {
    uint8_t this_node_allele_choice;
    uint8_t this_node_resolved_allele;
    bool parent_resolved = ((__builtin_popcount(possible_parent_state)) > 1);
    if (shared.child_idx_range[node_idx].child_size == 0) {
        //leaf node, assume no ambiguity
        auto this_state =
            two_bit_to_one_hot(shared.boundary1_major_allele[node_idx] & 0xf);
        //mutation vector for par allele different from this allele
        v4si different = {0, 0, 0, 0};
        //another mutation to this_state
        different[this_state] = 1;
        for (int par_nuc_idx = 0; par_nuc_idx < 4; par_nuc_idx++) {
            auto par_nuc = two_bit_to_one_hot(par_nuc_idx);
            //remain zero mutations if following parent
            if (!(possible_parent_state & par_nuc)) {
                back_mut.mutation_to[par_nuc_idx] = different;
            }
        }
        this_node_resolved_allele = this_state;
    } else {
        v4si first_two = __builtin_ia32_blendps(back_mut.mutation_to[0],
                                                back_mut.mutation_to[1], 2);
        v4si back_mutation_count = __builtin_ia32_blendps(
            back_mut.mutation_to[2], back_mut.mutation_to[3], 8);
        back_mutation_count =
            __builtin_ia32_blendps(first_two, back_mutation_count, 12);
        back_mut.back_mutations += back_mutation_count;
        this_node_allele_choice =
            update_back_mutation(shared.boundary1_major_allele[node_idx],
                                 parent_allele_not_matched, back_mut);
        if (parent_resolved) {
            this_node_resolved_allele =resolve_allele(this_node_allele_choice, possible_parent_state);
        } else {
            this_node_resolved_allele = 0xff;
        }
    }
    return std::make_pair(this_node_resolved_allele, this_node_allele_choice);
}

back_mut_info gather_back_mutation_count(size_t node_idx,
                                         bool need_back_mut_count,
                                         uint8_t possible_parent_state,
                                         long depend_idx, pos_info &shared) {
    // should need back mutation count if parent is not resolved
    bool parent_resolved = __builtin_popcount(possible_parent_state) == 1;
    assert(need_back_mut_count || parent_resolved);
    auto temp = get_child_allele_possibility(
        possible_parent_state, shared.boundary1_major_allele[node_idx]);
    uint8_t possible_this_major_allele = temp.second;
    bool this_node_unresolved =
        temp.first && __builtin_popcount(possible_this_major_allele) == 1;
    bool this_node_need_back_mutation =
        this_node_unresolved || need_back_mut_count;
    bool will_undecided =
        temp.first && ((!parent_resolved) || this_node_unresolved);
    auto child_end = shared.child_idx_range[node_idx].first_child_bfs_idx +
                     shared.child_idx_range[node_idx].child_size;
    back_mut_info back_mut;
    if (this_node_need_back_mutation) {
        back_mut.clear();
    }
    size_t this_undecided_idx;
    if (will_undecided) {
        this_undecided_idx = shared.undecided_nodes.size();
        shared.undecided_nodes.emplace_back(node_idx, depend_idx);
    }
    for (size_t child_idx =
             shared.child_idx_range[node_idx].first_child_bfs_idx;
         child_idx < child_end; child_idx++) {
        if (!this_node_need_back_mutation) {
            gather_back_mutation_count(child_idx, this_node_need_back_mutation,
                                       possible_this_major_allele, -2, shared);
        } else {
            back_mut += gather_back_mutation_count(
                child_idx, this_node_need_back_mutation,
                possible_this_major_allele,
                (will_undecided) ? this_undecided_idx : depend_idx, shared);
        }
    }
    if (this_node_need_back_mutation) {
        // Leaf node
        auto this_node_alleles = get_allele_choice(
            node_idx, possible_parent_state, shared, temp.first, back_mut);
        if (this_node_unresolved) {
            shared.undecided_nodes[this_undecided_idx].resolved_allele =
                this_node_alleles.first;
            shared.undecided_nodes[this_undecided_idx].allele_choice =
                this_node_alleles.second;
            return back_mut;
        }
    }
    set_mutation(one_hot_to_two_bit(possible_parent_state),
                 one_hot_to_two_bit(possible_this_major_allele), shared.base,
                 shared.mut_output[node_idx]);
    return back_mut;
}
void resolve_undecided_nodes(std::vector<undecided_node_info> &undecided_nodes,std::vector<tbb::concurrent_vector<MAT::Mutation>> &mut_output,const MAT::Mutation& base){
    for (auto& undecided : undecided_nodes) {
        if (undecided.resolved_allele==0xff) {
            assert(undecided.dependent_par_node_idx!=-2);
            auto parent_state=undecided_nodes[undecided.dependent_par_node_idx].resolved_allele;
            undecided.resolved_allele=resolve_allele(undecided.allele_choice, parent_state);
            set_mutation(parent_state, undecided.resolved_allele, base, mut_output[undecided.node_idx]);
        }
    }
}
void Fitch_Sankoff_Minimize_back_mutation(
    const std::vector<backward_pass_range> &child_idx_range,
    const Mutation_Annotated_Tree::Mutation &base,
    std::vector<tbb::concurrent_vector<Mutation_Annotated_Tree::Mutation>>
        &output,const mutated_t &mutated){
    std::vector<uint8_t> boundary1_major_alleles(child_idx_range.size());
    FS_backward_pass(child_idx_range, boundary1_major_alleles, mutated, base.get_ref_one_hot());
    std::vector<undecided_node_info> undecided_nodes;
    pos_info shared{child_idx_range,undecided_nodes,output,base,boundary1_major_alleles};
    gather_back_mutation_count(0,false,base.get_ref_one_hot(),-2,shared);
    resolve_undecided_nodes(undecided_nodes,output,base);
}

