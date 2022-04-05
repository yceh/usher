#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <signal.h>
#include <vector>
enum Segment_Type { A = 0, T = 1, G = 2, C = 3, N, O, R };
#ifndef NDEBUG
#define assert(exp, ...)                                                       \
    if (!(exp)) {                                                              \
        fprintf(stderr, __VA_ARGS__);                                          \
        raise(SIGTRAP);                                                        \
    }
#else
#define assert(exp, msg)
#endif

std::vector<uint8_t> ref_2bit;
typedef float v4f __attribute__((vector_size(16)));
class Segment_Element {
    uint32_t type_and_length;

  public:
    uint32_t start_position_inclusive;
    float down_blen;
    float up_blen;
    Segment_Element(Segment_Type t, uint32_t start_position_inclusive,
                    uint32_t segment_length = 0, float down_blen = 0,
                    float up_blen = 0)
        : type_and_length(segment_length | (t << 29) | ((up_blen != -1) << 28)),
          start_position_inclusive(start_position_inclusive),
          down_blen(down_blen), up_blen(up_blen) {}
    uint32_t get_segment_length() const {
        return type_and_length & (0x0fff'ffff);
    }
    Segment_Type get_seg_type() const {
        return Segment_Type(type_and_length >> 29);
    }
    bool is_from_root() const { return 1 & (type_and_length >> 28); }
    uint32_t get_end_position_exclusive() const {
        return start_position_inclusive + get_segment_length();
    }
    int get_allele() const {
        int allele = get_seg_type();
        if (allele == Segment_Type::R) {
            allele = ref_2bit[start_position_inclusive];
        }
    }
};

union Segment {
    Segment_Element head;
    v4f prob;
};
//#define check_merge_order
struct Position_Checker {
#ifdef check_merge_order
    uint32_t last_put_position_end_exclusive;
    Position_Checker() : last_put_position_end_exclusive(0) {}
#endif
    void operator()(uint32_t start_inclusive, uint32_t end_exclusive) {
#ifdef check_merge_order
        assert(start_inclusive == last_put_position_end_exclusive,
               "Not contiguous, prev %d, curr %d\n",
               last_put_position_end_exclusive, start_inclusive);
        last_put_position_end_exclusive = end_exclusive;
#endif
    }
};
typedef std::vector<Segment> Gene_List;
// f callable with start_position, end_position, head_first,
// head_second,prob_first, prob_second
template <typename iter_type> static void retire_segment(iter_type &iter) {
    if (iter->head.get_seg_type() == Segment_Type::O) {
        iter++;
    }
    iter++;
}
template <typename F>
void merge_gene_list(const Gene_List &first, const Gene_List &second, F &f) {
    auto first_iter = first.begin();
    auto second_iter = second.begin();
    Position_Checker position_checker;
    while (first_iter != first.end()) {
        // Assert the second iterator is actually relevent
        assert(second_iter->head.get_end_position_exclusive() >
                   first_iter->head.start_position_inclusive,
               "second irrelevent, first_pos :%d, second_end %d\n",
               first_iter->head.start_position_inclusive,
               second_iter->head.get_end_position_exclusive());
        while (true) {
            // Maybe left over from previous segment in the second genelist
            uint32_t start_position =
                std::min(first_iter->head.start_position_inclusive,
                         second_iter->head.start_position_inclusive);
            auto second_end_position =
                second_iter->head.get_end_position_exclusive();
            auto first_end_position =
                first_iter->head.get_end_position_exclusive();
            uint32_t end_position =
                std::min(second_end_position, first_end_position);
            position_checker(start_position, end_position);
            f(start_position, end_position, first_iter->head, second_iter->head,
              (first_iter + 1)->prob, (second_iter + 1)->prob);
            // second interval completely covered by the first interval
            if (second_end_position <= first_end_position) {
                retire_segment(second_iter);
            }
            // first interval have been covered by preceding segments in the
            // second interval
            if (second_iter->head.start_position_inclusive >=
                first_iter->head.get_end_position_exclusive()) {
                retire_segment(first_iter);
            }
        }
    }
    assert(second_iter == second.end(), "second gene list not consumed\n");
}
static float sum_vect(v4f in) { return (in[0] + in[1]) + (in[2] + in[3]); }
struct Mut_Parameter {
    // Rate matrix
    v4f src_to_dst[4];
    v4f dst_to_src[4]; // transposed
    v4f root_frequencies;
    std::vector<float> culmulative_rate;
    float get_cum_rate(uint32_t end, uint32_t start) const {
        return culmulative_rate[end] - culmulative_rate[start - 1];
    }
    v4f marginalize(float downLen, v4f downProb) const {
        v4f tot2 = src_to_dst[0] * downProb[0];
        for (int i = 1; i < 4; i++) {
            tot2 += src_to_dst[i] * downProb[i];
        }
        tot2 *= downLen;
        tot2 += downProb;
        return tot2;
    }
    v4f get_upDown_ref_prob(int idx, float rootLen, float downLen,
                            v4f downProb) const {
        v4f trans_prob = dst_to_src[idx] * rootLen;
        trans_prob[idx] += 1;
        trans_prob *= root_frequencies;

        return trans_prob * marginalize(downLen, downProb);
    }
};

struct Append_Prob_Functor {
    const Mut_Parameter &mut_parameter;
    const float blen;
    float LKcost;
    float totalFactor;
    Append_Prob_Functor(const Mut_Parameter &mut_parameter, const float blen)
        : mut_parameter(mut_parameter), blen(blen), LKcost(0), totalFactor(1) {}
    void operator()(uint32_t start_position, uint32_t end_position,
                    const Segment_Element &head_first,
                    const Segment_Element &head_second, const v4f &prob_first,
                    const v4f &prob_second) {
        auto first_type = head_first.get_seg_type();
        auto second_type = head_second.get_seg_type();
        // If either one is N, not contributing to likelihood
        if (first_type == Segment_Type::N || second_type == Segment_Type::N) {
            return;
        }
        // Both R, possibly end_position!=start_position+1
        if (first_type == Segment_Type::R && second_type == Segment_Type::R) {
            // STRANGE: Not delta: Since it is spliting a branch,
            // the likelihood of going here from ancestor([5]? up_blen)
            // and to decendant ([3]) should have been taken care of in the
            //  original calculation, so should only consider blen?

            // Following is literal translation
            LKcost += (blen + head_first.down_blen + head_first.down_blen) *
                      mut_parameter.get_cum_rate(start_position, end_position);
            return;
        }
        // The rest should only have one nucleotide
        assert((start_position + 1) == end_position, "Longer strech\n");
        // Two Os
        if (first_type == Segment_Type::O) {
            if (second_type == Segment_Type::O) {
                totalFactor *=
                    (sum_vect(prob_second > 0.5
                                  ? mut_parameter.marginalize(
                                        blen + head_first.down_blen, prob_first)
                                  : 0)) /
                    sum_vect(prob_first);
            } else {
                auto second_allele = head_second.get_allele();
                assert(second_allele < 4, "expecting second allele be ATGC\n");
                totalFactor *=
                    ((prob_first[second_allele] +
                      (blen + head_first.down_blen) *
                          sum_vect(mut_parameter.dst_to_src[second_allele] *
                                   prob_first)) /
                     sum_vect(prob_first));
            }
        } else {
            auto first_allele = head_first.get_allele();
            assert(first_allele < 4, "expecting first allele be ATGC\n");
            auto second_allele = head_second.get_allele();
            if (second_allele == Segment_Type::O) {
                if (head_second.is_from_root()) {
                    totalFactor *=
                        (sum_vect(mut_parameter.get_upDown_ref_prob(
                             first_allele, head_first.down_blen,
                             blen + head_first.up_blen, prob_second)) /
                         mut_parameter.root_frequencies[first_allele]);
                } else {
                    totalFactor *=
                        (sum_vect((prob_second > 0.5
                                       ? mut_parameter.src_to_dst[first_allele]
                                       : 0)) *
                         (blen + head_first.down_blen));
                }
            } else {
                assert(second_allele < 4, "expecting second allele be ATGC\n");
                if (first_allele == second_allele) {
                    LKcost +=
                        (blen + head_first.down_blen + head_first.down_blen) *
                        mut_parameter.get_cum_rate(start_position,
                                                   end_position);
                } else {
                    if (head_second.is_from_root()) {
                        auto prob = mut_parameter.dst_to_src[first_allele] *
                                    mut_parameter.dst_to_src[second_allele] *
                                    (blen + head_first.up_blen) *
                                    head_first.down_blen;
                        auto acc_prob =
                            sum_vect(prob) +
                            mut_parameter.root_frequencies[first_allele] *
                                mut_parameter
                                    .dst_to_src[second_allele][first_allele] *
                                (head_first.up_blen + blen)+
                            mut_parameter.root_frequencies[second_allele] *
                                mut_parameter
                                    .dst_to_src[first_allele][second_allele] *
                                (head_first.down_blen)
                            ;
                        totalFactor*=(acc_prob/mut_parameter.root_frequencies[first_allele]);
                    }else {
                        totalFactor*=mut_parameter.src_to_dst[first_allele][second_allele]*(blen+head_first.down_blen);
                    }
                }
            }
        }
        LKcost+=std::log(totalFactor);
        totalFactor=1;
    }
};

static float append_prob(const Mut_Parameter &mut_parameter, const float blen,
                        const Gene_List &first, const Gene_List &second) {
    Append_Prob_Functor functor(mut_parameter, blen);
    merge_gene_list(first, second, functor);
    return functor.LKcost;
}