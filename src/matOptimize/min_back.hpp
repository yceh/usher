#include <cstdint>
#include <array>
#include <cstring>
struct backward_info_per_allele {
  private:
    unsigned int mut_nuc_count[4]; //[par_nuc][mut_nuc]
  public:
    unsigned int back_mutations_count;
    unsigned int parsimony;

    unsigned int get_total_muts_to_particular_allele(uint8_t mut_nuc) const {
        return mut_nuc_count[mut_nuc];
    }

    unsigned int& get_total_muts_to_particular_allele(uint8_t mut_nuc) {
        return mut_nuc_count[mut_nuc];
    }

    void add_other_mutations(uint8_t excluded_mut_nuc,
                             const backward_info_per_allele &other) {
            for (int mut_nuc = 0; mut_nuc < 4; mut_nuc++) {
                get_total_muts_to_particular_allele(mut_nuc) +=
                    (mut_nuc == excluded_mut_nuc)
                        ? 0
                        : other.get_total_muts_to_particular_allele(mut_nuc);
            }
    }
    void clear() {
        /*for (int par_nuc = 0; par_nuc < 4; par_nuc++) {
            for (int mut_nuc = 0; mut_nuc < 4; mut_nuc++) {
                get_mut_count(mut_nuc, par_nuc)=0;
            }
        }*/
        memset(mut_nuc_count, 0,sizeof(unsigned int)*4);
        back_mutations_count=0;
    }
};
struct backward_info: public std::array<backward_info_per_allele, 4>{
    void clear(){
        for(int idx=0;idx<4;idx++){
            (*this)[idx].clear();
            (*this)[idx].parsimony=0;
        }
    }
    void clear_leaf(int nuc){
        for(int idx=0;idx<4;idx++){
            (*this)[idx].clear();
            (*this)[idx].parsimony=(idx==nuc)?0:2;
        }
    }
};
