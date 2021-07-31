#include <cstdint>
#include <array>
#include <cstring>
struct backward_info_per_allele {
  private:
    unsigned int mutations[4][4]; //[par_nuc][mut_nuc]
  public:
    unsigned int back_mutations_count;
    unsigned int parsimony;
    unsigned int get_mut_count(uint8_t mut_nuc, uint8_t par_nuc) const {
        return mutations[par_nuc][mut_nuc];
    }
    unsigned int &get_mut_count(uint8_t mut_nuc, uint8_t par_nuc) {
        return mutations[par_nuc][mut_nuc];
    }
    unsigned int get_total_muts_to_particular_allele(uint8_t mut_nuc) const {
        unsigned int out = 0;
        for (int par_nuc = 0; par_nuc < 4; par_nuc++) {
            out += get_mut_count(mut_nuc, par_nuc);
        }
        return out;
    }
    void add_other_mutations(uint8_t excluded_mut_nuc,
                             const backward_info_per_allele &other) {
        for (int par_nuc = 0; par_nuc < 4; par_nuc++) {
            // add the rest toghether
            for (int mut_nuc = 0; mut_nuc < 4; mut_nuc++) {
                get_mut_count(mut_nuc, par_nuc) +=
                    (mut_nuc == excluded_mut_nuc)
                        ? 0
                        : other.get_mut_count(mut_nuc, par_nuc);
            }
        }
    }
    void clear() {
        /*for (int par_nuc = 0; par_nuc < 4; par_nuc++) {
            for (int mut_nuc = 0; mut_nuc < 4; mut_nuc++) {
                get_mut_count(mut_nuc, par_nuc)=0;
            }
        }*/
        memset(mutations, 0,sizeof(unsigned int)*16);
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
