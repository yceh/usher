#ifndef LITE_HEADER
#define LITE_HEADER
#include "../check_samples.hpp"
#include "../mutation_annotated_tree.hpp"
namespace MAT = Mutation_Annotated_Tree;
struct Profitable_Moves{
    int score_change;
    MAT::Node* src;
    std::vector<MAT::Node*> src_to_LCA;
    std::vector<MAT::Node*> dst_to_LCA;
    MAT::Node* LCA;
    bool new_node;
    int radius_left;
    std::vector<MAT::Mutation> mutations_relative_to_parent;
    Profitable_Moves():score_change(0){}
    Profitable_Moves(int score_change,MAT::Node* LCA, MAT::Node* dst,const std::vector<MAT::Mutation>& mutations_relative_to_parent):score_change(score_change),dst_to_LCA({dst}),LCA(LCA),mutations_relative_to_parent(mutations_relative_to_parent){}
    MAT::Node* get_src()const{
        return src;
    }
    MAT::Node* get_dst()const{
        return dst_to_LCA.front();
    }
    template<typename F>
    void apply_nodes(F f){
        f(src);
        for(auto node: src_to_LCA){
            f(node);
        }
        for(auto node: dst_to_LCA){
            f(node);
        }
    } 
};
typedef Profitable_Moves* Profitable_Moves_ptr_t;
typedef std::vector<Profitable_Moves_ptr_t> output_t;
void find_place(MAT::Node* src,output_t &output,unsigned int radius_left);
void placement_prep(MAT::Tree *t);
void check_descendent(MAT::Tree& tree);
void get_mutations_relative_to_root(MAT::Node* node,std::vector<MAT::Mutation>& muts);
void compare_mutations(
    std::vector<MAT::Mutation> &original_mutations,
    std::vector<MAT::Mutation> &new_mutations);
void apply_moves(std::vector<Profitable_Moves_ptr_t> &all_moves, MAT::Tree &t,
                 std::vector<MAT::Node *> &bfs_ordered_nodes,
                 tbb::concurrent_vector<MAT::Node *> &to_filter
#ifdef CHECK_STATE_REASSIGN
                 ,
                Original_State_t& original_state
#endif
);
void individual_move(MAT::Node* src,MAT::Node* dst,MAT::Node* LCA,output_t& out);
void check_clean(const std::vector<MAT::Node *> &bfs_ordered_nodes);
void remove_nodes(MAT::Node* const root,size_t idx,MAT::Tree* tree);
size_t optimize_tree(std::vector<MAT::Node *> &bfs_ordered_nodes,
              tbb::concurrent_vector<MAT::Node *> &nodes_to_search,
              MAT::Tree &t,int radius,FILE* log
              #ifndef NDEBUG
              , Original_State_t origin_states
            #endif
              );
void save_final_tree(MAT::Tree &t, Original_State_t& origin_states,
                            const std::string &output_path);
#endif