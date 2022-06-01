//
//  vide_decoder.hpp
//  LDPCDecoders
//
//  Created by Jon on 5/10/22.
//

#ifndef vide_decoder_hpp
#define vide_decoder_hpp

#include "HGP_code.hpp"
#include <Eigen/Dense>
#include "dense_z2_solver.hpp"

class vide_decoder : public HGP_code {
public:
    vide_decoder(std::string filename, std::mt19937& engine, int h_bb, int h_cc);
    bool decode();
    void debug();
    
private:
    int h_bb;
    int h_cc;
    std::vector<int> L;
    std::vector<int> neighbors_in_R;
    int n_L;
    int n_R;
    std::vector<int> to_visit;
    int n_to_visit;
    std::vector<int> set_L;
    std::vector<int> set_R;
    
    Eigen::VectorXi restr_correction;
    Eigen::MatrixXi restr_bit_nbhd;
    dense_z2_solver rspace_check;
    dense_z2_solver correctable_check;
    
    void reset();
    void add_to_R(int chk, boost::multi_array<int,2>& chk_nbhd);
    void init(std::vector<int>& syndrome, boost::multi_array<int,2>& chk_nbhd);
    void visit_last(boost::multi_array<int,2>& chk_nbhd, boost::multi_array<int,2>& bit_nbhd);
    bool decode_helper(boost::multi_array<int,2>& bit_nbhd, boost::multi_array<int,2>& chk_nbhd, std::vector<int>& syndrome, std::vector<int>& errs);
    bool L_correctable(boost::multi_array<int,2>& bit_nbhd, std::vector<int>& syndrome);
};

#endif /* vide_decoder_hpp */
