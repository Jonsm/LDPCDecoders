//
//  HGP_code.hpp
//  LDPCDecoders
//
//  Created by Jon on 5/9/22.
//

#ifndef HGP_code_hpp
#define HGP_code_hpp

#include <random>
#include <string>
#include <vector>
#include "boost/multi_array.hpp"
#include <fstream>

class HGP_code {
public:
    HGP_code(std::string filename, std::mt19937& engine);
    void clear();
    void add_errors(std::vector<float>& error_p); //I,X,Y,Z
    void get_n_errors(std::vector<int>& errors_out);
    bool check_correctness();
    
protected:
    std::mt19937& engine;
    std::uniform_real_distribution<float> dis;
    
    int n;
    int m;
    int dv;
    int dc;
    std::string id;
    boost::multi_array<int,2> bit_nbhd_c; //A
    boost::multi_array<int,2> chk_nbhd_c; //B
    
    int n_chk_q;
    int n_q;
    int dc_q;
    
    std::vector<int> dv_q;
    boost::multi_array<int,2> x_bit_nbhd; //A*A first then B*B
    boost::multi_array<int,2> z_bit_nbhd; //A*A first then B*B
    boost::multi_array<int,2> x_chk_nbhd; //A*B
    boost::multi_array<int,2> z_chk_nbhd; //A*B
    
    std::vector<int> x_syndrome; //syndrome of x checks
    std::vector<int> z_syndrome; //syndrome of z checks
    std::vector<int> x_errs;
    std::vector<int> z_errs;
    std::vector<int> n_errs;
    
    std::pair<int,int> chk_coord(int chk);
    std::pair<int,int> bit_coord(int bit);
    void add_error(int b, boost::multi_array<int,2>& bit_nbhd, std::vector<int>& syndrome, std::vector<int>& errs);
    
private:
    std::string read_param(std::ifstream& file, std::string label);
    void read_matrix(std::ifstream& file, boost::multi_array<int,2>& matrix, std::string label);
    void read_file(std::string filename);
    void gen_bits();
    void gen_checks();
    void gen_q_code();
    bool check_correct_graph();
    bool check_correct_stabilizers();
};

#endif /* HGP_code_hpp */
