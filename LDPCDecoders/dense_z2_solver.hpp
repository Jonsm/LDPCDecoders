//
//  gaussian_z2.hpp
//  LDPCDecoders
//
//  Created by Jon on 5/24/22.
//

#ifndef gaussian_z2_hpp
#define gaussian_z2_hpp

#include <Eigen/Dense>
#include <vector>

class dense_z2_solver {
public:
    dense_z2_solver();
    dense_z2_solver(Eigen::MatrixXi& mat);
    dense_z2_solver(int w, int h);
    void reset(const Eigen::Block<Eigen::MatrixXi>& mat_new);
    bool solve(Eigen::VectorXi& soln, Eigen::VectorXi& rhs);
    
private:
    int w;
    int h;
    Eigen::VectorXi rhs_tmp;
    Eigen::MatrixXi mat;
    std::vector<int> set_pivots_col;
    std::vector<int> set_pivots_row;
    int n_pivots;
    std::vector<int> pivots_row;
    std::vector<int> pivots_col;
    Eigen::MatrixXi pivots_added_to_row;
    
    void init_pivots();
};

#endif /* gaussian_z2_hpp */
