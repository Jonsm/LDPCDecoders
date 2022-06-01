//
//  sparse_z2_solver.hpp
//  LDPCDecoders
//
//  Created by Jon on 6/1/22.
//

#ifndef sparse_z2_solver_hpp
#define sparse_z2_solver_hpp

#include <Eigen/Sparse>
#include <vector>

class sparse_z2_solver {
public:
    void compute(Eigen::SparseMatrix<int>& mat);
    bool solve(std::vector<int>& rhs, std::vector<int>& sol);
    
private:
    int rows = 0;
    int cols = 0;
    std::vector<int> set_piv_row;
    std::vector<int> set_piv_col;
    std::vector<int> pivs_col;
    std::vector<int> pivs_row;
    int n_piv;
    
    std::vector<std::vector<int>> rows_added_to;
    Eigen::VectorXi tmp_col;
    
    std::vector<int> rhs_tmp;
    
    void resize(Eigen::SparseMatrix<int>& mat);
    void compute_col(int col);
};

#endif /* sparse_z2_solver_hpp */
