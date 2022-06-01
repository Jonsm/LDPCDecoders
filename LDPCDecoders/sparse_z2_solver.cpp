//
//  sparse_z2_solver.cpp
//  LDPCDecoders
//
//  Created by Jon on 6/1/22.
//

#include "sparse_z2_solver.hpp"
#include <iostream>

using namespace Eigen;
using namespace std;

void sparse_z2_solver::resize(Eigen::SparseMatrix<int> &mat) {
    int new_rows = (int)mat.rows();
    int new_cols = (int)mat.cols();
    
    if (new_rows > rows) {
        set_piv_row.resize(new_rows);
        pivs_row.resize(new_rows);
        pivs_col.resize(new_rows);
        rows_added_to.resize(new_rows);
        tmp_col.resize(new_rows);
        rhs_tmp.resize(new_rows);
    }
    if (new_cols > cols) {
        set_piv_col.resize(new_cols);
    }
    
    fill(set_piv_row.begin(), set_piv_row.end(), -1);
    fill(set_piv_col.begin(), set_piv_col.end(), -1);
    
    n_piv = 0;
    rows = new_rows;
    cols = new_cols;
}

void sparse_z2_solver::compute_col(int col) {
    for (int i = 0; i < n_piv; i++) {
        for (int j = 0; j < rows_added_to[i].size(); j++) {
            tmp_col[rows_added_to[i][j]] ^= tmp_col[pivs_row[i]];
        }
    }
    
    int pivot_row = -1;
    for (int row = 0; row < rows; row++) {
        if (tmp_col[row] != 0 && set_piv_row[row] == -1) {
            pivot_row = row;
            break;
        }
    }
    
    if (pivot_row != -1) {
        rows_added_to[n_piv].clear();
        for (int row = 0; row < rows; row++) {
            if (tmp_col[row] != 0 && row != pivot_row) {
                rows_added_to[n_piv].push_back(row);
//                tmp_col[row] = 0;
            }
        }
        
        pivs_row[n_piv] = pivot_row;
        pivs_col[n_piv] = col;
        set_piv_row[pivot_row] = n_piv;
        set_piv_col[col] = n_piv;
        n_piv++;
    }
}

void sparse_z2_solver::compute(Eigen::SparseMatrix<int> &mat) {
    resize(mat);
    
    for (int col = 0; col < mat.outerSize(); col++) {
        tmp_col.head(rows) = mat.col(col);
        compute_col(col);
    }
}

bool sparse_z2_solver::solve(std::vector<int> &rhs, std::vector<int> &sol) {
    copy(rhs.begin(), rhs.end(), rhs_tmp.begin());
    
    for (int i = 0; i < n_piv; i++) {
        if (rhs_tmp[pivs_row[i]]) {
            for (int j = 0; j < rows_added_to[i].size(); j++) {
                rhs_tmp[rows_added_to[i][j]] ^= rhs_tmp[pivs_row[i]];
            }
        }
    }
    
    fill(sol.begin(), sol.end(), 0);
    
    for (int row = 0; row < rows; row++) {
        if (rhs_tmp[row]) {
            if (set_piv_row[row] != -1) {
                sol[pivs_col[set_piv_row[row]]] = 1;
            } else {
                return false;
            }
        }
    }
    
    return true;
}
