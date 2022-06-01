//
//  gaussian_z2.cpp
//  LDPCDecoders
//
//  Created by Jon on 5/24/22.
//

#include "dense_z2_solver.hpp"
#include <iostream>

using namespace Eigen;
using namespace std;

void dense_z2_solver::init_pivots() {
    for (int row = 0; row < h; row++) {
        for (int col = 0; col < w; col++) {
            if (n_pivots == min(w,h)) {
                break;
            }
            
            if (mat(row,col) && set_pivots_col[col] == -1) {
                set_pivots_col[col] = n_pivots;
                set_pivots_row[row] = n_pivots;
                pivots_row[n_pivots] = row;
                pivots_col[n_pivots] = col;
                n_pivots++;
                
                int rows_added_to=0;
                for (int row2 = 0; row2 < h; row2++) {
                    if (row2 != row && mat(row2,col)) {
                        mat.row(row2) += mat.row(row);
                        for (int col2 = 0; col2 < w; col2++) {
                            mat(row2,col2) %= 2;
                        }
                        
                        pivots_added_to_row(n_pivots-1,rows_added_to)=row2;
                        rows_added_to++;
                    }
                }
                
                break;
            }
        }
    }
}

dense_z2_solver::dense_z2_solver() :
w(0),
h(0)
{
    
}

dense_z2_solver::dense_z2_solver(MatrixXi& mat) :
mat(mat),
rhs_tmp(mat.rows()),
set_pivots_col(mat.cols()),
set_pivots_row(mat.rows()),
pivots_row(mat.cols()),
pivots_col(mat.cols()),
pivots_added_to_row(MatrixXi::Constant(mat.rows(),mat.rows(),-1)),
n_pivots(0)
{
    w = (int)mat.cols();
    h = (int)mat.rows();
    fill(set_pivots_col.begin(), set_pivots_col.end(), -1);
    fill(set_pivots_row.begin(), set_pivots_row.end(), -1);
    
    init_pivots();
}

dense_z2_solver::dense_z2_solver(int w, int h) :
w(w),
h(h),
rhs_tmp(h),
set_pivots_col(w),
set_pivots_row(h),
pivots_row(w),
pivots_col(w),
pivots_added_to_row(MatrixXi::Constant(h,h,-1)),
mat(w,h)
{
    fill(set_pivots_col.begin(), set_pivots_col.end(), -1);
    fill(set_pivots_row.begin(), set_pivots_row.end(), -1);
}


void dense_z2_solver::reset(const Block<MatrixXi>& mat_new) {
    int w_new = (int)mat_new.cols();
    int h_new = (int)mat_new.rows();
    
    if (w_new > w) {
        set_pivots_col.resize(w_new);
        pivots_row.resize(w_new);
        pivots_col.resize(w_new);
    }
    if (h_new > h) {
        rhs_tmp.resize(h_new);
        set_pivots_row.resize(h_new);
        pivots_added_to_row.resize(h_new,h_new);
    }
    if (w_new > w || h_new > h) {
        mat.resize(h_new,w_new);
    }
    
    w = w_new;
    h = h_new;
    
    n_pivots = 0;
    fill(set_pivots_col.begin(), set_pivots_col.begin()+w, -1);
    fill(set_pivots_row.begin(), set_pivots_row.begin()+w, -1);
    pivots_added_to_row.topLeftCorner(h,h) = MatrixXi::Constant(h,h,-1);
    mat.topLeftCorner(h,w) = mat_new;
    
    init_pivots();
}

bool dense_z2_solver::solve(Eigen::VectorXi &soln, Eigen::VectorXi &rhs) {
    rhs_tmp = rhs;
    fill(soln.begin(), soln.end(), 0);
    for (int i = 0; i < n_pivots; i++) {
        if (rhs_tmp[pivots_row[i]]) {
            int j = 0;
            while (pivots_added_to_row(i,j) != -1) {
                rhs_tmp[pivots_added_to_row(i,j)] ^= 1;
                j++;
            }
        }
    }
    
    for (int row = 0; row < h; row++) {
        if (rhs_tmp[row]) {
            if (set_pivots_row[row] != -1) {
                soln[pivots_col[set_pivots_row[row]]] = 1;
            } else {
                return false;
            }
        }
    }
    
    return true;
}
