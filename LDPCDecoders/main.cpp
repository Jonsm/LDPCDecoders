//
//  main.cpp
//  LDPCDecoders
//
//  Created by Jon on 5/9/22.
//

#include <iostream>
#include <string>
#include "vide_decoder.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "sparse_z2_solver.hpp"

using namespace std;

void test() {
    Eigen::SparseMatrix<int> s(6,5);
    s.insert(1,1) = 1;
    s.insert(2,2) = 1;
    s.insert(4,1) = 1;
    s.insert(4,2) = 1;
//    for (int k = 0; k < s.outerSize(); k++) {
//        for (Eigen::SparseMatrix<int>::InnerIterator it(s,k); it; ++it) {
//            cout << it.row() << ","<< it.col() << endl;
//        }
//    }
    
    sparse_z2_solver solver;
    solver.compute(s);
        
    cout << s << endl;
}

int main(int argc, const char * argv[]) {
    mt19937 engine;
    //    engine.seed(std::random_device{}());
    engine.seed(10);
    string filename = argv[1];
    //    vector<float> p_err {0.99,0.005,0.0,0.005};
    //    float h = .8;
    int h_bb = 1;
    int h_cc = 1;
    
    vide_decoder code(filename, engine, h_bb, h_cc);
    code.debug();

    
//    test();
}
