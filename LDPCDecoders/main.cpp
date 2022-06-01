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
//#include "gaussian_z2.hpp"
//#include "int_z2.hpp"
#include "sparse_z2_solver.hpp"

using namespace std;

//class int_z2 {
//public:
//    bool b;
//    int_z2() { b = false; }
//    int_z2(bool b) : b(b) {}
//    int_z2(int a) : b(a!=0) {}
//    int_z2(double a) : b(a!=0) {}
//    int_z2 operator* (const int_z2 m) const {return m.b & b;}
//    int_z2 operator/ (const int_z2 m) const {return b;}
//    int_z2 operator/ (const int m) const {return b;}
//    int_z2 operator+ (const int_z2 m) const {return m.b ^ b;}
//    int_z2 operator- (const int_z2 m) const {return m.b ^ b;}
//    int_z2 operator- (const int m) const {return m ^ b;}
//    int_z2 operator+= (const int_z2 m) {b ^= m.b; return b;}
//    int_z2 operator*= (const int_z2 m) {b &= m.b; return b;}
//    int_z2 operator/= (const int_z2 m) {return b;}
//    int_z2 operator-= (const int_z2 m) {b ^= m.b; return b;}
//    bool operator== (const int_z2 m) const {return m.b == b;}
//    bool operator== (const int m) const {return m == b;}
//    bool operator!= (const int_z2 m) const {return m.b != b;}
//    bool operator<= (const int_z2 m) const {return m.b <= b;}
//    friend ostream& operator<<(ostream& os, const int_z2& m);
//    operator int() const {return b;}
//};
//
//ostream& operator<<(ostream& os, const int_z2& m) { os << m.b; return os; }
//
//namespace Eigen {
//template<> struct NumTraits<int_z2>
//{
//    typedef int Real;
//    typedef int_z2 Nested;
//    typedef int_z2 Literal;
//    enum {
//        IsComplex = 0,
//        IsInteger = 1,
//        IsSigned = 0,
//        RequireInitialization = 0,
//        ReadCost = 1,
//        AddCost = 2,
//        MulCost = 2
//    };
//    static inline Real epsilon() { return 0; }
//    static inline Real dummy_precision() { return 0; }
//    static inline int digits10() { return 0; }
//
//};
//}

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
    
    //    vide_decoder code(filename, engine, h_bb, h_cc);
    //    code.debug();
    
    //
    //    Eigen::Matrix<mybool, 10, 10> m = Eigen::Matrix<mybool,10,10>::Constant(0);
//    Eigen::SparseMatrix<int_z2> s(4,5);
//    //    m(1,1)=1;
//    s.insert(0,0) = 1;
//    s.insert(1,1) = 1;
//    s.insert(1,3) = 1;
//    s.insert(2,3) = 1;
////    s.insert(3,4) = 1;
//    s.insert(1,2) = 1;
//    s.makeCompressed();
//    Eigen::SparseQR<Eigen::SparseMatrix<int_z2>, Eigen::NaturalOrdering<Eigen::SparseMatrix<int_z2>::StorageIndex> > solver;
//    solver.setPivotThreshold(0);
//    solver.compute(s);
//    cout << "success" << (solver.info()==Eigen::Success) << endl;
//    Eigen::Vector<int_z2,-1> rhs = Eigen::Vector<int_z2,-1>::Constant(4,0);
//    rhs(2) = 1;
//    Eigen::Vector<int_z2,-1> sol(5);
//    sol = solver.solve(rhs);
//    cout << "success" << (solver.info()==Eigen::Success) << endl;
////    cout << solver.error() << endl;
//
////    const mybool x1 = 1;
////    const mybool x2 = 0;
////    cout << (x1 != x2) << endl;
//
////    const mybool x = 1;
////    int i = x;
////    cout << i << endl;
//
//    cout << s << endl;
//    cout << rhs << endl;
//    cout << "=====" << endl;
//    cout << sol << endl;
    
    test();
}
