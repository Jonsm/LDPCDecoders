//
//  int_z2.hpp
//  LDPCDecoders
//
//  Created by Jon on 6/1/22.
//

#ifndef int_z2_hpp
#define int_z2_hpp

#include <Eigen/Core>
#include <iostream>

class int_z2 {
public:
    bool b;
    int_z2() { b = false; }
    int_z2(bool b) : b(b) {}
    int_z2(int a) : b(a!=0) {}
    int_z2(double a) : b(a!=0) {}
    int_z2 operator* (const int_z2 m) const {return m.b & b;}
//    int_z2 operator* (const int m) const {return m * b;}
    int_z2 operator/ (const int_z2 m) const {return b;}
//    int_z2 operator/ (const int m) const {return b;}
    int_z2 operator+ (const int_z2 m) const {return m.b ^ b;}
    int_z2 operator- (const int_z2 m) const {return m.b ^ b;}
    int_z2 operator- () const {return b;}
//    int_z2 operator- (const int m) const {return m ^ b;}
    int_z2 operator+= (const int_z2 m) {b ^= m.b; return b;}
    int_z2 operator*= (const int_z2 m) {b &= m.b; return b;}
    int_z2 operator/= (const int_z2 m) {return b;}
    int_z2 operator-= (const int_z2 m) {b ^= m.b; return b;}
    bool operator== (const int_z2 m) const {return m.b == b;}
//    bool operator== (const int m) const {return m == b;}
    bool operator!= (const int_z2 m) const {return m.b != b;}
    bool operator<= (const int_z2 m) const {return m.b <= b;}
    bool operator>= (const int_z2 m) const {return m.b >= b;}
    bool operator< (const int_z2 m) const {return m.b < b;}
    bool operator> (const int_z2 m) const {return m.b > b;}
    friend std::ostream& operator<<(std::ostream& os, const int_z2& m);
//    operator int() const {return b;}
};

namespace Eigen {
template<> struct NumTraits<int_z2>
{
    typedef int_z2 Real;
    typedef int_z2 Nested;
    typedef int_z2 Literal;
    enum {
        IsComplex = 0,
        IsInteger = 1,
        IsSigned = 0,
        RequireInitialization = 0,
        ReadCost = 1,
        AddCost = 2,
        MulCost = 2
    };
    static inline Real epsilon() { return 0; }
    static inline Real dummy_precision() { return 0; }
    static inline int digits10() { return 0; }
    
};
}

int_z2 sqrt(int_z2 val);
int_z2 abs(int_z2 val);

#endif /* int_z2_hpp */
