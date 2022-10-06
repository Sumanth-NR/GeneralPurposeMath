/**
 * @file mod_int.hpp
 * @author Sumanth N R
 * @note This is a part of the General Purpose Library
 */

#pragma once

#include <iostream>
using namespace std;

using Int = long long int;


namespace Sig {

class ModInt {
    Int val;
protected:
    void normalize() { val %= M; if (val < 0) val = M + val; }

public:
    static Int M;
    static Int MOD() { return M; }

    ModInt(): val(0) {}
    ModInt(const int &x) { val = x; normalize(); }
    ModInt(const Int& x) { val = x; normalize(); }
    ModInt(const ModInt &x) { val = x.val; }

    [[nodiscard]] Int getVal() const { return val; }

    void operator += (const ModInt& x) { val += x.getVal(); normalize(); }
    void operator -= (const ModInt& x) { val -= x.getVal(); normalize(); }
    void operator *= (const ModInt& x) { val *= x.getVal(); normalize(); }

    friend bool operator == (const ModInt &m1, const ModInt &m2) { return m1.val == m2.val; }
    friend bool operator != (const ModInt &m1, const ModInt &m2) { return m1.val != m2.val; }

    friend ModInt operator + (const ModInt &m1, const ModInt &m2) { return { m1.val + m2.val }; }
    friend ModInt operator - (const ModInt &m1, const ModInt &m2) { return { m1.val - m2.val }; }
    friend ModInt operator * (const ModInt &m1, const ModInt &m2) { return { m1.val * m2.val }; }

    static ModInt pow(const ModInt& a, const Int &k) {
        ModInt ans = 1, t = a;
        for (Int p = k; p; p /= 2, t *= t) if (p % 2 == 1) ans *= t;
        return ans;
    }

    static ModInt fermatInv(const ModInt &x) { return pow(x, M - 2); }

    void operator /= (const ModInt &x) { *this *= fermatInv(x); }
    friend ModInt operator / (const ModInt &m1, const ModInt &m2) { return { m1 * fermatInv(m2) }; }

    friend std::istream& operator >> (std::istream &i, ModInt &m) { i >> m.val; return i; }
    friend std::ostream& operator << (std::ostream &o, const ModInt &m) { o << m.val; return o; }

    ~ModInt() = default;
};
Int ModInt::M = 1;

}
