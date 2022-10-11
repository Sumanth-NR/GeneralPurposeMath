#pragma once

#include "mod_int.hpp"
#include <tuple>


namespace Sig {

/**
 * @brief Generic GCD Function
 * @tparam T Typename
 * @return tuple<T, T, T> {x, y, g} such that a*x + b*y = gcd(a, b)
 *
 * @note The function requires the following
 * @note 1. Constructor to explicitly convert 0 and 1 to T
 * @note 2. Overloaded !=, /, %=, -, * operators on T
 */
template<typename T>
tuple<T, T, T> extEuclid(T a, T b) {
    // The T(0) is helpful since C++ allows only 1 level of implicit conversion
    T x = T(1), y = T(0);
    T u = T(0), v = T(1), q, temp;

    while (b != T(0)) {
        q = a / b; a %= b; std::swap(a, b);
        temp = u; u = x - q*temp; x = temp;
        temp = v; v = y - q*temp; y = temp;
    }

    return { x, y, a };  // Returning x, y, and gcd(a, b)
}

/**
 * @brief Tonelli Shanks Algorithm
 * @tparam T Typename
 *
 * (Unimplemented)
 */
template<typename T>
ModInt tonelliShanks(const ModInt &b) {
    // (p-1) = 2^k * m
    const Int &p = ModInt::M;
    Int m = p - 1, k = 1;
    while(m%2 == 0)
        m /= 2, k++;
    return T(0);
}

}
