/**
 * @file simple_functions.hpp
 * @author Sumanth N R
 * @note This is a part of the General Purpose Library
 */

#pragma once
#include "mod_int.hpp"
#include <tuple>

namespace Sig {

/**
 * @brief Generic Extended Euclidian Algorithm
 * @tparam F Field
 * @return tuple<F, F, F> {x, y, g} such that a*x + b*y = gcd(a, b)
 *
 * @note The function requires the following
 * @note 1. Constructor to explicitly convert 0 and 1 to F
 * @note 2. Overloaded !=, /, %=, -, * operators on F
 */
template<typename F>
std::tuple<F, F, F> extEuclid(F a, F b) {
	// The F(0) is helpful since C++ allows only 1 level of implicit conversion
	F x = F(1), y = F(0);
	F u = F(0), v = F(1), q, temp;

	while (b != F(0)) {
		q = a / b;
		a %= b;
		std::swap(a, b);
		temp = u;
		u = x - q * temp;
		x = temp;
		temp = v;
		v = y - q * temp;
		y = temp;
	}

	return {x, y, a};// Returning x, y, and gcd(a, b)
}

/**
 * @brief Returns the GCD of a and b
 */
template<typename F>
F gcd(const F &a, const F &b) {
	return std::get<2>(extEuclid(a, b));
}

/**
 * @brief Returns a random int64 in the range [min, max)
 */
Int getRandom(Int const& min = 0, Int const& max = 0) {
	return (((Int)(unsigned int)rand() << 32) + (Int)(unsigned int)rand()) % (max - min) + min;
}

}