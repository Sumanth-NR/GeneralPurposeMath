/**
* @file quad_residue.hpp
* @author Sumanth N R
* @note This is a part of the General Purpose Library
*/

#pragma once
#include "mod_int.hpp"

namespace Sig {

/**
 * @brief Returns a random non-quadratic residue modulo p
 */
ModInt quadraticNonResidue() {
	Int p = ModInt::M;
	// for (Int i = 17; i < 16 + p; i++) {
	// if (ModInt::pow(i, (p-1)/2) == 1)
	// continue;
	// return i;
	// }
	// return -1;
	while (true) {
		Int i = getRandom(2, p - 1);
		if (ModInt::pow(i, (p - 1) / 2) == 1)
			continue;
		return i;
	}
}

/**
 * @brief Returns the least non-negative k such that x^(2^k) = 1 (mod p) \n
 * Helper function for Tonelli-Shanks algorithm
 */
Int __tonelliFind(ModInt x) {
	for (Int ans = 0; true; x = x * x, ans++)
		if (x == 1) return ans;
}

/**
 * @brief Converts x to x^[2^k] (mod p) \n
 */
ModInt powTwoPow(const ModInt &x, const Int &k) {
	ModInt ans = x;
	for (Int i = 0; i < k; i++)
		ans = ans * ans;
	return ans;
}

/**
 * @brief Tonelli Shanks Algorithm
 * @param a
 * @return Quadratic residue of a modulo p
 */
ModInt quadraticResidueTonelliShanks(const ModInt &a) {
	const Int &p = ModInt::M;
	/**
	 * Express p-1 as 2^t * m
	 */
	Int m = p - 1, t = 0;
	while (m % 2 == 0) {
		m /= 2;
		t++;
	}

	/**
	 * Let b = a^m \n
	 * Find smallest k such that b^[2^k] = 1 \n
	 * If k = t, => a^[(p-1)/2] = -1 (mod p) \n
	 *   Hence, a is a quadratic non-residue
	 */
	ModInt b = ModInt::pow(a, m);
	Int k = __tonelliFind(b);
	if (k == t)
		throw std::invalid_argument("b is a quadratic non-residue");

	/**
	 * Set x = a^[(m+1)/2] \n
	 * If k = 0, x^2 = a (mod p) and x is one of the two solutions
	 */
	ModInt x = ModInt::pow(a, (m + 1) / 2);
	if (k == 0)
		return x;

	/**
	 * Find a quadratic non-residue r \n
	 * We know that r^[(p-1)/2] = -1 (mod p) \n
	 */
	ModInt r = quadraticNonResidue();
	ModInt s = ModInt::pow(r, m);
	ModInt S = powTwoPow(s, t - k);
	while (k > 0) {
		b *= S;
		x *= powTwoPow(s, t - k - 1);
		k = __tonelliFind(b);
		S = powTwoPow(s, t - k);
	}

	return x;
}

}