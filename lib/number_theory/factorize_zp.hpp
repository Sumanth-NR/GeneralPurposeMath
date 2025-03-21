/**
* @file factorize_zp.hpp
* @author Sumanth N R
* @note This is a part of the General Purpose Library
*/

#pragma once
#include <polynomial.hpp>
#include "simple_functions.hpp"
#include <map>

namespace Sig {

// Some Helper Functions


/**
 * @brief Returns a random monic Polynomial of degree = deg
 */
Polynomial<ModInt> getRandomPolynomial(const Int &deg) {
	std::vector<ModInt> p(deg+1);
	for (Int i = 0; i < deg; i++)
		p[i] = getRandom(0, ModInt::M);
	p[deg] = 1;
	return {p };
}

/**
 * @brief Returns a^n mod m
 */
Polynomial<ModInt> modExp(const Polynomial<ModInt> &a, const Polynomial<ModInt> &m, Int n) {
	Polynomial<ModInt> ans = {1};
	Polynomial<ModInt> t = a % m;  // Space : O(degree)
	for (; n; n /= 2) {
		if (n % 2 == 1)
			ans = (ans * t) % m;
		t = (t * t) % m;
	}
	return ans;
}


// Phase 1 : Square Free Factorization

/**
 * @brief Returns the square free factorization of a polynomial
 */
Polynomial<ModInt> squareFreeFactorization(const Polynomial<ModInt> &f) {
	/**
	 * Express f(x) = g(x) * [h(x)]^p where h(x) is the largest degree
	 * polynomial such that h(x)^p divides f(x)
	 *
	 * Now, square free part of g(x) = g(x) / gcd(g(x), g'(x))
	 * We recursively solve for the squarefree part of h(x)
	 *
	 * To find h(x), notice f'(x) = g'(x) [h(x)]^p + p g(x) [h(x)]^(p-1) h'(x)
	 * Since p = 0, we have, f'(x) = g'(x) [h(x)]^p
	 * So, we keep checking f''(x), f'''(x), ...
	 *
	 */

	Polynomial<ModInt> h, g = f;
	while (not g.clean().isZero()) {
		h = g;
		g = Polynomial<ModInt>::differential(g).clean();
	}

	g = f / h;
	g = g / gcd(g, Polynomial<ModInt>::differential(g));

	if (h.degree() == 0) {
		return g;
	}

	// Here, h(x) is of the form H(x^p)
	for (Int i = 1; i * ModInt::M <= h.degree(); i += 1) {
		h.setCoef(i, h[i * ModInt::M]);
		h.setCoef(i * ModInt::M, 0);
	}
	h.clean();

	return g * squareFreeFactorization(h);
}


// Phase 2 : Distinct Degree Factorization
/**
 * @brief
 * @param f
 * @note Assumes Phase 1 is complete
 */
std::map<int, Polynomial<ModInt>> distinctDegreeFactorization(Polynomial<ModInt> g) {
	std::map<int, Polynomial<ModInt>> ans;
	Polynomial<ModInt> x = std::vector<ModInt>{0, 1}, t;
	for (Int i = 1, I = ModInt::M; not g.isConstant(); i++, I *= ModInt::M) {
		t = modExp(x, g, I) - x;
		t = gcd(g, t);
		if (not t.isZero())
			ans[i] = t;
	}
	return ans;
}


// Phase 3 : Uniform-degree irreducible factorization
/**
 * @brief
 * @param f
 * @return
 */
std::vector<Polynomial<ModInt>> factorizeUniqueDegree_i(const Polynomial<ModInt> &f, Int I) {
	if (f.degree() <= I) return { f };

	/**
	 * Pick a random g(x) in Zp[x] of degree less than f(x)
	 */

	Int pPowI = 1;
	for (Int i = 1; i <= I; i++) pPowI *= ModInt::M;
	Int b = (pPowI - 1) / 2;

	Polynomial<ModInt> h1(1), h2(1), t;
	while (h1.clean().isConstant() or h2.clean().isConstant()) {
		Polynomial<ModInt> g = getRandomPolynomial(getRandom(1, f.degree()));
		t = modExp(g, f, b) - std::vector<ModInt>{1};
		h1 = gcd(f, t);
		h2 = f / h1;
	}

	std::vector<Polynomial<ModInt>> ans = factorizeUniqueDegree_i(h1, I);
	for (auto &x : factorizeUniqueDegree_i(h2, I))
		ans.push_back(x);

	return ans;
}

}
