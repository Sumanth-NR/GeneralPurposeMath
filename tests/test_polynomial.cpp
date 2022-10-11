//
// Created by sig on 10/10/22.
//

#include <iostream>
#include "polynomial.hpp"
#include "number_theory.hpp"

using std::cout, std::cin;
using namespace Sig;

int main() {
	ModInt::M = 5;
	Polynomial<ModInt> p1({-1, 0, 1}), p2({-1, 0, 0, 1});
	Polynomial<ModInt>::prettyPrint = true;

	cout << "p1: " << p1 << '\n';
	cout << "p2: " << p2 << '\n';

	auto [a, b, g] = extEuclid(p1, p2);

	cout << a << '\n' << b << '\n' << g << '\n';
	cout << a*p1 << '\n';
	cout << b*p2 << '\n';

	cout << "(" << a << ") * (" << p1 << ") + (" << b << ") * (" << p2 << ")";

	cout << '\n' << '\n';
}
