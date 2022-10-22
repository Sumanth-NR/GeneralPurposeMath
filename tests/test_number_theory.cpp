#include <iostream>
#include <number_theory.hpp>

using namespace Sig;

#define MOD_INT_DEBUG 1

int main() {

#if MOD_INT_DEBUG == 0
	ModInt::M = 1e9 + 7;
	Int n = 2000045;
	try {
		auto x = quadraticResidueTonelliShanks(n);
		cout << x << '\n';
		cout << x * x << '\n';
	} catch (const std::exception &e) {
		cout << e.what() << '\n';
		cout << ModInt::pow(n, (ModInt::M - 1) / 2) << '\n';
	}

#elif MOD_INT_DEBUG == 1
	ModInt::M = 421;
	Int counter = 0;
	for (Int i = 1; i < ModInt::M; i++) {
		try {
			quadraticResidueTonelliShanks(i);
			std::cout << 1;
			counter++;
		} catch (const std::exception &e) {
			std::cout << 0;
		}
		if (i % 100 == 0)
			std::cout << '\n';
	}
	std::cout << '\n' << counter << '\n';

#endif

	return 0;
}
