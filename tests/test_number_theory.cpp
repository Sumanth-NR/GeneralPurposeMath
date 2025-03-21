#include <iostream>
#include <number_theory.hpp>

using namespace Sig;

#define NT_DEBUG 2

int main() {

#if NT_DEBUG == 0
	ModInt::M = 1e9 + 7;
	Int n = 20000;
	try {
		auto x = quadraticResidueTonelliShanks(n);
		std::cout << x << '\n';
		std::cout << x * x << '\n';
	} catch (const std::exception &e) {
		std::cout << e.what() << '\n';
		std::cout << ModInt::pow(n, (ModInt::M - 1) / 2).getVal() - ModInt::M << '\n';
	}
	n = 200004;
	try {
		auto x = quadraticResidueTonelliShanks(n);
		std::cout << x << '\n';
		std::cout << x * x << '\n';
	} catch (const std::exception &e) {
		std::cout << e.what() << '\n';
		std::cout << ModInt::pow(n, (ModInt::M - 1) / 2).getVal() - ModInt::M << '\n';
	}

#elif NT_DEBUG == 1
	ModInt::M = 421;  // prime
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

#elif NT_DEBUG == 2
	ModInt::M = 1e9 + 7;
	Int n = 5;
	try {
		auto x = quadraticResidueTonelliShanks(n);
		std::cout << x << '\n';
		std::cout << x * x << '\n';
	} catch (const std::exception &e) {
		std::cout << e.what() << '\n';
		std::cout << ModInt::pow(n, (ModInt::M - 1) / 2).getVal() - ModInt::M << '\n';
	}
#endif

	return 0;
}
