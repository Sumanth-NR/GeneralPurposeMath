#include <iostream>
#include "matrix.hpp"

#define DEBUG 3

using namespace Sig;

int main() {

#if DEBUG == 1
	Matrix<double> x({{1, 2}, {4, 5}, {6, 0}});
	Matrix<double> y({{2, 3, 6}, {6, 7, 9}});
	Matrix<double> prajwal({{1, 3, 2}, {7, 1, 11}, {9, 13, 1}});
	Matrix<double> I3 = Matrix<double>::eye(3);
	auto inv = prajwal.gaussianInverse();
	std::cout << inv << '\n';
	std::cout << prajwal * inv << '\n';

#elif DEBUG == 2
	// {{1}, {2}, {3}} is the answer
	Matrix<double> A({{1, 1, 1}, {3, 2, 1}, {9, 4, 1}});
	Matrix<double> B({{6}, {10}, {20}});
	auto X = A.solveGaussian(B);
	std::cout << X << '\n';

#elif DEBUG == 3


#endif

	return 0;
}
