#include <iostream>
#include "matrix.hpp"

#define DEBUG 6

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
	Matrix<double> A({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
	auto [L, U] = A.LU();
	std::cout << L << '\n';
	std::cout << U << '\n';
	std::cout << L * U << '\n';

#elif DEBUG == 4
	Matrix<double> A({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
	std::cout << A << '\n';
	A(1, 2) = 1;
	std::cout << A << '\n';

#elif DEBUG == 5
	Matrix<double> A({{1, 2, 3}, {2, 5, 6}, {3, 6, 9}});
	auto E = A.givensAllEigenValues();
	std::cout << E << '\n';

#elif DEBUG == 6
	Matrix<double> A({{3, 2, 4, 3},
	                  {2, 0, 2, 7},
	                  {4, 2, 3, 4},
	                  {3, 7, 4, 5}
	});
	auto E = A.givensAllEigenValues();
	std::cout << E << '\n';

#endif

	return 0;
}
