#include <iostream>
#include "matrix.hpp"

using namespace std;
using namespace Sig;

int main() {
	Matrix<int> m(5, 8);
	Matrix<double> x({{1, 2}, {4, 5}, {6, 0}});
	Matrix<double> y({{2, 3, 6}, {6, 7, 9}});
	cout << x*y;
}
