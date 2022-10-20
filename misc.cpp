#include <iostream>
#include <functional>
#include <sig.hpp>
using namespace std;
using namespace Sig;

using F = double;

void printAfter1(Matrix<F> x, std::function<Matrix<F>(const Matrix<F>&)> &f) {
	x = f(x);
	cout << x;
}

int main() {
	int y = 10;

	function<int(int)> intFunc = [&y] (int x) -> int {
		return x + y;
	};

	int n = 10;
	function<Matrix<F>(const Matrix<F>&)> func = [&n] (const Matrix<F> &x) {
		auto y = x;
		for (int i = 0; i < y.getNumRows(); i++) {
			for (int j = 0; j < y.getNumCols(); j++)
				y(i+1, j+1) += n;
		}
		return y;
	};

	Matrix<F> X = vector<F>{1, 2, 4};
	printAfter1(X, func);
}
