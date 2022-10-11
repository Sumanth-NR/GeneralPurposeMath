#include <iostream>
#include <sig.hpp>
using namespace std;
using namespace Sig;

template<typename F>
void func(const Matrix<F> &mt) {
    // Do shit
    auto [m, n] = mt.shape();
    cout << m << ' ' << n << '\n';
}

int main() {
    Matrix<int> m(5, 8);
    Matrix<double> x({{1, 2}, {4, 5}, {6, 0}});
    Matrix<double> y({{2, 3, 6}, {6, 7, 9}});
    cout << x*y;
}
