#include "sig.hpp"
#include <iostream>

using namespace std;
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

//    auto q1 = p1 / Polynomial<ModInt>({1, 1});
//    auto q2 = p2 / Polynomial<ModInt>({1, 1});
//    cout << q1 << '\n' << q2 << '\n';
//    auto [a1, b1, g1] = extEuclid(q1, q2);
//    cout << a1 << '\n';
//    cout << b1 << '\n';
//    cout << g1 << '\n';

    cout.flush();
    return 0;
}
