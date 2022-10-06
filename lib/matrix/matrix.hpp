#pragma once

#include <utility>
#include <tuple>
#include <vector>

using std::vector;
using std::pair;


namespace Sig {

template<typename F>
class Matrix {
    int m, n;
    vector<vector<F>> mat;
public:
    Matrix() { m = 0, n = 0; }
    Matrix(const int &rowSize, const int &colSize): m(rowSize), n(colSize) { mat.resize(m, vector<F>(n)); }
    Matrix(const vector<vector<F>> &matrix): mat(matrix) {}

    void cleanReShape(int rowSize, int colSize) {
        m = rowSize, n = colSize;
        mat.resize(m, vector<F>(n));
    }

    [[nodiscard]] pair<int, int> shape() const { return { m, n }; }

    friend Matrix<F> operator + (const Matrix<F> &a, const Matrix<F> &b) {
        if (a.shape().second != b.shape().first) {
        }
    }


};

}
