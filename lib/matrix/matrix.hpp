#pragma once

#include <utility>
#include <tuple>
#include <vector>

using std::vector;
using std::pair;


namespace Sig {

/**
 * @brief Generic Matrix Class
 * @tparam F Field
 */
template<typename F>
class Matrix {
	vector<vector<F>> mat;
	const static Matrix<F> trash;
protected:
	void swapRows(int r1, int r2) {
		F temp;
		for (int i = 0; i < mat[i].size(); i++) {
			temp = mat[r1][i];
			mat[r1][i] = mat[r2][i];
			mat[r2][i] = temp;
		}
	}
public:
	Matrix() { mat = {{0}}; }
	Matrix(const int &rowSize, const int &colSize) { cleanReShape(rowSize, colSize); }
	Matrix(const pair<int, int> &shape) { cleanReShape(shape.first, shape.second); }
	Matrix(const vector<vector<F>> &matrix): mat(matrix) {}

	static Matrix<F> eye(const int &n) {
		Matrix<F> ans(n, n);
		for (int i = 0; i < n; i++) ans[i][i] = 1;
		return ans;
	}

	void cleanReShape(const int &rowSize, const int &colSize) {
		if (rowSize <= 0 or colSize <= 0)
			throw invalid_argument("Row Size and Column Sizes must be positive");
		mat.resize(rowSize, vector<F>(colSize));
	}

	int getNumRows() const { return mat.size(); }
	int getNumCols() const { return mat[0].size(); }
	[[nodiscard]] pair<int, int> shape() const { return { getNumRows(), getNumCols() }; }

	friend Matrix<F> operator + (const Matrix<F> &a, const Matrix<F> &b) {
		if (a.shape() != b.shape())
			throw invalid_argument("The matrices must have the same shape for adding");
		Matrix<F> ans(a.shape());
		for (int i = 0; i < a.getNumRows(); i++)
			for (int j = 0; j < a.getNumCols(); j++)
				ans.mat[i][j] = a.mat[i][j] + b.mat[i][j];
		return ans;
	}

	friend Matrix<F> operator - (const Matrix<F> &a, const Matrix<F> &b) {
		if (a.shape() != b.shape())
			throw invalid_argument("The matrices must have the same shape for adding");
		Matrix<F> ans(a.shape());
		for (int i = 0; i < a.getNumRows(); i++)
			for (int j = 0; j < a.getNumCols(); j++)
				ans.mat[i][j] = a.mat[i][j] - b.mat[i][j];
		return ans;
	}

	friend Matrix<F> operator * (const Matrix<F> &a, const Matrix<F> &b) {
		if (a.shape().second != b.shape().first)
			throw invalid_argument("The matrices must have the same shape for adding");
		Matrix<F> ans(a.getNumRows(), b.getNumCols());
		for (int i = 0; i < a.getNumRows(); i++)
			for (int j = 0; j < b.getNumCols(); j++)
				for (int k = 0; k < a.getNumCols(); k++)
					ans.mat[i][j] += a.mat[i][k] * b.mat[k][j];
		return ans;
	}

	friend std::istream& operator >> (std::istream &i, const Matrix<F> &m) {
		for(int r = 0; r < m.getNumRows(); r++)
			for(int c = 0; c < m.getNumCols(); c++)
				i >> m.mat[r][c];
		return i;
	}

	friend std::ostream& operator << (std::ostream &o, const Matrix<F> &m) {
		for(int r = 0; r < m.getNumRows(); r++) {
			for(int c = 0; c < m.getNumCols(); c++)
				o << m.mat[r][c] << ' ';
			o << '\n';
		}
		return o;
	}

	/**
	 * Calculates and returns the Transformation Matrix which when applied
	 * will convert the current matrix into Echelon Form
	 * @param inplace Applies the changes in place
	 * @return Transformation Matrix<F> T
	 */
	Matrix<F> transUTFailed(bool inplace = false) {
		if (getNumRows() < getNumCols())
			throw invalid_argument("The matrix must have more rows than columns");
		vector<vector<F>> &A = inplace ? mat : * new vector<vector<F>>(mat);
		Matrix<F> T = eye(getNumRows());
		for (int i = 0; i < getNumRows(); i++) {
			for (int j = i; A[j][i] != F(0) and  j < getNumRows(); j++) {
				if (A.mat[j][i] != 0) {
					A.swapRows(i, j);
					T.swapRows(i, j);
					break;
				}
			}
			if (A[i][i] == F(0))
				continue;
			for (int j = i + 1; j < getNumRows(); j++) {

			}
		}
		if (not inplace) delete *A;
	}

	/**
	 * Uses elementary row operations to convert the matrix A into Upper Triangular
	 *   and simultaneously applies the same operations to the matrix B
	 * @param A Matrix A
	 * @param B Matrix B
	 */
	static void transUT(Matrix<F> &A, Matrix<F> &B) {
		if (B.getNumRows() != A.getNumRows())
			throw invalid_argument("The matrices must have the same number of rows");
		// Handles at row i and column c
		for (int i = 0, c = 0; i < A.getNumRows() and c < A.getNumCols(); c++, i++) {
			// Find the first non-zero element in the k-th column
			if (A.mat[i][c] == F(0)) {
				for (int j = i + 1; j < A.getNumRows(); j++) {
					if (A.mat[j][c] != F(0)) {
						A.swapRows(i, j);
						B.swapRows(i, j);
						break;
					}
				}
			}
			if (A.mat[i][c] == F(0)) {
				i--;
				continue;
			}
			// Make all the elements in the k-th column below the i-th row zero
			for (int j = i + 1; j < A.getNumRows(); j++) {
				for (int k = c; k < A.getNumCols(); k++)
					A.mat[j][k] = A.mat[i][c] * A.mat[j][k] - A.mat[j][c] * A.mat[i][k];
				for (int k = 0; k < B.getNumCols(); k++)
					B.mat[j][k] = A.mat[i][c] * B.mat[j][k] - B.mat[j][c] * B.mat[i][k];
			}
		}
	}

	/**
	 * Uses elementary row operations to convert the matrix A into Lower Triangular
	 *   and simultaneously applies the same operations to the matrix B
	 * @param A Matrix A
	 * @param B Matrix B
	 */
	static void transLT(Matrix<F> &A, Matrix<F> &B) {
		if (B.getNumRows() != A.getNumRows())
			throw invalid_argument("The matrices must have the same number of rows");
		// Handles at row i and column c
		for (int i = A.getNumRows() - 1, c = A.getNumCols() - 1; i >= 0 and c >= 0; i--, c--) {
			// Find the first non-zero element in the k-th column
			if (A.mat[i][c] == F(0)) {
				for (int j = i - 1; j >= 0; j--) {
					if (A.mat[j][c] != F(0)) {
						A.swapRows(i, j);
						B.swapRows(i, j);
						break;
					}
				}
			}
			if (A.mat[i][c] == F(0)) {
				i++;
				continue;
			}
			// Make all the elements in the k-th column above the i-th row zero
			for (int j = i - 1; j >= 0; j--) {
				for (int k = c; k >= 0; k--)
					A.mat[j][k] = A.mat[i][c] * A.mat[j][k] - A.mat[j][c] * A.mat[i][k];
				for (int k = 0; k < B.getNumCols(); k++)
					B.mat[j][k] = A.mat[i][c] * B.mat[j][k] - B.mat[j][c] * B.mat[i][k];
			}
		}
	}

	/**
	 * Uses elementary row operations to convert the matrix A into Diagonal Matrix
	 *   and simultaneously applies the same operations to the matrix B
	 */
	static void reduceDiag(Matrix<F> &A, Matrix<F> &B) {
		Matrix<F>::transUT(A, B);
		Matrix<F>::transLT(A, B);
	}
};

}
