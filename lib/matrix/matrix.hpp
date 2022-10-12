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
		for (int i = 0; i < getNumCols(); i++) {
			temp = mat[r1][i];
			mat[r1][i] = mat[r2][i];
			mat[r2][i] = temp;
		}
	}

	void verifySquare() const {
		if (getNumRows() != getNumCols())
			throw std::invalid_argument("Matrix is not square");
	}
public:
	Matrix() { mat = {{0}}; }
	Matrix(const int &rowSize, const int &colSize) { cleanReShape(rowSize, colSize); }
	Matrix(const pair<int, int> &shape) { cleanReShape(shape.first, shape.second); }
	Matrix(const vector<vector<F>> &matrix): mat(matrix) {}

	static Matrix<F> eye(const int &n) {
		Matrix<F> ans(n, n);
		for (int i = 0; i < n; i++) ans.mat[i][i] = 1;
		return ans;
	}

	void cleanReShape(const int &rowSize, const int &colSize) {
		if (rowSize <= 0 or colSize <= 0)
			throw std::invalid_argument("Row Size and Column Sizes must be positive");
		mat.resize(rowSize, vector<F>(colSize));
	}

	[[nodiscard]] int getNumRows() const { return mat.size(); }
	[[nodiscard]] int getNumCols() const { return mat[0].size(); }
	[[nodiscard]] pair<int, int> shape() const { return { getNumRows(), getNumCols() }; }

	friend Matrix<F> operator + (const Matrix<F> &a, const Matrix<F> &b) {
		if (a.shape() != b.shape())
			throw std::invalid_argument("The matrices must have the same shape for adding");
		Matrix<F> ans(a.shape());
		for (int i = 0; i < a.getNumRows(); i++)
			for (int j = 0; j < a.getNumCols(); j++)
				ans.mat[i][j] = a.mat[i][j] + b.mat[i][j];
		return ans;
	}

	friend Matrix<F> operator - (const Matrix<F> &a, const Matrix<F> &b) {
		if (a.shape() != b.shape())
			throw std::invalid_argument("The matrices must have the same shape for adding");
		Matrix<F> ans(a.shape());
		for (int i = 0; i < a.getNumRows(); i++)
			for (int j = 0; j < a.getNumCols(); j++)
				ans.mat[i][j] = a.mat[i][j] - b.mat[i][j];
		return ans;
	}

	friend Matrix<F> operator * (const Matrix<F> &a, const Matrix<F> &b) {
		if (a.shape().second != b.shape().first)
			throw std::invalid_argument("The matrix sizes not appropriate for multiplication");
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
			throw std::invalid_argument("The matrix must have more rows than columns");
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
	static void gaussUT(Matrix<F> &A, Matrix<F> &B) {
		if (B.getNumRows() != A.getNumRows())
			throw std::invalid_argument("The matrices must have the same number of rows");
		// Handles at row i1 and column j1
		for (int i1 = 0, j1 = 0; i1 < A.getNumRows() and j1 < A.getNumCols(); j1++, i1++) {
			// Find the first non-zero element in the j2-th column
			for (int i2 = i1; i2 < A.getNumRows(); i2++) {
				if (A.mat[i2][j1] != F(0)) {
					A.swapRows(i1, i2);
					B.swapRows(i1, i2);
					break;
				}
			}
			if (A.mat[i1][j1] == F(0)) {
				i1--;
				continue;
			}
			// Make all the elements in the j2-th column below the i1-th row zero
			for (int i2 = i1 + 1; i2 < A.getNumRows(); i2++) {
				for (int j2 = j1 + 1; j2 < A.getNumCols(); j2++)
					A.mat[i2][j2] -= A.mat[i2][j1] * A.mat[i1][j2] / A.mat[i1][j1];
				for (int j2 = 0; j2 < B.getNumCols(); j2++)
					B.mat[i2][j2] -= A.mat[i2][j1] * B.mat[i1][j2] / A.mat[i1][j1];
				A.mat[i2][j1] = F(0);
			}
		}
	}

	/**
	 * Uses elementary row operations to convert the matrix A into Lower Triangular
	 *   and simultaneously applies the same operations to the matrix B
	 * @param A Matrix A
	 * @param B Matrix B
	 */
	static void gaussLT(Matrix<F> &A, Matrix<F> &B) {
		if (B.getNumRows() != A.getNumRows())
			throw std::invalid_argument("The matrices must have the same number of rows");
		// Handles at row i1 and column j1
		for (int i1 = A.getNumRows() - 1, j1 = A.getNumCols() - 1; i1 >= 0 and j1 >= 0; i1--, j1--) {
			// Find the first non-zero element in the k-th column
			for (int i2 = i1; i2 >= 0; i2--) {
				if (A.mat[i2][j1] != F(0)) {
					A.swapRows(i1, i2);
					B.swapRows(i1, i2);
					break;
				}
			}
			if (A.mat[i1][j1] == F(0)) {
				i1++;
				continue;
			}
			// Make all the elements in the k-th column above the i-th row zero
			for (int i2 = i1 - 1; i2 >= 0; i2--) {
				for (int j2 = j1 - 1; j2 >= 0; j2--)
					A.mat[i2][j2] -= A.mat[i2][j1] * A.mat[i1][j2] / A.mat[i1][j1];
				for (int j2 = 0; j2 < B.getNumCols(); j2++)
					B.mat[i2][j2] -= A.mat[i2][j1] * B.mat[i1][j2] / A.mat[i1][j1];
				A.mat[i2][j1] = F(0);
			}
		}
	}

	/**
	 * Uses elementary row operations to convert the matrix A into Diagonal Matrix
	 *   and simultaneously applies the same operations to the matrix B
	 */
	static void reduceDiag(Matrix<F> &A, Matrix<F> &B) {
		A.verifySquare();
		Matrix<F>::gaussLT(A, B);
		Matrix<F>::gaussUT(A, B);
	}

	static void reduceIdentity(Matrix<F> &A, Matrix &B) {
		reduceDiag(A, B);
		if (A.mat[A.getNumRows()-1][A.getNumCols()-1] == F(0))
			throw std::invalid_argument("The Matrix is singular");
		for (int i = 0; i < A.getNumRows(); i++) {
			for (int j = 0; j < A.getNumCols(); j++)
				B.mat[i][j] /= A.mat[i][i];
			A.mat[i][i] = F(1);
		}
	}

	Matrix<F> gaussianInverse() {
		verifySquare();
		Matrix<F> B = eye(getNumRows());
		Matrix<F> A = *this;
		reduceIdentity(A, B);
		return B;
	}

	/**
	 * Solves X = A^-1B using Gaussian Elimination
	 * @param B Row matrix having the same number of
	 * @return A^-1 * B
	 */
	Matrix<F> solveGaussian(const Matrix<F> &B) {
		verifySquare();
		if (B.getNumRows() != getNumRows())
			throw std::invalid_argument("The matrices must have the same number of rows");
		Matrix<F> A1 = *this;
		Matrix<F> B1 = B;
		reduceIdentity(A1, B1);
		return B1;
	}

	/**
	 * LU Decomposition
	 * @return A pair of matrices (L, U) such that A = L * U
	 */
	std::pair<Matrix<F>, Matrix<F>> LU() {
		verifySquare();
		Matrix<F> L = eye(getNumRows());
		Matrix<F> U = *this;
		for (int i = 0; i < getNumRows(); i++) {
			for (int j = i + 1; j < getNumRows(); j++) {
				L.mat[j][i] = U.mat[j][i] / U.mat[i][i];
				for (int k = i; k < getNumRows(); k++)
					U.mat[j][k] -= L.mat[j][i] * U.mat[i][k];
			}
		}
		return { L, U };
	}
};

}
