#pragma once

#include <utility>
#include <tuple>
#include <vector>
#include <cmath>

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

	/**
	 * Transforms column c1 and column c2 into \n
	 * c1 = f1*c1 + f2*c2 \n
	 * c2 = -f2*c1 + f1*c2 \n
	 *
	 * This is particularly useful when performing Orthogonal Transformations
	 * where f1 = cos(theta) for some angle theta and f2 = sin(theta)
	 *
	 * @param c1 Column 1
	 * @param c2 Column 2
	 * @param f1 Scalar ( Usually cos(theta) )
	 * @param f2 Scalar ( Usually sin(theta) )
	 */
	void transformCols(int c1, int c2, F f1, F f2) {
		if (c1 < 0 || c1 >= getNumCols() || c2 < 0 || c2 >= getNumCols())
			throw std::invalid_argument("Invalid column index");
		F x1, x2;
		for (int i = 0; i < getNumRows(); i++) {
			x1 = mat[i][c1], x2 = mat[i][c2];
			mat[i][c1] = f1 * x1 + f2 * x2;
			mat[i][c2] = f1 * x2 - f2 * x1;
		}
	}
	/**
	 * Transforms row r1 and row r2 into \n
	 * r1 = f1*r1 + f2*r2 \n
	 * r2 = -f2*r1 + f1*r2 \n
	 *
	 * This is particularly useful when performing Orthogonal Transformations
	 * where f1 = cos(theta) for some angle theta and f2 = sin(theta)
	 *
	 * @param r1 Row 1
	 * @param r2 Row 2
	 * @param f1 Scalar ( Usually cos(theta) )
	 * @param f2 Scalar ( Usually sin(theta) )
	 */
	void transformRows(int r1, int r2, F f1, F f2) {
		if (r1 < 0 || r1 >= getNumCols() || r2 < 0 || r2 >= getNumCols())
			throw std::invalid_argument("Invalid column index");
		F x1, x2;
		for (int i = 0; i < getNumCols(); i++) {
			x1 = mat[r1][i], x2 = mat[r2][i];
			mat[r1][i] = f1 * x1 + f2 * x2;
			mat[r2][i] = f1 * x2 - f2 * x1;
		}
	}
public:
	Matrix() { mat = {{0}}; }
	Matrix(const int &rowSize, const int &colSize) { cleanReShape(rowSize, colSize); }
	Matrix(const pair<int, int> &shape) { cleanReShape(shape.first, shape.second); }
	Matrix(const vector<vector<F>> &matrix): mat(matrix) {}

	/**
	 * @brief Returns the reference to the element at (row, col)
	 * @param i Row index (1-indexed)
	 * @param j Column index (1-indexed)
	 * @return Reference to the element A[i][j] (1-indexed)
	 */
	F& operator () (const int &row, const int &col) {
		if (row <= 0 or row > getNumRows() or col <= 0 or col > getNumCols())
			throw std::out_of_range("Index out of range");
		return mat[row-1][col-1];
	}

	/**
	 * @return The transpose of the matrix
	 */
	Matrix<F> getTranspose() const {
		Matrix<F> ans(getNumCols(), getNumRows());
		for (int i = 0; i < getNumRows(); i++)
			for (int j = 0; j < getNumCols(); j++)
				ans.mat[j][i] = mat[i][j];
		return ans;
	}

	/**
	 * @return Returns an identity Matrix of size n
	 */
	static Matrix<F> eye(const int &n) {
		Matrix<F> ans(n, n);
		for (int i = 0; i < n; i++) ans.mat[i][i] = 1;
		return ans;
	}

	/**
	 * @brief Converts into Null Matrix of size rowSize * colSize
	 * @param rowSize
	 * @param colSize
	 */
	void cleanReShape(const int &rowSize, const int &colSize) {
		if (rowSize <= 0 or colSize <= 0)
			throw std::invalid_argument("Row Size and Column Sizes must be positive");
		mat.resize(rowSize, vector<F>(colSize));
	}

	[[nodiscard]] int getNumRows() const { return mat.size(); }
	[[nodiscard]] int getNumCols() const { return mat[0].size(); }
	[[nodiscard]] pair<int, int> shape() const { return { getNumRows(), getNumCols() }; }

	/**
	 * @brief Matrix Addition
	 * @param A Matrix A
	 * @param B Matrix B
	 * @return A * B
	 */
	friend Matrix<F> operator + (const Matrix<F> &A, const Matrix<F> &B) {
		if (A.shape() != B.shape())
			throw std::invalid_argument("The matrices must have the same shape for adding");
		Matrix<F> ans(A.shape());
		for (int i = 0; i < A.getNumRows(); i++)
			for (int j = 0; j < A.getNumCols(); j++)
				ans.mat[i][j] = A.mat[i][j] + B.mat[i][j];
		return ans;
	}

	/**
	 * @brief Matrix Subtraction
	 * @param A Matrix A
	 * @param B Matrix B
	 * @return A * B
	 */
	friend Matrix<F> operator - (const Matrix<F> &A, const Matrix<F> &B) {
		if (A.shape() != B.shape())
			throw std::invalid_argument("The matrices must have the same shape for adding");
		Matrix<F> ans(A.shape());
		for (int i = 0; i < A.getNumRows(); i++)
			for (int j = 0; j < A.getNumCols(); j++)
				ans.mat[i][j] = A.mat[i][j] - B.mat[i][j];
		return ans;
	}

	/**
	 * @brief Matrix multiplication
	 * @param A Matrix A
	 * @param B Matrix B
	 * @return A * B
	 */
	friend Matrix<F> operator * (const Matrix<F> &A, const Matrix<F> &B) {
		if (A.shape().second != B.shape().first)
			throw std::invalid_argument("The matrix sizes not appropriate for multiplication");
		Matrix<F> ans(A.getNumRows(), B.getNumCols());
		for (int i = 0; i < A.getNumRows(); i++)
			for (int j = 0; j < B.getNumCols(); j++)
				for (int k = 0; k < A.getNumCols(); k++)
					ans.mat[i][j] += A.mat[i][k] * B.mat[k][j];
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
	 * @param A Matrix A
	 * @param B Matrix B
	 */
	static void reduceDiag(Matrix<F> &A, Matrix<F> &B) {
		A.verifySquare();
		Matrix<F>::gaussLT(A, B);
		Matrix<F>::gaussUT(A, B);
	}

	/**
	 * Uses elementary row operations to convert the matrix A into Diagonal Matrix
	 * @param A Matrix
	 * @param B Matrix
	 */
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

	/**
	 * @return The inverse of the matrix A using Gauss-Jordan Elimination
	 */
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

	/**
	 * Given's Method to find all the eigenvalues of a matrix
	 *   (for symmetric matrices)\n
	 * Only considers the Upper Triangular part of the matrix
	 *   (if the matrix is not symmetric)\n
	 * @return std::vector<F> List of all the eigen-values of A
	 */
	[[maybe_unused]] Matrix<F> givensAllEigenValues() {
		verifySquare();
		const int n = getNumRows();
		Matrix<F> B = *this;
		Matrix<F> S = eye(n);

		// Converting A into tridigonal Matrix
		for (int i = 0; i < n; i++) {
			for (int j = i+2; j < n; j++) {
				// Make aij = 0

				// tan(theta) = aij / ai(i+1)
				F d = std::sqrt(B.mat[i][j] * B.mat[i][j] + B.mat[i][i+1] * B.mat[i][i+1]);
				F s = B.mat[i][j] / d;  // sin(theta)
				F c = B.mat[i][i+1] / d;  // cos(theta)

				// Orthogonal Transformations
				B.transformCols(i+1, j, c, s);
				B.transformRows(i+1, j, c, s);
				S.transformCols(i+1, j, c, s);
			}
		}

		/**
		 * Now, S is an Orthogonal Matrix, B is a tridigonal matrix and
		 * S' * A * S = B
		 */

		return B;
	}
};

}
