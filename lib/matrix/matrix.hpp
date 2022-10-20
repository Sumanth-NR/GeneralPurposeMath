/**
* @file matrix.hpp
* @author Sumanth N R
* @note This is a part of the General Purpose Library
*/

#pragma once

#include <utility>
#include <tuple>
#include <vector>
#include <functional>
#include <cmath>

using std::vector;
using std::pair;


namespace Sig {

/**
 * @brief Returns the absolute value
 * @tparam T typename
 */
template<typename T>
T abs(const T &a) { return a > T(0) ? a : -a; }

/**
 * @brief Generic Matrix Class
 * @tparam F Field
 */
template<typename F>
class Matrix {
	vector<vector<F>> mat;
protected:
	/**
	 * @brief Swaps rows ri and rj
	 */
	void swapRows(int r1, int r2) {
		F temp;
		for (int i = 0; i < getNumCols(); i++) {
			temp = mat[r1][i];
			mat[r1][i] = mat[r2][i];
			mat[r2][i] = temp;
		}
	}

	/**
	 * @brief Verifies if the matrix is a square matrix, throws an error otherwise
	 * @throw std::invalid_argument if the matrix is not square
	 */
	void verifySquare() const {
		if (getNumRows() != getNumCols())
			throw std::invalid_argument("Matrix is not square");
	}

	/**
	 * @brief Verifies if the matrix is a diagonally dominant matrix, throws an error otherwise
	 * @throw std::invalid_argument if the matrix is not diagonally dominant
	 */
	void verifyDiagonallyDominant() const {
		F norm1;
		for (int i = 0; i < getNumRows(); i++) {
			norm1 = 0;
			for (int j = 0; j < getNumCols(); j++)
				if (i != j)
					norm1 += abs(mat[i][j]);
			if (abs(mat[i][i]) < norm1)
				throw std::invalid_argument("Matrix is not diagonally dominant");
		}
	}

	/**
	 * @brief Returns the 1 norm of a column matrix
	 */
	static F norm1(const Matrix<F> &m) {
		if (m.getNumCols() != 1)
			throw std::invalid_argument("Matrix is not a column matrix");
		F ans = F(0);
		for (int i = 0; i < m.getNumRows(); i++)
			ans += abs(m.mat[i][0]);
		return ans;
	}

	/**
	 * @brief Returns the infinity norm of a row matrix
	 */
	static F normInf(const Matrix<F> &m) {
		if (m.getNumCols() != 1)
			throw std::invalid_argument("Matrix is not a column matrix");
		F ans = F(0);
		for (int i = 0; i < m.getNumCols(); i++)
			if (abs(m.mat[0][i]) > ans)
				ans = abs(m.mat[0][i]);
		return ans;
	}

	/**
	 * @brief Uses guassian row operations on A (simultaneously to B) to covert A into Upper Triangular
	 */
	static void gaussianReduceToUT(Matrix<F> &A, Matrix<F> &B) {
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
	 * @brief Uses guassian row operations on A (simultaneously to B) to covert A into Lower Triangular
	 */
	static void gaussianReduceToLT(Matrix<F> &A, Matrix<F> &B) {
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
	 * @brief Uses guassian row operations on A (simultaneously to B) to covert A into Diagonal Matrix
	 */
	static void gaussianReduceToDiag(Matrix<F> &A, Matrix<F> &B) {
		A.verifySquare();
		Matrix<F>::gaussianReduceToLT(A, B);
		Matrix<F>::gaussianReduceToUT(A, B);
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

	/**
	 * @brief General Iterative Scheme to perform X = f(X)
	 *
	 * Iteratively performs X -> AX + B until there is a convergence \n
	 * Throws an error if the number of iterations exceeds the maximum number of iterations
	 *
	 * @param X Initial Guess
	 * @param f Function to be iterated
	 * @param verbose If true, prints the number of iterations and the error at every 10th iteration
	 * @param eps Tolerance
	 * @param maxIter Maximum number of iterations
	 * @param norm Norm to be used to check for convergence
	 * @return The number of iterations required to converge
	 *
	 * @note If you want to print more frequently, change Matrix::ITER_PRINT_FREQ
	 */
	static int generalIterativeScheme(Matrix<F> &X, std::function<Matrix<F>(const Matrix<F>&)> f,
			bool verbose, F eps = F(1e-9), int maxIter = 1e5,
			std::function<F(const Matrix<F> &)> norm = Matrix<F>::norm1) {
		Matrix<F> X1;
		if (verbose)
			std::cout << "Beginning Iteration...\n";
		for (int i = 1; i <= maxIter; i++) {
			X1 = X;
			X = f(X1);
			F err = norm(X - X1);
			if (verbose and i % ITER_PRINT_FREQ == 0)
				std::cout << "Iteration " << i << " Error = " << err << '\n';
			if (norm(X - X1) < eps) {
				if (verbose)
					std::cout << "Final Error = " << err << '\n';
					std::cout << "Converged in " << i << " iterations\n";
				return i;
			}
		}
		throw std::runtime_error("Max Iterations Reached");
	}

public:
	static int ITER_PRINT_FREQ;

	/**
	 * @brief Creates a 1 by 1 matrix with value 0
	 */
	Matrix() { mat = {{0}}; }

	/**
	 * @brief Constructs a null matrix of size rowSize x colSize
	 */
	Matrix(const int &rowSize, const int &colSize) { cleanReShape(rowSize, colSize); }

	/**
	 * @brief Constructs a null matrix of size shape.first x shape.second
	 */
	explicit Matrix(const pair<int, int> &shape) { cleanReShape(shape.first, shape.second); }

	/**
	 * @brief Constructs a column matrix with the given values
	 */
	Matrix(const vector<F> &col) {
		cleanReShape(col.size(), 1);
		for (int i = 0; i < col.size(); i++)
			mat[i][0] = col[i];
	}

	/**
	 * @brief Constructs a matrix using the given vector of vectors
	 */
	Matrix(const vector<vector<F>> &matrix): mat(matrix) {}

	/**
	 * @brief Returns the reference to the element at (row, col)
	 * @param row Row index (1-indexed)
	 * @param col Column index (1-indexed)
	 * @return Reference to the element A[i][j] (1-indexed)
	 *
	 */
	F& operator () (const int &row, const int &col) {
		if (row <= 0 or row > getNumRows() or col <= 0 or col > getNumCols())
			throw std::out_of_range("Index out of range");
		return mat[row-1][col-1];
	}

	/**
	 * @brief Returns the transpose of the matrix
	 */
	Matrix<F> getTranspose() const {
		Matrix<F> ans(getNumCols(), getNumRows());
		for (int i = 0; i < getNumRows(); i++)
			for (int j = 0; j < getNumCols(); j++)
				ans.mat[j][i] = mat[i][j];
		return ans;
	}

	/**
	 * @brief Returns an identity Matrix of size n
	 */
	static Matrix<F> eye(const int &n) {
		Matrix<F> ans(n, n);
		for (int i = 0; i < n; i++) ans.mat[i][i] = 1;
		return ans;
	}

	/**
	 * @brief Returns a column matrix with all elements equal to 1
	 */
	static Matrix<F> ones(const int &n) {
		Matrix<F> ans(n, 1);
		for (int i = 0; i < n; i++)
			ans.mat[i][0] = 1;
		return ans;
	}

	/**
	 * @brief Converts the matrix into Null Matrix of size rowSize * colSize
	 */
	void cleanReShape(const int &rowSize, const int &colSize) {
		if (rowSize <= 0 or colSize <= 0)
			throw std::invalid_argument("Row Size and Column Sizes must be positive");
		mat.resize(rowSize, vector<F>(colSize));
	}

	/**
	 * @brief Returns (int) the number of rows of the matrix
	 */
	[[nodiscard]] int getNumRows() const { return mat.size(); }

	/**
	 * @brief Returns (int) the number of columns of the matrix
	 */
	[[nodiscard]] int getNumCols() const { return mat[0].size(); }

	/**
	 * @brief Returns (pair<int, int>) { numRows, numCols }
	 */
	[[nodiscard]] pair<int, int> shape() const { return { getNumRows(), getNumCols() }; }

	/**
	 * @brief Matrix Addition, Returns A + B
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
	 * @brief Matrix Subtraction, Returns A - B
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
	 * @brief Matrix Multiplication, Returns A * B
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

	/**
	 * @brief Stream Insertion
	 */
	friend std::istream& operator >> (std::istream &i, const Matrix<F> &m) {
		for(int r = 0; r < m.getNumRows(); r++)
			for(int c = 0; c < m.getNumCols(); c++)
				i >> m.mat[r][c];
		return i;
	}

	/**
	 * @brief Stream Extraction
	 */
	friend std::ostream& operator << (std::ostream &o, const Matrix<F> &m) {
		for(int r = 0; r < m.getNumRows(); r++) {
			for(int c = 0; c < m.getNumCols(); c++)
				o << m.mat[r][c] << ' ';
			o << '\n';
		}
		return o;
	}

	/**
	 * @brief Scalar Multiplication, Returns a * B
	 */
	void operator *= (const F &a) {
		for (int i = 0; i < getNumRows(); i++)
			for (int j = 0; j < getNumCols(); j++)
				mat[i][j] *= a;
	}

	/**
	 * @brief Scalar Division, performs element wise division
	 */
	void operator /= (const F &a) {
		for (int i = 0; i < getNumRows(); i++)
			for (int j = 0; j < getNumCols(); j++)
				mat[i][j] /= a;
	}

	/**
	 * @brief Uses guassian row operations to convert A to identity and B to inv(A)*B
	 * @note A is modified, if you want to preserve A, use A.gaussianInvTimes(B)
	 */
	static void guassianAInvB(Matrix<F> &A, Matrix &B) {
		gaussianReduceToDiag(A, B);
		if (A.mat[A.getNumRows()-1][A.getNumCols()-1] == F(0))
			throw std::invalid_argument("The Matrix is singular");
		for (int i = 0; i < A.getNumRows(); i++) {
			for (int j = 0; j < A.getNumCols(); j++)
				B.mat[i][j] /= A.mat[i][i];
			A.mat[i][i] = F(1);
		}
	}

	/**
	 * @brief Uses guassian row operations on A and I simultaneously and returns inv(A)
	 */
	Matrix<F> gaussianInverse() const {
		verifySquare();
		Matrix<F> B = eye(getNumRows());
		Matrix<F> A = *this;
		guassianAInvB(A, B);
		return B;
	}

	/**
	 * @brief Uses guassian row operations on A and B simultaneously and returns inv(A)*B
	 */
	Matrix<F> gaussianInvTimes(const Matrix<F> &B) const {
		verifySquare();
		if (B.getNumRows() != getNumRows())
			throw std::invalid_argument("The matrices must have the same number of rows");
		Matrix<F> A1 = *this;
		Matrix<F> B1 = B;
		guassianAInvB(A1, B1);
		return B1;
	}

	/**
	 * @brief LU Decomposition
	 * @return std::pair of matrices LU such that A = L * U (LU Decomposition)
	 */
	std::pair<Matrix<F>, Matrix<F>> LU() {
		verifySquare();
		/**
		 * Initially L = I and U = A
		 */
		Matrix<F> L = eye(getNumRows());
		Matrix<F> U = *this;
		/**
		 * Applying the transformations from top to bottom, mid to right j > i
		 */
		for (int i = 0; i < getNumRows(); i++) {
			for (int j = i + 1; j < getNumRows(); j++) {
				/**
				 * Making U[j][i] = 0 \n
				 * Set L[j][i] = U[j][i] / U[i][i] \n
				 * Then perform U[j] =  L[j][i] * U[i]
				 */
				L.mat[j][i] = U.mat[j][i] / U.mat[i][i];
				for (int k = i; k < getNumRows(); k++)
					U.mat[j][k] -= L.mat[j][i] * U.mat[i][k];
			}
		}
		/**
		 * After the transformations, L is a lower triangular matrix and U is an upper triangular matrix
		 */
		return { L, U };
	}

	/**
	 * @brief Returns the largest eigen value and eigen vector using Power Method
	 * @return std::pair<F, Matrix<F>> Largest Eigen Value and Eigen Vector
	 */
	std::pair<F, Matrix<F>> powerMethod(bool verbose = false, double eps = F(1e-9)) {
		verifySquare();
		const int n = getNumRows();
		Matrix<F> x = ones(n);

		if (verbose)
			std::cout << "Starting Power Method\n";

		const Matrix<F> &A = *this;

		/**
		 * Run the iterative scheme x1 = Ax and divide x1 by infinity norm \n
		 * This is done until the difference between 1-norm(x1) and 1-norm(x) is less than eps \n
		 */
		try {
			// iterateVector(x, *this, Matrix<F>{n, 1}, 1000, eps, true);
			generalIterativeScheme(x, [&A](const Matrix<F> &x) {
				Matrix<F> x1 = A * x;
				F a = F(0);
				for (int i = 0; i < x1.getNumRows(); i++)
					if (abs(x1.mat[i][0]) > abs(a))
						a = x1.mat[i][0];
				x1 /= a;
				return x1;
			}, verbose, eps);
		} catch (std::exception& e) {
			throw std::runtime_error("Matrix::powerMethod() : " + std::string(e.what()));
		}

		if (verbose)
			std::cout << "Terminating Power Method\n";

		/**
		 * After getting the eigen vector from the iterative scheme \n
		 * Eigen Value = a1 / a2 where \n
		 * a1 = Element with largest magnitude in A*X \n
		 * a2 = Element with largest magnitude in X \n
		 */
		auto x1 = A * x;
		F a1 = 0, a2 = 0;
		for (int i = 0; i < n; i++) {
			if (abs(x1.mat[i][0]) > abs(a1))
				a1 = x1.mat[i][0];
			if (abs(x.mat[i][0]) > abs(a2))
				a2 = x.mat[i][0];
		}

		return { a1/a2 , x };
	}

	/**
	 * @brief Solves Ax = b using Jacobi Iterative Scheme
	 * @return Returns column vector x
	 */
	Matrix<F> solveJacobi(const Matrix<F> &b, bool verbose = false, double eps = F(1e-9)) {
		verifySquare();
		verifyDiagonallyDominant();
		if (b.getNumCols() != 1)
			throw std::runtime_error("Matrix::solveGaussSeidel() : b must be a column vector");
		const int n = getNumRows();

		if (verbose)
			std::cout << "Starting Jacobi Iterative Scheme\n";

		Matrix<F> X = ones(n);
		const Matrix<F> &A = *this;

		try {
			generalIterativeScheme(X,  [&A, &b, &n](const Matrix<F> &x) {
				  Matrix<F> x1(x.shape());
				  for (int i = 0; i < n; i++) {
					  for (int j = 0; j < n; j++)
						  if (i != j)
							  x1.mat[i][0] -= A.mat[i][j] * x.mat[j][0];
					  x1.mat[i][0] += b.mat[i][0];
					  x1.mat[i][0] /= A.mat[i][i];
				  }
				  return x1;
			  }, verbose, eps);
		} catch (std::exception& e) {
			throw std::runtime_error("Matrix::solveGaussSeidel() : " + std::string(e.what()));
		}

		if (verbose)
			std::cout << "Terminating Jacobi Iterative Scheme\n";

		return X;
	}

	/**
	 * @brief Solves Solves Ax = b using Gauss-Seidel Iterative Scheme
	 * @return Returns column vector x
	 */
	Matrix<F> solveGaussSeidel(const Matrix<F> &b, bool verbose, double eps = F(1e-9)) {
		verifySquare();
		verifyDiagonallyDominant();
		if (b.getNumCols() != 1)
			throw std::runtime_error("Matrix::solveGaussSeidel() : b must be a column vector");
		const int n = getNumRows();

		if (verbose)
			std::cout << "Starting Gauss-Seidel Iterative Scheme\n";

		Matrix<F> X = ones(n);
		auto A = *this;

		try {
			generalIterativeScheme(X,  [&A, &b, &n](const Matrix<F> &x) {
				  Matrix<F> x1(x.shape());
				  for (int i = 0; i < n; i++) {
					  for (int j = 0; j < i; j++)
						  x1.mat[i][0] -= A.mat[i][j] * x1.mat[j][0];
					  for (int j = i+1; j < n; j++)
						  x1.mat[i][0] -= A.mat[i][j] * x.mat[j][0];
					  x1.mat[i][0] += b.mat[i][0];
					  x1.mat[i][0] /= A.mat[i][i];
				  }
				  return x1;
			  }, verbose, eps);
		} catch (std::exception& e) {
			throw std::runtime_error("Matrix::solveGaussSeidel() : " + std::string(e.what()));
		}

		if (verbose)
			std::cout << "Terminating Gauss-Seidel Iterative Scheme\n";

		return X;
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

		/**
		 * We need to convert A into a Tridigonal Matrix using Orthogonal Transformations \n
		 *  (i.e. B = S^T * A * S) for some orthogonal matrix S and a tridigonal matrix B \n
		 * Start from making a[i][j] = 0 from left to right and top to bottom \n
		 */
		for (int i = 0; i < n; i++) {
			for (int j = i+2; j < n; j++) {
				/**
				 * Make a[i][j] = 0 using a[i][i+1] without disturbing the other zeroes \n
				 * tan(theta) = a[i][j] / a[i][i+1] \n
				 * sin(theta) = a[i][j] / sqrt(a[i][j]^2 + a[i][i+1]^2) \n
				 * cos(theta) = a[i][i+1] / sqrt(a[i][j]^2 + a[i][i+1]^2) \n
				 */
				F d = std::sqrt(B.mat[i][j] * B.mat[i][j] + B.mat[i][i+1] * B.mat[i][i+1]);
				F s = B.mat[i][j] / d;  // sin(theta)
				F c = B.mat[i][i+1] / d;  // cos(theta)

				/**
				 * We use Orthogonal Transformations to make a[i][j] = 0 \n
				 * \n
				 * First, \n
				 * B.col[i+1] = cos * B.col[i+1] + sin * B.col[j] \n
				 * B.col[j] = -sin * B.col[i+1] + cos * B.col[j] \n
				 * are applied simultaneously \n
				 * \n
				 * Then, \n
				 * B.row[i+1] = cos * B.row[i+1] + sin * B.row[j] \n
				 * B.row[j] = -sin * B.row[i+1] + cos * B.row[j] \n
				 * are applied simultaneously \n
				 * \n
				 * Only Column Operations are done to S, \n
				 * S.col[i+1] = cos * S.col[i+1] + sin * S.col[j] \n
				 * S.col[j] = -sin * S.col[i+1] + cos * S.col[j] \n
				 * simultaneously \n
				 */
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

template <typename F>
int Matrix<F>::ITER_PRINT_FREQ = 10;

}
