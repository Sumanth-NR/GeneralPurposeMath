/**
 * @file polynomial.hpp
 * @author Sumanth N R
 * @note This is a part of the General Purpose Library
 */

#pragma once

#include <vector>
#include <utility>
#include <ostream>


namespace Sig {

template<typename T>
T max(const T &a, const T &b) { return a > b ? a : b; }

template<typename T>
class Polynomial {
    std::vector<T> cof;

public:

    static bool prettyPrint;

    // Constructors

    /**
     * @brief Default Constructor, Creates a Polynomial of degree 0, with value 0
     */
    Polynomial() {
        this->cof.resize(1);
    }

    /**
     * @brief Copy Constructor
     */
    Polynomial(const Polynomial<T> &p) {
        for (const T &i: p.cof)
            this->cof.push_back(i);
    }

    /**
     * @brief Converts a co-efficient vector where v[i] is the
     * co-efficient of x^i to a Polynomial
     */
    Polynomial(const std::vector<T> &v) {
        if (not v.empty()) for (const T &i: v) this->cof.push_back(i);
        else Polynomial<T>();
    }

    /**
     * @brief Creates a Polynomial of degree 0, with value x
     */
    Polynomial(const T &x) {
        cof.resize(1, x);
    }

    // Some helper functions

    /**
     * @brief Returns the coefficient vector
     */
    std::vector<T> getCoEfs() {
        return cof;
    }

    /**
     * @brief Returns degree of the Polynomial (Includes the leading 0 co-efficients if not cleaned)
     */
    [[nodiscard]] int degree() const {
        return cof.size() - 1;
    }

    /**
    * @brief Removes the initial co-efficients with 0 value
    */
    Polynomial<T>& clean() {
        while (cof.size() > 1 and cof.back() == T(0))
            cof.pop_back();
		return *this;
    }

    /**
     * @brief Returns the leading co-efficient
     */
    T lead() const {
        return cof.back();
    }

	/**
	 * @brief Checks if the Polynomial is zero
	 */
	[[nodiscard]] bool isZero() const {
		return { cof.size() == 1 and cof[0] == T(0) };
	}

	/**
	 * @brief Checks if the Polynomial constant with constant val
	 */
	[[nodiscard]] bool isConstant(const T &val = 1) const {
		return { cof.size() == 1 and cof[0] == T(val) };
	}

    /**
     * @brief Updates the co-efficient to the given value
     */
    void setCoef(int k, const T &newVal) {
        if (k < 0 or k > degree()) return;
        cof[k] = newVal;
    }

    /**
     * @brief Evaluates the input_helper at x and returns P(x)
     */
    virtual T eval(const T &x) {
        T ans = T(0), cur = T(1);
        for (int i = 0; i < cof.size(); i++, cur *= x)
            ans += cof[i] * cur;
        return ans;
    }

    /**
     * @brief Returns the differential of the polynomial
     */
    static Polynomial<T> differential(const Polynomial<T> &p) {
        std::vector<T> nP(p.degree()-1 + 1);
        for (int i = 1; i <= p.degree(); i++)
            nP[i-1] = p.cof[i] * T(i);
        return { nP };
    }

    /**
     * @brief Multiplies x^k to the polynomial
     */
    void multiplyX(int k) {
        if (k <= 0) return;
        for (int i = 0; i < k; i++)
            cof.push_back(T(0));
        for (int i = degree(); i >= k; i--)
            cof[i] = cof[i-k];
        for (int i = 0; i < k; i++)
            cof[i] = T(0);
    }

    /**
    * @brief Helper Function for Division
    * @return { q, r } where q and r are quotient and remainder of p1 / p2
    */
    static std::pair<Polynomial<T>, Polynomial<T>> divide(const Polynomial<T>& _p1, const Polynomial<T>& _p2) {
        // Find q such that p1 = q*p2 + r
        Polynomial<T> p = _p1, p2 = _p2; p.clean(); p2.clean();
        if (p.degree() < p2.degree()) return { Polynomial<T>(), p };
        Polynomial<T> q; q.cof.resize(p.degree() - p2.degree() + 1);
        while (p.degree() >= p2.degree() and p != Polynomial<T>(0)) {
            T val = p.lead() / p2.lead();
            q.cof[p.degree() - p2.degree()] = val;
            for (int i = p2.degree(), j = p.degree(); i >= 0; i--, j--)
                p.cof[j] -= val * p2[i];
            p.clean();
        }
        return { q, p };
    }


    // Operator Overloads

    /**
     * @brief Scalar Multiplication
     */
    Polynomial<T>& operator *= (const T &c) {
        for (T &i: this->cof) i *= c;
		return *this;
    }

	/**
	 * @brief Polynomial Multiplication
	 */
	Polynomial<T>& operator *= (const Polynomial<T> &p) {
		*this = *this * p;
		return *this;
	}

    /**
     * @brief Scalar Division
     */
	Polynomial<T>& operator /= (const T &c) {
        for (T &i: this->cof) i /= c;
		return *this;
    }

    /**
     * @brief Polynomial Modulo
     */
	Polynomial<T>& operator %= (const Polynomial<T> &p) {
        *this = divide(*this, p).second;
		return *this;
    }

    /**
     * @brief Returns the coefficient of x^i
     */
    T operator [] (const int &i) const {
        if (i < 0 or i > degree()) return T(0);
        return cof[i];
    }

    /**
     * @brief Generic Addition, returns p1 + p2
     */
    friend Polynomial<T> operator + (const Polynomial<T>& p1, const Polynomial<T>& p2) {
        Polynomial<T> ans;
        ans.cof.resize(max<int>(p1.degree(), p2.degree()) + 1);
        for (int i = 0; i <= ans.degree(); i++)
            ans.cof[i] = p1[i] + p2[i];
        ans.clean();
        return ans;
    }

    /**
     * @brief Generic Subtraction, returns p1 - p2
     */
    friend Polynomial<T> operator - (const Polynomial<T>& p1, const Polynomial<T>& p2) {
        Polynomial<T> ans;
        ans.cof.resize(max<int>(p1.degree(), p2.degree()) + 1);
        for (int i = 0; i <= ans.degree(); i++)
            ans.cof[i] = p1[i] - p2[i];
        ans.clean();
        return ans;
    }

    /**
     * @brief Generic Multiplication, returns p1 * p2
     */
    friend Polynomial<T> operator * (const Polynomial<T>& p1, const Polynomial<T>& p2) {
        Polynomial<T> ans;
        ans.cof.resize(p1.degree() + p2.degree() + 1);
        for (int i = 0; i <= p1.degree(); i++)
            for (int j = 0; j <= p2.degree(); j++)
                ans.cof[i+j] += p1.cof[i] * p2.cof[j];
        ans.clean();
        return ans;
    }

    /**
     * @brief Generic Division with ignored reminder, returns p1 / p2
     */
    friend Polynomial<T> operator / (const Polynomial<T>& p1, const Polynomial<T>& p2) {
        return Polynomial<T>::divide(p1, p2).first;
    }

    /**
     * @brief Generic Modulo, returns p1 % p2
     */
    friend Polynomial<T> operator % (const Polynomial<T>& p1, const Polynomial<T>& p2) {
        return Polynomial<T>::divide(p1, p2).second;
    }

    /**
     * @brief Returns true if all the coefficients are the same. False otherwise
     */
    friend bool operator == (const Polynomial<T> &p1, const Polynomial<T> &p2) {
        return p1.cof == p2.cof;
    }

    /**
     * @brief Returns false if all the coefficients are the same. True otherwise
     */
    friend bool operator != (const Polynomial<T> &p1, const Polynomial<T> &p2) {
        return p1.cof != p2.cof;
    }

    /**
     * @brief Prints the Coefficients
     * @param o ostream operator
     * @param p Polynomial
     * @return ostream operator o
     */
    friend std::ostream& operator << (std::ostream &o, const Polynomial<T> &p) {
        if (prettyPrint) {
            for (int i = p.degree(); i >= 0; i--) {
                if (i != p.degree() and p.cof[i] == T(0))
                    continue;
                o << (i == p.degree() ? "" : " + ");
                if (p.cof[i] != T(1) or i == 0)
                    o << p.cof[i];
                if (i == 0)
                    continue;
                o << "x";
                if (i == 1)
                    continue;
                o << "^" << i;
            }
        } else {
            for (int i = 0; i <= p.degree(); i++)
                o << p.cof[i] << ' ';
        }
        return o;
    }

    /**
     * @brief Destructor
     */
    ~Polynomial() = default;
};

template<typename T>
bool Polynomial<T>::prettyPrint = false;

}
