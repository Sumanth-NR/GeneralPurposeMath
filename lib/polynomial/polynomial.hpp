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
    * @tparam T Typename
    */
    Polynomial() {
        this->cof.resize(1);
    }

    /**
    * @brief Copy Constructor
    * @tparam T Typename
    * @param p
    */
    Polynomial(const Polynomial<T> &p) {
        for (const T &i: p.cof)
            this->cof.push_back(i);
    }

    /**
    * @tparam T Typename
    * @param v The co-efficients vector where v[i] is the co-efficient of x^i
    */
    Polynomial(const std::vector<T> &v) {
        if (not v.empty()) for (const T &i: v) this->cof.push_back(i);
        else Polynomial<T>();
    }

    /**
     * @param x Initial Value
     */
    Polynomial(const T &x) {
        cof.resize(1, x);
    }

    // Some helper functions

    /**
     * @return The coefficient vector
     */
    std::vector<T> getCoEfs() {
        return cof;
    }

    /**
    * @return Degree of the Polynomial (Includes the leading 0 co-efficients if not cleaned)
    */
    [[nodiscard]] int degree() const {
        return cof.size() - 1;
    }

    /**
    * @brief Removes the initial co-efficients with 0 value
    */
    void clean() {
        while (cof.size() > 1 and cof.back() == T(0))
            cof.pop_back();
    }

    /**
    * @return The leading co-efficient
    */
    T lead() const {
        return cof.back();
    }

    /**
     * @brief Updates a co-efficient to the given value
     * @param k Index
     * @param newVal new co-efficient
     */
    void setCoef(int k, const T &newVal) {
        if (k < 0 or k > degree()) return;
        cof[k] = newVal;
    }

    /**
     * @brief Evaluates the input_helper at x
     * @param x
     * @return P(x)
     */
    virtual T eval(const T &x) {
        T ans = T(0), cur = T(1);
        for (int i = 0; i < cof.size(); i++, cur *= x)
            ans += cof[i] * cur;
        return ans;
    }

    /**
     * @brief The differential of the polynomial
     * @return
     */
    static Polynomial<T> ddx(const Polynomial<T> &p) {
        std::vector<T> nP(p.degree()-1 + 1);
        for (int i = 1; i <= p.degree(); i++)
            nP[i-1] = p.cof[i] * T(i);
        return { nP };
    }

    /**
    * Multiplies x^k to the polynomial
    * @param k Power of x to multiply, k >= 0
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
    * @param c Scalar
    */
    void operator *= (const T &c) {
        for (T &i: this->cof) i *= c;
    }

    /**
    * @brief Scalar Division
    * @param c Scalar
    */
    void operator /= (const T &c) {
        for (T &i: this->cof) i /= c;
    }

    /**
     * @brief Modulo
     * @param p Polynomial
     */
    void operator %= (const Polynomial<T> &p) {
        *this = divide(*this, p).second;
    }

    /**
     * @return The coefficient of x^i
     */
    T operator [] (const int &i) const {
        if (i < 0 or i > degree()) return T(0);
        return cof[i];
    }

    /**
    * @brief Generic Addition
    * @return Polynomial p1 + p2
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
    * @brief Generic Subtraction
    * @return Polynomial p1 - p2
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
    * @brief Generic Multiplication
    * @return Polynomial p1 * p2
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
    * @brief Generic Division, Reminder Ignored
    * @return Quotient Polynomial when p1 is divided by p2
    */
    friend Polynomial<T> operator / (const Polynomial<T>& p1, const Polynomial<T>& p2) {
        return Polynomial<T>::divide(p1, p2).first;
    }

    /**
    * @brief Generic Modulo
    * @return Reminder Polynomial when p1 is divided by p2
    */
    friend Polynomial<T> operator % (const Polynomial<T>& p1, const Polynomial<T>& p2) {
        return Polynomial<T>::divide(p1, p2).second;
    }

    /**
     * @brief Generic Equality Check
     * @param p1 Polynomial 1
     * @param p2 Polynomial 2
     * @return true if p1 == p2, false otherwise
     */
    friend bool operator == (const Polynomial<T> &p1, const Polynomial<T> &p2) {
        return p1.cof == p2.cof;
    }

    /**
     * @brief Generic Inequality Check
     * @param p1 Polynomial 1
     * @param p2 Polynomial 2
     * @return true if p1 != p2, false otherwise
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
