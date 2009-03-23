/***************************************************************************
 *            taylor_variable.h
 *
 *  Copyright 2008  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file taylor_variable.h
 *  \brief Differential algebra variables with a sparse representation.
 */

#ifndef ARIADNE_TAYLOR_VARIABLE_H
#define ARIADNE_TAYLOR_VARIABLE_H

#include <map>

#include "macros.h"
#include "array.h"
#include "vector.h"
#include "multi_index.h"
#include "polynomial.h"

namespace Ariadne {

template<class T1, class T2> class Product;
template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Polynomial;

class TaylorVariable;
template<> class Vector<TaylorVariable>;

// Split the variable over two domains, subdividing along the independent variable j.
pair<TaylorVariable,TaylorVariable> split(const TaylorVariable& x, uint j);
pair< Vector<TaylorVariable>, Vector<TaylorVariable> > split(const Vector<TaylorVariable>& x, uint j);

// Scale the variabe by post-composing with an affine map taking the interval \a ivl to the unit interval
TaylorVariable unscale(const TaylorVariable& x, const Interval& ivl);
Vector<TaylorVariable> unscale(const Vector<TaylorVariable>& x, const Vector<Interval>& bx);

// Scale the variable by post-composing with an affine map taking the unit interval to \a ivl.
TaylorVariable scale(const TaylorVariable& x, const Interval& ivl);
Vector<TaylorVariable> scale(const Vector<TaylorVariable>& x, const Vector<Interval>& bx);

// Evaluate an array of Taylor variables on a vector.
Vector<Interval> evaluate(const Vector<TaylorVariable>& x, const Vector<Interval>& sy);
Interval evaluate(const TaylorVariable& x, const Vector<Interval>& sy);


// Compose an array of Taylor variables with another, after scaling by the interval vectors
Vector<TaylorVariable> compose(const Vector<TaylorVariable>& x, const Vector<Interval>& bx, const Vector<TaylorVariable>& y);

// Wrappers for univariate composition
TaylorVariable compose(const TaylorVariable& x, const Vector<Interval>& bx, const Vector<TaylorVariable>& y);
TaylorVariable compose(const TaylorVariable& x, const Interval& b, const TaylorVariable& y);


// Antidifferentiation operator
TaylorVariable antiderivative(const TaylorVariable& x, const Interval& dk, uint k);
Vector<TaylorVariable> antiderivative(const Vector<TaylorVariable>& x, const Interval& dk, uint k);

// Embed the variable in a space of higher dimension
TaylorVariable embed(const TaylorVariable& tv, uint as, uint b);
Vector<TaylorVariable> embed(const Vector<TaylorVariable>& tvs, uint as, uint b);

// Test if a variable refines another
bool refines(const TaylorVariable& tv1, const TaylorVariable& tv2);
bool refines(const Vector<TaylorVariable>& tv1, const Vector<TaylorVariable>& tv2);



/*! \brief A class representing a quantity depending on other quantities. 
 *
 * See also TaylorFunction, TaylorSet.
 */
class TaylorVariable
{
    typedef Polynomial<Float> PolynomialType;
    //typedef MapPolynomial PolynomialType;
    static const Float _zero;
    PolynomialType _polynomial;
    Float _error;
    double _sweep_threshold;
    uint _maximum_degree;
  private:
    static double _default_sweep_threshold;
    static uint _default_maximum_degree;
  public:
    static const double em;
    static const double ec;
  public:
    typedef Float ScalarType;
    //! \brief The type used to index the coefficients.
    typedef MultiIndex IndexType;
    //! \brief The type used for the coefficients.
    typedef Float ValueType;
    //! \brief An iterator through the (index,coefficient) pairs of the polynomial polynomial.
    typedef PolynomialType::iterator iterator;
    //! \brief A constant iterator through the (index,coefficient) pairs of the polynomial polynomial.
    typedef PolynomialType::const_iterator const_iterator;

    //@{
    /*! \name Constructors and destructors. */
    //! \brief Default constructor.
    TaylorVariable();
    //! \brief Construct a TaylorVariable in \a as arguments.
    TaylorVariable(uint as);
    //! \brief Construct from a map giving the polynomial polynomial and a constant giving the error.
    TaylorVariable(const std::map<MultiIndex,Float>& d, const Float& e);
    //! \brief Contruct a %TaylorVariable in \a as arguments of degree \a d
    //! from the raw data given by \a ptr, with error given by \a err.
    TaylorVariable(uint as, uint deg, const double* ptr, const double& err);
    //! \brief Construct a %TaylorVariable in \a as arguments of degree \a deg,
    //! with values given by the double-precision arguments. The last arguments
    //! is the error of the polynomial approximation.
    TaylorVariable(uint as, uint deg, double d0, ...);
    //! \brief Construct a %TaylorVariable in \a as arguments with \a nnz nonzero terms,
    //! given in a list of indices a[0],...,a[k-1] followed by the coefficient c.
    //! The last argument gives the error of the polynomial approximation.
    TaylorVariable(uint as, uint nnz, uint a0, ...);
    //@}

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to a constant, keeping the same number of arguments.
    TaylorVariable& operator=(const Float& c);
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorVariable& operator=(const Interval& c);
    //! \brief Test if the quantity is a better approximation than \a t throughout the domain.
    bool refines(const TaylorVariable& t);
    //@}

    //@{
    /*! \name Data access */
    //! \brief The polynomial polynomial.
    const PolynomialType& polynomial() const { return this->_polynomial; }
    //! \brief A reference to the polynomial polynomial.
    PolynomialType& polynomial() { return this->_polynomial; }
    //! \brief The error of the polynomial polynomial over the domain.
    const Float& error() const { return this->_error; }
    //! \brief A reference to the error of the polynomial polynomial over the domain.
    Float& error() { return this->_error; }
    //! \brief The constant term in the polynomial polynomial.
    const Float& value() const { return (*this)[MultiIndex::zero(this->argument_size())]; }
    //! \brief A reference to the constant term in the polynomial polynomial.
    Float& value() { return (*this)[MultiIndex::zero(this->argument_size())]; }
    //! \brief The coefficient of the gradient term \f$df/dx_j\f$.
    const Float& gradient(uint j) const { return (*this)[MultiIndex::unit(this->argument_size(),j)]; }
    //! \brief A reference to the coefficient of the gradient term \f$df/dx_j\f$.
    Float& gradient(uint j) { return (*this)[MultiIndex::unit(this->argument_size(),j)]; }

    //! \brief Set the error of the polynomial polynomial.
    void set_error(const Float& ne) {
        ARIADNE_ASSERT(ne>=0); this->_error=ne; }
    //! \brief Set the constant term in the polynomial polynomial.
    void set_value(const Float& c) {
        this->_polynomial[MultiIndex::zero(this->argument_size())]=c; }
    //! \brief Set the coefficient of the term \f$df/dx_j\f$.
    void set_gradient(uint j, const Float& c) {
        this->_polynomial[MultiIndex::unit(this->argument_size(),j)]=c; }

    //! \brief The coefficient of the term in $x^a$.
    const Float& operator[](const MultiIndex& a) const { return this->_polynomial[a]; }
    //! \brief A read/write reference to the coefficient of the term in $x^a$.
    Float& operator[](const MultiIndex& a) { return this->_polynomial[a]; }

    //! \brief The coefficient of the term \f$df/dx_j\f$.
    const Float& operator[](uint j) const {
        return (*this)[MultiIndex::unit(this->argument_size(),j)]; }
    //! \brief A reference to the coefficient of the term \f$df/dx_j\f$.
    Float& operator[](uint j) {
        return (*this)[MultiIndex::unit(this->argument_size(),j)]; }

    //! \brief An iterator to the first term in the polynomial polynomial.
    iterator begin() { return this->_polynomial.begin(); }
    //! \brief A constant iterator to the first term in the polynomial polynomial.
    const_iterator begin() const { return this->_polynomial.begin(); }
    //! \brief An iterator to the end of the polynomial polynomial.
    iterator end() { return this->_polynomial.end(); }
    //! \brief A constant iterator to the end of the polynomial polynomial.
    const_iterator end() const { return this->_polynomial.end(); }
    //! \brief An iterator to the term with index \a.
    iterator find(const MultiIndex& a) { return this->_polynomial.find(a); }
    //! \brief A constant iterator to the term with index \a.
    const_iterator find(const MultiIndex& a) const { return this->_polynomial.find(a); }

    //! \brief The number of variables in the argument of the quantity.
    uint argument_size() const { return this->_polynomial.argument_size(); }
    //! \brief The maximum degree of terms in the polynomial polynomial.
    uint degree() const { return (--this->_polynomial.end())->first.degree(); }
    //! \brief The number of nonzero terms in the polynomial polynomial.
    uint nnz() const { return this->_polynomial.size(); }
    //@}

    //@{
    /*! \name Named constructors. */
    //! \brief Construct the zero quantity in \a as independent variables.
    static TaylorVariable zero(uint as) {
        TaylorVariable r(as); r.set_value(0.0); return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorVariable constant(uint as, const Float& c) {
        TaylorVariable r(as); r.set_value(c); return r; }
    //! \brief Construct the quantity \f$c+x_j\f$ in \a as independent variables.
    static TaylorVariable variable(uint as, const Float& c, uint j) {
        TaylorVariable r(as); r.set_value(c); r.set_gradient(j,1.0); return r; }
    //! \brief Construct the quantity \f$c+\sum g_jx_j\f$.
    static TaylorVariable affine(const Float& c, const Vector<Float>& g) {
        TaylorVariable r(g.size()); r.set_value(c);
        for(uint j=0; j!=g.size(); ++j) { r.set_gradient(j,g[j]); } return r; }
    //! \brief Construct the quantity \f$c+\sum g_jx_j \pm e\f$.
    static TaylorVariable affine(const Float& x, const Vector<Float>& g, const Float& e) {
        TaylorVariable r(g.size()); r.set_value(x); r.set_error(e);
        for(uint j=0; j!=g.size(); ++j) { r.set_gradient(j,g[j]); } return r; }

    static Vector<TaylorVariable> zeroes(uint rs, uint as);
    static Vector<TaylorVariable> constants(uint as, const Vector<Float>& c);
    static Vector<TaylorVariable> variables(const Vector<Float>& x);
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    bool operator==(const TaylorVariable& sd) const {
        return this->_polynomial==sd._polynomial && this->_error == sd._error; }
    //! \brief Inequality operator.
    bool operator!=(const TaylorVariable& sd) const {
        return !(*this==sd); }
    //@}

    //@{
    /*! \name Function operations. */
    //! \brief The domain of the quantity, always given by \f$[-1,1]^{\mathrm{as}}\f$.
    Vector<Interval> domain() const;
    //! \brief An over-approximation to the range of the quantity.
    Interval range() const;
    //! \brief Evaluate the quantity at the point \a x.
    Interval evaluate(const Vector<Float>& x) const;
    //! \brief Evaluate the quantity over the interval of points \a x.
    Interval evaluate(const Vector<Interval>& x) const;
    //@}

    //@{
    /*! \name Simplification operations. */
    //! \brief Truncate to the default maximum degree of the quantity.
    TaylorVariable& truncate();
    //! \brief Truncate to degree \a deg.
    TaylorVariable& truncate(uint deg);
    //! \brief Truncate all terms with any coefficient higher than \a a.
    TaylorVariable& truncate(const MultiIndex& a);
    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    TaylorVariable& sweep();
    //! \brief Remove all terms whose coefficient has magnitude less than \a eps.
    TaylorVariable& sweep(double eps);
    //! \brief Remove all terms whose degree is higher than \a deg or
    //! whose coefficient has magnitude less than \a eps.
    TaylorVariable& clean(uint deg, double eps);
    //! \brief Remove all terms which have high degree or small magnitude.
    TaylorVariable& clean();
    //@}

    //@{
    /*! \name Accuracy parameters. */
    //! \brief .
    static void set_default_maximum_degree(uint md) { _default_maximum_degree=md; }
    //! \brief .
    static void set_default_sweep_threshold(double me) { ARIADNE_ASSERT(me>=0.0); _default_sweep_threshold=me; }
    //! \brief .
    static uint default_maximum_degree() { return _default_maximum_degree; }
    //! \brief .
    static double default_sweep_threshold() { return _default_sweep_threshold; }
    //@}

    //@{
    /*! \name Accuracy parameters. */
    //! \brief .
    void set_maximum_degree(uint md) { this->_maximum_degree=md; }
    //! \brief .
    void set_sweep_threshold(double me) { ARIADNE_ASSERT(me>=0.0); this->_sweep_threshold=me; }
    //! \brief .
    uint maximum_degree() const { return this->_maximum_degree; }
    //! \brief .
    double sweep_threshold() const { return this->_sweep_threshold; }
    //@}

    //@{
    /*! \name Arithmetic operations. */
    //! \brief Inplace addition of another variable.
    friend TaylorVariable& operator+=(TaylorVariable& x, const TaylorVariable& y);
    //! \brief Inplace subtraction of another variable.
    friend TaylorVariable& operator-=(TaylorVariable& x, const TaylorVariable& y);
    //! \brief Inplace addition of a product of two variables.
    friend TaylorVariable& operator+=(TaylorVariable& x, const Product<TaylorVariable,TaylorVariable>& y);
    //! \brief Inplace addition of an exact floating-point constant.
    friend TaylorVariable& operator+=(TaylorVariable& x, const Float& c);
    //! \brief Inplace addition of an interval constant.
    friend TaylorVariable& operator+=(TaylorVariable& x, const Interval& c);
    //! \brief Inplace subtraction of an exact floating-point constant.
    friend TaylorVariable& operator-=(TaylorVariable& x, const Float& c);
    //! \brief Inplace subtraction of an interval constant.
    friend TaylorVariable& operator-=(TaylorVariable& x, const Interval& c);
    //! \brief Inplace multiplication by an exact scalar.
    friend TaylorVariable& operator*=(TaylorVariable& x, const Float& c);
    //! \brief Inplace multiplication by an approximate scalar.
    friend TaylorVariable& operator*=(TaylorVariable& x, const Interval& c);
    //! \brief Inplace division by an exact scalar.
    friend TaylorVariable& operator/=(TaylorVariable& x, const Float& c);
    //! \brief Inplace division by an approximate scalar.
    friend TaylorVariable& operator/=(TaylorVariable& x, const Interval& c);

    //! \brief Unary plus.
    friend TaylorVariable operator+(const TaylorVariable& x);
    //! \brief Unary minus.
    friend TaylorVariable operator-(const TaylorVariable& x);
    //! \brief Addition.
    friend TaylorVariable operator+(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Subtraction.
    friend TaylorVariable operator-(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Multiplication.
    friend TaylorVariable operator*(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Division.
    friend TaylorVariable operator/(const TaylorVariable& x, const TaylorVariable& y);

    //! \brief Addition of a scakar.
    friend TaylorVariable operator+(const TaylorVariable& x, const Float& c);
    //! \brief Subtraction of a scakar.
    friend TaylorVariable operator-(const TaylorVariable& x, const Float& c);
    //! \brief Multiplication by a scakar.
    friend TaylorVariable operator*(const TaylorVariable& x, const Float& c);
    //! \brief Division by a scakar.
    friend TaylorVariable operator/(const TaylorVariable& x, const Float& c);
    //! \brief Addition of a scakar.
    friend TaylorVariable operator+(const TaylorVariable& x, const Interval& c);
    //! \brief Subtraction of a scakar.
    friend TaylorVariable operator-(const TaylorVariable& x, const Interval& c);
    //! \brief Multiplication by a scakar.
    friend TaylorVariable operator*(const TaylorVariable& x, const Interval& c);
    //! \brief Division by a scakar.
    friend TaylorVariable operator/(const TaylorVariable& x, const Interval& c);
    //! \brief Addition of a scakar.
    friend TaylorVariable operator+(const Float& c, const TaylorVariable& x);
    //! \brief Subtraction from a scakar.
    friend TaylorVariable operator-(const Float& c, const TaylorVariable& x);
    //! \brief Multiplication by a scakar.
    friend TaylorVariable operator*(const Float& c, const TaylorVariable& x);
    //! \brief Division through a scalar.
    friend TaylorVariable operator/(const Float& c, const TaylorVariable& x);
    //! \brief Addition of a scakar.
    friend TaylorVariable operator+(const Interval& c, const TaylorVariable& x);
    //! \brief Subtraction from a scakar.
    friend TaylorVariable operator-(const Interval& c, const TaylorVariable& x);
    //! \brief Multiplication by a scakar.
    friend TaylorVariable operator*(const Interval& c, const TaylorVariable& x);
    //! \brief Division through a scalar.
    friend TaylorVariable operator/(const Interval& c, const TaylorVariable& x);

    //! \brief Maximum. Throws an error if one variable is not greater than the other
    //! over the entire domain.
    friend TaylorVariable max(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Minimum. Throws an error if one variable is not greater than the other
    //! over the entire domain.
    friend TaylorVariable min(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Addition.
    friend TaylorVariable add(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Multiplication.
    friend TaylorVariable mul(const TaylorVariable& x, const TaylorVariable& y);
    //! \brief Absolute value. Throws an error if one variable is neither greater than
    //! zero over the entire domain, nor less than zero over the entire domain.
    friend TaylorVariable abs(const TaylorVariable& x);
    //! \brief Negation.
    friend TaylorVariable neg(const TaylorVariable& x);
    //! \brief Reciprocal.
    friend TaylorVariable rec(const TaylorVariable& x);
    //! \brief Square.
    friend TaylorVariable sqr(const TaylorVariable& x);
    //! \brief Power.
    friend TaylorVariable pow(const TaylorVariable& x, int n);
    //! \brief Square root.
    friend TaylorVariable sqrt(const TaylorVariable& x);
    //! \brief Natural exponent.
    friend TaylorVariable exp(const TaylorVariable& x);
    //! \brief Natural logarithm.
    friend TaylorVariable log(const TaylorVariable& x);
    //! \brief Sine.
    friend TaylorVariable sin(const TaylorVariable& x);
    //! \brief Cosine.
    friend TaylorVariable cos(const TaylorVariable& x);
    //! \brief Tangent.
    friend TaylorVariable tan(const TaylorVariable& x);
    //@}

    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const TaylorVariable& x);
    //@}

  public:
    std::string str() const;

    TaylorVariable& clobber();
    TaylorVariable& clobber(uint o);
    TaylorVariable& clobber(uint so, uint to);
};

TaylorVariable max(const TaylorVariable& x, const TaylorVariable& y);
TaylorVariable min(const TaylorVariable& x, const TaylorVariable& y);
TaylorVariable abs(const TaylorVariable& x);
TaylorVariable neg(const TaylorVariable& x);
TaylorVariable rec(const TaylorVariable& x);
TaylorVariable exp(const TaylorVariable& x);
TaylorVariable log(const TaylorVariable& x);
TaylorVariable sin(const TaylorVariable& x);
TaylorVariable cos(const TaylorVariable& x);
TaylorVariable tan(const TaylorVariable& x);


template<>
class Vector<TaylorVariable>
    : public ublas::vector<TaylorVariable>
{
  public:
    Vector() : ublas::vector<TaylorVariable>() { }
    Vector(uint rs) : ublas::vector<TaylorVariable>(rs) { for(uint i=0; i!=rs; ++i) { (*this)[i]=TaylorVariable(); } }
    Vector(uint rs, uint as) : ublas::vector<TaylorVariable>(rs) { for(uint i=0; i!=rs; ++i) { (*this)[i]=TaylorVariable(as); } }
    Vector(uint rs, const TaylorVariable& x) : ublas::vector<TaylorVariable>(rs) { for(uint i=0; i!=rs; ++i) { (*this)[i]=x; } }
    Vector(uint rs, uint as, uint deg, double d0, ...);
    template<class E> Vector(const ublas::vector_expression<E> &ve) : ublas::vector<TaylorVariable>(ve) { }
    uint result_size() const { return this->size(); }
    uint argument_size() const { ARIADNE_ASSERT(this->size()>0); return (*this)[0].argument_size(); }

    Vector<Float> value() const;
    Vector<Interval> evaluate(const Vector<Float>&) const;
    Vector<Interval> evaluate(const Vector<Interval>&) const;

    void check() const {
        for(uint i=0; i!=this->size(); ++i) {
            ARIADNE_ASSERT((*this)[0].argument_size()==(*this)[i].argument_size()); } }
};


inline Vector<TaylorVariable> TaylorVariable::zeroes(uint rs, uint as) {
    Vector<TaylorVariable> result(rs); for(uint i=0; i!=rs; ++i) {
        result[i]=TaylorVariable::zero(as); } return result; }
inline Vector<TaylorVariable> TaylorVariable::constants(uint as, const Vector<Float>& c) {
    Vector<TaylorVariable> result(c.size()); for(uint i=0; i!=c.size(); ++i) {
        result[i]=TaylorVariable::constant(as,c[i]); } return result; }
inline Vector<TaylorVariable> TaylorVariable::variables(const Vector<Float>& x) {
    Vector<TaylorVariable> result(x.size(),x.size()); for(uint i=0; i!=x.size(); ++i) {
        result[i]=TaylorVariable::variable(x.size(),x[i],i); } return result; }


} // namespace Ariadne

#endif // ARIADNE_TAYLOR_VARIABLE_H
