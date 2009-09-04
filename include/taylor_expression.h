/***************************************************************************
 *            taylor_expression.h
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

/*! \file taylor_expression.h
 *  \brief Approximate scalar functions on a bounded domain with a sparse representation.
 */

#ifndef ARIADNE_TAYLOR_EXPRESSION_H
#define ARIADNE_TAYLOR_EXPRESSION_H

#include <map>

#include "macros.h"
#include "array.h"
#include "vector.h"
#include "multi_index.h"
#include "taylor_model.h"

namespace Ariadne {

template<class T1, class T2> class Product;
template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Expansion;
template<class X> class Polynomial;
class TaylorModel;
class TaylorExpression;
class TaylorFunction;
class ScalarFunctionInterface;

// Remove the error term
TaylorExpression midpoint(const TaylorExpression& x);

// Restrict to a smaller domain. REQUIRED
TaylorExpression restrict(const TaylorExpression& x, const Vector<Interval>& d);
// Extend to a larger domain. REQUIRED
TaylorExpression extend(const TaylorExpression& x, const Vector<Interval>& d);

// Test if a variable refines another
bool refines(const TaylorExpression& tv1, const TaylorExpression& tv2);
// Test if two variables definitely represent different quantities
bool disjoint(const TaylorExpression& x1, const TaylorExpression& x2);
// Test if two variables definitely represent different quantities
TaylorExpression intersection(const TaylorExpression& x1, const TaylorExpression& x2);

// Evaluate an array of Taylor variables on a vector.
Interval evaluate(const TaylorExpression& x, const Vector<Interval>& sy);

// Set the value of the \a kth variable to c
TaylorExpression partial_evaluate(const TaylorExpression& x, uint k, const Interval& c);
// Restrict the \a kth variable to lie in the interval \a d.
TaylorExpression restrict(const TaylorExpression& x, uint k, const Interval& d);

// Compose with an expression.
TaylorExpression compose(const ScalarFunctionInterface& x, const Vector<TaylorExpression>& y);

// Split the variable over two domains, subdividing along the independent variable j.
pair<TaylorExpression,TaylorExpression> split(const TaylorExpression& x, uint j);


// Embed the variable in a space of higher dimension
TaylorExpression embed(const TaylorExpression& tv1, const Interval& d2);
TaylorExpression embed(const Vector<Interval>& d1, const TaylorExpression& tv2);

// Antidifferentiation operator
TaylorExpression antiderivative(const TaylorExpression& x, uint k);
TaylorExpression derivative(const TaylorExpression& x, uint k);

// Implicit function solver; solves f(x,h(x))=0 on dom1(f)
TaylorExpression implicit(const TaylorExpression& f);
// Implicit function solver solves f(g(x),h(x))=0 on dom(g)
TaylorExpression implicit(const ScalarFunctionInterface& f, const TaylorFunction& g);
// Implicit function solver solves f(x,h(x))=0 on d
TaylorExpression implicit(const ScalarFunctionInterface& f, const Vector<Interval>& d);

TaylorExpression crossing_time(const ScalarFunctionInterface& g, const FunctionInterface& f, const Vector<Interval>& d);


/*! \brief A TaylorExpression is a type of FunctionModel in which a the restriction of a scalar function \f$f:\R^n\rightarrow\R\f$ on a domain \f$D\f$ is approximated by polynomial \f$p\f$ with uniform error \f$e\f$.
 *
 * Formally, a TaylorExpression is a triple \f$(D,p,e)\f$ representing a set of continuous functions \f$\mathrm{T}(D,p,e)\f$ by
 * \f[ \mathrm{T}(D,p,e) = \{ f:\R^n\rightarrow \R \med \sup_{x\in D}|f(x)-p(x)| \leq e \} . \f]
 * Note that there is no need for the functions \f$f\f$ to be themselves polynomial, and that no information is given
 * about the values of \f$f\f$ outside of \f$D\f$. Information about the derivatives of \f$f\f$ is also unavailable.
 * However, integrals of \f$f\f$ can be computed.
 *
 * Internally, the polynomial \f$p\f$ is represented as the composition \f$p=m\circ s^{-1}\f$,
 * where \f$m:[-1,+1]^n\rightarrow\R\f$ and \f$s:[-1,+1]^n\rightarrow D\f$ is a scaling function,
 * \f$s_i(y_i)=(a_i+b_i)/2+(b_i-a_i)y_i/2\f$ where \f$D_i=[a_i,b_i]\f$ is the \f$i^\textrm{th}\f$ subinterval of \f$D\f$.
 *
 * When solving algebraic equations by iterative Newton-like methods, it is necessary to compute the derivatives of \f$f\f$.
 * For these applications, it suffices to compute the derivative of \f$p\f$, since only a uniform approximation to the solution is required.
 *
 * Finding exact bounds for the range of \f$p\f$ over \f$D\f$ is an NP-complete problem,
 * for but there are a number of techniques available.
 * 
 * \sa Expansion, TaylorModel, TaylorFunction, TaylorSet.
 */
class TaylorExpression
{
    typedef Vector<Interval> DomainType;
    typedef TaylorModel ModelType;
    typedef Expansion<Float> ExpansionType;
    typedef Float ErrorType;
    static const Float _zero;
    DomainType _domain;
    ModelType _model;
  public:
    //! \brief An iterator through the (index,coefficient) pairs of the expansion expansion.
    typedef ExpansionType::iterator iterator;
    //! \brief A constant iterator through the (index,coefficient) pairs of the expansion expansion.
    typedef ExpansionType::const_iterator const_iterator;

    //@{
    /*! \name Constructors and destructors. */
    //! \brief Default constructor.
    explicit TaylorExpression();
    //! \brief Construct a TaylorExpression over the domain \a d.
    explicit TaylorExpression(const DomainType& d);
    //! \brief Construct a TaylorExpression over the domain \a d, based on the scaled model \a m.
    explicit TaylorExpression(const DomainType& d, const TaylorModel& m);
    //! \brief Construct a TaylorExpression over the domain \a d, with scaled power series expansion \a f and error \a e.
    explicit TaylorExpression(const DomainType& d, const ExpansionType& f, const ErrorType& e=0);

    //! \brief Construct a TaylorExpression over the domain \a d from the expression \a f.
    explicit TaylorExpression(const DomainType& d, const ScalarFunctionInterface& f);
    //! \brief Construct a TaylorExpression over the domain \a d from the polynomial \a p.
    explicit TaylorExpression(const DomainType& d, const Polynomial<Float>& p);
    //! \brief Construct a TaylorExpression over the domain \a d from the interval polynomial \a p.
    explicit TaylorExpression(const DomainType& d, const Polynomial<Interval>& p);
    //@}

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to a constant, keeping the same number of arguments.
    TaylorExpression& operator=(const Float& c) { this->_model=c; return *this; }
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorExpression& operator=(const Interval& c) { this->_model=c; return *this; }
    //@}

    //@{
    /*! \name Named constructors. */
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorExpression constant(const DomainType& d, const Float& c);
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorExpression constant(const DomainType& d, const Interval& c);
    //! \brief Construct the quantity \f$x_j\f$ over the domain \a d.
    static TaylorExpression variable(const DomainType& d, unsigned int j);
    //! \brief Construct the quantity \f$c+\sum g_jx_j\f$ over the domain \a d.
    static TaylorExpression affine(const DomainType& d, const Float& c, const Vector<Float>& g);
    //! \brief Construct the quantity \f$c+\sum g_jx_j \pm e\f$ over domain \a d.
    static TaylorExpression affine(const DomainType& d, const Float& x, const Vector<Float>& g, const Float& e) ;

    //! \brief Return the vector of constants with values \a c over domain \a d.
    static Vector<TaylorExpression> constants(const DomainType& d, const Vector<Float>& c);
    //! \brief Return the vector of constants with interval values \a c over domain \a d.
    static Vector<TaylorExpression> constants(const DomainType& d, const Vector<Interval>& c);
    //! \brief Return the vector of variables with values \a x over domain \a d.
    static Vector<TaylorExpression> variables(const DomainType& d);
    //@}

    //@{
    /*! \name Data access */
    //! \brief The domain of the quantity.
    const DomainType& domain() const { return this->_domain; }
    //! \brief An over-approximation to the range of the quantity; not necessarily tight.
    const Interval codomain() const { return this->_model.range(); }
    //! \brief The scaled expansion over a unit box with error bound.
    const ModelType& model() const { return this->_model; }
    //! \brief The scaled expansion over a unit box.
    const ExpansionType& expansion() const { return this->_model._expansion; }
    //! \brief The error of the expansion over the domain.
    const ErrorType& error() const { return this->_model._error; }
    //! \brief A reference to the scaled expansion over a unit box with error bound.
    ModelType& model() { return this->_model; }
    //! \brief A reference to the expansion.
    ExpansionType& expansion() { return this->_model._expansion; }
    //! \brief A reference to the error of the expansion over the domain.
    ErrorType& error() { return this->_model._error; }
    //! \brief The centre of the expansion (the value of the constant term).
    Float centre() { return this->_model.value(); }
    //! \brief A reference to the constant term in the expansion.
    Float& value() { return this->_model.value(); }
    //! \brief The constant term in the expansion.
    const Float& value() const { return this->_model.value(); }
    //! \brief A polynomial representation.
    Polynomial<Interval> polynomial() const;

    //! \brief Set the error of the expansion.
    void set_error(const Float& ne) { this->_model.set_error(ne); }
    //! \brief Set the constant term in the expansion.
    void set_value(const Float& c) { this->_model.set_value(c); }

    //! \brief The coefficient of the term in $x^a$.
    const Float& operator[](const MultiIndex& a) const { return this->_model[a]; }
    //! \brief A read/write reference to the coefficient of the term in $x^a$.
    Float& operator[](const MultiIndex& a) { return this->_model[a]; }

    //! \brief An iterator to the first term in the expansion expansion.
    iterator begin() { return this->_model.begin(); }
    //! \brief A constant iterator to the first term in the expansion expansion.
    const_iterator begin() const { return this->_model.begin(); }
    //! \brief An iterator to the end of the expansion expansion.
    iterator end() { return this->_model.end(); }
    //! \brief A constant iterator to the end of the expansion expansion.
    const_iterator end() const { return this->_model.end(); }
    //! \brief An iterator to the term with index \a a.
    iterator find(const MultiIndex& a) { return this->_model.find(a); }
    //! \brief A constant iterator to the term with index \a a.
    const_iterator find(const MultiIndex& a) const { return this->_model.find(a); }

    //! \brief The number of variables in the argument of the quantity.
    uint argument_size() const { return this->_model.argument_size(); }
    //! \brief The maximum degree of terms in the expansion expansion.
    uint degree() const { return this->_model.degree(); }
    //! \brief The number of nonzero terms in the expansion expansion.
    uint number_of_nonzeros() const { return this->_model.number_of_nonzeros(); }
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    bool operator==(const TaylorExpression& tv) const;
    //! \brief Inequality operator.
    bool operator!=(const TaylorExpression& tv) const { return !(*this==tv); }
    //@}

    //@{
    /*! \name Function operations. */
    //! \brief An over-approximation to the range of the quantity.
    Interval range() const { return this->_model.range(); }
    //! \brief Evaluate the quantity at the point \a x.
    Interval evaluate(const Vector<Float>& x) const;
    //! \brief Evaluate the quantity over the interval of points \a x.
    Interval evaluate(const Vector<Interval>& x) const;
    //@}

    //@{
    /*! \name Simplification operations. */
    //! \brief Truncate to the default maximum degree of the quantity.
    TaylorExpression& truncate() { this->_model.truncate(); return *this; }
    //! \brief Truncate to degree \a deg.
    TaylorExpression& truncate(uint deg) { this->_model.truncate(deg); return *this; }
    //! \brief Truncate all terms with any coefficient higher than \a a.
    TaylorExpression& truncate(const MultiIndex& a) { this->_model.truncate(a); return *this; }
    //! \brief Truncate all terms with any coefficient higher than those given by \a b.
    TaylorExpression& truncate(const MultiIndexBound& b) { this->_model.truncate(b); return *this; }
    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    TaylorExpression& sweep() { this->_model.sweep(); return *this; }
    //! \brief Remove all terms whose coefficient has magnitude less than \a eps.
    TaylorExpression& sweep(double eps) { this->_model.sweep(eps); return *this; }
    //! \brief Remove all terms whose degree is higher than \a deg or
    //! whose coefficient has magnitude less than \a eps.
    TaylorExpression& clean(const TaylorModel::Accuracy& acc) { this->_model.clean(acc); return *this; }
    //! \brief Remove all terms which have high degree or small magnitude.
    TaylorExpression& clean() { this->_model.clean(); return *this; }
    //@}

    //@{
    /*! \name Accuracy parameters. */
    //! \brief Specify a bound on the terms which may be present in the expansion.
    void set_maximum_index(MultiIndexBound md) { this->_model.set_maximum_index(md); }
    //! \brief Specify the maximum degree \a md for terms which may be present in the expansion.
    void set_maximum_degree(uint md) { this->_model.set_maximum_degree(md); }
    //! \brief Specify the minimum absolute value \a me for coefficients of terms which may be present in the expansion.
    void set_sweep_threshold(double me) { ARIADNE_ASSERT(me>=0.0); this->_model.set_sweep_threshold(me); }
    //! \copydoc TaylorModel::maximum_index()
    MultiIndexBound maximum_index() const { return this->_model.maximum_index(); }
    //! \copydoc TaylorModel::maximum_degree()
    uint maximum_degree() const { return this->_model.maximum_degree(); }
    //! \copydoc TaylorModel::sweep_threshold()
    double sweep_threshold() const { return this->_model.sweep_threshold(); }
    //@}

    //@{
    /*! \name Non-arithmetic operations. */
    //! \brief Test if the quantity is a better approximation than \a t throughout the domain.
    friend bool refines(const TaylorExpression& x1, const TaylorExpression& x2);
    //! \brief Test if the quantities are disjoint.
    friend bool disjoint(const TaylorExpression& x1, const TaylorExpression& x2);
    //! \brief Test if the quantities are disjoint.
    friend TaylorExpression intersection(const TaylorExpression& x1, const TaylorExpression& x2);
    //! \brief Restrict to a subdomain.
    friend TaylorExpression restrict(const TaylorExpression& x, const DomainType& d);
    //! \brief Restrict the values of the \a k<sup>th</sup> variable to the subinterval \a d.
    friend TaylorExpression restrict(const TaylorExpression& x, uint k, const Interval& d);
    //! \brief Extend over a larger domain. Only possible if the larger domain is only larger where the smaller domain is a singleton.
    //! The extension is performed keeping \a x constant over the new coordinates.
    friend TaylorExpression extend(const TaylorExpression& x, const DomainType& d);
    //@}

    /*! \name Arithmetic operations. */
    //! \brief Inplace addition of another variable.
    friend TaylorExpression& operator+=(TaylorExpression& x, const TaylorExpression& y);
    //! \brief Inplace subtraction of another variable.
    friend TaylorExpression& operator-=(TaylorExpression& x, const TaylorExpression& y);
    //! \brief Inplace addition of a product of two variables.
    friend TaylorExpression& operator+=(TaylorExpression& x, const Product<TaylorExpression,TaylorExpression>& y);
    //! \brief Inplace addition of an exact floating-point constant.
    friend TaylorExpression& operator+=(TaylorExpression& x, const Float& c);
    //! \brief Inplace addition of an interval constant.
    friend TaylorExpression& operator+=(TaylorExpression& x, const Interval& c);
    //! \brief Inplace subtraction of an exact floating-point constant.
    friend TaylorExpression& operator-=(TaylorExpression& x, const Float& c);
    //! \brief Inplace subtraction of an interval constant.
    friend TaylorExpression& operator-=(TaylorExpression& x, const Interval& c);
    //! \brief Inplace multiplication by an exact scalar.
    friend TaylorExpression& operator*=(TaylorExpression& x, const Float& c);
    //! \brief Inplace multiplication by an approximate scalar.
    friend TaylorExpression& operator*=(TaylorExpression& x, const Interval& c);
    //! \brief Inplace division by an exact scalar.
    friend TaylorExpression& operator/=(TaylorExpression& x, const Float& c);
    //! \brief Inplace division by an approximate scalar.
    friend TaylorExpression& operator/=(TaylorExpression& x, const Interval& c);

    //! \brief Unary plus.
    friend TaylorExpression operator+(const TaylorExpression& x);
    //! \brief Unary minus.
    friend TaylorExpression operator-(const TaylorExpression& x);
    //! \brief Addition.
    friend TaylorExpression operator+(const TaylorExpression& x, const TaylorExpression& y);
    //! \brief Subtraction.
    friend TaylorExpression operator-(const TaylorExpression& x, const TaylorExpression& y);
    //! \brief Multiplication.
    friend TaylorExpression operator*(const TaylorExpression& x, const TaylorExpression& y);
    //! \brief Division.
    friend TaylorExpression operator/(const TaylorExpression& x, const TaylorExpression& y);

    //! \brief Addition of a scakar.
    friend TaylorExpression operator+(const TaylorExpression& x, const Float& c);
    //! \brief Subtraction of a scakar.
    friend TaylorExpression operator-(const TaylorExpression& x, const Float& c);
    //! \brief Multiplication by a scakar.
    friend TaylorExpression operator*(const TaylorExpression& x, const Float& c);
    //! \brief Division by a scakar.
    friend TaylorExpression operator/(const TaylorExpression& x, const Float& c);
    //! \brief Addition of a scakar.
    friend TaylorExpression operator+(const TaylorExpression& x, const Interval& c);
    //! \brief Subtraction of a scakar.
    friend TaylorExpression operator-(const TaylorExpression& x, const Interval& c);
    //! \brief Multiplication by a scakar.
    friend TaylorExpression operator*(const TaylorExpression& x, const Interval& c);
    //! \brief Division by a scakar.
    friend TaylorExpression operator/(const TaylorExpression& x, const Interval& c);
    //! \brief Addition of a scakar.
    friend TaylorExpression operator+(const Float& c, const TaylorExpression& x);
    //! \brief Subtraction from a scakar.
    friend TaylorExpression operator-(const Float& c, const TaylorExpression& x);
    //! \brief Multiplication by a scakar.
    friend TaylorExpression operator*(const Float& c, const TaylorExpression& x);
    //! \brief Division through a scalar.
    friend TaylorExpression operator/(const Float& c, const TaylorExpression& x);
    //! \brief Addition of a scakar.
    friend TaylorExpression operator+(const Interval& c, const TaylorExpression& x);
    //! \brief Subtraction from a scakar.
    friend TaylorExpression operator-(const Interval& c, const TaylorExpression& x);
    //! \brief Multiplication by a scakar.
    friend TaylorExpression operator*(const Interval& c, const TaylorExpression& x);
    //! \brief Division through a scalar.
    friend TaylorExpression operator/(const Interval& c, const TaylorExpression& x);

    //! \brief Maximum. Throws an error if one variable is not greater than the other
    //! over the entire domain.
    friend TaylorExpression max(const TaylorExpression& x, const TaylorExpression& y);
    //! \brief Minimum. Throws an error if one variable is not greater than the other
    //! over the entire domain.
    friend TaylorExpression min(const TaylorExpression& x, const TaylorExpression& y);
    //! \brief Addition.
    friend TaylorExpression add(const TaylorExpression& x, const TaylorExpression& y);
    //! \brief Multiplication.
    friend TaylorExpression mul(const TaylorExpression& x, const TaylorExpression& y);
    //! \brief Absolute value. Throws an error if one variable is neither greater than
    //! zero over the entire domain, nor less than zero over the entire domain.
    friend TaylorExpression abs(const TaylorExpression& x);
    //! \brief Negation.
    friend TaylorExpression neg(const TaylorExpression& x);
    //! \brief Reciprocal.
    friend TaylorExpression rec(const TaylorExpression& x);
    //! \brief Square.
    friend TaylorExpression sqr(const TaylorExpression& x);
    //! \brief Power.
    friend TaylorExpression pow(const TaylorExpression& x, int n);
    //! \brief Square root.
    friend TaylorExpression sqrt(const TaylorExpression& x);
    //! \brief Natural exponent.
    friend TaylorExpression exp(const TaylorExpression& x);
    //! \brief Natural logarithm.
    friend TaylorExpression log(const TaylorExpression& x);
    //! \brief Sine.
    friend TaylorExpression sin(const TaylorExpression& x);
    //! \brief Cosine.
    friend TaylorExpression cos(const TaylorExpression& x);
    //! \brief Tangent.
    friend TaylorExpression tan(const TaylorExpression& x);
    //! \brief Inverse sine.
    friend TaylorExpression asin(const TaylorExpression& x);
    //! \brief Inverse cosine.
    friend TaylorExpression acos(const TaylorExpression& x);
    //! \brief Inverse tangent.
    friend TaylorExpression atan(const TaylorExpression& x);
    //@}

    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const TaylorExpression& x);
    //@}

  public:
    TaylorExpression& clobber() { this->_model.clobber(); return *this; }
    TaylorExpression& clobber(uint o) { this->_model.clobber(o); return *this; }
    TaylorExpression& clobber(uint so, uint to) { this->_model.clobber(so,to); return *this; }
};


inline tribool operator>(const TaylorExpression& x, const Float& c) {
    Interval r=x.range(); if(r.lower()>c) { return true; } else if(r.upper()<=c) { return false; } else { return indeterminate; } }
inline tribool operator<(const TaylorExpression& x, const Float& c) {
    Interval r=x.range(); if(r.lower()<c) { return true; } else if(r.upper()>=c) { return false; } else { return indeterminate; } }

inline tribool operator>(const TaylorExpression& x, const TaylorExpression& y) { return (x-y)>0; }
inline tribool operator<(const TaylorExpression& x, const TaylorExpression& y) { return (x-y)<0; }

inline TaylorExpression& operator+=(TaylorExpression& x, const Float& c) {
    x._model+=c; return x; }
inline TaylorExpression& operator-=(TaylorExpression& x, const Float& c) {
    x._model-=c; return x; }
inline TaylorExpression& operator*=(TaylorExpression& x, const Float& c) {
    x._model*=c; return x; }
inline TaylorExpression& operator/=(TaylorExpression& x, const Float& c) {
    x._model*=c; return x; }

inline TaylorExpression& operator+=(TaylorExpression& x, const Interval& c) {
    x._model+=c; return x; }
inline TaylorExpression& operator-=(TaylorExpression& x, const Interval& c) {
    x._model-=c; return x; }
inline TaylorExpression& operator*=(TaylorExpression& x, const Interval& c) {
    x._model*=c; return x; }
inline TaylorExpression& operator/=(TaylorExpression& x, const Interval& c) {
    x._model*=c; return x; }


inline TaylorExpression operator+(const TaylorExpression& x) {
    return TaylorExpression(x._domain,x._model); }
inline TaylorExpression operator-(const TaylorExpression& x) {
    return TaylorExpression(x._domain,-x._model); }

inline TaylorExpression operator+(const TaylorExpression& x, const Float& c) {
    return TaylorExpression(x._domain,x._model+c); }
inline TaylorExpression operator-(const TaylorExpression& x, const Float& c) {
    return TaylorExpression(x._domain,x._model-c); }
inline TaylorExpression operator*(const TaylorExpression& x, const Float& c) {
    return TaylorExpression(x._domain,x._model*c); }
inline TaylorExpression operator/(const TaylorExpression& x, const Float& c) {
    return TaylorExpression(x._domain,x._model/c); }
inline TaylorExpression operator+(const TaylorExpression& x, const Interval& c) {
    return TaylorExpression(x._domain,x._model+c); }
inline TaylorExpression operator-(const TaylorExpression& x, const Interval& c) {
    return TaylorExpression(x._domain,x._model-c); }
inline TaylorExpression operator*(const TaylorExpression& x, const Interval& c) {
    return TaylorExpression(x._domain,x._model*c); }
inline TaylorExpression operator/(const TaylorExpression& x, const Interval& c) {
    return TaylorExpression(x._domain,x._model/c); }
inline TaylorExpression operator+(const Float& c, const TaylorExpression& x) {
    return TaylorExpression(x._domain,c+x._model); }
inline TaylorExpression operator-(const Float& c, const TaylorExpression& x) {
    return TaylorExpression(x._domain,c-x._model); }
inline TaylorExpression operator*(const Float& c, const TaylorExpression& x) {
    return TaylorExpression(x._domain,c*x._model); }
inline TaylorExpression operator/(const Float& c, const TaylorExpression& x) {
    return TaylorExpression(x._domain,c/x._model); }
inline TaylorExpression operator+(const Interval& c, const TaylorExpression& x) {
    return TaylorExpression(x._domain,c+x._model); }
inline TaylorExpression operator-(const Interval& c, const TaylorExpression& x) {
    return TaylorExpression(x._domain,c-x._model); }
inline TaylorExpression operator*(const Interval& c, const TaylorExpression& x) {
    return TaylorExpression(x._domain,c*x._model); }
inline TaylorExpression operator/(const Interval& c, const TaylorExpression& x) {
    return TaylorExpression(x._domain,c/x._model); }


inline TaylorExpression abs(const TaylorExpression& x) {
    return TaylorExpression(x._domain,abs(x._model)); }
inline TaylorExpression neg(const TaylorExpression& x) {
    return TaylorExpression(x._domain,abs(x._model)); }
inline TaylorExpression rec(const TaylorExpression& x) {
    return TaylorExpression(x._domain,rec(x._model)); }
inline TaylorExpression sqr(const TaylorExpression& x) {
    return TaylorExpression(x._domain,sqr(x._model)); }
inline TaylorExpression pow(const TaylorExpression& x, int n) {
    return TaylorExpression(x._domain,pow(x._model,n)); }
inline TaylorExpression sqrt(const TaylorExpression& x) {
    return TaylorExpression(x._domain,sqrt(x._model)); }
inline TaylorExpression exp(const TaylorExpression& x) {
    return TaylorExpression(x._domain,exp(x._model)); }
inline TaylorExpression log(const TaylorExpression& x) {
    return TaylorExpression(x._domain,log(x._model)); }
inline TaylorExpression sin(const TaylorExpression& x) {
    return TaylorExpression(x._domain,sin(x._model)); }
inline TaylorExpression cos(const TaylorExpression& x) {
    return TaylorExpression(x._domain,cos(x._model)); }
inline TaylorExpression tan(const TaylorExpression& x) {
    return TaylorExpression(x._domain,tan(x._model)); }
inline TaylorExpression asin(const TaylorExpression& x) {
    return TaylorExpression(x._domain,asin(x._model)); }
inline TaylorExpression acos(const TaylorExpression& x) {
    return TaylorExpression(x._domain,acos(x._model)); }
inline TaylorExpression atan(const TaylorExpression& x) {
    return TaylorExpression(x._domain,atan(x._model)); }



inline Interval evaluate(const TaylorExpression& tv, const Vector<Interval>& x) {
    return tv.evaluate(x); }

inline TaylorExpression antiderivative(const TaylorExpression& x, uint k) {
    Interval sf=rad_ivl(x.domain()[k]);
    return TaylorExpression(x.domain(),antiderivative(x.model(),k)*sf); }

inline TaylorExpression derivative(const TaylorExpression& x, uint k) {
    Interval sf=1/rad_ivl(x.domain()[k]);
    return TaylorExpression(x.domain(),derivative(x.model(),k)*sf); }

inline TaylorExpression embed(const TaylorExpression& tv1, const Interval& dom2) {
    return TaylorExpression(join(tv1.domain(),dom2),embed(tv1.model(),1u)); }
inline TaylorExpression embed(const TaylorExpression& tv1, const Vector<Interval>& dom2) {
    return TaylorExpression(join(tv1.domain(),dom2),embed(tv1.model(),dom2.size())); }
inline TaylorExpression embed(const Vector<Interval>& dom1, const TaylorExpression& tv2) {
    return TaylorExpression(join(dom1,tv2.domain()),embed(dom1.size(),tv2.model())); }





} // namespace Ariadne

#endif // ARIADNE_TAYLOR_EXPRESSION_H
