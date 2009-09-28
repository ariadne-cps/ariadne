/***************************************************************************
 *            taylor_function.h
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

/*! \file taylor_function.h
 *  \brief Over-approximations of functions based on Taylor expansions.
 */

#ifndef ARIADNE_TAYLOR_FUNCTION_H
#define ARIADNE_TAYLOR_FUNCTION_H

#include <iosfwd>
#include "numeric.h"
#include "vector.h"
#include "taylor_model.h"

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Polynomial;

class ScalarFunction;
class VectorFunction;
class MultiIndex;
class TaylorModel;
class ScalarTaylorFunction;
class VectorTaylorFunction;

// Remove the error term
ScalarTaylorFunction midpoint(const ScalarTaylorFunction& x);

// Restrict to a smaller domain. REQUIRED
ScalarTaylorFunction restrict(const ScalarTaylorFunction& x, const Vector<Interval>& d);
// Extend to a larger domain. REQUIRED
ScalarTaylorFunction extend(const ScalarTaylorFunction& x, const Vector<Interval>& d);

// Test if a variable refines another
bool refines(const ScalarTaylorFunction& tv1, const ScalarTaylorFunction& tv2);
// Test if two variables definitely represent different quantities
bool disjoint(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2);
// Test if two variables definitely represent different quantities
ScalarTaylorFunction intersection(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2);

// Evaluate an array of Taylor variables on a vector.
Interval evaluate(const ScalarTaylorFunction& x, const Vector<Interval>& sy);

// Set the value of the \a kth variable to c
ScalarTaylorFunction partial_evaluate(const ScalarTaylorFunction& x, uint k, const Interval& c);
// Restrict the \a kth variable to lie in the interval \a d.
ScalarTaylorFunction restrict(const ScalarTaylorFunction& x, uint k, const Interval& d);

// Compose with an expression.
ScalarTaylorFunction compose(const ScalarFunction& x, const VectorTaylorFunction& y);

// Split the variable over two domains, subdividing along the independent variable j.
pair<ScalarTaylorFunction,ScalarTaylorFunction> split(const ScalarTaylorFunction& x, uint j);


// Embed the variable in a space of higher dimension
ScalarTaylorFunction embed(const ScalarTaylorFunction& tv1, const Interval& d2);
ScalarTaylorFunction embed(const Vector<Interval>& d1, const ScalarTaylorFunction& tv2);

// Antidifferentiation operator
ScalarTaylorFunction antiderivative(const ScalarTaylorFunction& x, uint k);
ScalarTaylorFunction derivative(const ScalarTaylorFunction& x, uint k);

// Implicit function solver; solves f(x,h(x))=0 on dom1(f)
ScalarTaylorFunction implicit(const ScalarTaylorFunction& f);
// Implicit function solver solves f(g(x),h(x))=0 on dom(g)
ScalarTaylorFunction implicit(const ScalarFunction& f, const VectorTaylorFunction& g);
// Implicit function solver solves f(x,h(x))=0 on d
ScalarTaylorFunction implicit(const ScalarFunction& f, const Vector<Interval>& d);

ScalarTaylorFunction crossing_time(const ScalarFunction& g, const VectorFunction& f, const Vector<Interval>& d);

class VectorTaylorFunctionElementReference;

/*! \brief A ScalarTaylorFunction is a type of FunctionModel in which a the restriction of a scalar function \f$f:\R^n\rightarrow\R\f$ on a domain \f$D\f$ is approximated by polynomial \f$p\f$ with uniform error \f$e\f$.
 *
 * Formally, a ScalarTaylorFunction is a triple \f$(D,p,e)\f$ representing a set of continuous functions \f$\mathrm{T}(D,p,e)\f$ by
 * \f[ \mathrm{T}(D,p,e) = \{ f:\R^n\rightarrow \R \mid \sup_{x\in D}|f(x)-p(x)| \leq e \} . \f]
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
 * \sa Expansion, TaylorModel, VectorTaylorFunction, TaylorSet.
 */
class ScalarTaylorFunction
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
    explicit ScalarTaylorFunction();
    //! \brief Construct a ScalarTaylorFunction over the domain \a d.
    explicit ScalarTaylorFunction(const DomainType& d);
    //! \brief Construct a ScalarTaylorFunction over the domain \a d, based on the scaled model \a m.
    explicit ScalarTaylorFunction(const DomainType& d, const TaylorModel& m);
    //! \brief Construct a ScalarTaylorFunction over the domain \a d, with scaled power series expansion \a f and error \a e.
    explicit ScalarTaylorFunction(const DomainType& d, const ExpansionType& f, const ErrorType& e=0);

    //! \brief Construct a ScalarTaylorFunction over the domain \a d from the expression \a f.
    explicit ScalarTaylorFunction(const DomainType& d, const ScalarFunction& f);
    //! \brief Construct a ScalarTaylorFunction over the domain \a d from the polynomial \a p.
    explicit ScalarTaylorFunction(const DomainType& d, const Polynomial<Float>& p);
    //! \brief Construct a ScalarTaylorFunction over the domain \a d from the interval polynomial \a p.
    explicit ScalarTaylorFunction(const DomainType& d, const Polynomial<Interval>& p);
    //@}

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to a constant, keeping the same number of arguments.
    ScalarTaylorFunction& operator=(const Float& c) { this->_model=c; return *this; }
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    ScalarTaylorFunction& operator=(const Interval& c) { this->_model=c; return *this; }
    //@}

    //@{
    /*! \name Named constructors. */
    //! \brief Construct a constant quantity in \a as independent variables.
    static ScalarTaylorFunction constant(const DomainType& d, const Float& c);
    //! \brief Construct a constant quantity in \a as independent variables.
    static ScalarTaylorFunction constant(const DomainType& d, const Interval& c);
    //! \brief Construct the quantity \f$x_j\f$ over the domain \a d.
    static ScalarTaylorFunction variable(const DomainType& d, unsigned int j);
    //! \brief Construct the quantity \f$c+\sum g_jx_j\f$ over the domain \a d.
    static ScalarTaylorFunction affine(const DomainType& d, const Float& c, const Vector<Float>& g);
    //! \brief Construct the quantity \f$c+\sum g_jx_j \pm e\f$ over domain \a d.
    static ScalarTaylorFunction affine(const DomainType& d, const Float& x, const Vector<Float>& g, const Float& e) ;

    //! \brief Return the vector of constants with values \a c over domain \a d.
    static Vector<ScalarTaylorFunction> constants(const DomainType& d, const Vector<Float>& c);
    //! \brief Return the vector of constants with interval values \a c over domain \a d.
    static Vector<ScalarTaylorFunction> constants(const DomainType& d, const Vector<Interval>& c);
    //! \brief Return the vector of variables with values \a x over domain \a d.
    static Vector<ScalarTaylorFunction> variables(const DomainType& d);
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
    bool operator==(const ScalarTaylorFunction& tv) const;
    //! \brief Inequality operator.
    bool operator!=(const ScalarTaylorFunction& tv) const { return !(*this==tv); }
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
    ScalarTaylorFunction& truncate() { this->_model.truncate(); return *this; }
    //! \brief Truncate to degree \a deg.
    ScalarTaylorFunction& truncate(uint deg) { this->_model.truncate(deg); return *this; }
    //! \brief Truncate all terms with any coefficient higher than \a a.
    ScalarTaylorFunction& truncate(const MultiIndex& a) { this->_model.truncate(a); return *this; }
    //! \brief Truncate all terms with any coefficient higher than those given by \a b.
    ScalarTaylorFunction& truncate(const MultiIndexBound& b) { this->_model.truncate(b); return *this; }
    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    ScalarTaylorFunction& sweep() { this->_model.sweep(); return *this; }
    //! \brief Remove all terms whose coefficient has magnitude less than \a eps.
    ScalarTaylorFunction& sweep(double eps) { this->_model.sweep(eps); return *this; }
    //! \brief Remove all terms whose degree is higher than \a deg or
    //! whose coefficient has magnitude less than \a eps.
    ScalarTaylorFunction& clean(const TaylorModel::Accuracy& acc) { this->_model.clean(acc); return *this; }
    //! \brief Remove all terms which have high degree or small magnitude.
    ScalarTaylorFunction& clean() { this->_model.clean(); return *this; }
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
    friend bool refines(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2);
    //! \brief Test if the quantities are disjoint.
    friend bool disjoint(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2);
    //! \brief Test if the quantities are disjoint.
    friend ScalarTaylorFunction intersection(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2);
    //! \brief Restrict to a subdomain.
    friend ScalarTaylorFunction restrict(const ScalarTaylorFunction& x, const DomainType& d);
    //! \brief Restrict the values of the \a k<sup>th</sup> variable to the subinterval \a d.
    friend ScalarTaylorFunction restrict(const ScalarTaylorFunction& x, uint k, const Interval& d);
    //! \brief Extend over a larger domain. Only possible if the larger domain is only larger where the smaller domain is a singleton.
    //! The extension is performed keeping \a x constant over the new coordinates.
    friend ScalarTaylorFunction extend(const ScalarTaylorFunction& x, const DomainType& d);
    //@}

    /*! \name Arithmetic operations. */
    //! \brief Inplace addition of another variable.
    friend ScalarTaylorFunction& operator+=(ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
    //! \brief Inplace subtraction of another variable.
    friend ScalarTaylorFunction& operator-=(ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
    //! \brief Inplace addition of a product of two variables.
    friend ScalarTaylorFunction& operator+=(ScalarTaylorFunction& x, const Product<ScalarTaylorFunction,ScalarTaylorFunction>& y);
    //! \brief Inplace addition of an exact floating-point constant.
    friend ScalarTaylorFunction& operator+=(ScalarTaylorFunction& x, const Float& c);
    //! \brief Inplace addition of an interval constant.
    friend ScalarTaylorFunction& operator+=(ScalarTaylorFunction& x, const Interval& c);
    //! \brief Inplace subtraction of an exact floating-point constant.
    friend ScalarTaylorFunction& operator-=(ScalarTaylorFunction& x, const Float& c);
    //! \brief Inplace subtraction of an interval constant.
    friend ScalarTaylorFunction& operator-=(ScalarTaylorFunction& x, const Interval& c);
    //! \brief Inplace multiplication by an exact scalar.
    friend ScalarTaylorFunction& operator*=(ScalarTaylorFunction& x, const Float& c);
    //! \brief Inplace multiplication by an approximate scalar.
    friend ScalarTaylorFunction& operator*=(ScalarTaylorFunction& x, const Interval& c);
    //! \brief Inplace division by an exact scalar.
    friend ScalarTaylorFunction& operator/=(ScalarTaylorFunction& x, const Float& c);
    //! \brief Inplace division by an approximate scalar.
    friend ScalarTaylorFunction& operator/=(ScalarTaylorFunction& x, const Interval& c);

    //! \brief Unary plus.
    friend ScalarTaylorFunction operator+(const ScalarTaylorFunction& x);
    //! \brief Unary minus.
    friend ScalarTaylorFunction operator-(const ScalarTaylorFunction& x);
    //! \brief Addition.
    friend ScalarTaylorFunction operator+(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
    //! \brief Subtraction.
    friend ScalarTaylorFunction operator-(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
    //! \brief Multiplication.
    friend ScalarTaylorFunction operator*(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
    //! \brief Division.
    friend ScalarTaylorFunction operator/(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);

    //! \brief Addition of a scakar.
    friend ScalarTaylorFunction operator+(const ScalarTaylorFunction& x, const Float& c);
    //! \brief Subtraction of a scakar.
    friend ScalarTaylorFunction operator-(const ScalarTaylorFunction& x, const Float& c);
    //! \brief Multiplication by a scakar.
    friend ScalarTaylorFunction operator*(const ScalarTaylorFunction& x, const Float& c);
    //! \brief Division by a scakar.
    friend ScalarTaylorFunction operator/(const ScalarTaylorFunction& x, const Float& c);
    //! \brief Addition of a scakar.
    friend ScalarTaylorFunction operator+(const ScalarTaylorFunction& x, const Interval& c);
    //! \brief Subtraction of a scakar.
    friend ScalarTaylorFunction operator-(const ScalarTaylorFunction& x, const Interval& c);
    //! \brief Multiplication by a scakar.
    friend ScalarTaylorFunction operator*(const ScalarTaylorFunction& x, const Interval& c);
    //! \brief Division by a scakar.
    friend ScalarTaylorFunction operator/(const ScalarTaylorFunction& x, const Interval& c);
    //! \brief Addition of a scakar.
    friend ScalarTaylorFunction operator+(const Float& c, const ScalarTaylorFunction& x);
    //! \brief Subtraction from a scakar.
    friend ScalarTaylorFunction operator-(const Float& c, const ScalarTaylorFunction& x);
    //! \brief Multiplication by a scakar.
    friend ScalarTaylorFunction operator*(const Float& c, const ScalarTaylorFunction& x);
    //! \brief Division through a scalar.
    friend ScalarTaylorFunction operator/(const Float& c, const ScalarTaylorFunction& x);
    //! \brief Addition of a scakar.
    friend ScalarTaylorFunction operator+(const Interval& c, const ScalarTaylorFunction& x);
    //! \brief Subtraction from a scakar.
    friend ScalarTaylorFunction operator-(const Interval& c, const ScalarTaylorFunction& x);
    //! \brief Multiplication by a scakar.
    friend ScalarTaylorFunction operator*(const Interval& c, const ScalarTaylorFunction& x);
    //! \brief Division through a scalar.
    friend ScalarTaylorFunction operator/(const Interval& c, const ScalarTaylorFunction& x);

    //! \brief Maximum. Throws an error if one variable is not greater than the other
    //! over the entire domain.
    friend ScalarTaylorFunction max(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
    //! \brief Minimum. Throws an error if one variable is not greater than the other
    //! over the entire domain.
    friend ScalarTaylorFunction min(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
    //! \brief Addition.
    friend ScalarTaylorFunction add(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
    //! \brief Multiplication.
    friend ScalarTaylorFunction mul(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
    //! \brief Absolute value. Throws an error if one variable is neither greater than
    //! zero over the entire domain, nor less than zero over the entire domain.
    friend ScalarTaylorFunction abs(const ScalarTaylorFunction& x);
    //! \brief Negation.
    friend ScalarTaylorFunction neg(const ScalarTaylorFunction& x);
    //! \brief Reciprocal.
    friend ScalarTaylorFunction rec(const ScalarTaylorFunction& x);
    //! \brief Square.
    friend ScalarTaylorFunction sqr(const ScalarTaylorFunction& x);
    //! \brief Power.
    friend ScalarTaylorFunction pow(const ScalarTaylorFunction& x, int n);
    //! \brief Square root.
    friend ScalarTaylorFunction sqrt(const ScalarTaylorFunction& x);
    //! \brief Natural exponent.
    friend ScalarTaylorFunction exp(const ScalarTaylorFunction& x);
    //! \brief Natural logarithm.
    friend ScalarTaylorFunction log(const ScalarTaylorFunction& x);
    //! \brief Sine.
    friend ScalarTaylorFunction sin(const ScalarTaylorFunction& x);
    //! \brief Cosine.
    friend ScalarTaylorFunction cos(const ScalarTaylorFunction& x);
    //! \brief Tangent.
    friend ScalarTaylorFunction tan(const ScalarTaylorFunction& x);
    //! \brief Inverse sine.
    friend ScalarTaylorFunction asin(const ScalarTaylorFunction& x);
    //! \brief Inverse cosine.
    friend ScalarTaylorFunction acos(const ScalarTaylorFunction& x);
    //! \brief Inverse tangent.
    friend ScalarTaylorFunction atan(const ScalarTaylorFunction& x);
    //@}

    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const ScalarTaylorFunction& x);
    //@}

  public:
    ScalarTaylorFunction& clobber() { this->_model.clobber(); return *this; }
    ScalarTaylorFunction& clobber(uint o) { this->_model.clobber(o); return *this; }
    ScalarTaylorFunction& clobber(uint so, uint to) { this->_model.clobber(so,to); return *this; }
};


inline tribool operator>(const ScalarTaylorFunction& x, const Float& c) {
    Interval r=x.range(); if(r.lower()>c) { return true; } else if(r.upper()<=c) { return false; } else { return indeterminate; } }
inline tribool operator<(const ScalarTaylorFunction& x, const Float& c) {
    Interval r=x.range(); if(r.lower()<c) { return true; } else if(r.upper()>=c) { return false; } else { return indeterminate; } }

inline tribool operator>(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y) { return (x-y)>0; }
inline tribool operator<(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y) { return (x-y)<0; }

inline ScalarTaylorFunction& operator+=(ScalarTaylorFunction& x, const Float& c) {
    x._model+=c; return x; }
inline ScalarTaylorFunction& operator-=(ScalarTaylorFunction& x, const Float& c) {
    x._model-=c; return x; }
inline ScalarTaylorFunction& operator*=(ScalarTaylorFunction& x, const Float& c) {
    x._model*=c; return x; }
inline ScalarTaylorFunction& operator/=(ScalarTaylorFunction& x, const Float& c) {
    x._model*=c; return x; }

inline ScalarTaylorFunction& operator+=(ScalarTaylorFunction& x, const Interval& c) {
    x._model+=c; return x; }
inline ScalarTaylorFunction& operator-=(ScalarTaylorFunction& x, const Interval& c) {
    x._model-=c; return x; }
inline ScalarTaylorFunction& operator*=(ScalarTaylorFunction& x, const Interval& c) {
    x._model*=c; return x; }
inline ScalarTaylorFunction& operator/=(ScalarTaylorFunction& x, const Interval& c) {
    x._model*=c; return x; }


inline ScalarTaylorFunction operator+(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,x._model); }
inline ScalarTaylorFunction operator-(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,-x._model); }

inline ScalarTaylorFunction operator+(const ScalarTaylorFunction& x, const Float& c) {
    return ScalarTaylorFunction(x._domain,x._model+c); }
inline ScalarTaylorFunction operator-(const ScalarTaylorFunction& x, const Float& c) {
    return ScalarTaylorFunction(x._domain,x._model-c); }
inline ScalarTaylorFunction operator*(const ScalarTaylorFunction& x, const Float& c) {
    return ScalarTaylorFunction(x._domain,x._model*c); }
inline ScalarTaylorFunction operator/(const ScalarTaylorFunction& x, const Float& c) {
    return ScalarTaylorFunction(x._domain,x._model/c); }
inline ScalarTaylorFunction operator+(const ScalarTaylorFunction& x, const Interval& c) {
    return ScalarTaylorFunction(x._domain,x._model+c); }
inline ScalarTaylorFunction operator-(const ScalarTaylorFunction& x, const Interval& c) {
    return ScalarTaylorFunction(x._domain,x._model-c); }
inline ScalarTaylorFunction operator*(const ScalarTaylorFunction& x, const Interval& c) {
    return ScalarTaylorFunction(x._domain,x._model*c); }
inline ScalarTaylorFunction operator/(const ScalarTaylorFunction& x, const Interval& c) {
    return ScalarTaylorFunction(x._domain,x._model/c); }
inline ScalarTaylorFunction operator+(const Float& c, const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,c+x._model); }
inline ScalarTaylorFunction operator-(const Float& c, const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,c-x._model); }
inline ScalarTaylorFunction operator*(const Float& c, const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,c*x._model); }
inline ScalarTaylorFunction operator/(const Float& c, const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,c/x._model); }
inline ScalarTaylorFunction operator+(const Interval& c, const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,c+x._model); }
inline ScalarTaylorFunction operator-(const Interval& c, const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,c-x._model); }
inline ScalarTaylorFunction operator*(const Interval& c, const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,c*x._model); }
inline ScalarTaylorFunction operator/(const Interval& c, const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,c/x._model); }

ScalarTaylorFunction operator+(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
ScalarTaylorFunction operator-(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
ScalarTaylorFunction operator*(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
ScalarTaylorFunction operator/(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);

inline ScalarTaylorFunction abs(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,abs(x._model)); }
inline ScalarTaylorFunction neg(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,abs(x._model)); }
inline ScalarTaylorFunction rec(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,rec(x._model)); }
inline ScalarTaylorFunction sqr(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,sqr(x._model)); }
inline ScalarTaylorFunction pow(const ScalarTaylorFunction& x, int n) {
    return ScalarTaylorFunction(x._domain,pow(x._model,n)); }
inline ScalarTaylorFunction sqrt(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,sqrt(x._model)); }
inline ScalarTaylorFunction exp(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,exp(x._model)); }
inline ScalarTaylorFunction log(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,log(x._model)); }
inline ScalarTaylorFunction sin(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,sin(x._model)); }
inline ScalarTaylorFunction cos(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,cos(x._model)); }
inline ScalarTaylorFunction tan(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,tan(x._model)); }
inline ScalarTaylorFunction asin(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,asin(x._model)); }
inline ScalarTaylorFunction acos(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,acos(x._model)); }
inline ScalarTaylorFunction atan(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,atan(x._model)); }



inline Interval evaluate(const ScalarTaylorFunction& tv, const Vector<Interval>& x) {
    return tv.evaluate(x); }

inline ScalarTaylorFunction antiderivative(const ScalarTaylorFunction& x, uint k) {
    Interval sf=rad_ivl(x.domain()[k]);
    return ScalarTaylorFunction(x.domain(),antiderivative(x.model(),k)*sf); }

inline ScalarTaylorFunction derivative(const ScalarTaylorFunction& x, uint k) {
    Interval sf=1/rad_ivl(x.domain()[k]);
    return ScalarTaylorFunction(x.domain(),derivative(x.model(),k)*sf); }

inline ScalarTaylorFunction embed(const ScalarTaylorFunction& tv1, const Interval& dom2) {
    return ScalarTaylorFunction(join(tv1.domain(),dom2),embed(tv1.model(),1u)); }
inline ScalarTaylorFunction embed(const ScalarTaylorFunction& tv1, const Vector<Interval>& dom2) {
    return ScalarTaylorFunction(join(tv1.domain(),dom2),embed(tv1.model(),dom2.size())); }
inline ScalarTaylorFunction embed(const Vector<Interval>& dom1, const ScalarTaylorFunction& tv2) {
    return ScalarTaylorFunction(join(dom1,tv2.domain()),embed(dom1.size(),tv2.model())); }




Vector<Interval> evaluate(const VectorTaylorFunction& f, const Vector<Interval>& x);
VectorTaylorFunction partial_evaluate(const VectorTaylorFunction& f, uint k, const Interval& c);
VectorTaylorFunction embed(const VectorTaylorFunction& tv1, const Vector<Interval>& d2);
VectorTaylorFunction embed(const VectorTaylorFunction& tv1, const Interval& d2);
VectorTaylorFunction embed(const Vector<Interval>& d1, const VectorTaylorFunction& tv2);
VectorTaylorFunction restrict(const VectorTaylorFunction&, const Vector<Interval>& bx);
bool refines(const VectorTaylorFunction&, const VectorTaylorFunction&);
bool disjoint(const VectorTaylorFunction&, const VectorTaylorFunction&);
VectorTaylorFunction intersection(const VectorTaylorFunction&, const VectorTaylorFunction&);
ScalarTaylorFunction compose(const ScalarFunction&, const VectorTaylorFunction&);
ScalarTaylorFunction compose(const ScalarTaylorFunction&, const VectorTaylorFunction&);
VectorTaylorFunction compose(const VectorTaylorFunction&, const VectorTaylorFunction&);
VectorTaylorFunction compose(const VectorFunction&, const VectorTaylorFunction&);
VectorTaylorFunction antiderivative(const VectorTaylorFunction&, uint);
VectorTaylorFunction implicit(const VectorTaylorFunction&);
ScalarTaylorFunction implicit(const ScalarFunction&, const VectorTaylorFunction&);
VectorTaylorFunction flow(const VectorTaylorFunction& vf, const Vector<Interval>& d, const Interval& h, uint o);
VectorTaylorFunction flow(const VectorTaylorFunction& vf, const Vector<Interval>& d, const Float& h, uint o);
VectorTaylorFunction flow(const VectorFunction& vf, const Vector<Interval>& d, const Float& h, uint o);
VectorTaylorFunction parameterised_flow(const VectorTaylorFunction& vf, const Vector<Interval>& d, const Float& h, uint o);

ScalarTaylorFunction unchecked_compose(const ScalarTaylorFunction&, const VectorTaylorFunction&);
VectorTaylorFunction unchecked_compose(const VectorTaylorFunction&, const VectorTaylorFunction&);
VectorTaylorFunction unchecked_implicit(const VectorTaylorFunction&);
VectorTaylorFunction unchecked_flow(const VectorTaylorFunction& vf, const Vector<Interval>& d, const Interval& h, uint o);
VectorTaylorFunction unchecked_flow(const VectorTaylorFunction& vf, const Vector<Interval>& d, const Float& h, uint o);


/*! \brief A taylor_model with multivalued output using the TaylorModel class.
 *
 *  See also TaylorModel, ScalarTaylorFunction, VectorTaylorFunction.
 */
class VectorTaylorFunction {
    friend class VectorTaylorFunctionElementReference;

    typedef Float R;
    typedef Interval I;
  public:
    /*! \brief Default constructor constructs a Taylor model of order zero with no arguments and no result variables. */
    VectorTaylorFunction();

    /*! \brief Construct the zero vector function over an unspecified domain. */
    explicit VectorTaylorFunction(unsigned int result_size);

    /*! \brief Construct from a result size and a domain. */
    VectorTaylorFunction(unsigned int result_size, const Vector<Interval>& domain);

    /*! \brief Construct a vector function all of whose components are the same. */
    VectorTaylorFunction(unsigned int result_size, const ScalarTaylorFunction& scalar_function);

    /*! \brief Construct from a domain and the expansion. */
    VectorTaylorFunction(const Vector<Interval>& domain,
                   const Vector< Expansion<Float> >& expansion);

    /*! \brief Construct from a domain, and expansion and errors. */
    VectorTaylorFunction(const Vector<Interval>& domain,
                   const Vector< Expansion<Float> >& expansion,
                   const Vector<Float>& error);

    /*! \brief Construct from a domain and the models. */
    explicit VectorTaylorFunction(const Vector<Interval>& domain, const Vector<TaylorModel>& variables);

    /*! \brief Construct from a domain and a function. */
    VectorTaylorFunction(const Vector<Interval>& domain,
                   const VectorFunction& function);

    /*! \brief Construct from a domain, a function, and accuracy paramters. */
    VectorTaylorFunction(const Vector<Interval>& domain,
                   const VectorFunction& function,
                   shared_ptr<TaylorModel::Accuracy> accuracy_ptr);

    /*! \brief Construct from a domain and a polynomial. */
    VectorTaylorFunction(const Vector<Interval>& domain,
                   const Vector< Polynomial<Float> >& polynomial);

    /*! \brief Construct from a domain and a n interval polynomial. */
    VectorTaylorFunction(const Vector<Interval>& domain,
                   const Vector< Polynomial<Interval> >& polynomial);

    /*! \brief Construct from a vector of scalar Taylor functions. */
    explicit VectorTaylorFunction(const Vector<ScalarTaylorFunction>& components);

    /*! \brief Construct from a vector expression. */
    template<class E> explicit VectorTaylorFunction(const boost::numeric::ublas::vector_expression<E>& ve);

    /*! \brief Equality operator. */
    bool operator==(const VectorTaylorFunction& p) const;
    /*! \brief Inequality operator. */
    bool operator!=(const VectorTaylorFunction& p) const;

    // Data access
    /*! \brief The accuracy parameter used to control approximation of the Taylor function. */
    shared_ptr<TaylorModel::Accuracy> accuracy_ptr() const;
    /*! \brief Set the accuracy parameter used to control approximation of the Taylor function. */
    void set_accuracy(shared_ptr<TaylorModel::Accuracy> acc);
    /*! \brief The data used to define the domain of the Taylor model. */
    const Vector<Interval>& domain() const;
    /*! \brief A rough bound for the range of the function. */
    const Vector<Interval> codomain() const;
    /*! \brief The centre of the Taylor model. */
    const Vector<Float> centre() const;
    /*! \brief The range of the Taylor model. */
    const Vector<Interval> range() const;
    /*! \brief The data used to define the centre of the Taylor model. */
    const Vector<TaylorModel>& models() const;

    /*! \brief The size of the argument. */
    uint argument_size() const;
    /*! \brief The size of the result. */
    uint result_size() const;

    /*! \brief Get the \a ith Taylor variable */
    ScalarTaylorFunction get(uint i) const;
    /*! \brief Set the \a ith Taylor variable */
    void set(uint i, const ScalarTaylorFunction& te);
    /*! \brief The \a ith Taylor variable */
    ScalarTaylorFunction operator[](uint i) const;
    /*! \brief The \a ith Taylor variable */
    VectorTaylorFunctionElementReference operator[](uint i);
    /*! \brief Evaluate the Taylor model at the point \a x. */
    Vector<Interval> evaluate(const Vector<Interval>& x) const;
    /*! \brief Evaluate the Taylor model at the point \a x. */
    Vector<Interval> evaluate(const Vector<Float>& x) const;
    /*! \brief Compute an approximation to Jacobian derivative of the Taylor model sat the point \a x. */
    Matrix<Interval> jacobian(const Vector<Interval>& x) const;

    /*! \brief Truncate to a model of lower order and/or smoothness. */
    VectorTaylorFunction& truncate(ushort degree);
    /*! \brief Set the error to zero. */
    VectorTaylorFunction& clobber();

    /*! \brief The constant Taylor model with range \a r and argument domain \a d. */
    static VectorTaylorFunction constant(const Vector<Interval>& d, const Vector<Interval>& r);
    /*! \brief The constant Taylor model with result \a c and argument domain \a d. */
    static VectorTaylorFunction constant(const Vector<Interval>& d, const Vector<Float>& c);
    /*! \brief The identity Taylor model on domain \a d. */
    static VectorTaylorFunction identity(const Vector<Interval>& d);

    /*! \brief Convert to an interval polynomial. */
    Vector< Polynomial<Interval> > polynomial() const;

    /*! \brief Truncate terms higher than \a bd. */
    VectorTaylorFunction& truncate(const MultiIndexBound& bd);

    /*! \brief Write to an output stream. */
    std::ostream& write(std::ostream& os) const;

    /*! \brief Inplace addition. */
    friend VectorTaylorFunction& operator+=(VectorTaylorFunction& f, const VectorTaylorFunction& g);
    /*! \brief Inplace subtraction. */
    friend VectorTaylorFunction& operator-=(VectorTaylorFunction& f, const VectorTaylorFunction& g);
    /*! \brief Inplace addition. */
    friend VectorTaylorFunction& operator+=(VectorTaylorFunction& f, const Vector<Interval>& c);
    /*! \brief Inplace subtraction. */
    friend VectorTaylorFunction& operator-=(VectorTaylorFunction& f, const Vector<Interval>& c);
    /*! \brief Inplace scalar multiplication. */
    friend VectorTaylorFunction& operator*=(VectorTaylorFunction& f, const Float& c);
    /*! \brief Inplace scalar division. */
    friend VectorTaylorFunction& operator/=(VectorTaylorFunction& f, const Float& c);

    /*! \brief Negation. */
    friend VectorTaylorFunction operator-(const VectorTaylorFunction& f);
    /*! \brief Addition. */
    friend VectorTaylorFunction operator+(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2);
    /*! \brief Subtraction. */
    friend VectorTaylorFunction operator-(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2);

    /*! \brief Addition of a constant. */
    friend VectorTaylorFunction operator+(const VectorTaylorFunction& f, const Vector<Float>& c);
    /*! \brief Subtraction of a constant. */
    friend VectorTaylorFunction operator-(const VectorTaylorFunction& f, const Vector<Float>& c);
    /*! \brief Multiplication by a scalar. */
    friend VectorTaylorFunction operator*(const Float& c, const VectorTaylorFunction& f);
    /*! \brief Multiplication by a scalar. */
    friend VectorTaylorFunction operator*(const VectorTaylorFunction& f, const Float& c);
    /*! \brief Division by a scalar. */
    friend VectorTaylorFunction operator/(const VectorTaylorFunction& f, const Float& c);
    /*! \brief Addition of a constant. */
    friend VectorTaylorFunction operator+(const VectorTaylorFunction& f, const Vector<Interval>& c);
    /*! \brief Subtraction of a constant. */
    friend VectorTaylorFunction operator-(const VectorTaylorFunction& f, const Vector<Interval>& c);
    /*! \brief Multiplication by a scalar. */
    friend VectorTaylorFunction operator*(const Interval& c, const VectorTaylorFunction& f);
    /*! \brief Multiplication by a scalar. */
    friend VectorTaylorFunction operator*(const VectorTaylorFunction& f, const Interval& c);
    /*! \brief Division by a scalar. */
    friend VectorTaylorFunction operator/(const VectorTaylorFunction& f, const Interval& c);
    /*! \brief Multiplication by a matrix. */
    friend VectorTaylorFunction operator*(const Matrix<Float>& A, const VectorTaylorFunction& f);
    /*! \brief Multiplication by a matrix. */
    friend VectorTaylorFunction operator*(const Matrix<Interval>& A, const VectorTaylorFunction& f);

    //! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
    friend VectorTaylorFunction compose(const VectorFunction& f, const VectorTaylorFunction& g);
    //! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
    friend ScalarTaylorFunction compose(const ScalarTaylorFunction& f, const VectorTaylorFunction& g);
    //! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
    friend VectorTaylorFunction compose(const VectorTaylorFunction& f, const VectorTaylorFunction& g);
    //! \brief Antiderivative of \a f with respect to variable \a k.
    friend VectorTaylorFunction antiderivative(const VectorTaylorFunction& f, uint k);
    //! \brief The flow of the vector field \a vf defined over a space domain \a d over a time interval \a t.
    friend VectorTaylorFunction flow(const VectorTaylorFunction& vf, const Vector<Interval>& d, const Interval& t, uint o);
    //! \brief Compute the implicit function of \a f satisfying \f$f(c,h(c))=0\f$,
    //! where \f$c\f$ is the centre of the domain of \f$f\f$.
    friend VectorTaylorFunction implicit(const VectorTaylorFunction& f);
    //! \brief Compute the inverse function of \a f based at the centre of the domain. */
    friend VectorTaylorFunction inverse(const VectorTaylorFunction& f);
    //! \brief Compute the inverse function of \a f based at \f$f(c)\f$. */
    friend VectorTaylorFunction inverse(const VectorTaylorFunction& f, const Vector<Float>& c);
    //! \brief Compute the function \f$(f,g)(x)=(f(x),g(x))\f$.
    friend VectorTaylorFunction join(const VectorTaylorFunction& f, const VectorTaylorFunction& g);
    friend VectorTaylorFunction join(const VectorTaylorFunction& f, const ScalarTaylorFunction& g);
    friend VectorTaylorFunction join(const ScalarTaylorFunction& f, const ScalarTaylorFunction& g);
    friend VectorTaylorFunction join(const ScalarTaylorFunction& f, const VectorTaylorFunction& g);
    //! \brief Compute the function \f$(f\oplus g)(x,y)=(f(x),g(y))\f$.
    friend VectorTaylorFunction combine(const VectorTaylorFunction& f, const VectorTaylorFunction& g);
    friend VectorTaylorFunction combine(const VectorTaylorFunction& f, const ScalarTaylorFunction& g);
    friend VectorTaylorFunction combine(const ScalarTaylorFunction& f, const VectorTaylorFunction& g);
    friend VectorTaylorFunction combine(const ScalarTaylorFunction& f, const ScalarTaylorFunction& g);
    //! \brief Restrict the function \a f to a subdomain \a d.
    friend VectorTaylorFunction restrict(const VectorTaylorFunction& f, const Vector<Interval>& d);
    //! \brief Tests if a function \a f refines another function \a g.
    //! To be a refinement, the domain of \a f must contain the domain of \a g.
    friend bool refines(const VectorTaylorFunction& f, const VectorTaylorFunction& g);

    // For compatibility wit Vector.
    uint size() const { return this->result_size(); }
  private:
    array< array<Interval> > _powers(const Vector<Interval>&) const;
    void _compute_jacobian() const;
    void _set_argument_size(uint n);
    uint _compute_maximum_component_size() const;
    void _resize(uint rs, uint as, ushort d, ushort s);

  private:
    /* Domain of definition. */
    Vector<Interval> _domain;
    Vector<TaylorModel> _models;
};

VectorTaylorFunction join(const VectorTaylorFunction& f, const VectorTaylorFunction& g);
VectorTaylorFunction join(const VectorTaylorFunction& f, const ScalarTaylorFunction& g);
VectorTaylorFunction join(const ScalarTaylorFunction& f, const ScalarTaylorFunction& g);
VectorTaylorFunction join(const ScalarTaylorFunction& f, const VectorTaylorFunction& g);
VectorTaylorFunction combine(const VectorTaylorFunction& f, const VectorTaylorFunction& g);
VectorTaylorFunction combine(const VectorTaylorFunction& f, const ScalarTaylorFunction& g);
VectorTaylorFunction combine(const ScalarTaylorFunction& f, const VectorTaylorFunction& g);
VectorTaylorFunction combine(const ScalarTaylorFunction& f, const ScalarTaylorFunction& g);

std::ostream& operator<<(std::ostream&, const VectorTaylorFunction&);


template<class E> VectorTaylorFunction::VectorTaylorFunction(const boost::numeric::ublas::vector_expression<E>& ve)            : _domain(), _models(ve().size())
{
    if(ve().size()!=0) { this->_domain=ve()[0].domain(); }
    for(uint i=0; i!=ve().size(); ++i) { this->set(i,ve()[i]); }
}

class VectorTaylorFunctionElementReference
{
    friend class ScalarTaylorFunction;
    friend class VectorTaylorFunction;
 public:
    VectorTaylorFunctionElementReference(VectorTaylorFunction& c, uint i) : _c(&c), _i(i) { }
    operator ScalarTaylorFunction () const { return this->_c->get(this->_i); }
    void operator=(const VectorTaylorFunctionElementReference& x) { this->_c->set(this->_i,x._c->get(x._i)); }
    void operator=(const ScalarTaylorFunction& x) { this->_c->set(this->_i,x); }
    void set_error(double e) { this->_c->_models[this->_i].set_error(e); }
  private:
    VectorTaylorFunction* _c; uint _i;
};


} // namespace Ariadne

#endif // ARIADNE_TAYLOR_FUNCTION_H
