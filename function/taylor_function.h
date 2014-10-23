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
#include "utility/container.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "function/taylor_model.h"

#include "function/function_interface.h"
#include "function/function_mixin.h"
#include "function/function_model.h"

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Polynomial;

template<class X> class ScalarFunction;
typedef EffectiveScalarFunction EffectiveScalarFunction;
typedef ValidatedScalarFunction ValidatedScalarFunction;
template<class X> class VectorFunction;
typedef EffectiveVectorFunction EffectiveVectorFunction;
typedef ValidatedVectorFunction ValidatedVectorFunction;

class MultiIndex;
template<class X> class TaylorModel;
class ScalarTaylorFunction;
class VectorTaylorFunction;
class TaylorFunctionFactory;

inline ApproximateNumber med_apprx(ExactInterval const& ivl) {
    return ApproximateNumber(half_exact(add_approx(ivl.lower_raw(),ivl.upper_raw())));
}

inline ApproximateNumber rad_apprx(ExactInterval const& ivl) {
    return ApproximateNumber(half_exact(sub_approx(ivl.upper_raw(),ivl.lower_raw())));
}

inline ValidatedNumber med_val(ExactInterval const& ivl) {
    return half(ivl.lower()+ivl.upper());
}

inline ValidatedNumber rad_val(ExactInterval const& ivl) {
    return half(ivl.upper()-ivl.lower());
}

template<template<class> class T> inline T<ApproximateNumber> unscale(const T<ApproximateNumber>& x, const ExactInterval& d) {
    ApproximateNumber c(med_apprx(d));
    ApproximateNumber r(rad_apprx(d));
    return (x-c)/r;
}

template<template<class> class T> inline T<ValidatedNumber> unscale(const T<ValidatedNumber>& x, const ExactInterval& d) {
    ValidatedNumber c(med_val(d));
    ValidatedNumber r(rad_val(d));
    return (x-c)/r;
}


inline ApproximateNumber unscale(const ApproximateNumber& x, const ExactInterval& d) {
    ApproximateNumber c(med_apprx(d));
    ApproximateNumber r(rad_apprx(d));
    return (x-c)/r;
}

inline ValidatedNumber unscale(const ValidatedNumber& x, const ExactInterval& d) {
    ValidatedNumber c(med_val(d));
    ValidatedNumber r(rad_val(d));
    return (x-c)/r;
}

template<class X> Vector<X> unscale(const Vector<X>& x, const ExactBox& d) {
    Vector<X> r(x);
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=unscale(x[i],d[i]);
    }
    return r;
}

// Remove the error term
ScalarTaylorFunction midpoint(const ScalarTaylorFunction& x);

// Set the value of the \a kth variable to c
ScalarTaylorFunction partial_evaluate(const ScalarTaylorFunction& f, uint k, const RawFloat& c);
ScalarTaylorFunction partial_evaluate(const ScalarTaylorFunction& f, uint k, const ValidatedNumber& c);
// Evaluate a scalar Taylor function on a vector.
ValidatedNumber evaluate(const ScalarTaylorFunction& x, const Vector<ValidatedNumber>& sy);

// Restrict the \a kth variable to lie in the interval \a d.
ScalarTaylorFunction restrict(const ScalarTaylorFunction& x, uint k, const ExactInterval& d);
// Restrict to a smaller domain. REQUIRED
ScalarTaylorFunction restrict(const ScalarTaylorFunction& x, const ExactBox& d);
// Extend to a larger domain. REQUIRED
ScalarTaylorFunction extend(const ScalarTaylorFunction& x, const ExactBox& d);

// Compose with an function.
ScalarTaylorFunction compose(const ValidatedScalarFunction& x, const VectorTaylorFunction& y);

// Substitute \a h into the \a k<sup>th</sup> argument of \a f.
ScalarTaylorFunction substitute(const ScalarTaylorFunction& f, uint k, const ScalarTaylorFunction& h);

// Test if a variable refines another
bool refines(const ScalarTaylorFunction& tv1, const ScalarTaylorFunction& tv2);
// Test if two variables definitely represent different quantities
bool disjoint(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2);
// Test if two variables definitely represent different quantities
ScalarTaylorFunction intersection(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2);


// Split the variable over two domains, subdividing along the independent variable j.
Pair<ScalarTaylorFunction,ScalarTaylorFunction> split(const ScalarTaylorFunction& x, uint j);


// Embed the variable in a space of higher dimension
ScalarTaylorFunction embed(const ScalarTaylorFunction& tv1, const ExactInterval& d2);
ScalarTaylorFunction embed(const ExactBox& d1, const ScalarTaylorFunction& tv2);
ScalarTaylorFunction embed(const ExactBox& d1, const ScalarTaylorFunction& tv2, const ExactBox& d3);

// Antidifferentiation operator
ScalarTaylorFunction antiderivative(const ScalarTaylorFunction& x, uint k, ExactNumber c);
ScalarTaylorFunction antiderivative(const ScalarTaylorFunction& x, uint k);
ScalarTaylorFunction derivative(const ScalarTaylorFunction& x, uint k);



class VectorTaylorFunctionElementReference;

/*! \ingroup FunctionModelSubModule
 *  \brief A ScalarTaylorFunction is a type of FunctionModel in which a the restriction of a scalar function \f$f:\R^n\rightarrow\R\f$ on a domain \f$D\f$ is approximated by polynomial \f$p\f$ with uniform error \f$e\f$.
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
 * \sa Expansion, TaylorModel, VectorTaylorFunction, TaylorConstrainedImageSet.
 */
class ScalarTaylorFunction
    : public ScalarFunctionModelMixin<ScalarTaylorFunction, ValidatedTag>
{
  public:
    typedef ExactBox DomainType;
    typedef ExactInterval RangeType;
    typedef ValidatedNumber NumericType;
    typedef TaylorModel<ValidatedNumber> ModelType;
    typedef Expansion<CoefficientType> ExpansionType;
    typedef Ariadne::ErrorType ErrorType;
  private:
    static const CoefficientType _zero;
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
    //explicit ScalarTaylorFunction(const DomainType& d);
    explicit ScalarTaylorFunction(const DomainType& d, Sweeper swp);
    //! \brief Construct a ScalarTaylorFunction over the domain \a d, based on the scaled model \a m.
    explicit ScalarTaylorFunction(const DomainType& d, const TaylorModel<ValidatedNumber>& m);
    explicit ScalarTaylorFunction(const DomainType& d, const Expansion<CoefficientType>& p, const ErrorType& e, const Sweeper& swp);
    explicit ScalarTaylorFunction(const DomainType& d, const Expansion<RawFloat>& p, const double& e, const Sweeper& swp);

    explicit ScalarTaylorFunction(const ScalarFunctionModel<ValidatedTag>& f);
    ScalarTaylorFunction& operator=(const ScalarFunctionModel<ValidatedTag>& f);

    //! \brief Construct a ScalarTaylorFunction over the domain \a d from the function \a f.
    explicit ScalarTaylorFunction(const DomainType& d, const ValidatedScalarFunction& f, Sweeper swp);
    //@}

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    ScalarTaylorFunction& operator=(const ValidatedNumber& c) { this->_model=c; return *this; }
    template<class X, typename std::enable_if<std::is_same<X,int>::value,int>::type=0>
        ScalarTaylorFunction& operator=(const X& c) { return (*this)=static_cast<ValidatedNumber>(c); }
    //@}

    //@{
    /*! \name Named constructors. */
    //! \brief Construct the identity function over the one-dimensional domain \a d.
    static ScalarTaylorFunction identity(const ExactInterval& d, Sweeper swp);
    //! \brief Construct a zero function over domain \a d.
    static ScalarTaylorFunction zero(const DomainType& d, Sweeper swp);
    //! \brief Construct a constant quantity in \a as independent variables.
    static ScalarTaylorFunction constant(const DomainType& d, const ValidatedNumber& c, Sweeper swp);
    template<class X> static ScalarTaylorFunction constant(const DomainType& d, const X& c, Sweeper swp) {
        return constant(d,numeric_cast<ValidatedNumber>(c),swp); }
    //! \brief Construct the quantity \f$x_j\f$ over the domain \a d.
    static ScalarTaylorFunction coordinate(const DomainType& d, unsigned int j, Sweeper swp);
    //! \brief Construct the quantity \f$x_j\f$ over the domain \a d (Deprecated).
    static ScalarTaylorFunction variable(const DomainType& d, unsigned int j, Sweeper swp);
    //! \brief Construct the quantity \f$c+\sum g_jx_j\f$ over the domain \a d.
    static ScalarTaylorFunction affine(const DomainType& d, const CoefficientType& c, const Vector<CoefficientType>& g, Sweeper swp);
    //! \brief Construct the quantity \f$c+\sum g_jx_j \pm e\f$ over domain \a d.
    static ScalarTaylorFunction affine(const DomainType& d, const CoefficientType& x, const Vector<CoefficientType>& g, const ErrorType& e, Sweeper swp) ;

    //! \brief Return the vector of constants with values \a c over domain \a d.
    static Vector<ScalarTaylorFunction> constants(const DomainType& d, const Vector<ExactNumber>& c, Sweeper swp);
    //! \brief Return the vector of constants with interval values \a c over domain \a d.
    static Vector<ScalarTaylorFunction> constants(const DomainType& d, const Vector<ValidatedNumber>& c, Sweeper swp);
    //! \brief Return the vector of variables with values \a x over domain \a d.
    static Vector<ScalarTaylorFunction> variables(const DomainType& d, Sweeper swp);
    //! \brief Return the vector of variables in the range \a imin to \a imax with values \a x over domain \a d.
    static Vector<ScalarTaylorFunction> variables(const DomainType& d, uint imin, uint imax, Sweeper swp);
    //@}

    //@{
    /*! \name Data access */
    /*! \brief The accuracy parameter used to control approximation of the Taylor function. */
    Sweeper sweeper() const;
    //! \brief The domain of the quantity.
    const DomainType& domain() const { return this->_domain; }
    //! \brief An over-approximation to the range of the quantity; not necessarily tight.
    const ExactInterval codomain() const { UpperInterval rng=this->_model.range(); return ExactInterval(rng.lower_raw(),rng.upper_raw()); }
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
    CoefficientType centre() { return this->_model.value(); }
    //! \brief A reference to the constant term in the expansion.
    CoefficientType& value() { return this->_model.value(); }
    //! \brief The constant term in the expansion.
    const CoefficientType& value() const { return this->_model.value(); }
    //! \brief The gradient at the centre of the domain.
    const CoefficientType gradient_value(Nat i) const {
        return make_exact(this->_model.gradient(i)/this->_domain[i].radius()); }
    //! \brief A polynomial representation.
    Polynomial<ValidatedNumber> polynomial() const;
    //! \brief A multivalued function equal to the model on the domain.
    ValidatedScalarFunction function() const;

    //! \brief Set the error of the expansion.
    void set_error(const ErrorType& ne) { this->_model.set_error(ne); }
    //! \brief Set the constant term in the expansion.
    void set_value(const CoefficientType& c) { this->_model.set_value(c); }

    //! \brief The coefficient of the term in $x^a$.
    const CoefficientType& operator[](const MultiIndex& a) const { return this->_model[a]; }
    //! \brief A read/write reference to the coefficient of the term in $x^a$.
    CoefficientType& operator[](const MultiIndex& a) { return this->_model[a]; }

    //! \brief An iterator to the first term in the expansion expansion.
    iterator begin() { return this->_model.begin(); }
    //! \brief A constant iterator to the first term in the expansion expansion.
    const_iterator begin() const { return this->_model.begin(); }
    //! \brief An iterator to the end of the expansion expansion.
    iterator end() { return this->_model.end(); }
    //! \brief A constant iterator to the end of the expansion expansion.
    const_iterator end() const { return this->_model.end(); }

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
    UpperInterval range() const { return this->_model.range(); }
    //! \brief Evaluate the quantity at the point \a x.
    ApproximateNumber evaluate(const Vector<ApproximateNumber>& x) const;
    //! \brief Evaluate the quantity over the interval of points \a x.
    ValidatedNumber evaluate(const Vector<ValidatedNumber>& x) const;
    //! \brief Evaluate the quantity over the interval of points \a x.
    ValidatedNumber operator()(const Vector<ValidatedNumber>& x) const;
    using ScalarFunctionMixin< ScalarTaylorFunction, ValidatedTag>::evaluate;

    //@}

    //@{
    /*! \name Simplification operations. */
   //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    ScalarTaylorFunction& sweep() { this->_model.sweep(); return *this; }
    //! \brief Remove all terms whose degree is higher than \a deg or
    //! whose coefficient has magnitude less than \a eps.
    ScalarTaylorFunction& sweep(const SweeperInterface& swp) { this->_model.sweep(swp); return *this; }
    //@}

    //@{
    /*! \name Accuracy parameters. */
    //! \copydoc TaylorModel<ValidatedNumber>::set_sweeper()
    void set_sweeper(const Sweeper& swp) { this->_model.set_sweeper(swp); }
    //@}

    //@{
    /*! \name Non-arithmetic operations. */
    //! \brief Restrict to a subdomain.
    void restrict(const DomainType& d);
    //! \brief Test if the quantity is a better approximation than \a t throughout the domain.
    friend bool refines(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2);
    //! \brief Test if the quantities are disjoint.
    friend bool disjoint(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2);
    //! \brief Test if the quantities are disjoint.
    friend ScalarTaylorFunction intersection(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2);
    //! \brief Restrict to a subdomain.
    friend ScalarTaylorFunction restrict(const ScalarTaylorFunction& x, const DomainType& d);
    //! \brief Restrict the values of the \a k<sup>th</sup> variable to the subinterval \a d.
    friend ScalarTaylorFunction restrict(const ScalarTaylorFunction& x, uint k, const ExactInterval& d);
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
    //! \brief Inplace addition of an interval constant.
    friend ScalarTaylorFunction& operator+=(ScalarTaylorFunction& x, const ValidatedNumber& c);
    //! \brief Inplace subtraction of an interval constant.
    friend ScalarTaylorFunction& operator-=(ScalarTaylorFunction& x, const ValidatedNumber& c);
    //! \brief Inplace multiplication by an approximate scalar.
    friend ScalarTaylorFunction& operator*=(ScalarTaylorFunction& x, const ValidatedNumber& c);
    //! \brief Inplace division by an approximate scalar.
    friend ScalarTaylorFunction& operator/=(ScalarTaylorFunction& x, const ValidatedNumber& c);

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
    std::ostream& write(std::ostream& os) const;
    /*! \brief Write a full representation to an output stream. */
    std::ostream& repr(std::ostream& os) const;
    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const ScalarTaylorFunction& x);
    //@}

  public:
    Void clobber() { this->_model.clobber(); }
  private:
    friend class TaylorFunctionFactory;
    friend class ScalarFunctionMixin<ScalarTaylorFunction, ValidatedTag>;
    friend class ScalarFunctionModelMixin<ScalarTaylorFunction, ValidatedTag>;
    template<class T> void _compute(T& r, const Vector<T>& a) const {
        typedef typename T::NumericType R;
        r=Ariadne::horner_evaluate(this->_model.expansion(),Ariadne::unscale(a,this->_domain))
            + convert_error<R>(this->_model.error());
    }
    ScalarTaylorFunction* _derivative(uint j) const;
    ScalarTaylorFunction* _clone() const;
    ScalarTaylorFunction* _create() const;
    VectorFunctionModelInterface<ValidatedTag>* _create_identity() const;
    VectorFunctionModelInterface<ValidatedTag>* _create_vector(uint i) const;
};

template<> struct Arithmetic<RawFloat,ScalarTaylorFunction> { typedef ScalarTaylorFunction ResultType; };
template<> struct Arithmetic<ValidatedNumber,ScalarTaylorFunction> { typedef ScalarTaylorFunction ResultType; };
template<> struct Arithmetic<EffectiveNumber,ScalarTaylorFunction> { typedef ScalarTaylorFunction ResultType; };
template<> struct Arithmetic<ScalarTaylorFunction,RawFloat> { typedef ScalarTaylorFunction ResultType; };
template<> struct Arithmetic<ScalarTaylorFunction,ValidatedNumber> { typedef ScalarTaylorFunction ResultType; };
template<> struct Arithmetic<ScalarTaylorFunction,EffectiveNumber> { typedef ScalarTaylorFunction ResultType; };
template<> struct Arithmetic<ScalarTaylorFunction,ScalarTaylorFunction> { typedef ScalarTaylorFunction ResultType; };

inline Tribool operator>(const ScalarTaylorFunction& x, const RawFloat& c) {
    UpperInterval r=x.range(); if(r.lower_raw()>c) { return true; } else if(r.upper_raw()<=c) { return false; } else { return indeterminate; } }
inline Tribool operator<(const ScalarTaylorFunction& x, const RawFloat& c) {
    UpperInterval r=x.range(); if(r.lower_raw()<c) { return true; } else if(r.upper_raw()>=c) { return false; } else { return indeterminate; } }

inline Tribool operator>(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y) { return (x-y)>0; }
inline Tribool operator<(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y) { return (x-y)<0; }

inline ScalarTaylorFunction operator+(const EffectiveScalarFunction& f1, const ScalarTaylorFunction& tf2) {
    return ScalarTaylorFunction(tf2.domain(),f1,tf2.sweeper())+tf2; }
inline ScalarTaylorFunction operator-(const EffectiveScalarFunction& f1, const ScalarTaylorFunction& tf2) {
    return ScalarTaylorFunction(tf2.domain(),f1,tf2.sweeper())-tf2; }
inline ScalarTaylorFunction operator*(const EffectiveScalarFunction& f1, const ScalarTaylorFunction& tf2) {
    return ScalarTaylorFunction(tf2.domain(),f1,tf2.sweeper())*tf2; }
inline ScalarTaylorFunction operator/(const EffectiveScalarFunction& f1, const ScalarTaylorFunction& tf2) {
    return ScalarTaylorFunction(tf2.domain(),f1,tf2.sweeper())/tf2; }
inline ScalarTaylorFunction operator+(const ScalarTaylorFunction& tf1, const EffectiveScalarFunction& f2) {
    return tf1+ScalarTaylorFunction(tf1.domain(),f2,tf1.sweeper()); }
inline ScalarTaylorFunction operator-(const ScalarTaylorFunction& tf1, const EffectiveScalarFunction& f2) {
    return tf1-ScalarTaylorFunction(tf1.domain(),f2,tf1.sweeper()); }
inline ScalarTaylorFunction operator*(const ScalarTaylorFunction& tf1, const EffectiveScalarFunction& f2) {
    return tf1*ScalarTaylorFunction(tf1.domain(),f2,tf1.sweeper()); }
inline ScalarTaylorFunction operator/(const ScalarTaylorFunction& tf1, const EffectiveScalarFunction& f2) {
    return tf1/ScalarTaylorFunction(tf1.domain(),f2,tf1.sweeper()); }

ScalarTaylorFunction& operator+=(ScalarTaylorFunction& f, const ValidatedNumber& c);
ScalarTaylorFunction& operator-=(ScalarTaylorFunction& f, const ValidatedNumber& c);
ScalarTaylorFunction& operator*=(ScalarTaylorFunction& f, const ValidatedNumber& c);
ScalarTaylorFunction& operator/=(ScalarTaylorFunction& f, const ValidatedNumber& c);
ScalarTaylorFunction operator+(const ScalarTaylorFunction& x);
ScalarTaylorFunction operator-(const ScalarTaylorFunction& x);
ScalarTaylorFunction operator+(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
ScalarTaylorFunction operator-(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
ScalarTaylorFunction operator*(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
ScalarTaylorFunction operator/(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);

template<class X> inline EnableIfNumeric<X,ScalarTaylorFunction>&
operator+=(ScalarTaylorFunction& f, const X& c) { f.model()+=ValidatedNumber(c); return f; }
template<class X> inline EnableIfNumeric<X,ScalarTaylorFunction>&
operator-=(ScalarTaylorFunction& f, const X& c) { f.model()-=ValidatedNumber(c); return f; }
template<class X> inline EnableIfNumeric<X,ScalarTaylorFunction>&
operator*=(ScalarTaylorFunction& f, const X& c) { f.model()*=ValidatedNumber(c); return f; }
template<class X> inline EnableIfNumeric<X,ScalarTaylorFunction>&
operator/=(ScalarTaylorFunction& f, const X& c) { f.model()/=ValidatedNumber(c); return f; }

template<class X> inline EnableIfNumeric<X,ScalarTaylorFunction>
operator+(const ScalarTaylorFunction& f, const X& c) { ScalarTaylorFunction r(f); r+=ValidatedNumber(c); return r; }
template<class X> inline EnableIfNumeric<X,ScalarTaylorFunction>
operator-(const ScalarTaylorFunction& f, const X& c) { ScalarTaylorFunction r(f); r+=neg(ValidatedNumber(c)); return r; }
template<class X> inline EnableIfNumeric<X,ScalarTaylorFunction>
operator*(const ScalarTaylorFunction& f, const X& c) { ScalarTaylorFunction r(f); r*=ValidatedNumber(c); return r; }
template<class X> inline EnableIfNumeric<X,ScalarTaylorFunction>
operator/(const ScalarTaylorFunction& f, const X& c) { ScalarTaylorFunction r(f); r*=(rec(ValidatedNumber(c))); return r; }
template<class X> inline EnableIfNumeric<X,ScalarTaylorFunction>
operator+(const X& c,const ScalarTaylorFunction& f) { ScalarTaylorFunction r(f); r+=ValidatedNumber(c); return r; }
template<class X> inline EnableIfNumeric<X,ScalarTaylorFunction>
operator-(const X& c,const ScalarTaylorFunction& f) { ScalarTaylorFunction r(neg(f)); r+=ValidatedNumber(c); return r; }
template<class X> inline EnableIfNumeric<X,ScalarTaylorFunction>
operator*(const X& c,const ScalarTaylorFunction& f) { ScalarTaylorFunction r(f); r*=ValidatedNumber(c); return r; }
template<class X> inline EnableIfNumeric<X,ScalarTaylorFunction>
operator/(const X& c,const ScalarTaylorFunction& f) { ScalarTaylorFunction r(rec(f)); r*=ValidatedNumber(c); return r; }

inline ScalarTaylorFunction abs(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,abs(x._model)); }
inline ScalarTaylorFunction neg(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,-x._model); }
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



inline ScalarTaylorFunction antiderivative(const ScalarTaylorFunction& x, uint k) {
    ValidatedNumber sf=rad_val(x.domain()[k]);
    return ScalarTaylorFunction(x.domain(),antiderivative(x.model(),k)*sf); }

inline ScalarTaylorFunction derivative(const ScalarTaylorFunction& x, uint k) {
    ValidatedNumber sf=rec(rad_val(x.domain()[k]));
    return ScalarTaylorFunction(x.domain(),derivative(x.model(),k)*sf); }

inline ScalarTaylorFunction embed(const ScalarTaylorFunction& tv1, const ExactInterval& dom2) {
    return ScalarTaylorFunction(join(tv1.domain(),dom2),embed(tv1.model(),1u)); }
inline ScalarTaylorFunction embed(const ScalarTaylorFunction& tv1, const ExactBox& dom2) {
    return ScalarTaylorFunction(join(tv1.domain(),dom2),embed(tv1.model(),dom2.size())); }
inline ScalarTaylorFunction embed(const ExactBox& dom1, const ScalarTaylorFunction& tv2) {
    return ScalarTaylorFunction(join(dom1,tv2.domain()),embed(dom1.size(),tv2.model())); }
inline ScalarTaylorFunction embed(const ExactBox& dom1, const ScalarTaylorFunction& tv2,const ExactBox& dom3) {
    return ScalarTaylorFunction(join(join(dom1,tv2.domain()),dom3),embed(embed(dom1.size(),tv2.model()),dom3.size())); }



Vector<ValidatedNumber> evaluate(const VectorTaylorFunction& f, const Vector<ValidatedNumber>& x);
VectorTaylorFunction partial_evaluate(const VectorTaylorFunction& f, uint k, const ExactNumber& c);
VectorTaylorFunction partial_evaluate(const VectorTaylorFunction& f, uint k, const ValidatedNumber& c);
VectorTaylorFunction embed(const VectorTaylorFunction& tv1, const ExactBox& d2);
VectorTaylorFunction embed(const VectorTaylorFunction& tv1, const ExactInterval& d2);
VectorTaylorFunction embed(const ExactBox& d1, const VectorTaylorFunction& tv2);
VectorTaylorFunction embed(const ExactBox& d1, const VectorTaylorFunction& tv2,const ExactBox& d3);
VectorTaylorFunction restrict(const VectorTaylorFunction&, const ExactBox& bx);
bool refines(const VectorTaylorFunction&, const VectorTaylorFunction&);
bool disjoint(const VectorTaylorFunction&, const VectorTaylorFunction&);
VectorTaylorFunction intersection(const VectorTaylorFunction&, const VectorTaylorFunction&);
ScalarTaylorFunction compose(const EffectiveScalarFunction&, const VectorTaylorFunction&);
ScalarTaylorFunction compose(const ValidatedScalarFunction&, const VectorTaylorFunction&);
ScalarTaylorFunction compose(const ScalarTaylorFunction&, const VectorTaylorFunction&);
VectorTaylorFunction compose(const VectorTaylorFunction&, const VectorTaylorFunction&);
VectorTaylorFunction compose(const ValidatedVectorFunction&, const VectorTaylorFunction&);
VectorTaylorFunction compose(const EffectiveVectorFunction&, const VectorTaylorFunction&);
VectorTaylorFunction antiderivative(const VectorTaylorFunction&, uint);
VectorTaylorFunction antiderivative(const VectorTaylorFunction&, uint, ExactNumber);
VectorTaylorFunction derivative(const VectorTaylorFunction&, uint);
NormType norm(const ScalarTaylorFunction& f);
NormType distance(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2);
NormType distance(const VectorTaylorFunction& f1, const EffectiveVectorFunction& f2);


ValidatedNumber unchecked_evaluate(const ScalarTaylorFunction&, const Vector<ValidatedNumber>&);
Vector<ValidatedNumber> unchecked_evaluate(const VectorTaylorFunction&, const Vector<ValidatedNumber>&);
ScalarTaylorFunction unchecked_compose(const ScalarTaylorFunction&, const VectorTaylorFunction&);
VectorTaylorFunction unchecked_compose(const VectorTaylorFunction&, const VectorTaylorFunction&);


/*! \ingroup FunctionModelSubModule
 *  \brief A Taylor function model with multivalued codomain built from the TaylorModel class.
 *
 *  See also TaylorModel, ScalarTaylorFunction, VectorTaylorFunction.
 */
class VectorTaylorFunction
    : public VectorFunctionModelMixin<VectorTaylorFunction,ValidatedTag>
{
    friend class VectorTaylorFunctionElementReference;
  public:
    typedef ExactBox DomainType;

    /*! \brief Default constructor constructs a Taylor model of order zero with no arguments and no result variables. */
    VectorTaylorFunction();

    /*! \brief Construct the zero vector function over an unspecified domain. */
    explicit VectorTaylorFunction(unsigned int result_size);

    /*! \brief Construct from a result size and a domain. */
    VectorTaylorFunction(unsigned int result_size, const ExactBox& domain, Sweeper swp);

    /*! \brief Construct a vector function all of whose components are the same. */
    VectorTaylorFunction(unsigned int result_size, const ScalarTaylorFunction& scalar_function);

    /*! \brief Construct from a domain and the expansion. */
    VectorTaylorFunction(const ExactBox& domain,
                         const Vector< Expansion<CoefficientType> >& expansion,
                         Sweeper swp);

    /*! \brief Construct from a domain, and expansion and errors. */
    VectorTaylorFunction(const ExactBox& domain,
                         const Vector< Expansion<CoefficientType> >& expansion,
                         const Vector<ErrorType>& error,
                         Sweeper swp);

    /*! \brief Construct from a domain, and expansion and errors. */
    VectorTaylorFunction(const ExactBox& domain,
                         const Vector< Expansion<RawFloat> >& expansion,
                         const Vector<RawFloat>& error,
                         Sweeper swp);

    /*! \brief Construct from a domain, and expansion and errors. */
    VectorTaylorFunction(const ExactBox& domain,
                         const Vector< Expansion<RawFloat> >& expansion,
                         Sweeper swp);

    /*! \brief Construct from a domain and the models. */
    explicit VectorTaylorFunction(const ExactBox& domain, const Vector< TaylorModel<ValidatedNumber> >& variables);

    /*! \brief Construct from a domain, a function, and a sweeper determining the accuracy. */
    VectorTaylorFunction(const ExactBox& domain,
                         const ValidatedVectorFunction& function,
                         const Sweeper& sweeper);

    /*! \brief Construct from a vector of scalar Taylor functions. */
    explicit VectorTaylorFunction(const Vector<ScalarTaylorFunction>& components);

    /*! \brief Construct from a list of scalar Taylor functions. */
    explicit VectorTaylorFunction(const List<ScalarTaylorFunction>& components);

    /*! \brief Construct from an initializer list of scalar Taylor functions. */
    VectorTaylorFunction(std::initializer_list<ScalarTaylorFunction> components);

    /*! \brief Construct from a vector expression. */
    template<class E> explicit VectorTaylorFunction(const VectorExpression<E>& ve);

    explicit VectorTaylorFunction (const VectorFunctionModel<ValidatedTag>& f);
    VectorTaylorFunction& operator=(const VectorFunctionModel<ValidatedTag>& f);

    /*! \brief Equality operator. */
    bool operator==(const VectorTaylorFunction& p) const;
    /*! \brief Inequality operator. */
    bool operator!=(const VectorTaylorFunction& p) const;

    // Data access
    /*! \brief The sweeper used to control approximation of the Taylor function. */
    Sweeper sweeper() const;
    /*! \brief Set the sweeper used to control approximation of the Taylor function. */
    void set_sweeper(Sweeper swp);
    /*! \brief The data used to define the domain of the Taylor model. */
    const ExactBox& domain() const;
    /*! \brief A rough bound for the range of the function. */
    const ExactBox codomain() const;
    /*! \brief The centre of the Taylor model. */
    const Vector<CoefficientType> centre() const;
    /*! \brief The range of the Taylor model. */
    const UpperBox range() const;
    /*! \brief The data used to define the Taylor models. */
    const Vector< TaylorModel<ValidatedNumber> >& models() const;
    /*! \brief The data used to define the centre of the Taylor models. */
    const Vector< Expansion<CoefficientType> > expansions() const;

    /*! \brief The \a i<sup>th</sup> Taylor model used to define the function. */
    const TaylorModel<ValidatedNumber>& model(uint i) const;
    /*! \brief The \a i<sup>th</sup> Taylor model used to define the function. */
    TaylorModel<ValidatedNumber>& model(uint i);

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
    Vector<ValidatedNumber> evaluate(const Vector<ValidatedNumber>& x) const;

    using VectorFunctionMixin< VectorTaylorFunction, ValidatedTag>::evaluate;

    /*! \brief Evaluate the Taylor model at the point \a x. */
    Vector<ValidatedNumber> operator()(const Vector<ValidatedNumber>& x) const;
    Vector<ApproximateNumber> operator()(const Vector<ApproximateNumber>& x) const;
    /*! \brief Evaluate the Taylor model at the point \a x. */
    Vector<ApproximateNumber> evaluate(const Vector<ApproximateNumber>& x) const;
    /*! \brief Compute an approximation to Jacobian derivative of the Taylor model sat the point \a x. */
    Matrix<ValidatedNumber> jacobian(const Vector<ValidatedNumber>& x) const;

    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    VectorTaylorFunction& sweep();
    //! \brief Remove all terms as specified by \a sweeper.
    VectorTaylorFunction& sweep(const SweeperInterface& sweeper);
    /*! \brief Set the error to zero. */
    Void clobber();

    /*! \brief The constant Taylor model with range \a r and argument domain \a d. */
    static VectorTaylorFunction constant(const ExactBox& d, const Vector<ValidatedNumber>& r, Sweeper swp);
    /*! \brief The constant Taylor model with result \a c and argument domain \a d. */
    static VectorTaylorFunction constant(const ExactBox& d, const Vector<ExactNumber>& c, Sweeper swp);
    /*! \brief The identity Taylor model on domain \a d. */
    static VectorTaylorFunction identity(const ExactBox& d, Sweeper swp);
    //! \brief Return the vector of variables in the range with values \a x over domain \a d.
    static VectorTaylorFunction projection(const ExactBox& d, uint imin, uint imax, Sweeper swp);

    /*! \brief Convert to an interval polynomial. */
    Vector< Polynomial<ValidatedNumber> > polynomials() const;
    /*! \brief The vector of roundoff/truncation errors of each component. */
    Vector< ErrorType > const errors() const;
    /*! \brief The maximum roundoff/truncation error of the components. */
    ErrorType const error() const;
    //! \brief A multivalued function equal to the model on the domain.
    ValidatedVectorFunction function() const;

    /*! \brief Truncate terms higher than \a bd. */
    VectorTaylorFunction& truncate(const MultiIndexBound& bd);
    /*! \brief Restrict to a subdomain. */
    Void restrict(const ExactBox& d);
    //! \brief Adjoin a scalar function.
    Void adjoin(const ScalarTaylorFunction& sf);

    /*! \brief Write to an output stream. */
    std::ostream& write(std::ostream& os) const;

    /*! \brief Write a full representation to an output stream. */
    std::ostream& repr(std::ostream& os) const;

    /*! \brief Inplace addition. */
    friend VectorTaylorFunction& operator+=(VectorTaylorFunction& f, const VectorTaylorFunction& g);
    /*! \brief Inplace subtraction. */
    friend VectorTaylorFunction& operator-=(VectorTaylorFunction& f, const VectorTaylorFunction& g);
    /*! \brief Inplace addition. */
    friend VectorTaylorFunction& operator+=(VectorTaylorFunction& f, const Vector<ValidatedNumber>& c);
    /*! \brief Inplace subtraction. */
    friend VectorTaylorFunction& operator-=(VectorTaylorFunction& f, const Vector<ValidatedNumber>& c);
    /*! \brief Inplace scalar multiplication. */
    friend VectorTaylorFunction& operator*=(VectorTaylorFunction& f, const ValidatedNumber& c);
    /*! \brief Inplace scalar division. */
    friend VectorTaylorFunction& operator/=(VectorTaylorFunction& f, const ValidatedNumber& c);

    /*! \brief Negation. */
    friend VectorTaylorFunction operator-(const VectorTaylorFunction& f);
    /*! \brief Addition. */
    friend VectorTaylorFunction operator+(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2);
    /*! \brief Subtraction. */
    friend VectorTaylorFunction operator-(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2);
    /*! \brief Multiplication. */
    friend VectorTaylorFunction operator*(const VectorTaylorFunction& f1, const ScalarTaylorFunction& f2);
    /*! \brief Multiplication. */
    friend VectorTaylorFunction operator*(const ScalarTaylorFunction& f1, const VectorTaylorFunction& f2);

    /*! \brief Addition of a constant. */
    friend VectorTaylorFunction operator+(const VectorTaylorFunction& f, const Vector<ValidatedNumber>& c);
    /*! \brief Subtraction of a constant. */
    friend VectorTaylorFunction operator-(const VectorTaylorFunction& f, const Vector<ValidatedNumber>& c);
    /*! \brief Multiplication by a scalar. */
    friend VectorTaylorFunction operator*(const ValidatedNumber& c, const VectorTaylorFunction& f);
    /*! \brief Multiplication by a scalar. */
    friend VectorTaylorFunction operator*(const VectorTaylorFunction& f, const ValidatedNumber& c);
    /*! \brief Division by a scalar. */
    friend VectorTaylorFunction operator/(const VectorTaylorFunction& f, const ValidatedNumber& c);
    /*! \brief Multiplication by a matrix. */
    friend VectorTaylorFunction operator*(const Matrix<ExactNumber>& A, const VectorTaylorFunction& f);
    /*! \brief Multiplication by a matrix. */
    friend VectorTaylorFunction operator*(const Matrix<ValidatedNumber>& A, const VectorTaylorFunction& f);

    //! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
    friend VectorTaylorFunction compose(const EffectiveVectorFunction& f, const VectorTaylorFunction& g);
    //! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
    friend ScalarTaylorFunction compose(const ScalarTaylorFunction& f, const VectorTaylorFunction& g);
    //! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
    friend VectorTaylorFunction compose(const VectorTaylorFunction& f, const VectorTaylorFunction& g);
    //! \brief Antiderivative of \a f with respect to variable \a k.
    friend VectorTaylorFunction antiderivative(const VectorTaylorFunction& f, uint k);
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
    friend VectorTaylorFunction restrict(const VectorTaylorFunction& f, const ExactBox& d);
    //! \brief Tests if a function \a f refines another function \a g.
    //! To be a refinement, the domain of \a f must contain the domain of \a g.
    friend bool refines(const VectorTaylorFunction& f, const VectorTaylorFunction& g);

    // For compatibility wit Vector.
    uint size() const { return this->result_size(); }
  private:
    Array< Array<ValidatedNumber> > _powers(const Vector<ValidatedNumber>&) const;
    void _compute_jacobian() const;
    void _set_argument_size(uint n);
    uint _compute_maximum_component_size() const;
    void _resize(uint rs, uint as, ushort d, ushort s);
    virtual ScalarTaylorFunction* _get(uint i) const { return new ScalarTaylorFunction(this->_domain,this->_models[i]); }
    virtual VectorTaylorFunction* _clone() const;
    virtual VectorTaylorFunction* _create() const;
    virtual VectorTaylorFunction* _create_identity() const;
    virtual ScalarTaylorFunction* _create_zero() const;
  private:
    friend class VectorFunctionMixin<VectorTaylorFunction,ValidatedTag>;
    friend class TaylorFunctionFactory;
    template<class X> void _compute(Vector<X>& r, const Vector<X>& a) const;
  private:
    /* Domain of definition. */
    ExactBox _domain;
    Vector< TaylorModel<ValidatedNumber> > _models;
};

VectorTaylorFunction operator-(const VectorTaylorFunction& f1, const EffectiveVectorFunction& f2);

// Set the value of the \a kth variable to c
VectorTaylorFunction partial_evaluate(const VectorTaylorFunction& f, uint k, const ValidatedNumber& c);
// Evaluate a scalar Taylor function on a vector.
Vector<ValidatedNumber> evaluate(const VectorTaylorFunction& f, const Vector<ValidatedNumber>& c);

// Restrict the \a kth variable to lie in the interval \a d.
VectorTaylorFunction restrict(const VectorTaylorFunction& f, uint k, const ExactInterval& d);
// Restrict to a smaller domain. REQUIRED
VectorTaylorFunction restrict(const VectorTaylorFunction& f, const ExactBox& d);
// Extend to a larger domain. REQUIRED
VectorTaylorFunction extend(const VectorTaylorFunction& f, const ExactBox& d);

// The argument size of the result is the same as that of \a e, and must be either the same as that of \a f, or one less.
VectorTaylorFunction compose(const VectorTaylorFunction& f, const VectorTaylorFunction& e);
// Compose a vector function with a Taylor function.
VectorTaylorFunction compose(const EffectiveVectorFunction& f, const VectorTaylorFunction& e);

// Substitute \a h into the \a k<sup>th</sup> argument of \a f.
VectorTaylorFunction substitute(const VectorTaylorFunction& f, uint k, const ScalarTaylorFunction& h);

// Split the domain into halves along the \a j<sup>th</sup> coordinate.
Pair<VectorTaylorFunction,VectorTaylorFunction> split(const VectorTaylorFunction& x, uint j);

VectorTaylorFunction join(const VectorTaylorFunction& f, const VectorTaylorFunction& g);
VectorTaylorFunction join(const VectorTaylorFunction& f, const ScalarTaylorFunction& g);
VectorTaylorFunction join(const ScalarTaylorFunction& f, const ScalarTaylorFunction& g);
VectorTaylorFunction join(const ScalarTaylorFunction& f, const VectorTaylorFunction& g);
VectorTaylorFunction combine(const VectorTaylorFunction& f, const VectorTaylorFunction& g);
VectorTaylorFunction combine(const VectorTaylorFunction& f, const ScalarTaylorFunction& g);
VectorTaylorFunction combine(const ScalarTaylorFunction& f, const VectorTaylorFunction& g);
VectorTaylorFunction combine(const ScalarTaylorFunction& f, const ScalarTaylorFunction& g);

NormType norm(const VectorTaylorFunction& f);

std::ostream& operator<<(std::ostream&, const VectorTaylorFunction&);

// Conversion operatations
Polynomial<ValidatedNumber> polynomial(const ScalarTaylorFunction& tfn);
Vector< Polynomial<ValidatedNumber> > polynomial(const VectorTaylorFunction& tfn);
List< Polynomial<ValidatedNumber> > polynomials(const List<ScalarTaylorFunction>& tfns);

// Sanitised output
std::ostream& operator<<(std::ostream&, const Representation<ScalarTaylorFunction>&);
std::ostream& operator<<(std::ostream&, const Representation<VectorTaylorFunction>&);
template<class F> struct ModelsRepresentation { const F* pointer; double threshold; };
template<class F> ModelsRepresentation<F> model_repr(const F& f, double swpt) { ModelsRepresentation<F> r={&f,swpt}; return r; }
std::ostream& operator<<(std::ostream&,const ModelsRepresentation<ScalarTaylorFunction>&);
std::ostream& operator<<(std::ostream&,const ModelsRepresentation< List<ScalarTaylorFunction> >&);
std::ostream& operator<<(std::ostream&,const ModelsRepresentation<VectorTaylorFunction>&);
template<class F> struct PolynomialRepresentation { const F* pointer; double threshold; List<String> names; };
template<class F> PolynomialRepresentation<F> polynomial_repr(const F& f, double swpt) { PolynomialRepresentation<F> r={&f,swpt}; return r; }
template<class F> PolynomialRepresentation<F> polynomial_repr(const F& f, double swpt, const List<String>& names) { PolynomialRepresentation<F> r={&f,swpt,names}; return r; }
std::ostream& operator<<(std::ostream&,const PolynomialRepresentation<ScalarTaylorFunction>&);
std::ostream& operator<<(std::ostream&,const PolynomialRepresentation< List<ScalarTaylorFunction> >&);
std::ostream& operator<<(std::ostream&,const PolynomialRepresentation<VectorTaylorFunction>&);


template<class E> VectorTaylorFunction::VectorTaylorFunction(const VectorExpression<E>& ve) : _domain(), _models(ve().size())
{
    if(ve().size()!=0) { this->_domain=ve().zero_element().domain(); }
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
    const TaylorModel<ValidatedNumber>& model() const { return this->_c->_models[this->_i]; }
    ErrorType error() const { return this->_c->_models[this->_i].error(); }
    void set_error(const ErrorType& e) { this->_c->_models[this->_i].set_error(e); }
    void sweep() { this->_c->_models[this->_i].sweep(); }
    template<class X> X evaluate(const Vector<X>& x) const { return this->_c->get(this->_i).evaluate(x); }
    template<class X> X operator()(const Vector<X>& x) const { return this->_c->get(this->_i).operator()(x); }
    friend std::ostream& operator<<(std::ostream& os, const VectorTaylorFunctionElementReference& t) { return os<<ScalarTaylorFunction(t); }
  private:
    VectorTaylorFunction* _c; uint _i;
};


class TaylorFunctionFactory
    : public FunctionModelFactoryInterface<ValidatedTag>
{
    Sweeper _sweeper;
  public:
    explicit TaylorFunctionFactory(Sweeper sweeper) : _sweeper(sweeper) { }
    Sweeper sweeper() const { return this->_sweeper; }
    TaylorFunctionFactory* clone() const { return new TaylorFunctionFactory(this->_sweeper); }
    Void write(OutputStream& os) const { os << "TaylorFunctionFactory( sweeper=" << this->_sweeper << " )"; }
    ScalarTaylorFunction create(const ExactBox& domain, const ValidatedScalarFunctionInterface& function) const;
    VectorTaylorFunction create(const ExactBox& domain, const ValidatedVectorFunctionInterface& function) const;
    ScalarTaylorFunction create_zero(const ExactBox& domain) const;
    ScalarTaylorFunction create_constant(const ExactBox& domain, ValidatedNumber c) const;
    ScalarTaylorFunction create_coordinate(const ExactBox& domain, uint k) const;
    VectorTaylorFunction create_zero(uint i, const ExactBox& domain) const;
    ScalarTaylorFunction create_identity(const ExactInterval& domain) const;
    VectorTaylorFunction create_identity(const ExactBox& domain) const;
  private:
    ScalarTaylorFunction* _create(const ExactBox& domain, const ValidatedScalarFunctionInterface& function) const;
    VectorTaylorFunction* _create(const ExactBox& domain, const ValidatedVectorFunctionInterface& function) const;
};



} // namespace Ariadne

#endif // ARIADNE_TAYLOR_FUNCTION_H
