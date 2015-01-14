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
#include "utility/exceptions.h"
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
typedef ScalarFunction<Validated> ValidatedScalarFunction;
template<class X> class VectorFunction;
typedef VectorFunction<Validated> ValidatedVectorFunction;

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
    for(SizeType i=0; i!=r.size(); ++i) {
        r[i]=unscale(x[i],d[i]);
    }
    return r;
}



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
    typedef ValidatedTaylorModel ModelType;
    typedef ModelType::CodomainType CodomainType;
    typedef ModelType::RangeType RangeType;
    typedef ModelType::ExpansionType ExpansionType;
    typedef ModelType::CoefficientType CoefficientType;
    typedef ModelType::ErrorType ErrorType;
    typedef ModelType::NumericType NumericType;
    typedef ModelType::NumericType::Paradigm Paradigm;
    typedef ScalarFunction<Paradigm> FunctionType;
  private:
    static const CoefficientType _zero;
    DomainType _domain;
    ModelType _model;
  public:
    //! \brief An Iterator through the (index,coefficient) pairs of the expansion expansion.
    typedef ExpansionType::Iterator Iterator;
    //! \brief A constant Iterator through the (index,coefficient) pairs of the expansion expansion.
    typedef ExpansionType::ConstIterator ConstIterator;

    //@{
    /*! \name Constructors and destructors. */
    //! \brief Default constructor.
    explicit ScalarTaylorFunction();
    //! \brief Construct a ScalarTaylorFunction over the domain \a d.
    //explicit ScalarTaylorFunction(const DomainType& d);
    explicit ScalarTaylorFunction(const DomainType& d, Sweeper swp);
    //! \brief Construct a ScalarTaylorFunction over the domain \a d, based on the scaled model \a m.
    explicit ScalarTaylorFunction(const DomainType& d, const ModelType& m);

    explicit ScalarTaylorFunction(const DomainType& d, const Expansion<CoefficientType>& p, const ErrorType& e, const Sweeper& swp);
    explicit ScalarTaylorFunction(const DomainType& d, const Expansion<RawFloat>& p, const RawFloat& e, const Sweeper& swp);

    explicit ScalarTaylorFunction(const ScalarFunctionModel<ValidatedTag>& f);
    ScalarTaylorFunction& operator=(const ScalarFunctionModel<ValidatedTag>& f);

    //! \brief Construct a ScalarTaylorFunction over the domain \a d from the function \a f.
    explicit ScalarTaylorFunction(const DomainType& d, const ValidatedScalarFunction& f, Sweeper swp);
    //@}

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    ScalarTaylorFunction& operator=(const NumericType& c) { this->_model=c; return *this; }
    //@}

    //@{
    /*! \name Named constructors. */
    //! \brief Construct a zero function over domain \a d.
    static ScalarTaylorFunction zero(const DomainType& d, Sweeper swp);
    //! \brief Construct a constant quantity in \a as independent variables.
    static ScalarTaylorFunction constant(const DomainType& d, const NumericType& c, Sweeper swp);
    //! \brief Construct the quantity \f$x_j\f$ over the domain \a d.
    static ScalarTaylorFunction coordinate(const DomainType& d, SizeType j, Sweeper swp);

    //! \brief Construct the quantity \f$c+\sum g_jx_j\f$ over the domain \a d. // DEPRECATED
    static ScalarTaylorFunction affine(const DomainType& d, const CoefficientType& c, const Vector<CoefficientType>& g, Sweeper swp);
    //! \brief Construct the quantity \f$c+\sum g_jx_j \pm e\f$ over domain \a d. // DEPRECATED
    static ScalarTaylorFunction affine(const DomainType& d, const CoefficientType& x, const Vector<CoefficientType>& g, const ErrorType& e, Sweeper swp) ;

    //! \brief Return the vector of constants with interval values \a c over domain \a d.
    static Vector<ScalarTaylorFunction> constants(const DomainType& d, const Vector<NumericType>& c, Sweeper swp);
    //! \brief Return the vector of variables with values \a x over domain \a d.
    static Vector<ScalarTaylorFunction> coordinates(const DomainType& d, Sweeper swp);
    //! \brief Return the vector of variables in the range \a imin to \a imax with values \a x over domain \a d.
    static Vector<ScalarTaylorFunction> coordinates(const DomainType& d, SizeType imin, SizeType imax, Sweeper swp);
    //@}

    //@{
    /*! \name Prototype constructors. */
    ScalarTaylorFunction create_zero() const;
    ScalarTaylorFunction create_constant(NumericType const& c) const;
    ScalarTaylorFunction create_coordinate(SizeType j) const;
    //@}

    //@{
    /*! \name Data access */
    //! \brief The domain of the quantity.
    const DomainType& domain() const { return this->_domain; }
    //! \brief The scaled expansion over a unit box with error bound.
    const ModelType& model() const { return this->_model; }
    //! \brief A reference to the scaled expansion over a unit box with error bound.
    ModelType& model() { return this->_model; }

    //! \brief An over-approximation to the range of the quantity; not necessarily tight.
    const ExactInterval codomain() const { this->_model.codomain(); }
    //! \brief The scaled expansion over a unit box.
    const ExpansionType& expansion() const { return this->_model.expansion(); }
    //! \brief The error of the expansion over the domain.
    const ErrorType& error() const { return this->_model.error(); }
    /*! \brief The accuracy parameter used to control approximation of the Taylor function. */
    Sweeper sweeper() const { return this->_model.sweeper(); }
    //! \brief A reference to the expansion.
    ExpansionType& expansion() { return this->_model.expansion(); }
    //! \brief A reference to the error of the expansion over the domain.
    ErrorType& error() { return this->_model.error(); }

    //! \brief A reference to the constant term in the expansion.
    CoefficientType& value() { return this->_model.value(); }
    //! \brief The constant term in the expansion.
    const CoefficientType& value() const { return this->_model.value(); }
    //! \brief The gradient at the centre of the domain.
    const CoefficientType gradient_value(SizeType i) const { return make_exact(this->_model.gradient(i)/this->_domain[i].radius()); }

    //! \brief A polynomial representation.
    Polynomial<ValidatedFloat> polynomial() const;
    //! \brief A multivalued function equal to the model on the domain.
    ValidatedScalarFunction function() const;

    //! \brief Set the error of the expansion.
    Void set_error(const ErrorType& ne) { this->_model.set_error(ne); }
    //! \brief Set the constant term in the expansion.
    Void set_value(const CoefficientType& c) { this->_model.set_value(c); }

    //! \brief The coefficient of the term in $x^a$.
    const CoefficientType& operator[](const MultiIndex& a) const { return this->_model[a]; }
    //! \brief A read/write reference to the coefficient of the term in $x^a$.
    CoefficientType& operator[](const MultiIndex& a) { return this->_model[a]; }

    //! \brief The number of variables in the argument of the quantity.
    SizeType argument_size() const { return this->_model.argument_size(); }
    //! \brief The maximum degree of terms in the expansion expansion.
    DegreeType degree() const { return this->_model.degree(); }
    //! \brief The number of nonzero terms in the expansion expansion.
    SizeType number_of_nonzeros() const { return this->_model.number_of_nonzeros(); }
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    Bool operator==(const ScalarTaylorFunction& tv) const;
    //! \brief Inequality operator.
    Bool operator!=(const ScalarTaylorFunction& tv) const { return !(*this==tv); }
    //@}

    //@{
    /*! \name Function operations. */
    //! \brief An over-approximation to the range of the function.
    UpperInterval range() const { return this->_model.range(); }
    //! \brief Evaluate the function at the point \a x.
    ValidatedNumber operator()(const Vector<ValidatedNumber>& x) const;
    ValidatedNumber operator()(const Vector<ExactNumber>& x) const;
    ApproximateNumber operator()(const Vector<ApproximateNumber>& x) const;

    /*! \brief Compute an approximation to gradient derivative of the function at the point \a x. */
    Vector<NumericType> gradient(const Vector<NumericType>& x) const;
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
    //! \copydoc TaylorModel<ValidatedFloat>::set_sweeper()
    Void set_sweeper(const Sweeper& swp) { this->_model.set_sweeper(swp); }
    //@}

    //@{
    /*! \name Non-arithmetic operations. */
    //! \brief Restrict to a subdomain.
    Void restrict(const DomainType& d);


    //@}

    /*! \name Arithmetic operations. */

/*
    //! \brief Inplace addition of another variable.
    friend ScalarTaylorFunction& operator+=(ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
    //! \brief Inplace subtraction of another variable.
    friend ScalarTaylorFunction& operator-=(ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
    //! \brief Inplace addition of a product of two variables.
    friend ScalarTaylorFunction& operator+=(ScalarTaylorFunction& x, const Product<ScalarTaylorFunction,ScalarTaylorFunction>& y);
    //! \brief Inplace addition of an interval constant.
    friend ScalarTaylorFunction& operator+=(ScalarTaylorFunction& x, const NumericType& c);
    //! \brief Inplace subtraction of an interval constant.
    friend ScalarTaylorFunction& operator-=(ScalarTaylorFunction& x, const NumericType& c);
    //! \brief Inplace multiplication by an approximate scalar.
    friend ScalarTaylorFunction& operator*=(ScalarTaylorFunction& x, const NumericType& c);
    //! \brief Inplace division by an approximate scalar.
    friend ScalarTaylorFunction& operator/=(ScalarTaylorFunction& x, const NumericType& c);

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
    friend ScalarTaylorFunction pow(const ScalarTaylorFunction& x, Int n);
    //! \brief Square root.
    friend ScalarTaylorFunction sqrt(const ScalarTaylorFunction& x);
    //! \brief SizeTypeural exponent.
    friend ScalarTaylorFunction exp(const ScalarTaylorFunction& x);
    //! \brief SizeTypeural logarithm.
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
*/
    //@}

    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    OutputStream& write(OutputStream& os) const;
    /*! \brief Write a full representation to an output stream. */
    OutputStream& repr(OutputStream& os) const;
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const ScalarTaylorFunction& x);
    //@}

  public:
    Void clobber() { this->_model.clobber(); }
  private:
    friend class TaylorFunctionFactory;
    friend class ScalarFunctionMixin<ScalarTaylorFunction, ValidatedTag>;
    friend class ScalarFunctionModelMixin<ScalarTaylorFunction, ValidatedTag>;
    template<class T> Void _compute(T& r, const Vector<T>& a) const {
        typedef typename T::NumericType R;
        r=Ariadne::horner_evaluate(this->_model.expansion(),Ariadne::unscale(a,this->_domain))
            + convert_error<R>(this->_model.error());
    }
    ScalarTaylorFunction* _derivative(SizeType j) const;
    ScalarTaylorFunction* _clone() const;
    ScalarTaylorFunction* _create() const;
    VectorFunctionModelInterface<ValidatedTag>* _create_identity() const;
    VectorFunctionModelInterface<ValidatedTag>* _create_vector(SizeType i) const;
};

//! \brief Restrict to a subdomain.
ScalarTaylorFunction restriction(const ScalarTaylorFunction& x, const ExactBox& d);
//! \brief Extend over a larger domain. Only possible if the larger domain is only larger where the smaller domain is a singleton.
//! The extension is performed keeping \a x constant over the new coordinates. // DEPRECATED
ScalarTaylorFunction extension(const ScalarTaylorFunction& x, const ExactBox& d);

//! \brief Test if the quantity is a better approximation than \a t throughout the domain.
Bool refines(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2);
//! \brief Test if the function models are inconsistent with representing the same exact function.
Bool inconsistent(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2);
//! \brief Compute an over-approximation to the common refinement of \a x1 and \a x2.
ScalarTaylorFunction refinement(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2);


// Normed algebra
NormType norm(const ScalarTaylorFunction& f);

// Differential algebra
ScalarTaylorFunction antiderivative(const ScalarTaylorFunction& x, SizeType k, ExactFloat c);
ScalarTaylorFunction antiderivative(const ScalarTaylorFunction& x, SizeType k);
ScalarTaylorFunction derivative(const ScalarTaylorFunction& x, SizeType k);

// Function algebra
ScalarTaylorFunction embed(const ExactBox& dom1, const ScalarTaylorFunction& tv2,const ExactBox& dom3);
// Set the value of the \a kth variable to c
ScalarTaylorFunction partial_evaluate(const ScalarTaylorFunction& f, SizeType k, const ScalarTaylorFunction::NumericType& c);
// Evaluate a scalar Taylor function on a vector.
ScalarTaylorFunction::NumericType unchecked_evaluate(const ScalarTaylorFunction&, const Vector<ScalarTaylorFunction::NumericType>&);

// Compose with an function.
ScalarTaylorFunction compose(const ValidatedScalarFunction& x, const VectorTaylorFunction& y);
ScalarTaylorFunction unchecked_compose(const ScalarTaylorFunction&, const VectorTaylorFunction&);

// Split the variable over two domains, subdividing along the independent variable j.
Pair<ScalarTaylorFunction,ScalarTaylorFunction> split(const ScalarTaylorFunction& x, SizeType j);


ScalarTaylorFunction& operator+=(ScalarTaylorFunction& f, const ScalarTaylorFunction::NumericType& c);
ScalarTaylorFunction& operator-=(ScalarTaylorFunction& f, const ScalarTaylorFunction::NumericType& c);
ScalarTaylorFunction& operator*=(ScalarTaylorFunction& f, const ScalarTaylorFunction::NumericType& c);
ScalarTaylorFunction& operator/=(ScalarTaylorFunction& f, const ScalarTaylorFunction::NumericType& c);
ScalarTaylorFunction& operator+=(ScalarTaylorFunction& f1, const ScalarTaylorFunction& f2);
ScalarTaylorFunction& operator-=(ScalarTaylorFunction& f1, const ScalarTaylorFunction& f2);

ScalarTaylorFunction operator+(const ScalarTaylorFunction& f);
ScalarTaylorFunction operator-(const ScalarTaylorFunction& f);
ScalarTaylorFunction operator+(const ScalarTaylorFunction& f1, const ScalarTaylorFunction& f2);
ScalarTaylorFunction operator-(const ScalarTaylorFunction& f1, const ScalarTaylorFunction& f2);
ScalarTaylorFunction operator*(const ScalarTaylorFunction& f1, const ScalarTaylorFunction& f2);
ScalarTaylorFunction operator/(const ScalarTaylorFunction& f1, const ScalarTaylorFunction& f2);
ScalarTaylorFunction operator+(const ScalarTaylorFunction& f, const ScalarTaylorFunction::NumericType& c);
ScalarTaylorFunction operator-(const ScalarTaylorFunction& f, const ScalarTaylorFunction::NumericType& c);
ScalarTaylorFunction operator*(const ScalarTaylorFunction& f, const ScalarTaylorFunction::NumericType& c);
ScalarTaylorFunction operator/(const ScalarTaylorFunction& f, const ScalarTaylorFunction::NumericType& c);
ScalarTaylorFunction operator+(const ScalarTaylorFunction::NumericType& c, const ScalarTaylorFunction& f);
ScalarTaylorFunction operator-(const ScalarTaylorFunction::NumericType& c, const ScalarTaylorFunction& f);
ScalarTaylorFunction operator*(const ScalarTaylorFunction::NumericType& c, const ScalarTaylorFunction& f);
ScalarTaylorFunction operator/(const ScalarTaylorFunction::NumericType& c, const ScalarTaylorFunction& f);

ScalarTaylorFunction operator+(const ScalarTaylorFunction::FunctionType& f1, const ScalarTaylorFunction& tf2);
ScalarTaylorFunction operator-(const ScalarTaylorFunction::FunctionType& f1, const ScalarTaylorFunction& tf2);
ScalarTaylorFunction operator*(const ScalarTaylorFunction::FunctionType& f1, const ScalarTaylorFunction& tf2);
ScalarTaylorFunction operator/(const ScalarTaylorFunction::FunctionType& f1, const ScalarTaylorFunction& tf2);
ScalarTaylorFunction operator+(const ScalarTaylorFunction& tf1, const ScalarTaylorFunction::FunctionType& f2);
ScalarTaylorFunction operator-(const ScalarTaylorFunction& tf1, const ScalarTaylorFunction::FunctionType& f2);
ScalarTaylorFunction operator*(const ScalarTaylorFunction& tf1, const ScalarTaylorFunction::FunctionType& f2);
ScalarTaylorFunction operator/(const ScalarTaylorFunction& tf1, const ScalarTaylorFunction::FunctionType& f2);

ScalarTaylorFunction abs(const ScalarTaylorFunction& x);
ScalarTaylorFunction neg(const ScalarTaylorFunction& x);
ScalarTaylorFunction rec(const ScalarTaylorFunction& x);
ScalarTaylorFunction sqr(const ScalarTaylorFunction& x);
ScalarTaylorFunction pow(const ScalarTaylorFunction& x, Int n);
ScalarTaylorFunction sqrt(const ScalarTaylorFunction& x);
ScalarTaylorFunction exp(const ScalarTaylorFunction& x);
ScalarTaylorFunction log(const ScalarTaylorFunction& x);
ScalarTaylorFunction sin(const ScalarTaylorFunction& x);
ScalarTaylorFunction cos(const ScalarTaylorFunction& x);
ScalarTaylorFunction tan(const ScalarTaylorFunction& x);
ScalarTaylorFunction asin(const ScalarTaylorFunction& x);
ScalarTaylorFunction acos(const ScalarTaylorFunction& x);
ScalarTaylorFunction atan(const ScalarTaylorFunction& x);

// Remove the error term
ScalarTaylorFunction midpoint(const ScalarTaylorFunction& x);






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
    typedef ValidatedTaylorModel ModelType;
    typedef Box<ModelType::CodomainType> CodomainType;
    typedef Box<ModelType::RangeType> RangeType;
    typedef ModelType::ExpansionType ExpansionType;
    typedef ModelType::CoefficientType CoefficientType;
    typedef ModelType::ErrorType ErrorType;
    typedef ModelType::NumericType NumericType;

    /*! \brief Default constructor constructs a Taylor model of order zero with no arguments and no result variables. */
    VectorTaylorFunction();

    /*! \brief Construct the zero vector function over an unspecified domain. */
    explicit VectorTaylorFunction(SizeType result_size);

    /*! \brief Construct from a result size and a domain. */
    VectorTaylorFunction(SizeType result_size, const ExactBox& domain, Sweeper swp);

    /*! \brief Construct a vector function all of whose components are the same. */
    VectorTaylorFunction(SizeType result_size, const ScalarTaylorFunction& scalar_function);

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
    explicit VectorTaylorFunction(const ExactBox& domain, const Vector< ModelType >& variables);

    /*! \brief Construct from a domain, a function, and a sweeper determining the accuracy. */
    VectorTaylorFunction(const ExactBox& domain,
                         const ValidatedVectorFunction& function,
                         const Sweeper& sweeper);

    /*! \brief Construct from a vector of scalar Taylor functions. */
    explicit VectorTaylorFunction(const Vector<ScalarTaylorFunction>& components);

    /*! \brief Construct from a list of scalar Taylor functions. */
    explicit VectorTaylorFunction(const List<ScalarTaylorFunction>& components);

    /*! \brief Construct from an initializer list of scalar Taylor functions. */
    VectorTaylorFunction(InitializerList<ScalarTaylorFunction> components);

    /*! \brief Construct from a vector expression. */
    template<class E> explicit VectorTaylorFunction(const VectorExpression<E>& ve);

    explicit VectorTaylorFunction (const VectorFunctionModel<ValidatedTag>& f);
    VectorTaylorFunction& operator=(const VectorFunctionModel<ValidatedTag>& f);

    /*! \brief Equality operator. */
    Bool operator==(const VectorTaylorFunction& p) const;
    /*! \brief Inequality operator. */
    Bool operator!=(const VectorTaylorFunction& p) const;

    // Data access
    /*! \brief The sweeper used to control approximation of the Taylor function. */
    Sweeper sweeper() const;
    /*! \brief Set the sweeper used to control approximation of the Taylor function. */
    Void set_sweeper(Sweeper swp);
    /*! \brief The data used to define the domain of the Taylor model. */
    const ExactBox& domain() const;
    /*! \brief A rough bound for the range of the function. */
    const ExactBox codomain() const;
    /*! \brief The centre of the Taylor model. */
    const Vector<CoefficientType> centre() const;
    /*! \brief The range of the Taylor model. */
    const UpperBox range() const;
    /*! \brief The data used to define the Taylor models. */
    const Vector<ModelType>& models() const;
    Vector<ModelType>& models();
    /*! \brief The data used to define the centre of the Taylor models. */
    const Vector<Expansion<ExactFloat>> expansions() const;

    /*! \brief The \a i<sup>th</sup> Taylor model used to define the function. */
    const ModelType& model(SizeType i) const;
    /*! \brief The \a i<sup>th</sup> Taylor model used to define the function. */
    ModelType& model(SizeType i);

    /*! \brief The size of the argument. */
    SizeType argument_size() const;
    /*! \brief The size of the result. */
    SizeType result_size() const;

    // For compatibility wit Vector.
    SizeType size() const { return this->result_size(); }
    /*! \brief Get the \a ith Taylor variable */
    ScalarTaylorFunction get(SizeType i) const;
    /*! \brief Set the \a ith Taylor variable */
    Void set(SizeType i, const ScalarTaylorFunction& te);
    /*! \brief The \a ith Taylor variable */
    ScalarTaylorFunction operator[](SizeType i) const;
    /*! \brief The \a ith Taylor variable */
    VectorTaylorFunctionElementReference operator[](SizeType i);


    /*! \brief Evaluate the Taylor model at the point \a x. */
    Vector<ValidatedNumber> operator()(const Vector<ValidatedNumber>& x) const;
    Vector<ApproximateNumber> operator()(const Vector<ApproximateNumber>& x) const;
    Vector<ValidatedNumber> operator()(const Vector<ExactNumber>& x) const;
    /*! \brief Compute an approximation to Jacobian derivative of the Taylor model sat the point \a x. */
    Matrix<NumericType> jacobian(const Vector<NumericType>& x) const;

    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    VectorTaylorFunction& sweep();
    //! \brief Remove all terms as specified by \a sweeper.
    VectorTaylorFunction& sweep(const SweeperInterface& sweeper);
    /*! \brief Set the error to zero. */
    Void clobber();

    /*! \brief The constant Taylor model with range \a r and argument domain \a d. */
    static VectorTaylorFunction constant(const ExactBox& d, const Vector<NumericType>& r, Sweeper swp);
    /*! \brief The identity Taylor model on domain \a d. */
    static VectorTaylorFunction identity(const ExactBox& d, Sweeper swp);
    //! \brief Return the vector of variables in the range with values \a x over domain \a d.
    static VectorTaylorFunction projection(const ExactBox& d, SizeType imin, SizeType imax, Sweeper swp);

    /*! \brief Convert to an interval polynomial. */
    Vector<Polynomial<ValidatedFloat>> polynomials() const;
    /*! \brief The vector of roundoff/truncation errors of each component. */
    Vector<ErrorType> const errors() const;
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
    OutputStream& write(OutputStream& os) const;

    /*! \brief Write a full representation to an output stream. */
    OutputStream& repr(OutputStream& os) const;

  private:
    Array< Array<ValidatedNumber> > _powers(const Vector<ValidatedNumber>&) const;
    Void _compute_jacobian() const;
    Void _set_argument_size(SizeType n);
    SizeType _compute_maximum_component_size() const;
    Void _resize(SizeType rs, SizeType as, ushort d, ushort s);
    virtual ScalarTaylorFunction* _get(SizeType i) const { return new ScalarTaylorFunction(this->_domain,this->_models[i]); }
    virtual VectorTaylorFunction* _clone() const;
    virtual VectorTaylorFunction* _create() const;
    virtual VectorTaylorFunction* _create_identity() const;
    virtual ScalarTaylorFunction* _create_zero() const;
  private:
    friend class VectorFunctionMixin<VectorTaylorFunction,ValidatedTag>;
    friend class TaylorFunctionFactory;
    template<class X> Void _compute(Vector<X>& r, const Vector<X>& a) const;
  private:
    /* Domain of definition. */
    ExactBox _domain;
    Vector< ModelType > _models;
};


    /*! \brief Inplace addition. */
    VectorTaylorFunction& operator+=(VectorTaylorFunction& f, const VectorTaylorFunction& g);
    /*! \brief Inplace subtraction. */
    VectorTaylorFunction& operator-=(VectorTaylorFunction& f, const VectorTaylorFunction& g);
    /*! \brief Inplace addition. */
    VectorTaylorFunction& operator+=(VectorTaylorFunction& f, const Vector<VectorTaylorFunction::NumericType>& c);
    /*! \brief Inplace subtraction. */
    VectorTaylorFunction& operator-=(VectorTaylorFunction& f, const Vector<VectorTaylorFunction::NumericType>& c);
    /*! \brief Inplace scalar multiplication. */
    VectorTaylorFunction& operator*=(VectorTaylorFunction& f, const VectorTaylorFunction::NumericType& c);
    /*! \brief Inplace scalar division. */
    VectorTaylorFunction& operator/=(VectorTaylorFunction& f, const VectorTaylorFunction::NumericType& c);

    /*! \brief Negation. */
    VectorTaylorFunction operator-(const VectorTaylorFunction& f);
    /*! \brief Addition. */
    VectorTaylorFunction operator+(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2);
    /*! \brief Subtraction. */
    VectorTaylorFunction operator-(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2);
    /*! \brief Multiplication. */
    VectorTaylorFunction operator*(const ScalarTaylorFunction& f1, const VectorTaylorFunction& f2);
    /*! \brief Multiplication. */
    VectorTaylorFunction operator*(const VectorTaylorFunction& f1, const ScalarTaylorFunction& f2);
    /*! \brief Division. */
    VectorTaylorFunction operator/(const VectorTaylorFunction& f1, const ScalarTaylorFunction& f2);

    /*! \brief Addition of a constant. */
    VectorTaylorFunction operator+(const VectorTaylorFunction& f, const Vector<VectorTaylorFunction::NumericType>& c);
    /*! \brief Subtraction of a constant. */
    VectorTaylorFunction operator-(const VectorTaylorFunction& f, const Vector<VectorTaylorFunction::NumericType>& c);
    /*! \brief Multiplication by a scalar. */
    VectorTaylorFunction operator*(const VectorTaylorFunction::NumericType& c, const VectorTaylorFunction& f);
    /*! \brief Multiplication by a scalar. */
    VectorTaylorFunction operator*(const VectorTaylorFunction& f, const VectorTaylorFunction::NumericType& c);
    /*! \brief Division by a scalar. */
    VectorTaylorFunction operator/(const VectorTaylorFunction& f, const VectorTaylorFunction::NumericType& c);
    /*! \brief Multiplication by a matrix. */
    VectorTaylorFunction operator*(const Matrix<ExactNumber>& A, const VectorTaylorFunction& f);
    /*! \brief Multiplication by a matrix. */
    VectorTaylorFunction operator*(const Matrix<VectorTaylorFunction::NumericType>& A, const VectorTaylorFunction& f);

    VectorTaylorFunction operator+(const ValidatedVectorFunction& f1, const VectorTaylorFunction& tf2);
    VectorTaylorFunction operator-(const ValidatedVectorFunction& f1, const VectorTaylorFunction& tf2);
    VectorTaylorFunction operator*(const ValidatedScalarFunction& f1, const VectorTaylorFunction& tf2);
    VectorTaylorFunction operator*(const ValidatedVectorFunction& f1, const ScalarTaylorFunction& tf2);
    VectorTaylorFunction operator/(const ValidatedVectorFunction& f1, const ScalarTaylorFunction& tf2);
    VectorTaylorFunction operator+(const VectorTaylorFunction& tf1, const ValidatedVectorFunction& f2);
    VectorTaylorFunction operator-(const VectorTaylorFunction& tf1, const ValidatedVectorFunction& f2);
    VectorTaylorFunction operator*(const ScalarTaylorFunction& tf1, const ValidatedVectorFunction& f2);
    VectorTaylorFunction operator*(const VectorTaylorFunction& tf1, const ValidatedScalarFunction& f2);
    VectorTaylorFunction operator/(const VectorTaylorFunction& tf1, const ValidatedScalarFunction& f2);

    //! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
    ScalarTaylorFunction compose(const ValidatedScalarFunction& f, const VectorTaylorFunction& g);
    //! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
    VectorTaylorFunction compose(const ValidatedVectorFunction& f, const VectorTaylorFunction& g);
    //! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
    ScalarTaylorFunction compose(const ScalarTaylorFunction& f, const VectorTaylorFunction& g);
    //! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
    VectorTaylorFunction compose(const VectorTaylorFunction& f, const VectorTaylorFunction& g);

    //! \brief Weak derivative of \a f with respect to variable \a k.
    VectorTaylorFunction derivative(const VectorTaylorFunction& f, SizeType k);
    //! \brief Antiderivative of \a f with respect to variable \a k.
    VectorTaylorFunction antiderivative(const VectorTaylorFunction& f, SizeType k);
    //! \brief Antiderivative of \a f with respect to variable \a k, taking value \c 0 when \a x[k]=c.
    VectorTaylorFunction antiderivative(const VectorTaylorFunction& f, SizeType k, ExactFloat c);

    NormType norm(const VectorTaylorFunction& f);
    NormType distance(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2);
    NormType distance(const VectorTaylorFunction& f1, const ValidatedVectorFunction& f2);

    //! \brief Restrict the function \a f to a subdomain \a d.
    VectorTaylorFunction restriction(const VectorTaylorFunction& f, const ExactBox& d);
    //! \brief Restrict the function \a f to a larger domain \a d.
    VectorTaylorFunction extension(const VectorTaylorFunction& f, const ExactBox& d);

    VectorTaylorFunction embed(const ExactBox& d1, const VectorTaylorFunction& tv2,const ExactBox& d3);

    //! \brief Tests if a function \a f refines another function \a g.
    //! To be a refinement, the domain of \a f must contain the domain of \a g.
    Bool refines(const VectorTaylorFunction& f, const VectorTaylorFunction& g);
    Bool inconsistent(const VectorTaylorFunction&, const VectorTaylorFunction&);
    VectorTaylorFunction refinement(const VectorTaylorFunction&, const VectorTaylorFunction&);

    //! \brief Compute the function \f$(f \oplus g)(x)=(f(x),g(x))\f$.
    VectorTaylorFunction join(const VectorTaylorFunction& f, const VectorTaylorFunction& g);
    VectorTaylorFunction join(const VectorTaylorFunction& f, const ScalarTaylorFunction& g);
    VectorTaylorFunction join(const ScalarTaylorFunction& f, const ScalarTaylorFunction& g);
    VectorTaylorFunction join(const ScalarTaylorFunction& f, const VectorTaylorFunction& g);
    //! \brief Compute the function \f$(f\otimes g)(x,y)=(f(x),g(y))\f$.
    VectorTaylorFunction combine(const VectorTaylorFunction& f, const VectorTaylorFunction& g);
    VectorTaylorFunction combine(const VectorTaylorFunction& f, const ScalarTaylorFunction& g);
    VectorTaylorFunction combine(const ScalarTaylorFunction& f, const VectorTaylorFunction& g);
    VectorTaylorFunction combine(const ScalarTaylorFunction& f, const ScalarTaylorFunction& g);

    VectorTaylorFunction partial_evaluate(const VectorTaylorFunction& f, SizeType k, const VectorTaylorFunction::NumericType& c);

    Vector<VectorTaylorFunction::NumericType> unchecked_evaluate(const VectorTaylorFunction&, const Vector<VectorTaylorFunction::NumericType>&);
    ScalarTaylorFunction unchecked_compose(const ScalarTaylorFunction&, const VectorTaylorFunction&);
    VectorTaylorFunction unchecked_compose(const VectorTaylorFunction&, const VectorTaylorFunction&);

    // Split the domain into halves along the \a j<sup>th</sup> coordinate.
    Pair<VectorTaylorFunction,VectorTaylorFunction> split(const VectorTaylorFunction& x, SizeType j);

    OutputStream& operator<<(OutputStream&, const VectorTaylorFunction&);




// Conversion operatations
Polynomial<ValidatedFloat> polynomial(const ScalarTaylorFunction& tfn);
Vector< Polynomial<ValidatedFloat> > polynomials(const VectorTaylorFunction& tfn);


// Sanitised output
OutputStream& operator<<(OutputStream&, const Representation<ScalarTaylorFunction>&);
OutputStream& operator<<(OutputStream&, const Representation<VectorTaylorFunction>&);
template<class F> struct ModelRepresentation { const F* pointer; double threshold; };
template<class F> ModelRepresentation<F> model_representation(const F& f, double swpt) {
    ModelRepresentation<F> r={&f,swpt}; return r; }
OutputStream& operator<<(OutputStream&,const ModelRepresentation<ScalarTaylorFunction>&);
OutputStream& operator<<(OutputStream&,const ModelRepresentation<VectorTaylorFunction>&);
template<class F> struct PolynomialRepresentation { const F* pointer; double threshold; List<String> names; };
template<class F> PolynomialRepresentation<F> polynomial_representation(const F& f, double swpt) {
    PolynomialRepresentation<F> r={&f,swpt}; return r; }
template<class F> PolynomialRepresentation<F> polynomial_representation(const F& f, double swpt, const List<String>& names) {
    PolynomialRepresentation<F> r={&f,swpt,names}; return r; }
OutputStream& operator<<(OutputStream&,const PolynomialRepresentation<ScalarTaylorFunction>&);
OutputStream& operator<<(OutputStream&,const PolynomialRepresentation<VectorTaylorFunction>&);


template<class E> VectorTaylorFunction::VectorTaylorFunction(const VectorExpression<E>& ve)
    : _domain(), _models(ve().size(),ve().zero_element().model())
{
    if(ve().size()!=0) { this->_domain=ve().zero_element().domain(); }
    for(SizeType i=0; i!=ve().size(); ++i) { this->set(i,ve()[i]); }
}

class VectorTaylorFunctionElementReference
{
    friend class ScalarTaylorFunction;
    friend class VectorTaylorFunction;
    typedef ValidatedTaylorModel ModelType;
 public:
    VectorTaylorFunctionElementReference(VectorTaylorFunction& c, SizeType i) : _c(&c), _i(i) { }
    operator ScalarTaylorFunction () const { return this->_c->get(this->_i); }
    Void operator=(const VectorTaylorFunctionElementReference& x) { this->_c->set(this->_i,x._c->get(x._i)); }
    Void operator=(const ScalarTaylorFunction& x) { this->_c->set(this->_i,x); }
    ScalarTaylorFunction element() const { return this->_c->get(this->_i); }
    ExactBox const& domain() const { return this->_c->domain(); }
    const ModelType& model() const { return this->_c->_models[this->_i]; }
    ErrorType error() const { return this->_c->_models[this->_i].error(); }
    Void set_error(const ErrorType& e) { this->_c->_models[this->_i].set_error(e); }
    Void sweep() { this->_c->_models[this->_i].sweep(); }
    template<class X> X operator()(const Vector<X>& x) const { return this->_c->get(this->_i).operator()(x); }
    friend OutputStream& operator<<(OutputStream& os, const VectorTaylorFunctionElementReference& t) { return os<<ScalarTaylorFunction(t); }
  private:
    VectorTaylorFunction* _c; SizeType _i;
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
    ScalarTaylorFunction create_coordinate(const ExactBox& domain, SizeType k) const;
    VectorTaylorFunction create_zero(SizeType i, const ExactBox& domain) const;
    ScalarTaylorFunction create_identity(const ExactInterval& domain) const;
    VectorTaylorFunction create_identity(const ExactBox& domain) const;
  private:
    ScalarTaylorFunction* _create(const ExactBox& domain, const ValidatedScalarFunctionInterface& function) const;
    VectorTaylorFunction* _create(const ExactBox& domain, const ValidatedVectorFunctionInterface& function) const;
};



} // namespace Ariadne

#endif // ARIADNE_TAYLOR_FUNCTION_H
