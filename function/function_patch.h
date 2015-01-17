/***************************************************************************
 *            function_patch.h
 *
 *  Copyright 2008-15  Pieter Collins
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

/*! \file function_patch.h
 *  \brief Over-approximations of functions on box domains.
 */

#ifndef ARIADNE_FUNCTION_PATCH_H
#define ARIADNE_FUNCTION_PATCH_H

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

template<class T> using NumericType = typename T::NumericType;
template<class T> using FunctionType = typename T::FunctionType;
template<class T> using ScalarFunctionType = typename T::ScalarFunctionType;
template<class T> using VectorFunctionType = typename T::VectorFunctionType;

template<class M> class FunctionPatch;
template<class M> using ScalarFunctionPatch = FunctionPatch<M>;
template<class M> class VectorFunctionPatch;
template<class M> class VectorFunctionPatchElementReference;

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
template<class M> class FunctionPatch
    : public ScalarFunctionModelMixin<FunctionPatch<M>, typename M::Paradigm>
{
    typedef typename M::Paradigm P;
  public:
    typedef ExactBox DomainType;
    typedef M ModelType;
    typedef typename ModelType::CodomainType CodomainType;
    typedef typename ModelType::RangeType RangeType;
    typedef typename ModelType::ExpansionType ExpansionType;
    typedef typename ModelType::CoefficientType CoefficientType;
    typedef typename ModelType::ErrorType ErrorType;
    typedef typename ModelType::NumericType NumericType;
    typedef typename ModelType::Paradigm Paradigm;
    typedef ScalarFunction<Paradigm> FunctionType;
  private:
    static const CoefficientType _zero;
    DomainType _domain;
    ModelType _model;
  public:

    //@{
    /*! \name Constructors and destructors. */
    //! \brief Default constructor.
    explicit FunctionPatch();
    //! \brief Construct a FunctionPatch<M> over the domain \a d.
    //explicit FunctionPatch(const DomainType& d);
    explicit FunctionPatch(const DomainType& d, Sweeper swp);
    //! \brief Construct a FunctionPatch<M> over the domain \a d, based on the scaled model \a m.
    explicit FunctionPatch(const DomainType& d, const ModelType& m);

    explicit FunctionPatch(const ExactBox& d, const Expansion<ExactFloat>& p, const ErrorFloat& e, const Sweeper& swp);
    explicit FunctionPatch(const ExactBox& d, const Expansion<RawFloat>& p, const RawFloat& e, const Sweeper& swp);

    explicit FunctionPatch(const ScalarFunctionModel<ValidatedTag>& f);
    FunctionPatch& operator=(const ScalarFunctionModel<ValidatedTag>& f);

    //! \brief Construct a FunctionPatch over the domain \a d from the function \a f.
    explicit FunctionPatch(const DomainType& d, const ScalarFunctionType<M>& f, Sweeper swp);
    //@}

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    FunctionPatch<M>& operator=(const NumericType& c) { this->_model=c; return *this; }
    //@}

    //@{
    /*! \name Named constructors. */
    //! \brief Construct a zero function over domain \a d.
    static FunctionPatch<M> zero(const DomainType& d, Sweeper swp);
    //! \brief Construct a constant quantity in \a as independent variables.
    static FunctionPatch<M> constant(const DomainType& d, const NumericType& c, Sweeper swp);
    //! \brief Construct the quantity \f$x_j\f$ over the domain \a d.
    static FunctionPatch<M> coordinate(const DomainType& d, SizeType j, Sweeper swp);

    //! \brief Construct the quantity \f$c+\sum g_jx_j\f$ over the domain \a d. // DEPRECATED
    static FunctionPatch<M> affine(const DomainType& d, const CoefficientType& c, const Vector<CoefficientType>& g, Sweeper swp);
    //! \brief Construct the quantity \f$c+\sum g_jx_j \pm e\f$ over domain \a d. // DEPRECATED
    static FunctionPatch<M> affine(const DomainType& d, const CoefficientType& x, const Vector<CoefficientType>& g, const ErrorType& e, Sweeper swp) ;

    //! \brief Return the vector of constants with interval values \a c over domain \a d.
    static Vector<FunctionPatch<M>> constants(const DomainType& d, const Vector<NumericType>& c, Sweeper swp);
    //! \brief Return the vector of variables with values \a x over domain \a d.
    static Vector<FunctionPatch<M>> coordinates(const DomainType& d, Sweeper swp);
    //! \brief Return the vector of variables in the range \a imin to \a imax with values \a x over domain \a d.
    static Vector<FunctionPatch<M>> coordinates(const DomainType& d, SizeType imin, SizeType imax, Sweeper swp);
    //@}

    //@{
    /*! \name Prototype constructors. */
    FunctionPatch<M> create_zero() const;
    FunctionPatch<M> create_constant(NumericType const& c) const;
    FunctionPatch<M> create_coordinate(SizeType j) const;
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

    //! \brief The constant term in the expansion.
    const CoefficientType& value() const { return this->_model.value(); }
    //! \brief The gradient at the centre of the domain.
    const CoefficientType gradient_value(SizeType i) const { return make_exact(this->_model.gradient_value(i)/this->_domain[i].radius()); }

    //! \brief A polynomial representation.
    Polynomial<ValidatedFloat> polynomial() const;
    //! \brief A multivalued function equal to the model on the domain.
    ScalarFunctionType<M> function() const;

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
    Bool operator==(const FunctionPatch<M>& tv) const;
    //! \brief Inequality operator.
    Bool operator!=(const FunctionPatch<M>& tv) const { return !(*this==tv); }
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
    Covector<NumericType> gradient(const Vector<NumericType>& x) const;
    //@}

    //@{
    /*! \name Simplification operations. */
   //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    FunctionPatch<M>& sweep() { this->_model.sweep(); return *this; }
    //! \brief Remove all terms whose degree is higher than \a deg or
    //! whose coefficient has magnitude less than \a eps.
    FunctionPatch<M>& sweep(const SweeperInterface& swp) { this->_model.sweep(swp); return *this; }
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
    //@}

    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    OutputStream& write(OutputStream& os) const;
    /*! \brief Write a full representation to an output stream. */
    OutputStream& repr(OutputStream& os) const;
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const FunctionPatch<M>& x) {
        return x.write(os); }
    //@}

  public:
    Void clobber() { this->_model.clobber(); }
  private:
    friend class TaylorFunctionFactory;
    friend class ScalarFunctionMixin<FunctionPatch<M>, ValidatedTag>;
    friend class ScalarFunctionModelMixin<FunctionPatch<M>, ValidatedTag>;
    template<class T> Void _compute(T& r, const Vector<T>& a) const {
        typedef typename T::NumericType R;
        r=Ariadne::horner_evaluate(this->_model.expansion(),Ariadne::unscale(a,this->_domain))
            + convert_error<R>(this->_model.error());
    }
    FunctionPatch<M>* _derivative(SizeType j) const;
    FunctionPatch<M>* _clone() const;
    FunctionPatch<M>* _create() const;
    VectorFunctionModelInterface<Paradigm>* _create_identity() const;
    VectorFunctionModelInterface<Paradigm>* _create_vector(SizeType i) const;
};

//! \brief Restrict to a subdomain.
template<class M> FunctionPatch<M> restriction(const FunctionPatch<M>& x, const ExactBox& d);
//! \brief Extend over a larger domain. Only possible if the larger domain is only larger where the smaller domain is a singleton.
//! The extension is performed keeping \a x constant over the new coordinates. // DEPRECATED
template<class M> FunctionPatch<M> extension(const FunctionPatch<M>& x, const ExactBox& d);

//! \brief Test if the quantity is a better approximation than \a t throughout the domain.
template<class M> Bool refines(const FunctionPatch<M>& x1, const FunctionPatch<M>& x2);
//! \brief Test if the function models are inconsistent with representing the same exact function.
template<class M> Bool inconsistent(const FunctionPatch<M>& x1, const FunctionPatch<M>& x2);
//! \brief Compute an over-approximation to the common refinement of \a x1 and \a x2.
template<class M> FunctionPatch<M> refinement(const FunctionPatch<M>& x1, const FunctionPatch<M>& x2);


// Normed algebra
template<class M> NormType norm(const FunctionPatch<M>& f);

// Differential algebra
template<class M> FunctionPatch<M> antiderivative(const FunctionPatch<M>& x, SizeType k, ExactFloat c);
template<class M> FunctionPatch<M> antiderivative(const FunctionPatch<M>& x, SizeType k);
template<class M> FunctionPatch<M> derivative(const FunctionPatch<M>& x, SizeType k);

// Function algebra
template<class M> FunctionPatch<M> embed(const ExactBox& dom1, const FunctionPatch<M>& tv2,const ExactBox& dom3);
// Set the value of the \a kth variable to c
template<class M> FunctionPatch<M> partial_evaluate(const FunctionPatch<M>& f, SizeType k, const NumericType<M>& c);
// Evaluate a scalar Taylor function on a vector.
template<class M> NumericType<M> unchecked_evaluate(const FunctionPatch<M>&, const Vector<NumericType<M>>&);

// Compose with an function.
template<class M> FunctionPatch<M> compose(const FunctionType<M>& x, const VectorFunctionPatch<M>& y);
template<class M> FunctionPatch<M> unchecked_compose(const FunctionPatch<M>&, const VectorFunctionPatch<M>&);

// Split the variable over two domains, subdividing along the independent variable j.
template<class M> Pair<FunctionPatch<M>,FunctionPatch<M>> split(const FunctionPatch<M>& x, SizeType j);


template<class M> FunctionPatch<M>& operator+=(FunctionPatch<M>& f, const NumericType<M>& c);
template<class M> FunctionPatch<M>& operator-=(FunctionPatch<M>& f, const NumericType<M>& c);
template<class M> FunctionPatch<M>& operator*=(FunctionPatch<M>& f, const NumericType<M>& c);
template<class M> FunctionPatch<M>& operator/=(FunctionPatch<M>& f, const NumericType<M>& c);
template<class M> FunctionPatch<M>& operator+=(FunctionPatch<M>& f1, const FunctionPatch<M>& f2);
template<class M> FunctionPatch<M>& operator-=(FunctionPatch<M>& f1, const FunctionPatch<M>& f2);

template<class M> FunctionPatch<M> operator+(const FunctionPatch<M>& f);
template<class M> FunctionPatch<M> operator-(const FunctionPatch<M>& f);
template<class M> FunctionPatch<M> operator+(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2);
template<class M> FunctionPatch<M> operator-(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2);
template<class M> FunctionPatch<M> operator*(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2);
template<class M> FunctionPatch<M> operator/(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2);
template<class M> FunctionPatch<M> operator+(const FunctionPatch<M>& f, const NumericType<M>& c);
template<class M> FunctionPatch<M> operator-(const FunctionPatch<M>& f, const NumericType<M>& c);
template<class M> FunctionPatch<M> operator*(const FunctionPatch<M>& f, const NumericType<M>& c);
template<class M> FunctionPatch<M> operator/(const FunctionPatch<M>& f, const NumericType<M>& c);
template<class M> FunctionPatch<M> operator+(const NumericType<M>& c, const FunctionPatch<M>& f);
template<class M> FunctionPatch<M> operator-(const NumericType<M>& c, const FunctionPatch<M>& f);
template<class M> FunctionPatch<M> operator*(const NumericType<M>& c, const FunctionPatch<M>& f);
template<class M> FunctionPatch<M> operator/(const NumericType<M>& c, const FunctionPatch<M>& f);

template<class M> FunctionPatch<M> operator+(const typename FunctionPatch<M>::FunctionType& f1, const FunctionPatch<M>& tf2);
template<class M> FunctionPatch<M> operator-(const typename FunctionPatch<M>::FunctionType& f1, const FunctionPatch<M>& tf2);
template<class M> FunctionPatch<M> operator*(const typename FunctionPatch<M>::FunctionType& f1, const FunctionPatch<M>& tf2);
template<class M> FunctionPatch<M> operator/(const typename FunctionPatch<M>::FunctionType& f1, const FunctionPatch<M>& tf2);
template<class M> FunctionPatch<M> operator+(const FunctionPatch<M>& tf1, const typename FunctionPatch<M>::FunctionType& f2);
template<class M> FunctionPatch<M> operator-(const FunctionPatch<M>& tf1, const typename FunctionPatch<M>::FunctionType& f2);
template<class M> FunctionPatch<M> operator*(const FunctionPatch<M>& tf1, const typename FunctionPatch<M>::FunctionType& f2);
template<class M> FunctionPatch<M> operator/(const FunctionPatch<M>& tf1, const typename FunctionPatch<M>::FunctionType& f2);

template<class M> FunctionPatch<M> abs(const FunctionPatch<M>& x);
template<class M> FunctionPatch<M> neg(const FunctionPatch<M>& x);
template<class M> FunctionPatch<M> rec(const FunctionPatch<M>& x);
template<class M> FunctionPatch<M> sqr(const FunctionPatch<M>& x);
template<class M> FunctionPatch<M> pow(const FunctionPatch<M>& x, Int n);
template<class M> FunctionPatch<M> sqrt(const FunctionPatch<M>& x);
template<class M> FunctionPatch<M> exp(const FunctionPatch<M>& x);
template<class M> FunctionPatch<M> log(const FunctionPatch<M>& x);
template<class M> FunctionPatch<M> sin(const FunctionPatch<M>& x);
template<class M> FunctionPatch<M> cos(const FunctionPatch<M>& x);
template<class M> FunctionPatch<M> tan(const FunctionPatch<M>& x);
template<class M> FunctionPatch<M> asin(const FunctionPatch<M>& x);
template<class M> FunctionPatch<M> acos(const FunctionPatch<M>& x);
template<class M> FunctionPatch<M> atan(const FunctionPatch<M>& x);

// Remove the error term
template<class M> FunctionPatch<M> midpoint(const FunctionPatch<M>& x);

// Code
template<class M> FunctionPatch<M>::FunctionPatch(const ScalarFunctionModel<ValidatedTag>& f) {
     ARIADNE_ASSERT_MSG(dynamic_cast<const FunctionPatch<M>*>(f._ptr.operator->())," f="<<f);
     *this=dynamic_cast<const FunctionPatch<M>&>(*f._ptr);
}





/*! \ingroup FunctionModelSubModule
 *  \brief A Taylor function model with multivalued codomain built from the TaylorModel class.
 *
 *  See also TaylorModel, FunctionPatch<M>, VectorTaylorFunction.
 */
template<class M> class VectorFunctionPatch
    : public VectorFunctionModelMixin<VectorFunctionPatch<M>,typename M::Paradigm>
{
    friend class VectorFunctionPatchElementReference<M>;
    typedef typename M::Paradigm P;
  public:
    typedef ExactBox DomainType;
    typedef M ModelType;
    typedef Box<typename ModelType::CodomainType> CodomainType;
    typedef Box<typename ModelType::RangeType> RangeType;
    typedef typename ModelType::ExpansionType ExpansionType;
    typedef typename ModelType::CoefficientType CoefficientType;
    typedef typename ModelType::ErrorType ErrorType;
    typedef typename ModelType::NumericType NumericType;
    typedef P Paradigm;

    /*! \brief Default constructor constructs a Taylor model of order zero with no arguments and no result variables. */
    VectorFunctionPatch<M>();

    /*! \brief Construct the zero vector function over an unspecified domain. */
    explicit VectorFunctionPatch<M>(SizeType result_size);

    /*! \brief Construct from a result size and a domain. */
    VectorFunctionPatch<M>(SizeType result_size, const ExactBox& domain, Sweeper swp);

    /*! \brief Construct a vector function all of whose components are the same. */
    VectorFunctionPatch<M>(SizeType result_size, const FunctionPatch<M>& scalar_function);

    /*! \brief Construct from a domain and the expansion. */
    VectorFunctionPatch<M>(const ExactBox& domain,
                         const Vector<Expansion<ExactFloat>>& expansion,
                         Sweeper swp);

    /*! \brief Construct from a domain, and expansion and errors. */
    VectorFunctionPatch<M>(const ExactBox& domain,
                         const Vector<Expansion<ExactFloat>>& expansion,
                         const Vector<ErrorFloat>& error,
                         Sweeper swp);

    /*! \brief Construct from a domain, and expansion and errors. */
    VectorFunctionPatch<M>(const ExactBox& domain,
                         const Vector<Expansion<RawFloat>>& expansion,
                         const Vector<RawFloat>& error,
                         Sweeper swp);

    /*! \brief Construct from a domain, and expansion and errors. */
    VectorFunctionPatch<M>(const ExactBox& domain,
                         const Vector<Expansion<RawFloat>>& expansion,
                         Sweeper swp);

    /*! \brief Construct from a domain and the models. */
    explicit VectorFunctionPatch<M>(const ExactBox& domain, const Vector< ModelType >& variables);

    /*! \brief Construct from a domain, a function, and a sweeper determining the accuracy. */
    VectorFunctionPatch<M>(const ExactBox& domain,
                         const VectorFunctionType<M>& function,
                         const Sweeper& sweeper);

    /*! \brief Construct from a vector of scalar Taylor functions. */
    explicit VectorFunctionPatch<M>(const Vector<FunctionPatch<M>>& components);

    /*! \brief Construct from a list of scalar Taylor functions. */
    explicit VectorFunctionPatch<M>(const List<FunctionPatch<M>>& components);

    /*! \brief Construct from an initializer list of scalar Taylor functions. */
    VectorFunctionPatch<M>(InitializerList<FunctionPatch<M>> components);

    /*! \brief Construct from a vector expression. */
    template<class E> explicit VectorFunctionPatch<M>(const VectorExpression<E>& ve);

    explicit VectorFunctionPatch<M> (const VectorFunctionModel<ValidatedTag>& f);
    VectorFunctionPatch<M>& operator=(const VectorFunctionModel<ValidatedTag>& f);

    /*! \brief Equality operator. */
    Bool operator==(const VectorFunctionPatch<M>& p) const;
    /*! \brief Inequality operator. */
    Bool operator!=(const VectorFunctionPatch<M>& p) const;

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
    FunctionPatch<M> get(SizeType i) const;
    /*! \brief Set the \a ith Taylor variable */
    Void set(SizeType i, const FunctionPatch<M>& te);
    /*! \brief The \a ith Taylor variable */
    FunctionPatch<M> operator[](SizeType i) const;
    /*! \brief The \a ith Taylor variable */
    VectorFunctionPatchElementReference<M> operator[](SizeType i);


    /*! \brief Evaluate the Taylor model at the point \a x. */
    Vector<ValidatedNumber> operator()(const Vector<ValidatedNumber>& x) const;
    Vector<ApproximateNumber> operator()(const Vector<ApproximateNumber>& x) const;
    Vector<ValidatedNumber> operator()(const Vector<ExactNumber>& x) const;
    /*! \brief Compute an approximation to Jacobian derivative of the Taylor model sat the point \a x. */
    Matrix<NumericType> jacobian(const Vector<NumericType>& x) const;

    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    VectorFunctionPatch<M>& sweep();
    //! \brief Remove all terms as specified by \a sweeper.
    VectorFunctionPatch<M>& sweep(const SweeperInterface& sweeper);
    /*! \brief Set the error to zero. */
    Void clobber();

    /*! \brief The constant Taylor model with range \a r and argument domain \a d. */
    static VectorFunctionPatch<M> constant(const ExactBox& d, const Vector<NumericType>& r, Sweeper swp);
    /*! \brief The identity Taylor model on domain \a d. */
    static VectorFunctionPatch<M> identity(const ExactBox& d, Sweeper swp);
    //! \brief Return the vector of variables in the range with values \a x over domain \a d.
    static VectorFunctionPatch<M> projection(const ExactBox& d, SizeType imin, SizeType imax, Sweeper swp);

    /*! \brief Convert to an interval polynomial. */
    Vector<Polynomial<ValidatedFloat>> polynomials() const;
    /*! \brief The vector of roundoff/truncation errors of each component. */
    Vector<ErrorType> const errors() const;
    /*! \brief The maximum roundoff/truncation error of the components. */
    ErrorType const error() const;
    //! \brief A multivalued function equal to the model on the domain.
    VectorFunctionType<M> function() const;

    /*! \brief Truncate terms higher than \a bd. */
    VectorFunctionPatch<M>& truncate(const MultiIndexBound& bd);
    /*! \brief Restrict to a subdomain. */
    Void restrict(const ExactBox& d);
    //! \brief Adjoin a scalar function.
    Void adjoin(const FunctionPatch<M>& sf);

    /*! \brief Write to an output stream. */
    OutputStream& write(OutputStream& os) const;

    /*! \brief Write a full representation to an output stream. */
    OutputStream& repr(OutputStream& os) const;

  private:
    Void _compute_jacobian() const;
    Void _set_argument_size(SizeType n);
    SizeType _compute_maximum_component_size() const;
    Void _resize(SizeType rs, SizeType as, ushort d, ushort s);
    virtual ScalarFunctionPatch<M>* _get(SizeType i) const { return new FunctionPatch<M>(this->_domain,this->_models[i]); }
    virtual VectorFunctionPatch<M>* _clone() const;
    virtual VectorFunctionPatch<M>* _create() const;
    virtual VectorFunctionPatch<M>* _create_identity() const;
    virtual ScalarFunctionPatch<M>* _create_zero() const;
  private:
    friend class VectorFunctionMixin<VectorFunctionPatch<M>,ValidatedTag>;
    friend class TaylorFunctionFactory;
    template<class X> Void _compute(Vector<X>& r, const Vector<X>& a) const;
  private:
    /* Domain of definition. */
    ExactBox _domain;
    Vector< ModelType > _models;
};


/*! \brief Inplace addition. */
template<class M> VectorFunctionPatch<M>& operator+=(VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g);
/*! \brief Inplace subtraction. */
template<class M> VectorFunctionPatch<M>& operator-=(VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g);
/*! \brief Inplace addition. */
template<class M> VectorFunctionPatch<M>& operator+=(VectorFunctionPatch<M>& f, const Vector<NumericType<M>>& c);
/*! \brief Inplace subtraction. */
template<class M> VectorFunctionPatch<M>& operator-=(VectorFunctionPatch<M>& f, const Vector<NumericType<M>>& c);
/*! \brief Inplace scalar multiplication. */
template<class M> VectorFunctionPatch<M>& operator*=(VectorFunctionPatch<M>& f, const NumericType<M>& c);
/*! \brief Inplace scalar division. */
template<class M> VectorFunctionPatch<M>& operator/=(VectorFunctionPatch<M>& f, const NumericType<M>& c);

/*! \brief Negation. */
template<class M> VectorFunctionPatch<M> operator-(const VectorFunctionPatch<M>& f);
/*! \brief Addition. */
template<class M> VectorFunctionPatch<M> operator+(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2);
/*! \brief Subtraction. */
template<class M> VectorFunctionPatch<M> operator-(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2);
/*! \brief Multiplication. */
template<class M> VectorFunctionPatch<M> operator*(const FunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2);
/*! \brief Multiplication. */
template<class M> VectorFunctionPatch<M> operator*(const VectorFunctionPatch<M>& f1, const FunctionPatch<M>& f2);
/*! \brief Division. */
template<class M> VectorFunctionPatch<M> operator/(const VectorFunctionPatch<M>& f1, const FunctionPatch<M>& f2);

/*! \brief Addition of a constant. */
template<class M> VectorFunctionPatch<M> operator+(const VectorFunctionPatch<M>& f, const Vector<NumericType<M>>& c);
/*! \brief Subtraction of a constant. */
template<class M> VectorFunctionPatch<M> operator-(const VectorFunctionPatch<M>& f, const Vector<NumericType<M>>& c);
/*! \brief Multiplication by a scalar. */
template<class M> VectorFunctionPatch<M> operator*(const NumericType<M>& c, const VectorFunctionPatch<M>& f);
/*! \brief Multiplication by a scalar. */
template<class M> VectorFunctionPatch<M> operator*(const VectorFunctionPatch<M>& f, const NumericType<M>& c);
/*! \brief Division by a scalar. */
template<class M> VectorFunctionPatch<M> operator/(const VectorFunctionPatch<M>& f, const NumericType<M>& c);
/*! \brief Multiplication by a matrix. */
template<class M> VectorFunctionPatch<M> operator*(const Matrix<ExactNumber>& A, const VectorFunctionPatch<M>& f);
/*! \brief Multiplication by a matrix. */
template<class M> VectorFunctionPatch<M> operator*(const Matrix<NumericType<M>>& A, const VectorFunctionPatch<M>& f);

template<class M> VectorFunctionPatch<M> operator+(const VectorFunctionType<M>& f1, const VectorFunctionPatch<M>& tf2);
template<class M> VectorFunctionPatch<M> operator-(const VectorFunctionType<M>& f1, const VectorFunctionPatch<M>& tf2);
template<class M> VectorFunctionPatch<M> operator*(const ScalarFunctionType<M>& f1, const VectorFunctionPatch<M>& tf2);
template<class M> VectorFunctionPatch<M> operator*(const VectorFunctionType<M>& f1, const FunctionPatch<M>& tf2);
template<class M> VectorFunctionPatch<M> operator/(const VectorFunctionType<M>& f1, const FunctionPatch<M>& tf2);
template<class M> VectorFunctionPatch<M> operator+(const VectorFunctionPatch<M>& tf1, const VectorFunctionType<M>& f2);
template<class M> VectorFunctionPatch<M> operator-(const VectorFunctionPatch<M>& tf1, const VectorFunctionType<M>& f2);
template<class M> VectorFunctionPatch<M> operator*(const FunctionPatch<M>& tf1, const VectorFunctionType<M>& f2);
template<class M> VectorFunctionPatch<M> operator*(const VectorFunctionPatch<M>& tf1, const ScalarFunctionType<M>& f2);
template<class M> VectorFunctionPatch<M> operator/(const VectorFunctionPatch<M>& tf1, const ScalarFunctionType<M>& f2);

//! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
template<class M> FunctionPatch<M> compose(const FunctionType<M>& f, const VectorFunctionPatch<M>& g);
//! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
template<class M> VectorFunctionPatch<M> compose(const VectorFunctionType<M>& f, const VectorFunctionPatch<M>& g);
//! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
template<class M> FunctionPatch<M> compose(const FunctionPatch<M>& f, const VectorFunctionPatch<M>& g);
//! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
template<class M> VectorFunctionPatch<M> compose(const VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g);

//! \brief Weak derivative of \a f with respect to variable \a k.
template<class M> VectorFunctionPatch<M> derivative(const VectorFunctionPatch<M>& f, SizeType k);
//! \brief Antiderivative of \a f with respect to variable \a k.
template<class M> VectorFunctionPatch<M> antiderivative(const VectorFunctionPatch<M>& f, SizeType k);
//! \brief Antiderivative of \a f with respect to variable \a k, taking value \c 0 when \a x[k]=c.
template<class M> VectorFunctionPatch<M> antiderivative(const VectorFunctionPatch<M>& f, SizeType k, ExactFloat c);

template<class M> NormType norm(const VectorFunctionPatch<M>& f);
template<class M> NormType distance(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2);
template<class M> NormType distance(const VectorFunctionPatch<M>& f1, const VectorFunctionType<M>& f2);

//! \brief Restrict the function \a f to a subdomain \a d.
template<class M> VectorFunctionPatch<M> restriction(const VectorFunctionPatch<M>& f, const ExactBox& d);
//! \brief Restrict the function \a f to a larger domain \a d.
template<class M> VectorFunctionPatch<M> extension(const VectorFunctionPatch<M>& f, const ExactBox& d);

template<class M> VectorFunctionPatch<M> embed(const ExactBox& d1, const VectorFunctionPatch<M>& tv2,const ExactBox& d3);

//! \brief Tests if a function \a f refines another function \a g.
//! To be a refinement, the domain of \a f must contain the domain of \a g.
template<class M> Bool refines(const VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g);
template<class M> Bool inconsistent(const VectorFunctionPatch<M>&, const VectorFunctionPatch<M>&);
template<class M> VectorFunctionPatch<M> refinement(const VectorFunctionPatch<M>&, const VectorFunctionPatch<M>&);

//! \brief Compute the function \f$(f \oplus g)(x)=(f(x),g(x))\f$.
template<class M> VectorFunctionPatch<M> join(const VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g);
template<class M> VectorFunctionPatch<M> join(const VectorFunctionPatch<M>& f, const FunctionPatch<M>& g);
template<class M> VectorFunctionPatch<M> join(const FunctionPatch<M>& f, const FunctionPatch<M>& g);
template<class M> VectorFunctionPatch<M> join(const FunctionPatch<M>& f, const VectorFunctionPatch<M>& g);
//! \brief Compute the function \f$(f\otimes g)(x,y)=(f(x),g(y))\f$.
template<class M> VectorFunctionPatch<M> combine(const VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g);
template<class M> VectorFunctionPatch<M> combine(const VectorFunctionPatch<M>& f, const FunctionPatch<M>& g);
template<class M> VectorFunctionPatch<M> combine(const FunctionPatch<M>& f, const VectorFunctionPatch<M>& g);
template<class M> VectorFunctionPatch<M> combine(const FunctionPatch<M>& f, const FunctionPatch<M>& g);

template<class M> VectorFunctionPatch<M> partial_evaluate(const VectorFunctionPatch<M>& f, SizeType k, const NumericType<M>& c);

template<class M> Vector<NumericType<M>> unchecked_evaluate(const VectorFunctionPatch<M>&, const Vector<NumericType<M>>&);
template<class M> FunctionPatch<M> unchecked_compose(const FunctionPatch<M>&, const VectorFunctionPatch<M>&);
template<class M> VectorFunctionPatch<M> unchecked_compose(const VectorFunctionPatch<M>&, const VectorFunctionPatch<M>&);

// Split the domain into halves along the \a j<sup>th</sup> coordinate.
template<class M> Pair<VectorFunctionPatch<M>,VectorFunctionPatch<M>> split(const VectorFunctionPatch<M>& x, SizeType j);

template<class M> OutputStream& operator<<(OutputStream&, const VectorFunctionPatch<M>&);

// Conversion operatations
template<class M> Polynomial<ValidatedFloat> polynomial(const FunctionPatch<M>& tfn);
template<class M> Vector< Polynomial<ValidatedFloat> > polynomials(const VectorFunctionPatch<M>& tfn);


// Sanitised output
template<class M> OutputStream& operator<<(OutputStream&, const Representation<FunctionPatch<M>>&);
template<class M> OutputStream& operator<<(OutputStream&, const Representation<VectorFunctionPatch<M>>&);
template<class F> struct ModelRepresentation { const F* pointer; double threshold; };
template<class F> ModelRepresentation<F> model_representation(const F& f, double swpt) {
    ModelRepresentation<F> r={&f,swpt}; return r; }
template<class M> OutputStream& operator<<(OutputStream&,const ModelRepresentation<FunctionPatch<M>>&);
template<class M> OutputStream& operator<<(OutputStream&,const ModelRepresentation<VectorFunctionPatch<M>>&);
template<class F> struct PolynomialRepresentation { const F* pointer; double threshold; List<String> names; };
template<class F> PolynomialRepresentation<F> polynomial_representation(const F& f, double swpt) {
    PolynomialRepresentation<F> r={&f,swpt}; return r; }
template<class F> PolynomialRepresentation<F> polynomial_representation(const F& f, double swpt, const List<String>& names) {
    PolynomialRepresentation<F> r={&f,swpt,names}; return r; }
template<class M> OutputStream& operator<<(OutputStream&,const PolynomialRepresentation<FunctionPatch<M>>&);
template<class M> OutputStream& operator<<(OutputStream&,const PolynomialRepresentation<VectorFunctionPatch<M>>&);


template<class M> template<class E> VectorFunctionPatch<M>::VectorFunctionPatch(const VectorExpression<E>& ve)
    : _domain(), _models(ve().size(),ve().zero_element().model())
{
    if(ve().size()!=0) { this->_domain=ve().zero_element().domain(); }
    for(SizeType i=0; i!=ve().size(); ++i) { this->set(i,ve()[i]); }
}

template<class M> class VectorFunctionPatchElementReference
{
    typedef M ModelType;
    typedef typename M::NumericType NumericType;
    typedef VectorFunctionPatchElementReference<M> SelfType;
    typedef FunctionPatch<M> FunctionType;
 public:
    VectorFunctionPatchElementReference(VectorFunctionPatch<M>& c, SizeType i) : _c(&c), _i(i) { }
    operator ScalarFunctionPatch<M> () const { return this->_c->get(this->_i); }
    Void operator=(const VectorFunctionPatchElementReference<M>& x) { this->_c->set(this->_i,x._c->get(x._i)); }
    Void operator=(const ScalarFunctionPatch<M>& x) { this->_c->set(this->_i,x); }
    ScalarFunctionPatch<M> element() const { return this->_c->get(this->_i); }
    ExactBox const& domain() const { return this->_c->domain(); }
    const ModelType& model() const { return this->_c->_models[this->_i]; }
    ErrorType error() const { return this->_c->_models[this->_i].error(); }
    Void set_error(const ErrorType& e) { this->_c->_models[this->_i].set_error(e); }
    Void sweep() { this->_c->_models[this->_i].sweep(); }
    template<class X> X operator()(const Vector<X>& x) const { return this->_c->get(this->_i).operator()(x); }
    friend OutputStream& operator<<(OutputStream& os, const VectorFunctionPatchElementReference<M>& f) { return os<<f.element(); }
  public:
    friend FunctionType operator+(SelfType e1, SelfType e2) { return e1.element()+e2.element(); }
    friend FunctionType operator+(NumericType c, SelfType e) { return c+e.element(); }
    friend FunctionType operator+(SelfType e, NumericType c) { return e.element()+c; }
    friend FunctionType operator+(FunctionType f, SelfType e) { return f+e.element(); }
    friend FunctionType operator+(SelfType e, FunctionType f) { return e.element()+f; }
    friend FunctionType operator*(SelfType e1, SelfType e2) { return e1.element()*e2.element(); }
    friend FunctionType operator*(NumericType c, SelfType e) { return c*e.element(); }
    friend FunctionType operator*(SelfType e, NumericType c) { return e.element()*c; }
    friend FunctionType operator*(FunctionType f, SelfType e) { return f*e.element(); }
    friend FunctionType operator*(SelfType e, FunctionType f) { return e.element()*f; }
  private:
    VectorFunctionPatch<M>* _c; SizeType _i;
};

} // namespace Ariadne

#include "function_patch.tcc"

#endif // ARIADNE_FUNCTION_PATCH_H
