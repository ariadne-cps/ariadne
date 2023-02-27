/***************************************************************************
 *            function/taylor_model.hpp
 *
 *  Copyright  2008-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file function/taylor_model.hpp
 *  \brief ApproximateTag functions on a is_bounded domain with a sparse representation.
 */

#ifndef ARIADNE_TAYLOR_MODEL_HPP
#define ARIADNE_TAYLOR_MODEL_HPP

#include <map>

#include "utility/macros.hpp"
#include "utility/declarations.hpp"
#include "utility/array.hpp"
#include "utility/pointer.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/expansion.hpp"
#include "algebra/sweeper.hpp"
#include "algebra/operations.hpp"
#include "algebra/evaluate.hpp"
#include "function/domain.hpp"
#include "function/scaling.hpp"
#include "function/polynomial.hpp"

namespace Ariadne {

template<class FLT> class ZeroError {
  public:
    template<class... PRS> requires Constructible<Error<FLT>,PRS...> ZeroError(PRS...) { }
    FLT raw() const;
    typename FLT::PrecisionType precision() const;
    operator Error<FLT> () const;
};

template<class FLT> class UnknownError {
  public:
//    template<class... PRS> requires Constructible<Error<FLT>,PRS...> UnknownError(PRS...) { }
    template<class... PRS> UnknownError(PRS...) { }
    FLT raw() const;
    typename FLT::PrecisionType precision() const;
    operator PositiveApproximation<FLT> () const;
};

template<class T1, class T2> struct Product;

template<class P, class FLT> class TaylorModel;

//! \relates TaylorModel
//! \name Template shorthands and type synonyms for Taylor models
//!@{
template<class FLT> using ValidatedTaylorModel = TaylorModel<ValidatedTag,FLT>; //!< <p/>
template<class FLT> using ValidatedIntervalTaylorModel = TaylorModel<ValidatedTag,UpperInterval<FLT>>; //!< <p/>
template<class FLT> using ApproximateTaylorModel = TaylorModel<ApproximateTag,FLT>; //!< <p/>
using ValidatedTaylorModelDP = TaylorModel<ValidatedTag,FloatDP>; //!< <p/>
using ValidatedTaylorModelMP = TaylorModel<ValidatedTag,FloatMP>; //!< <p/>
using ValidatedIntervalTaylorModelDP = TaylorModel<ValidatedTag,FloatDPUpperInterval>; //!< <p/>
using ValidatedIntervalTaylorModelMP = TaylorModel<ValidatedTag,FloatMPUpperInterval>; //!< <p/>
using ApproximateTaylorModelDP = TaylorModel<ApproximateTag,FloatDP>; //!< <p/>
using ApproximateTaylorModelMP = TaylorModel<ApproximateTag,FloatMP>; //!< <p/>
//!@}

template<class P, class FLT> struct IsScalar< TaylorModel<P,FLT> > { static const Bool value = true; };

class IntersectionException;

class IntersectionException : public std::runtime_error {
  public:
    IntersectionException(const StringType& what) : std::runtime_error(what) { }
};

template<class P, class FLT> struct ModelNumericTraits;
template<class FLT> struct ModelNumericTraits<ValidatedTag,Interval<UpperBound<FLT>>>
    : public FunctionModelTraits<ValidatedTag,FLT>
{
    typedef Interval<UpperBound<FLT>> CoefficientType;
    typedef Interval<UpperBound<FLT>> NumericType;
//    typedef Bounds<FLT> NumericType;
    typedef FLT RawFloatType;
};
template<class FLT> struct ModelNumericTraits<ValidatedTag,Bounds<FLT>>
    : public FunctionModelTraits<ValidatedTag,FLT>
{
    typedef Bounds<FLT> CoefficientType;
    typedef Bounds<FLT> NumericType;
    typedef FLT RawFloatType;
};
template<class FLT> struct ModelNumericTraits<ValidatedTag,FLT>
    : public FunctionModelTraits<ValidatedTag,FLT>
{
    typedef FLT CoefficientType;
};
template<class FLT> struct ModelNumericTraits<ApproximateTag,FLT>
    : public FunctionModelTraits<ApproximateTag,FLT>
{
    typedef Approximation<FLT> CoefficientType;
};

template<class P, class FLT> using ModelNumericType = typename ModelNumericTraits<P,FLT>::NumericType;


template<class P, class FLT> struct ModelNumericTypedef;
template<class FLT> struct ModelNumericTypedef<ValidatedTag,UpperInterval<FLT>> { typedef UpperInterval<FLT> Type; };
template<class FLT> struct ModelNumericTypedef<ValidatedTag,FLT> { typedef FloatBounds<typename FLT::PrecisionType> Type; };
template<class FLT> struct ModelNumericTypedef<ApproximateTag,FLT> { typedef FloatApproximation<typename FLT::PrecisionType> Type; };

template<class P, class FLT> struct AlgebraOperations<TaylorModel<P,FLT>>
    : NormedAlgebraOperations<TaylorModel<P,FLT>>
{
    typedef ModelNumericType<P,FLT> X;
    typedef TaylorModel<P,FLT> ModelType;
    typedef ModelNumericType<P,FLT> NumericType;
    using NormedAlgebraOperations<TaylorModel<P,FLT>>::apply;
    static TaylorModel<P,FLT> apply(Nul,TaylorModel<P,FLT> const& tm);
    static TaylorModel<P,FLT> apply(Pos,TaylorModel<P,FLT> tm);
    static TaylorModel<P,FLT> apply(Neg,TaylorModel<P,FLT> tm);
    static TaylorModel<P,FLT> apply(Add,TaylorModel<P,FLT> tm, X const& c);
    static TaylorModel<P,FLT> apply(Mul,TaylorModel<P,FLT> tm, X const& c);
    static TaylorModel<P,FLT> apply(Add,TaylorModel<P,FLT> const& tm1, TaylorModel<P,FLT> const& tm2);
    static TaylorModel<P,FLT> apply(Sub,TaylorModel<P,FLT> const& tm1, TaylorModel<P,FLT> const& tm2);
    static TaylorModel<P,FLT> apply(Mul,TaylorModel<P,FLT> const& tm1, TaylorModel<P,FLT> const& tm2);
    static TaylorModel<P,FLT> apply(Min,TaylorModel<P,FLT> const& tm1, TaylorModel<P,FLT> const& tm2);
    static TaylorModel<P,FLT> apply(Max,TaylorModel<P,FLT> const& tm1, TaylorModel<P,FLT> const& tm2);
    static TaylorModel<P,FLT> apply(Max,TaylorModel<P,FLT> const& tm, X const& c);
    static TaylorModel<P,FLT> apply(Min,TaylorModel<P,FLT> const& tm, X const& c);
    static TaylorModel<P,FLT> apply(Max,X const& c, TaylorModel<P,FLT> const& tm);
    static TaylorModel<P,FLT> apply(Min,X const& c, TaylorModel<P,FLT> const& tm);
    static TaylorModel<P,FLT> apply(Abs,TaylorModel<P,FLT> const& tm);
};

template<class P, class FLT> class TaylorModel;

//! \ingroup FunctionModelSubModule
//! \brief A class representing polynomial approximation to a function, scaled to the unit box, with a uniform error bound.
//!  \tparam P	The information paradigm tag, which can be either ValidatedTag, indicating that the approximation has a known (uniform) error bound, or ApproximateTag, indicating that no error bound is available. See the \ref InformationSubModule for more details.
//!  \tparam FLT  The (floating-point) number type used internally for the coefficients of the polynomial, such as \ref FloatDP or \ref FloatMP.
//!  \tparam FLTE  The (floating-point) number type used internally for the error bound.
//!
//! \see Expansion, ValidatedScalarMultivariateTaylorFunctionModelDP, ValidatedVectorMultivariateTaylorFunctionModelDP, TaylorConstrainedImageSet.
template<class P, class FLT>
class TaylorModel
    : public DispatchElementaryAlgebraOperations<TaylorModel<P,FLT>,typename ModelNumericTraits<P,FLT>::NumericType>
    , public DispatchConcreteGenericAlgebraNumberOperations<TaylorModel<P,FLT>,typename ModelNumericTraits<P,FLT>::NumericType,Number<P>>
{
    typedef typename FLT::PrecisionType PR;
    typedef typename FLT::PrecisionType PRE;
    typedef typename ModelNumericTraits<P,FLT>::NumericType X;
  public:
    typedef PR PrecisionType;
    typedef PRE ErrorPrecisionType;

    typedef typename ModelNumericTraits<P,FLT>::RawFloatType RawFloatType;
    typedef typename ModelNumericTraits<P,FLT>::CoefficientType CoefficientType;
    typedef typename ModelNumericTraits<P,FLT>::ValueType ValueType;
    typedef typename ModelNumericTraits<P,FLT>::ErrorType ErrorType;
    typedef typename ModelNumericTraits<P,FLT>::ErrorValueType ErrorValueType;
    typedef typename ModelNumericTraits<P,FLT>::NormType NormType;
    typedef ReverseLexicographicIndexLess ComparisonType;
    typedef SortedExpansion<MultiIndex,CoefficientType,ComparisonType> ExpansionType;
    typedef Polynomial<MultiIndex,CoefficientType> PolynomialType;
    typedef Sweeper<RawFloatType> SweeperType;


    typedef typename ModelNumericTraits<P,FLT>::RangeType RangeType;
    typedef IntervalDomainType CodomainType;

    template<class X> using Argument = Vector<X>;
    template<class X> using Result = Scalar<X>;

    //! \brief The computational paradigm.
    typedef ValidatedTag Paradigm;
    //! \brief The properties needed to define the TaylorModel calculus.
    typedef Sweeper<RawFloatType> PropertiesType;

    //! \brief The type used for algebraic operations.
    typedef typename ModelNumericTraits<P,FLT>::NumericType NumericType;
    typedef Number<P> GenericNumericType;

    typedef Interval<FloatUpperBound<PR>> IntervalNumericType;
    typedef FloatBounds<PR> ValidatedNumericType;
    typedef FloatApproximation<PR> ApproximateNumericType;

//    typedef decltype(std::declval<NumericType>() < std::declval<NumericType>()) ComparisonType;

    typedef ValidatedScalarMultivariateFunction FunctionType;
    typedef ValidatedScalarMultivariateFunction ScalarFunctionType;
    typedef ValidatedVectorMultivariateFunction VectorFunctionType;

    //! \brief The type used to index the coefficients.
    typedef MultiIndex IndexType;

    //! \brief An Iterator through the (index,coefficient) pairs of the expansion.
    typedef typename ExpansionType::Iterator Iterator;
    //! \brief A constant Iterator through the (index,coefficient) pairs of the expansion.
    typedef typename ExpansionType::ConstIterator ConstIterator;
  private:
    PolynomialType _polynomial;
    ErrorType _error;
    mutable SweeperType _sweeper;
  public:
    //! \name Constructors and destructors.
    //!@{
    //
    //! \brief Default constructor.
    TaylorModel();
    //! \brief Construct a TaylorModel in \a as arguments with the given accuracy control.
    TaylorModel(SizeType as, SweeperType swp);
    //! \brief Construct from a map giving the expansion, a constant giving the error, and an accuracy parameter.
    TaylorModel(const Expansion<MultiIndex,ExactDouble>& f, const ExactDouble& e, SweeperType swp);
    TaylorModel(const Expansion<MultiIndex,CoefficientType>& f, const ErrorType& e, SweeperType swp);
    TaylorModel(const Expansion<MultiIndex,RawFloatType>& f, const RawFloatType& e, SweeperType swp);
    //! \brief Fast swap with another Taylor model.
    Void swap(TaylorModel<P,FLT>& tm);
    //! \brief The zero element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<P,FLT> create() const;
    //! \brief The zero element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<P,FLT> create_zero() const;
    //! \brief A constant element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<P,FLT> create_constant(NumericType c) const;
    TaylorModel<P,FLT> create_constant(GenericNumericType) const;
    //! \brief The \a j<sup>th</sup> coordinate element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<P,FLT> create_coordinate(SizeType j) const;
    //! \brief Set to zero.
    TaylorModel<P,FLT> create_ball(ErrorType e) const;
    //! \brief Set to zero.
    Void clear();
    //!@}

    //! \name Assignment to constant values.
    //!@{
    //
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorModel<P,FLT>& operator=(const NumericType& c);
    TaylorModel<P,FLT>& operator=(const GenericNumericType& c);
    //!@}

    //! \name Named constructors.
    //!@{
    //
    //! \brief Construct the zero quantity in \a as independent variables.
    static TaylorModel<P,FLT> zero(SizeType as, SweeperType swp) {
        TaylorModel<P,FLT> r(as,swp); return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<P,FLT> constant(SizeType as, const NumericType& c, SweeperType swp);
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<P,FLT> constant(SizeType as, const GenericNumericType& c, SweeperType swp);
    //! \brief Construct the quantity with expansion \f$x_j\f$ in \a as independent variables.
    static TaylorModel<P,FLT> coordinate(SizeType as, SizeType j, SweeperType swp) {
        TaylorModel<P,FLT> r(as,swp); r.set_gradient(j,1); return r; }
    //! \brief Construct a constant quantity in \a as independent variables with value zero and uniform error \a e
    static TaylorModel<P,FLT> error(SizeType as, ErrorType e, SweeperType swp) {
        TaylorModel<P,FLT> r(as,swp); r.set_error(e); return r; }
    //! \brief Construct a constant quantity in \a as independent variables with value zero and uniform error \a e
    static TaylorModel<P,FLT> ball(SizeType as, ErrorType e, SweeperType swp) {
        TaylorModel<P,FLT> r(as,swp); r.set_error(e); return r; }
    //! \brief Construct a constant quantity in \a as independent variables with value zero and uniform error \a e
    static TaylorModel<P,FLT> unit_ball(SizeType as, SweeperType swp) {
        TaylorModel<P,FLT> r(as,swp); r.set_error(1u); return r; }

    //! \brief Construct the quantity which scales the interval \a codom onto the unit interval.
    static TaylorModel<P,FLT> scaling(SizeType as, SizeType j, const IntervalDomainType& codom, SweeperType swp);

    //! \brief Return the vector of zero variables of size \a rs in \a as arguments.
    static Vector<TaylorModel<P,FLT>> zeros(SizeType rs, SizeType as, SweeperType swp);
    //! \brief Return the vector of constants with values \a c in \a as arguments.
    static Vector<TaylorModel<P,FLT>> constants(SizeType as, const Vector<NumericType>& c, SweeperType swp);
    //! \brief Return the vector of variables on the unit domain.
    static Vector<TaylorModel<P,FLT>> coordinates(SizeType as, SweeperType swp);

    //! \brief Return the vector scaling the box \a codom onto the unit box.
    static Vector<TaylorModel<P,FLT>> scalings(const BoxDomainType& codom, SweeperType swp);
    //!@}

    //! \name Comparison operators.
    //!@{
    //
    //! \brief Equality operator. Tests equality of representation, including error term.
    friend Bool same(const TaylorModel<P,FLT>& tm1, const TaylorModel<P,FLT>& tm2) {
        return same(tm1.expansion(), tm2.expansion()) && same(tm1._error, tm2._error); }

    decltype(auto) operator<(const TaylorModel<P,FLT>& sd) const {
        return (sd-*this)>0; }
    //! \brief Comparison with another Taylor model.
    decltype(auto) operator>(const TaylorModel<P,FLT>& sd) const {
        return (*this-sd)>0; }

    //! \brief Comparison with a scalar.
    decltype(auto) operator<(Int c) const {
        return this->range().upper_bound()<c; }
    //! \brief Comparison with a scalar.
    decltype(auto) operator>(Int c) const {
        return this->range().lower_bound()>c; }
    //!@}

    //! \name Data access
    //!@{
    //
    //! \brief The expansion.
    const ExpansionType& expansion() const { return this->_polynomial.expansion(); }
    //! \brief A reference to the expansion.
    ExpansionType& expansion() { return this->_polynomial.expansion(); }
    //! \brief The polynomial approximation.
    const PolynomialType& polynomial() const { return this->_polynomial; }
    //! \brief A reference to the polynomial approximation.
    PolynomialType& polynomial() { return this->_polynomial; }
    //! \brief The error of the expansion over the domain.
    const ErrorType& error() const { return this->_error; }
    //! \brief A reference to the error of the expansion over the domain.
    ErrorType& error() { return this->_error; }

    //! \brief The constant term in the expansion.
    const ValueType value() const { return this->average(); }
//    //! \brief The coefficient of the gradient term \f$df/dx_j\f$.
//    const ValueType gradient_value(SizeType j) const ;//{ return (*this)[MultiIndex::unit(this->argument_size(),j)]; }
//    Covector<ValueType> gradient_value() const { Covector<ValueType> r(this->argument_size());
//        for(SizeType j=0; j!=this->argument_size(); ++j) { r[j]=this->gradient_value(j); } return r; }

    //! \brief The constant term in the expansion.
    ValueType average() const;
    //! \brief The radius of the smallest ball about a constant function containing the model.
    NormType radius() const;
    //! \brief An over-approximation to the supremum norm.
    NormType norm() const;
    //! \brief An over-approximation to the supremum norm.
    friend typename TaylorModel<P,FLT>::NormType norm(const TaylorModel<P,FLT>& tm) { return tm.norm(); }
    friend typename TaylorModel<P,FLT>::NormType mag(const TaylorModel<P,FLT>& tm) { return tm.norm(); }


    //! \brief A value \c e such that analytic functions are evaluated to a tolerance of \c e. Equal to the sweep threshold.
    RawFloatType tolerance() const;


    //! \brief Set the constant term in the expansion.
    Void set_value(const CoefficientType& c) {
        this->expansion().set(MultiIndex::zero(this->argument_size()),c); }
    //! \brief Set the coefficient of the term \f$df/dx_j\f$.
    Void set_gradient(SizeType j, const CoefficientType& c) {
        this->expansion().set(MultiIndex::unit(this->argument_size(),j),c); }
    Void set_gradient(SizeType j,const Dyadic& c) {
        this->set_gradient(j,CoefficientType(c,this->precision())); }
     //! \brief Set the error of the expansion.
    Void set_error(const ErrorType& ne) { ARIADNE_ASSERT(ne.raw()>=0.0_x); this->_error=ne; }
    Void set_error(Nat m) { this->_error=m; }

    //! \brief The coefficient of the term in $x^a$.
    const CoefficientType& operator[](const MultiIndex& a) const { return this->expansion()[a]; }
    //! \brief A read/write reference to the coefficient of the term in $x^a$.
    CoefficientType& operator[](const MultiIndex& a) { return this->expansion().at(a); }

    //! \brief An Iterator to the first term in the expansion.
    Iterator begin() { return this->expansion().begin(); }
    //! \brief A constant Iterator to the first term in the expansion.
    ConstIterator begin() const { return this->expansion().begin(); }
    //! \brief An Iterator to the end of the expansion.
    Iterator end() { return this->expansion().end(); }
    //! \brief A constant Iterator to the end of the expansion.
    ConstIterator end() const { return this->expansion().end(); }

    //! \brief The number of variables in the argument of the quantity.
    SizeType argument_size() const { return this->expansion().argument_size(); }
    //! \brief The maximum degree of terms in the expansion.
    DegreeType degree() const;
    //! \brief The number of nonzero terms in the expansion.
    SizeType number_of_terms() const { return this->expansion().number_of_terms(); }
    SizeType number_of_nonzeros() const { return this->expansion().number_of_nonzeros(); }
    //!@}

    //! \name Function evaluation.
    //!@{
    //
    //! \brief The domain of the quantity, always given by \f$[-1,1]^{\mathrm{as}}\f$.
    UnitBox domain() const;
    //! \brief The codomain of the quantity.
    IntervalDomainType codomain() const;
    //! \brief An over-approximation to the range of the quantity.
    RangeType range() const;

    //! \brief Evaluate the quantity over the interval of points \a x.
    friend ArithmeticType<CoefficientType,ApproximateNumericType> evaluate(const TaylorModel<P,FLT>& f, const Vector<ApproximateNumericType>& x) {
        return TaylorModel<P,FLT>::_evaluate(f,x); }
    friend ArithmeticType<CoefficientType,ValidatedNumericType> evaluate(const TaylorModel<P,FLT>& f, const Vector<ValidatedNumericType>& x) {
        return TaylorModel<P,FLT>::_evaluate(f,x); }
    template<class X=NumericType> requires (not Same<X,ValidatedNumericType>) friend ArithmeticType<CoefficientType,ValidatedNumericType> evaluate(const TaylorModel<P,FLT>& f, const Vector<X>& x) {
        return TaylorModel<P,FLT>::_evaluate(f,x); }
    //! \brief Evaluate the gradient over the interval of points \a x.
    friend Covector<NumericType> gradient(const TaylorModel<P,FLT>& f, const Vector<NumericType>& x) {
        return TaylorModel<P,FLT>::_gradient(f,x); }
    //! \brief Substite \a c for the \a k th variable.
    friend TaylorModel<P,FLT> partial_evaluate(const TaylorModel<P,FLT>& x, SizeType k, NumericType c) {
        return TaylorModel<P,FLT>::_partial_evaluate(x,k,c); }
    //! \relates TaylorModel<P,FLT Compose a vector of Taylor models with another.
    friend TaylorModel<P,FLT> compose(const Unscaling& uf, const TaylorModel<P,FLT>& tg) {
        return TaylorModel<P,FLT>::_compose(uf,tg); }
    //! \brief Evaluate the quantity over the interval of points \a x.
    friend TaylorModel<P,FLT> compose(const TaylorModel<P,FLT>& f, const Vector<TaylorModel<P,FLT>>& g) {
        return TaylorModel<P,FLT>::_compose(f,g); }
    //! \brief Scale the variable by post-composing with an affine map taking the interval ivl to the unit interval
    friend TaylorModel<P,FLT> compose(const TaylorModel<P,FLT>& tf, const VectorUnscaling& u);

    friend TaylorModel<P,FLT> evaluate(const TaylorModel<P,FLT>& f, const Vector<TaylorModel<P,FLT>>& g) { return compose(f,g); }
    template<class A> ArithmeticType<CoefficientType,A> operator() (Vector<A> const&) const;
    //!@}

    //! \name Inplace modifications.
    //!@{
    //
    //! \brief Scales the model by a function mapping \a dom into the unit interval.
    Void unscale(IntervalDomainType const& dom);
    //! \brief Compute the antiderivative (in place).
    Void antidifferentiate(SizeType k);
    //! \brief Compute the weak derivative (in place).
    Void differentiate(SizeType k);
    //! \brief Estimate the gradient at a point.
    friend Covector<X> gradient(const TaylorModel<P,FLT>& x, const Vector<X>& v);
//    Covector<X> gradient(Vector<X> const& v) const;

    friend TaylorModel<P,FLT> unscale(TaylorModel<P,FLT> tm, IntervalDomainType const& dom) {
        tm.unscale(dom); return tm; }
    friend TaylorModel<P,FLT> antiderivative(TaylorModel<P,FLT> tm, SizeType k) {
        tm.antidifferentiate(k); return tm; }
    friend TaylorModel<P,FLT> derivative(TaylorModel<P,FLT> tm, SizeType k) {
        tm.differentiate(k); return tm; }
    //!@}

    //! \name Operations on the domain.
    //!@{
    //
    //! \brief Split the Taylor model \a tm, subdividing along the independent variable \a k, taking the lower/middle/upper \a part.
    friend TaylorModel<P,FLT> split(const TaylorModel<P,FLT>& tm, SizeType k, SplitPart part) {
        return TaylorModel<P,FLT>::_split(tm,k,part); }
    //! \brief Split the Taylor model \a tm, subdividing along the independent variable \a k.
    friend Pair<TaylorModel<P,FLT>,TaylorModel<P,FLT>> split(const TaylorModel<P,FLT>& tm, SizeType k) {
        return make_pair(split(tm,k,SplitPart::LOWER),split(tm,k,SplitPart::UPPER)); }
    //! \relates TaylorModel<P,FLT> \brief Embed the model in a space of higher dimension
    friend TaylorModel<P,FLT> embed(SizeType as1, const TaylorModel<P,FLT>& tm2, SizeType as3) {
        return TaylorModel<P,FLT>::_embed(as1,tm2,as3); }
    //!@}

    //! \name Paradigm-based operations.
    //!@{
    //
    //! \brief Test if two models are compatible with each other, in the sense that they could both represent the same function.
    //! Equivalent to checking that the two sets of possible functions intersect.
    //! If the result is \a true, then the models are definitely consistent; may return \a false even if the models are inconsistent.
    friend Bool consistent(const TaylorModel<P,FLT>& tm1, const TaylorModel<P,FLT>& tm2) {
        return TaylorModel<P,FLT>::_consistent(tm1,tm2); }
    //! \brief Test if one model is incompatible with another.
    //! Equivalent to checking that the two sets of possible functions are disjoint.
    //! If the result is \a true, then the models are definitely inconsistent; may return \a false even if the models are inconsistent.
    friend Bool inconsistent(const TaylorModel<P,FLT>& tm1, const TaylorModel<P,FLT>& tm2) {
        return TaylorModel<P,FLT>::_inconsistent(tm1,tm2); }
    //! \brief Test if one model refines another, in the sense that it represents a subset of possible functions.
    //! If the result is \a true, then the first model definitely refines the second;
    //! may return \a false even if the first model actually does refine the second.
    friend Bool refines(const TaylorModel<P,FLT>& tm1, const TaylorModel<P,FLT>& tm2) {
        return TaylorModel<P,FLT>::_refines(tm1,tm2); }
    //! \brief An over-approximation to the common refinement of two Taylor models.
    //! Since the intersection of represented sets of functions cannot be represented
    //! exactly in the class of TaylorModels, truncation errors as well as roundoff errors
    //! may be present. In the absence of roundoff errors, the result is a subset of both
    //! arguments, and is guaranteed to contain any function contained in both arguments.
    friend TaylorModel<P,FLT> refinement(const TaylorModel<P,FLT>& tm1, const TaylorModel<P,FLT>& tm2) {
        return TaylorModel<P,FLT>::_refinement(tm1,tm2); }
    //!@}

    //! \name Simplification operations.
    //!@{
    //
    //! \brief Embed the model in a space of higher dimension, placing the error in the final variable.
    friend TaylorModel<P,FLT> embed_error(const TaylorModel<P,FLT>& tm) {
        return TaylorModel<P,FLT>::_embed_error(tm); }
    //! \relates TaylorModel<P,FLT> \brief Abstract away the given variables.
    //! For example, the model tm(x0,x1,x2,x3) becomes tm'(x0,x1)=tm(x0,[-1,+1],x1,[-1,+1]) on discarding x1 and x3.
    friend TaylorModel<P,FLT>  discard_variables(const TaylorModel<P,FLT>& tm, const Array<SizeType>& variables) {
        return TaylorModel<P,FLT>::_discard_variables(tm,variables); }

    //! \brief Remove all terms based on the \a swp conditions.
    TaylorModel<P,FLT>& sweep(const SweeperType& swp);
    TaylorModel<P,FLT>& simplify(const PropertiesType& prp);

    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    TaylorModel<P,FLT>& sweep();
    TaylorModel<P,FLT>& simplify();
    //! \brief Sorts keys.
    TaylorModel<P,FLT>& sort();
    //! \brief Remove terms with the same keys. Assumes sorted.
    TaylorModel<P,FLT>& unique();
    //! \brief Sort the terms in index order and combine terms with the same index.
    TaylorModel<P,FLT>& cleanup();
    //! \brief Set the error to zero.
    //! WARNING: This method does not preserve rigour of the model approximation.
    TaylorModel<P,FLT>& clobber();
    //!@}

    //! \name Accuracy parameters.
    //!@{
    //
    //! \brief Specify a policy to use to remove low-impact terms.
    Void set_sweeper(SweeperType swp) { this->_sweeper=swp; }
    Void set_properties(PropertiesType prp) { this->_sweeper=prp; }
    //! \brief A shared pointer to an object using for removing low-impact terms.
    SweeperType sweeper() const { return this->_sweeper; }
    //! \brief The precision of the coefficients.
    PropertiesType properties() const { return this->sweeper(); }
    //! \brief The precision of the coefficients.
    PrecisionType precision() const { return this->sweeper().precision(); }
    //!@}

    //! \name Order operators.
    //!@{
    //
    //! \brief The pointwise maximum.
    template<class FF> friend TaylorModel<P,FF> max(const TaylorModel<P,FF>& x, const TaylorModel<P,FF>& y);
    //! \brief The pointwise minimum.
    template<class FF> friend TaylorModel<P,FF> min(const TaylorModel<P,FF>& x, const TaylorModel<P,FF>& y);
    //! \brief The pointwise absolute value.
    //! \details If the range of \a x definitely does not include 0, returns +x or -x. Otherwise, uses a uniform polynomial approximation to abs.
    template<class FF> friend TaylorModel<P,FF> abs(const TaylorModel<P,FF>& x);
    //!@}

    //! \name Stream input/output operators.
    //!@{
    //
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const TaylorModel<P,FLT>& x) { return x.str(os); }
    //!@}

  public:
    OutputStream& str(OutputStream&) const;
    OutputStream& repr(OutputStream&) const;
  public: // FIXME: Should be private
    Void _set_error(const RawFloat<PR>& ne) { ARIADNE_ASSERT(ne>=0); this->_error=ErrorType(ne); }
    Void _append(MultiIndex const& a, CoefficientType const& v) { this->expansion().append(a,v); }
    static Bool _consistent(const TaylorModel<P,FLT>& tm1, const TaylorModel<P,FLT>& tm2);
    static Bool _inconsistent(const TaylorModel<P,FLT>& tm1, const TaylorModel<P,FLT>& tm2);
    static Bool _refines(const TaylorModel<P,FLT>& tm1, const TaylorModel<P,FLT>& tm2);
    static TaylorModel<P,FLT> _refinement(const TaylorModel<P,FLT>& tm1, const TaylorModel<P,FLT>& tm2);
    static TaylorModel<P,FLT> _antiderivative(const TaylorModel<P,FLT>& tm, SizeType k);
    static TaylorModel<P,FLT> _weak_derivative(const TaylorModel<P,FLT>& tm, SizeType k);
    static TaylorModel<P,FLT> _embed_error(const TaylorModel<P,FLT>& tm);
    static TaylorModel<P,FLT>  _discard_variables(const TaylorModel<P,FLT>&, const Array<SizeType>& variables);
    static TaylorModel<P,FLT> _split(const TaylorModel<P,FLT>& tm, SizeType k, SplitPart part);
    static TaylorModel<P,FLT> _embed(SizeType as1, const TaylorModel<P,FLT>& tm2, SizeType as3);
    static TaylorModel<P,FLT> _compose(TaylorModel<P,FLT> const& tf, const Vector<TaylorModel<P,FLT>>& tg);
    static TaylorModel<P,FLT> _compose(Unscaling const& uf, TaylorModel<P,FLT> const& tg);
    static TaylorModel<P,FLT> _compose(TaylorModel<P,FLT> const& tf, VectorUnscaling const& u, const Vector<TaylorModel<P,FLT>>& tg);
    static TaylorModel<P,FLT> _partial_evaluate(const TaylorModel<P,FLT>& x, SizeType k, NumericType c);
    static Covector<X> _gradient(const TaylorModel<P,FLT>& x, Vector<X> const& v);
    static ArithmeticType<CoefficientType,IntervalNumericType> _evaluate(const TaylorModel<P,FLT>& x, Vector<IntervalNumericType> const& v);
    static ArithmeticType<CoefficientType,ValidatedNumericType> _evaluate(const TaylorModel<P,FLT>& x, Vector<ValidatedNumericType> const& v);
    static ArithmeticType<CoefficientType,ApproximateNumericType> _evaluate(const TaylorModel<P,FLT>& x, Vector<ApproximateNumericType> const& v);
};

// FIXME: Needed to dispatch gradient of ScaledFunctionPatch
template<class P, class FLT> template<class A> auto TaylorModel<P,FLT>::operator() (Vector<A> const& x) const -> ArithmeticType<CoefficientType,A> {
    return horner_evaluate(this->expansion(),x)+FloatBounds<typename FLT::PrecisionType>(-this->error(),+this->error());
}



template<class FLT> Bool inconsistent(const Vector<TaylorModel<ValidatedTag,FLT>>& tm1, const Vector<TaylorModel<ValidatedTag,FLT>>& tm2) {
    ARIADNE_ASSERT(tm1.size()==tm2.size());
    for(SizeType i=0; i!=tm1.size(); ++i) { if(inconsistent(tm1[i],tm2[i])) { return true; } }
    return false;
}

template<class FLT> Bool refines(const Vector<TaylorModel<ValidatedTag,FLT>>& tm1, const Vector<TaylorModel<ValidatedTag,FLT>>& tm2) {
    ARIADNE_ASSERT(tm1.size()==tm2.size());
    for(SizeType i=0; i!=tm1.size(); ++i) { if(not refines(tm1[i],tm2[i])) { return false; } }
    return true;
}

template<class FLT> Vector<TaylorModel<ValidatedTag,FLT>> refinement(const Vector<TaylorModel<ValidatedTag,FLT>>& tm1, const Vector<TaylorModel<ValidatedTag,FLT>>& tm2) {
    ARIADNE_ASSERT(tm1.size()==tm2.size());
    Vector<TaylorModel<ValidatedTag,FLT>> r(tm1.size());
    for(SizeType i=0; i!=tm1.size(); ++i) { r[i]=refinement(tm1[i],tm2[i]); }
    return r;
}





template<class P, class FLT, class X> decltype(auto) evaluate(const Vector<TaylorModel<P,FLT>>& tf, const Vector<X>& x)
{
    typedef decltype(evaluate(tf[0],x)) R;
    return Vector<R>(tf.size(),[&](SizeType i){return evaluate(tf[i],x);});
}

template<class FLT> Vector<TaylorModel<ValidatedTag,FLT>> partial_evaluate(const Vector<TaylorModel<ValidatedTag,FLT>>& tf, SizeType k, const FloatBounds<PrecisionType<FLT>>& c) {
    Vector<TaylorModel<ValidatedTag,FLT>> r(tf.size(),ValidatedTaylorModel<FLT>::zero(tf.zero_element().argument_size()-1u,tf.zero_element().sweeper()));
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=partial_evaluate(tf[i],k,c); }
    return r;
}

template<class FLT> Vector<TaylorModel<ValidatedTag,FLT>> compose(const Vector<TaylorModel<ValidatedTag,FLT>>& tf, const Vector<TaylorModel<ValidatedTag,FLT>>& tg) {
    Vector<TaylorModel<ValidatedTag,FLT>> r(tf.size(),tg.zero_element());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=compose(tf[i],tg); }
    return r;
}

template<class P, class FLT> Vector<TaylorModel<P,FLT>> unscale(const Vector<TaylorModel<P,FLT>>& x, const Vector<IntervalDomainType>& dom) {
    Vector<TaylorModel<P,FLT>> r(x.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=unscale(x[i],dom[i]); }
    return r;
}

template<class P, class FLT> Vector<TaylorModel<P,FLT>> compose(const VectorUnscaling& u, const Vector<TaylorModel<P,FLT>>& g) {
    ARIADNE_ASSERT_MSG(u.size()==g.size(),"u="<<u.domain()<<", g="<<g);
    Vector<TaylorModel<P,FLT>> r(u.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=compose(u[i],g[i]); }
    return r;
}

template<class P, class FLT> Vector<TaylorModel<P,FLT>> compose(const Vector<TaylorModel<P,FLT>>& x,
                                                            const VectorUnscaling& u,
                                                            const Vector<TaylorModel<P,FLT>>& y) {
    return compose(x,compose(u,y));
}

template<class P, class FLT> Vector<TaylorModel<P,FLT>> antiderivative(const Vector<TaylorModel<P,FLT>>& x, SizeType k) {
    Vector<TaylorModel<P,FLT>> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) { r[i]=antiderivative(x[i],k); }
    return r;
}

template<class P, class FLT> Vector<TaylorModel<P,FLT>> derivative(const Vector<TaylorModel<P,FLT>>& x, SizeType k) {
    Vector<TaylorModel<P,FLT>> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) { r[i]=derivative(x[i],k); }
    return r;
}

template<class P, class FLT> Vector<TaylorModel<P,FLT>> combine(const Vector<TaylorModel<P,FLT>>& x1, const Vector<TaylorModel<P,FLT>>& x2) {
    return join(embed(0u,x1,x2.zero_element().argument_size()),embed(x1.zero_element().argument_size(),x2,0u));
}

template<class FLT> Vector<TaylorModel<ValidatedTag,FLT>> embed(SizeType as1, const Vector<TaylorModel<ValidatedTag,FLT>>& x2, SizeType as3) {
    Vector<TaylorModel<ValidatedTag,FLT>> r(x2.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=embed(as1,x2[i],as3); }
    return r;
}

template<class P, class FLT> Vector<TaylorModel<P,FLT>> split(const Vector<TaylorModel<P,FLT>>& x, SizeType j, SplitPart h) {
    Vector<TaylorModel<P,FLT>> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) { r[i]=split(x[i],j,h); }
    return r;
}

template<class P, class FLT> Vector<typename TaylorModel<P,FLT>::RangeType> ranges(const Vector<TaylorModel<P,FLT>>& f) {
    Vector<typename TaylorModel<P,FLT>::RangeType> r(f.size()); for(SizeType i=0; i!=f.size(); ++i) { r[i]=f[i].range(); } return r;
}

template<class P, class FLT> Vector<typename TaylorModel<P,FLT>::ErrorType> errors(const Vector<TaylorModel<P,FLT>>& h) {
    Vector<typename TaylorModel<P,FLT>::ErrorType> e(h.size()); for(SizeType i=0; i!=h.size(); ++i) { e[i]=h[i].error(); } return e; }

template<class P, class FLT> Vector<typename TaylorModel<P,FLT>::NormType> norms(const Vector<TaylorModel<P,FLT>>& h) {
    Vector<typename TaylorModel<P,FLT>::NormType> r(h.size()); for(SizeType i=0; i!=h.size(); ++i) { r[i]=norm(h[i]); } return r; }

template<class P, class FLT> typename TaylorModel<P,FLT>::NormType norm(const Vector<TaylorModel<P,FLT>>& h) {
    typename TaylorModel<P,FLT>::NormType r=0u; for(SizeType i=0; i!=h.size(); ++i) { r=max(r,norm(h[i])); } return r;
}

template<class FLT> Matrix<Bounds<FLT>> jacobian(const Vector<TaylorModel<ValidatedTag,FLT>>& x, const Vector<Bounds<FLT>>& y);
template<class FLT> Matrix<Bounds<FLT>> jacobian(const Vector<TaylorModel<ValidatedTag,FLT>>& x, const Vector<Bounds<FLT>>& y, Array<SizeType>& p);
template<class FLT> Matrix<FLT> jacobian_value(const Vector<TaylorModel<ValidatedTag,FLT>>& x);
template<class FLT> Matrix<FLT> jacobian_value(const Vector<TaylorModel<ValidatedTag,FLT>>& x, const Array<SizeType>& p);
template<class FLT> Matrix<UpperInterval<FLT>> jacobian_range(const Vector<TaylorModel<ValidatedTag,FLT>>& x);
template<class FLT> Matrix<UpperInterval<FLT>> jacobian_range(const Vector<TaylorModel<ValidatedTag,FLT>>& x, const Array<SizeType>& p);

template<class FLT> TaylorModel<ValidatedTag,FLT> value_coefficients(TaylorModel<ValidatedTag,Bounds<FLT>> const& tm);
template<class FLT> TaylorModel<ValidatedTag,Bounds<FLT>> exact_coefficients(TaylorModel<ValidatedTag,Bounds<FLT>> const& tm);


} // namespace Ariadne

#endif // ARIADNE_TAYLOR_MODEL_HPP
