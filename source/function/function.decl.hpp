/***************************************************************************
 *            function.decl.hpp
 *
 *  Copyright  2016-20  Pieter Collins
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

/*! \file function.decl.hpp
 *  \brief Declarations of function types.
 */

#ifndef ARIADNE_FUNCTION_DECL_HPP
#define ARIADNE_FUNCTION_DECL_HPP

#include "../geometry/interval.decl.hpp"
#include "../geometry/box.decl.hpp"

namespace Ariadne {

template<class X> class Vector;

// Domain declarations
typedef Interval<FloatDPValue> IntervalDomainType;
typedef Box<Interval<FloatDPValue>> BoxDomainType;

// Function declarations
template<class P, class D, class C> class Function;

//@{
//! \name Type shorthands for function classes
template<class P, class D> using ScalarFunction = Function<P,D,IntervalDomainType>; //!< . \ingroup FunctionModule
template<class P, class D> using VectorFunction = Function<P,D,BoxDomainType>; //!< . \ingroup FunctionModule
template<class P, class C> using UnivariateFunction = Function<P,IntervalDomainType,C>; //!< . \ingroup FunctionModule
template<class P, class C> using MultivariateFunction = Function<P,BoxDomainType,C>; //!< . \ingroup FunctionModule

template<class P> using ScalarUnivariateFunction = Function<P,IntervalDomainType,IntervalDomainType>; //!< . \ingroup FunctionModule
template<class P> using VectorUnivariateFunction = Function<P,IntervalDomainType,BoxDomainType>; //!< . \ingroup FunctionModule
template<class P> using ScalarMultivariateFunction = Function<P,BoxDomainType,IntervalDomainType>; //!< . \ingroup FunctionModule
template<class P> using VectorMultivariateFunction = Function<P,BoxDomainType,BoxDomainType>; //!< . \ingroup FunctionModule

template<class D, class C> using ApproximateFunction = Function<ApproximateTag,D,C>; //!< . \ingroup FunctionModule
template<class D, class C> using ValidatedFunction = Function<ValidatedTag,D,C>; //!< . \ingroup FunctionModule
template<class D, class C> using EffectiveFunction = Function<EffectiveTag,D,C>; //!< . \ingroup FunctionModule
//@}

//@{
//! \name Type synonyms for function classes
using ApproximateScalarUnivariateFunction = ScalarUnivariateFunction<ApproximateTag>; //!< . \ingroup FunctionModule
using ValidatedScalarUnivariateFunction = ScalarUnivariateFunction<ValidatedTag>; //!< . \ingroup FunctionModule
using EffectiveScalarUnivariateFunction = ScalarUnivariateFunction<EffectiveTag>; //!< . \ingroup FunctionModule

using ApproximateScalarMultivariateFunction = ScalarMultivariateFunction<ApproximateTag>; //!< . \ingroup FunctionModule
using ValidatedScalarMultivariateFunction = ScalarMultivariateFunction<ValidatedTag>; //!< . \ingroup FunctionModule
using EffectiveScalarMultivariateFunction = ScalarMultivariateFunction<EffectiveTag>; //!< . \ingroup FunctionModule

using ApproximateVectorUnivariateFunction = VectorUnivariateFunction<ApproximateTag>; //!< . \ingroup FunctionModule
using ValidatedVectorUnivariateFunction = VectorUnivariateFunction<ValidatedTag>; //!< . \ingroup FunctionModule
using EffectiveVectorUnivariateFunction = VectorUnivariateFunction<EffectiveTag>; //!< . \ingroup FunctionModule

using ApproximateVectorMultivariateFunction = VectorMultivariateFunction<ApproximateTag>; //!< . \ingroup FunctionModule
using ValidatedVectorMultivariateFunction = VectorMultivariateFunction<ValidatedTag>; //!< . \ingroup FunctionModule
using EffectiveVectorMultivariateFunction = VectorMultivariateFunction<EffectiveTag>; //!< . \ingroup FunctionModule

using RealScalarUnivariateFunction = EffectiveScalarUnivariateFunction; //!< DEPRECATED \ingroup FunctionModule
using RealScalarMultivariateFunction = EffectiveScalarMultivariateFunction; //!< .
using RealVectorUnivariateFunction = EffectiveVectorUnivariateFunction; //!< .
using RealVectorMultivariateFunction = EffectiveVectorMultivariateFunction; //!< .
//@}

//! \ingroup FunctionModule
//! \brief An interface for generic functions.
//! \see Function
template<class P, class D, class C> class FunctionInterface;

//@{
//! \relates FunctionInterface
//! \name Type shorthands
template<class P, class D> using ScalarFunctionInterface = FunctionInterface<P,D,IntervalDomainType>; //!< .
template<class P, class D> using VectorFunctionInterface = FunctionInterface<P,D,BoxDomainType>; //!< .
template<class P, class C> using UnivariateFunctionInterface = FunctionInterface<P,IntervalDomainType,C>; //!< .
template<class P, class C> using MultivariateFunctionInterface = FunctionInterface<P,BoxDomainType,C>; //!< .
template<class P> using ScalarUnivariateFunctionInterface = FunctionInterface<P,IntervalDomainType,IntervalDomainType>; //!< .
template<class P> using VectorUnivariateFunctionInterface = FunctionInterface<P,IntervalDomainType,BoxDomainType>; //!< .
template<class P> using ScalarMultivariateFunctionInterface = FunctionInterface<P,BoxDomainType,IntervalDomainType>; //!< .
template<class P> using VectorMultivariateFunctionInterface = FunctionInterface<P,BoxDomainType,BoxDomainType>; //!< .
//@}

typedef ScalarUnivariateFunctionInterface<ApproximateTag> ApproximateScalarUnivariateFunctionInterface;
typedef ScalarUnivariateFunctionInterface<ValidatedTag> ValidatedScalarUnivariateFunctionInterface;
typedef ScalarUnivariateFunctionInterface<EffectiveTag> EffectiveScalarUnivariateFunctionInterface;

typedef VectorUnivariateFunctionInterface<ApproximateTag> ApproximateVectorUnivariateFunctionInterface;
typedef VectorUnivariateFunctionInterface<ValidatedTag> ValidatedVectorUnivariateFunctionInterface;
typedef VectorUnivariateFunctionInterface<EffectiveTag> EffectiveVectorUnivariateFunctionInterface;

typedef ScalarMultivariateFunctionInterface<ApproximateTag> ApproximateScalarMultivariateFunctionInterface;
typedef ScalarMultivariateFunctionInterface<ValidatedTag> ValidatedScalarMultivariateFunctionInterface;
typedef ScalarMultivariateFunctionInterface<EffectiveTag> EffectiveScalarMultivariateFunctionInterface;

typedef VectorMultivariateFunctionInterface<ApproximateTag> ApproximateVectorMultivariateFunctionInterface;
typedef VectorMultivariateFunctionInterface<ValidatedTag> ValidatedVectorMultivariateFunctionInterface;
typedef VectorMultivariateFunctionInterface<EffectiveTag> EffectiveVectorMultivariateFunctionInterface;

using ValidatedScalarUnivariateFunctionPatch = Function<ValidatedTag,IntervalDomainType,IntervalDomainType>;
using ValidatedVectorUnivariateFunctionPatch = Function<ValidatedTag,IntervalDomainType,BoxDomainType>;
using ValidatedScalarMultivariateFunctionPatch = Function<ValidatedTag,BoxDomainType,IntervalDomainType>;
using ValidatedVectorMultivariateFunctionPatch = Function<ValidatedTag,BoxDomainType,BoxDomainType>;


// Function models declarations


template<class P, class D, class C, class PR, class PRE=PR> class FunctionModelInterface;
template<class P, class D, class PR, class PRE=PR> using ScalarFunctionModelInterface = FunctionModelInterface<P,D,IntervalDomainType,PR,PRE>;
template<class P, class D, class PR, class PRE=PR> using VectorFunctionModelInterface = FunctionModelInterface<P,D,BoxDomainType,PR,PRE>;
template<class P, class C, class PR, class PRE=PR> using UnivariateFunctionModelInterface = FunctionModelInterface<P,IntervalDomainType,C,PR,PRE>;
template<class P, class C, class PR, class PRE=PR> using MultivariateFunctionModelInterface = FunctionModelInterface<P,BoxDomainType,C,PR,PRE>;
template<class P, class PR, class PRE=PR> using ScalarUnivariateFunctionModelInterface
    = FunctionModelInterface<P,IntervalDomainType,IntervalDomainType,PR,PRE>;
template<class P, class PR, class PRE=PR> using VectorUnivariateFunctionModelInterface
    = FunctionModelInterface<P,IntervalDomainType,BoxDomainType,PR,PRE>;
template<class P, class PR, class PRE=PR> using ScalarMultivariateFunctionModelInterface
    = FunctionModelInterface<P,BoxDomainType,IntervalDomainType,PR,PRE>;
template<class P, class PR, class PRE=PR> using VectorMultivariateFunctionModelInterface
    = FunctionModelInterface<P,BoxDomainType,BoxDomainType,PR,PRE>;

template<class P, class D, class C, class PR, class PRE=PR> class FunctionModel;

//@{
//! \relates FunctionModel
//! \name Type shorthands
template<class P, class D, class PR, class PRE=PR> using ScalarFunctionModel = FunctionModel<P,D,IntervalDomainType,PR,PRE>; //!< .
template<class P, class D, class PR, class PRE=PR> using VectorFunctionModel = FunctionModel<P,D,BoxDomainType,PR,PRE>; //!< .
template<class P, class C, class PR, class PRE=PR> using UnivariateFunctionModel = FunctionModel<P,IntervalDomainType,C,PR,PRE>; //!< .
template<class P, class C, class PR, class PRE=PR> using MultivariateFunctionModel = FunctionModel<P,BoxDomainType,C,PR,PRE>; //!< .
template<class P, class PR, class PRE=PR> using ScalarUnivariateFunctionModel = FunctionModel<P,IntervalDomainType,IntervalDomainType,PR,PRE>; //!< .
template<class P, class PR, class PRE=PR> using VectorUnivariateFunctionModel = FunctionModel<P,IntervalDomainType,BoxDomainType,PR,PRE>; //!< .
template<class P, class PR, class PRE=PR> using ScalarMultivariateFunctionModel = FunctionModel<P,BoxDomainType,IntervalDomainType,PR,PRE>; //!< .
template<class P, class PR, class PRE=PR> using VectorMultivariateFunctionModel = FunctionModel<P,BoxDomainType,BoxDomainType,PR,PRE>; //!< .

template<class PR, class PRE=PR> using ValidatedScalarUnivariateFunctionModel = ScalarUnivariateFunctionModel<ValidatedTag,PR,PRE>; //!< .
template<class PR, class PRE=PR> using ValidatedVectorUnivariateFunctionModel = VectorUnivariateFunctionModel<ValidatedTag,PR,PRE>; //!< .
template<class PR, class PRE=PR> using ValidatedScalarMultivariateFunctionModel = ScalarMultivariateFunctionModel<ValidatedTag,PR,PRE>; //!< .
template<class PR, class PRE=PR> using ValidatedVectorMultivariateFunctionModel = VectorMultivariateFunctionModel<ValidatedTag,PR,PRE>; //!< .

template<class PR, class PRE=PR> using ApproximateScalarUnivariateFunctionModel = ScalarUnivariateFunctionModel<ApproximateTag,PR,PRE>; //!< .
template<class PR, class PRE=PR> using ApproximateVectorUnivariateFunctionModel = VectorUnivariateFunctionModel<ApproximateTag,PR,PRE>; //!< .
template<class PR, class PRE=PR> using ApproximateScalarMultivariateFunctionModel = ScalarMultivariateFunctionModel<ApproximateTag,PR,PRE>; //!< .
template<class PR, class PRE=PR> using ApproximateVectorMultivariateFunctionModel = VectorMultivariateFunctionModel<ApproximateTag,PR,PRE>; //!< .

template<class P> using ScalarUnivariateFunctionModelDP = ScalarUnivariateFunctionModel<P,DoublePrecision>; //!< .
template<class P> using VectorUnivariateFunctionModelDP = VectorUnivariateFunctionModel<P,DoublePrecision>; //!< .
template<class P> using ScalarMultivariateFunctionModelDP = ScalarMultivariateFunctionModel<P,DoublePrecision>; //!< .
template<class P> using VectorMultivariateFunctionModelDP = VectorMultivariateFunctionModel<P,DoublePrecision>; //!< .

template<class P> using ScalarUnivariateFunctionModelMP = ScalarUnivariateFunctionModel<P,MultiplePrecision>; //!< .
template<class P> using VectorUnivariateFunctionModelMP = VectorUnivariateFunctionModel<P,MultiplePrecision>; //!< .
template<class P> using ScalarMultivariateFunctionModelMP = ScalarMultivariateFunctionModel<P,MultiplePrecision>; //!< .
template<class P> using VectorMultivariateFunctionModelMP = VectorMultivariateFunctionModel<P,MultiplePrecision>; //!< .
//@}

//@{
//! \relates FunctionModel
//! \name Type synonyms
using ValidatedScalarUnivariateFunctionModelDP = ScalarUnivariateFunctionModel<ValidatedTag,DoublePrecision>; //!< .
using ValidatedVectorUnivariateFunctionModelDP = VectorUnivariateFunctionModel<ValidatedTag,DoublePrecision>; //!< .
using ValidatedScalarMultivariateFunctionModelDP = ScalarMultivariateFunctionModel<ValidatedTag,DoublePrecision>; //!< .
using ValidatedVectorMultivariateFunctionModelDP = VectorMultivariateFunctionModel<ValidatedTag,DoublePrecision>; //!< .

using ApproximateScalarUnivariateFunctionModelDP = ScalarUnivariateFunctionModel<ApproximateTag,DoublePrecision>; //!< .
using ApproximateVectorUnivariateFunctionModelDP = VectorUnivariateFunctionModel<ApproximateTag,DoublePrecision>; //!< .
using ApproximateScalarMultivariateFunctionModelDP = ScalarMultivariateFunctionModel<ApproximateTag,DoublePrecision>; //!< .
using ApproximateVectorMultivariateFunctionModelDP = VectorMultivariateFunctionModel<ApproximateTag,DoublePrecision>; //!< .

using ValidatedScalarUnivariateFunctionModelMP = ScalarUnivariateFunctionModel<ValidatedTag,MultiplePrecision>; //!< .
using ValidatedVectorUnivariateFunctionModelMP = VectorUnivariateFunctionModel<ValidatedTag,MultiplePrecision>; //!< .
using ValidatedScalarMultivariateFunctionModelMP = ScalarMultivariateFunctionModel<ValidatedTag,MultiplePrecision>; //!< .
using ValidatedVectorMultivariateFunctionModelMP = VectorMultivariateFunctionModel<ValidatedTag,MultiplePrecision>; //!< .

using ApproximateScalarUnivariateFunctionModelMP = ScalarUnivariateFunctionModel<ApproximateTag,MultiplePrecision>; //!< .
using ApproximateVectorUnivariateFunctionModelMP = VectorUnivariateFunctionModel<ApproximateTag,MultiplePrecision>; //!< .
using ApproximateScalarMultivariateFunctionModelMP = ScalarMultivariateFunctionModel<ApproximateTag,MultiplePrecision>; //!< .
using ApproximateVectorMultivariateFunctionModelMP = VectorMultivariateFunctionModel<ApproximateTag,MultiplePrecision>; //!< .
//@}

template<class P, class PR, class PRE=PR> struct FunctionModelTraits;

template<class F> class UnknownError;

template<class PR, class PRE> struct FunctionModelTraits<ValidatedTag,PR,PRE> {
    typedef RawFloat<PR> F; typedef RawFloat<PRE> FE;
    typedef Value<F> ValueType; typedef Error<FE> ErrorType;
    typedef PositiveUpperBound<F> NormType; typedef Interval<UpperBound<F>> RangeType;
    typedef Bounds<F> NumericType; typedef ValidatedNumber GenericNumericType;
    typedef F RawFloatType;
};
template<class PR> struct FunctionModelTraits<ApproximateTag,PR> {
    typedef RawFloat<PR> F;
    typedef Approximation<F> ValueType; typedef UnknownError<F> ErrorType;
    typedef PositiveApproximation<F> NormType; typedef Interval<Approximation<F>> RangeType;
    typedef Approximation<F> NumericType; typedef ApproximateNumber GenericNumericType;
    typedef F RawFloatType;
};

template<class P, class PR, class PRE=PR> using CanonicalNumericType = typename FunctionModelTraits<P,PR,PRE>::NumericType;
template<class P, class PR> using CanonicalCoefficientType = typename FunctionModelTraits<P,PR>::CoefficientType;
template<class P, class PRE> using CanonicalErrorType = typename FunctionModelTraits<P,PRE,PRE>::ErrorType;

template<class P> using CanonicalNumeric64Type = typename FunctionModelTraits<P,DoublePrecision>::NumericType;
template<class P> using CanonicalCoefficient64Type = typename FunctionModelTraits<P,DoublePrecision>::CoefficientType;
template<class P> using CanonicalError64Type = typename FunctionModelTraits<P,DoublePrecision,DoublePrecision>::ErrorType;

template<class X> using PrecisionType = typename X::PrecisionType;
template<class X> using ErrorPrecisionType = typename X::ErrorPrecisionType;


template<class P> using ScalarUnivariateFunctionModelDPInterface = ScalarUnivariateFunctionModelInterface<P,DoublePrecision>;
template<class P> using VectorUnivariateFunctionModelDPInterface = VectorUnivariateFunctionModelInterface<P,DoublePrecision>;
template<class P> using ScalarMultivariateFunctionModelDPInterface = ScalarMultivariateFunctionModelInterface<P,DoublePrecision>;
template<class P> using VectorMultivariateFunctionModelDPInterface = VectorMultivariateFunctionModelInterface<P,DoublePrecision>;

using ValidatedScalarUnivariateFunctionModelDPInterface = ScalarUnivariateFunctionModelInterface<ValidatedTag,DoublePrecision>;
using ValidatedVectorUnivariateFunctionModelDPInterface = VectorUnivariateFunctionModelInterface<ValidatedTag,DoublePrecision>;
using ValidatedScalarMultivariateFunctionModelDPInterface = ScalarMultivariateFunctionModelInterface<ValidatedTag,DoublePrecision>;
using ValidatedVectorMultivariateFunctionModelDPInterface = VectorMultivariateFunctionModelInterface<ValidatedTag,DoublePrecision>;

using ApproximateScalarUnivariateFunctionModelDPInterface = ScalarUnivariateFunctionModelInterface<ApproximateTag,DoublePrecision>;
using ApproximateVectorUnivariateFunctionModelDPInterface = VectorUnivariateFunctionModelInterface<ApproximateTag,DoublePrecision>;
using ApproximateScalarMultivariateFunctionModelDPInterface = ScalarMultivariateFunctionModelInterface<ApproximateTag,DoublePrecision>;
using ApproximateVectorMultivariateFunctionModelDPInterface = VectorMultivariateFunctionModelInterface<ApproximateTag,DoublePrecision>;


template<class P, class PR, class PRE=PR> class FunctionModelFactoryInterface;
typedef FunctionModelFactoryInterface<ValidatedTag,DoublePrecision> ValidatedFunctionModelDPFactoryInterface;
template<class P, class PR, class PRE=PR> class FunctionModelFactory;
typedef FunctionModelFactory<ValidatedTag,DoublePrecision> ValidatedFunctionModelDPFactory;
template<class FMF, class D> class FunctionModelCreator;

class UniIndex;
class MultiIndex;

template<class I, class X> class Monomial;
template<class I, class X> class Polynomial;

//@{
//! \name Type shorthands and synonyms for Monomial classes.
template<class X> using UnivariateMonomial = Monomial<UniIndex,X>; //!< .
template<class X> using MultivariateMonomial = Monomial<MultiIndex,X>; //!< .

using FloatDPApproximationMultivariateMonomial = MultivariateMonomial<FloatDPApproximation>; //!< .
using FloatDPBoundsMultivariateMonomial = MultivariateMonomial<FloatDPBounds>; //!< .
using FloatMPApproximationMultivariateMonomial = MultivariateMonomial<FloatMPApproximation>; //!< .
using FloatMPBoundsMultivariateMonomial = MultivariateMonomial<FloatMPBounds>; //!< .
//@}

//@{
//! \name Type shorthands and synonyms for Polynomial classes
template<class X> using UnivariatePolynomial = Polynomial<UniIndex,X>; //!< . \ingroup FunctionModule
template<class X> using MultivariatePolynomial = Polynomial<MultiIndex,X>; //!< . \ingroup FunctionModule

using FloatDPApproximationMultivariatePolynomial = MultivariatePolynomial<FloatDPApproximation>; //!< . \ingroup FunctionModule
using FloatDPBoundsMultivariatePolynomial = MultivariatePolynomial<FloatDPBounds>; //!< . \ingroup FunctionModule
using FloatMPApproximationMultivariatePolynomial = MultivariatePolynomial<FloatMPApproximation>; //!< . \ingroup FunctionModule
using FloatMPBoundsMultivariatePolynomial = MultivariatePolynomial<FloatMPBounds>; //!< . \ingroup FunctionModule

using FloatDPApproximationMultivariatePolynomialVector = Vector<MultivariatePolynomial<FloatDPApproximation>>; //!< . \ingroup FunctionModule
using FloatDPBoundsMultivariatePolynomialVector = Vector<MultivariatePolynomial<FloatDPBounds>>; //!< . \ingroup FunctionModule
using FloatMPApproximationMultivariatePolynomialVector = Vector<MultivariatePolynomial<FloatMPApproximation>>; //!< . \ingroup FunctionModule
using FloatMPBoundsMultivariatePolynomialVector = Vector<MultivariatePolynomial<FloatMPBounds>>; //!< . \ingroup FunctionModule
//@}


} // namespace Ariadne

#endif
