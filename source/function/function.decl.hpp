/***************************************************************************
 *            function.decl.hpp
 *
 *  Copyright 2016-17  Pieter Collins
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

// Domain declarations
using IntervalDomainType = ExactIntervalType;
using BoxDomainType = ExactBoxType;

// Function declarations
template<class P, class D, class C> class Function;
template<class P, class D> using ScalarFunction = Function<P,D,IntervalDomainType>;
template<class P, class D> using VectorFunction = Function<P,D,BoxDomainType>;
template<class P, class C> using UnivariateFunction = Function<P,IntervalDomainType,C>;
template<class P, class C> using MultivariateFunction = Function<P,BoxDomainType,C>;

template<class P> using ScalarUnivariateFunction = Function<P,IntervalDomainType,IntervalDomainType>;
template<class P> using VectorUnivariateFunction = Function<P,IntervalDomainType,BoxDomainType>;
template<class P> using ScalarMultivariateFunction = Function<P,BoxDomainType,IntervalDomainType>;
template<class P> using VectorMultivariateFunction = Function<P,BoxDomainType,BoxDomainType>;

typedef ScalarUnivariateFunction<ApproximateTag> ApproximateScalarUnivariateFunction;
typedef ScalarUnivariateFunction<ValidatedTag> ValidatedScalarUnivariateFunction;
typedef ScalarUnivariateFunction<EffectiveTag> EffectiveScalarUnivariateFunction;
typedef EffectiveScalarUnivariateFunction RealScalarUnivariateFunction;

typedef ScalarMultivariateFunction<ApproximateTag> ApproximateScalarMultivariateFunction;
typedef ScalarMultivariateFunction<ValidatedTag> ValidatedScalarMultivariateFunction;
typedef ScalarMultivariateFunction<EffectiveTag> EffectiveScalarMultivariateFunction;
typedef EffectiveScalarMultivariateFunction RealScalarMultivariateFunction;

typedef VectorUnivariateFunction<ApproximateTag> ApproximateVectorUnivariateFunction;
typedef VectorUnivariateFunction<ValidatedTag> ValidatedVectorUnivariateFunction;
typedef VectorUnivariateFunction<EffectiveTag> EffectiveVectorUnivariateFunction;
typedef EffectiveVectorUnivariateFunction RealVectorUnivariateFunction;

typedef VectorMultivariateFunction<ApproximateTag> ApproximateVectorMultivariateFunction;
typedef VectorMultivariateFunction<ValidatedTag> ValidatedVectorMultivariateFunction;
typedef VectorMultivariateFunction<EffectiveTag> EffectiveVectorMultivariateFunction;
typedef EffectiveVectorMultivariateFunction RealVectorMultivariateFunction;


// Function interface declarations
template<class P, class D, class C> class FunctionInterface;
template<class P, class D> using ScalarFunctionInterface = FunctionInterface<P,D,IntervalDomainType>;
template<class P, class D> using VectorFunctionInterface = FunctionInterface<P,D,BoxDomainType>;
template<class P> using ScalarMultivariateFunctionInterface = FunctionInterface<P,BoxDomainType,IntervalDomainType>;
template<class P> using VectorMultivariateFunctionInterface = FunctionInterface<P,BoxDomainType,BoxDomainType>;

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

template<class P, class D, class C, class PR, class PRE=PR> class FunctionModel;
template<class P, class D, class PR, class PRE=PR> using ScalarFunctionModel = FunctionModel<P,D,IntervalDomainType,PR,PRE>;
template<class P, class D, class PR, class PRE=PR> using VectorFunctionModel = FunctionModel<P,D,BoxDomainType,PR,PRE>;
template<class P, class PR, class PRE=PR> using ScalarUnivariateFunctionModel = FunctionModel<P,IntervalDomainType,IntervalDomainType,PR,PRE>;
template<class P, class PR, class PRE=PR> using VectorUnivariateFunctionModel = FunctionModel<P,IntervalDomainType,BoxDomainType,PR,PRE>;
template<class P, class PR, class PRE=PR> using ScalarMultivariateFunctionModel = FunctionModel<P,BoxDomainType,IntervalDomainType,PR,PRE>;
template<class P, class PR, class PRE=PR> using VectorMultivariateFunctionModel = FunctionModel<P,BoxDomainType,BoxDomainType,PR,PRE>;

template<class P, class PR, class PRE=PR> struct FunctionModelTraits;


template<class PR> struct FunctionModelTraits<ApproximateTag,PR> {
    static_assert(IsSame<PR,DP>::value or IsSame<PR,MP>::value,"");
    typedef FloatApproximation<PR> CoefficientType;
    typedef Void ErrorType;
    typedef FloatApproximation<PR> NumericType;
};

template<class PR, class PRE> struct FunctionModelTraits<ValidatedTag,PR,PRE> {
    static_assert(IsSame<PR,DP>::value or IsSame<PR,MP>::value,"");
    typedef FloatValue<PR> CoefficientType;
    typedef FloatError<PRE> ErrorType;
    typedef FloatBounds<PR> NumericType;
};

template<class P, class PR, class PRE=PR> using CanonicalNumericType = typename FunctionModelTraits<P,PR,PRE>::NumericType;
template<class P, class PR> using CanonicalCoefficientType = typename FunctionModelTraits<P,PR>::CoefficientType;
template<class P, class PRE> using CanonicalErrorType = typename FunctionModelTraits<P,PRE,PRE>::ErrorType;

template<class P> using CanonicalNumeric64Type = typename FunctionModelTraits<P,DoublePrecision>::NumericType;
template<class P> using CanonicalCoefficient64Type = typename FunctionModelTraits<P,DoublePrecision>::CoefficientType;
template<class P> using CanonicalError64Type = typename FunctionModelTraits<P,DoublePrecision,DoublePrecision>::ErrorType;

template<class X> using PrecisionType = typename X::PrecisionType;
template<class X> using ErrorPrecisionType = typename X::ErrorPrecisionType;

template<class P> using ScalarMultivariateFunctionModelDP = ScalarFunctionModel<P,BoxDomainType,DoublePrecision>;
template<class P> using VectorMultivariateFunctionModelDP = VectorFunctionModel<P,BoxDomainType,DoublePrecision>;

template<class P> using ScalarMultivariateFunctionModelDPInterface = ScalarFunctionModelInterface<P,BoxDomainType,DoublePrecision>;
template<class P> using VectorMultivariateFunctionModelDPInterface = VectorFunctionModelInterface<P,BoxDomainType,DoublePrecision>;

using ValidatedScalarMultivariateFunctionModelDPInterface = ScalarFunctionModelInterface<ValidatedTag,BoxDomainType,DoublePrecision>;
using ValidatedVectorMultivariateFunctionModelDPInterface = VectorFunctionModelInterface<ValidatedTag,BoxDomainType,DoublePrecision>;

using ValidatedScalarMultivariateFunctionModelDP = ScalarFunctionModel<ValidatedTag,BoxDomainType,DoublePrecision>;
using ValidatedVectorMultivariateFunctionModelDP = VectorFunctionModel<ValidatedTag,BoxDomainType,DoublePrecision>;

using ApproximateScalarMultivariateFunctionModelDPInterface = ScalarFunctionModelInterface<ApproximateTag,BoxDomainType,DoublePrecision>;
using ApproximateVectorMultivariateFunctionModelDPInterface = VectorFunctionModelInterface<ApproximateTag,BoxDomainType,DoublePrecision>;

using ApproximateScalarMultivariateFunctionModelDP = ScalarFunctionModel<ApproximateTag,BoxDomainType,DoublePrecision>;
using ApproximateVectorMultivariateFunctionModelDP = VectorFunctionModel<ApproximateTag,BoxDomainType,DoublePrecision>;

template<class P, class PR=DoublePrecision, class PRE=PR> class FunctionModelFactoryInterface;
typedef FunctionModelFactoryInterface<ValidatedTag,DoublePrecision> ValidatedFunctionModelDPFactoryInterface;
template<class P, class PR=DoublePrecision, class PRE=PR> class FunctionModelFactory;
typedef FunctionModelFactory<ValidatedTag,DoublePrecision> ValidatedFunctionModelDPFactory;
template<class FMF, class D> class FunctionModelCreator;



} // namespace Ariadne

#endif
