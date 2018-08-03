/***************************************************************************
 *            function.decl.hpp
 *
 *  Copyright 2016-17  Pieter Collins
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
template<class P, class D=BoxDomainType> using ScalarFunction = Function<P,D,IntervalDomainType>;
template<class P, class D=BoxDomainType> using VectorFunction = Function<P,D,BoxDomainType>;
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

typedef ScalarFunction<ApproximateTag> ApproximateScalarFunction;
typedef ScalarFunction<ValidatedTag> ValidatedScalarFunction;
typedef ScalarFunction<EffectiveTag> EffectiveScalarFunction;
typedef EffectiveScalarFunction RealScalarFunction;

typedef VectorFunction<ApproximateTag> ApproximateVectorFunction;
typedef VectorFunction<ValidatedTag> ValidatedVectorFunction;
typedef VectorFunction<EffectiveTag> EffectiveVectorFunction;
typedef EffectiveVectorFunction RealVectorFunction;

typedef VectorUnivariateFunction<ApproximateTag> ApproximateVectorUnivariateFunction;
typedef VectorUnivariateFunction<ValidatedTag> ValidatedVectorUnivariateFunction;
typedef VectorUnivariateFunction<EffectiveTag> EffectiveVectorUnivariateFunction;

// Function interface declarations
template<class P, class D, class C> class FunctionInterface;
template<class P, class D=BoxDomainType> using ScalarFunctionInterface = FunctionInterface<P,D,IntervalDomainType>;
template<class P, class D=BoxDomainType> using VectorFunctionInterface = FunctionInterface<P,D,BoxDomainType>;

typedef ScalarFunctionInterface<ApproximateTag> ApproximateScalarFunctionInterface;
typedef ScalarFunctionInterface<ValidatedTag> ValidatedScalarFunctionInterface;
typedef ScalarFunctionInterface<EffectiveTag> EffectiveScalarFunctionInterface;

typedef VectorFunctionInterface<ApproximateTag> ApproximateVectorFunctionInterface;
typedef VectorFunctionInterface<ValidatedTag> ValidatedVectorFunctionInterface;
typedef VectorFunctionInterface<EffectiveTag> EffectiveVectorFunctionInterface;



// Function models declarations


template<class P, class D, class C, class PR, class PRE=PR> class FunctionModelInterface;
template<class P, class D, class PR, class PRE=PR> using ScalarFunctionModelInterface = FunctionModelInterface<P,D,IntervalDomainType,PR,PRE>;
template<class P, class D, class PR, class PRE=PR> using VectorFunctionModelInterface = FunctionModelInterface<P,D,BoxDomainType,PR,PRE>;

template<class P, class D, class C, class PR, class PRE=PR> class FunctionModel;
template<class P, class D, class PR, class PRE=PR> using ScalarFunctionModel = FunctionModel<P,D,IntervalDomainType,PR,PRE>;
template<class P, class D, class PR, class PRE=PR> using VectorFunctionModel = FunctionModel<P,D,BoxDomainType,PR,PRE>;

template<class P, class PR, class PRE=PR> struct FunctionModelTraits;

template<class P, class D=BoxDomainType> using ScalarFunctionModelDPInterface = ScalarFunctionModelInterface<P,D,DoublePrecision>;
template<class P, class D=BoxDomainType> using VectorFunctionModelDPInterface = VectorFunctionModelInterface<P,D,DoublePrecision>;
template<class P, class D=BoxDomainType> using ScalarFunctionModelDP = ScalarFunctionModel<P,D,DoublePrecision>;
template<class P, class D=BoxDomainType> using VectorFunctionModelDP = VectorFunctionModel<P,D,DoublePrecision>;

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

using ValidatedScalarFunctionModelDPInterface = ScalarFunctionModelInterface<ValidatedTag,BoxDomainType,DoublePrecision>;
using ValidatedVectorFunctionModelDPInterface = VectorFunctionModelInterface<ValidatedTag,BoxDomainType,DoublePrecision>;

using ValidatedScalarFunctionModelDP = ScalarFunctionModel<ValidatedTag,BoxDomainType,DoublePrecision>;
using ValidatedVectorFunctionModelDP = VectorFunctionModel<ValidatedTag,BoxDomainType,DoublePrecision>;

using ApproximateScalarFunctionModelDPInterface = ScalarFunctionModelInterface<ApproximateTag,BoxDomainType,DoublePrecision>;
using ApproximateVectorFunctionModelDPInterface = VectorFunctionModelInterface<ApproximateTag,BoxDomainType,DoublePrecision>;

using ApproximateScalarFunctionModelDP = ScalarFunctionModel<ApproximateTag,BoxDomainType,DoublePrecision>;
using ApproximateVectorFunctionModelDP = VectorFunctionModel<ApproximateTag,BoxDomainType,DoublePrecision>;

template<class P, class PR=DoublePrecision, class PRE=PR> class FunctionModelFactoryInterface;
typedef FunctionModelFactoryInterface<ValidatedTag,DoublePrecision> ValidatedFunctionModelDPFactoryInterface;
template<class P, class PR=DoublePrecision, class PRE=PR> class FunctionModelFactory;
typedef FunctionModelFactory<ValidatedTag,DoublePrecision> ValidatedFunctionModelDPFactory;
template<class FMF, class D> class FunctionModelCreator;



} // namespace Ariadne

#endif
