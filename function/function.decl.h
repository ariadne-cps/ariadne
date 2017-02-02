/***************************************************************************
 *            function.decl.h
 *
 *  Copyright 2016  Pieter Collins
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

/*! \file function.decl.h
 *  \brief Declarations of function types.
 */

#ifndef ARIADNE_FUNCTION_DECL_H
#define ARIADNE_FUNCTION_DECL_H

#include "geometry/interval.decl.h"
#include "geometry/box.decl.h"

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


template<class P, class PR, class PRE=PR> class ScalarFunctionModelInterface;
template<class P, class PR, class PRE=PR> class VectorFunctionModelInterface;

template<class P, class PR, class PRE=PR> class ScalarFunctionModel;
template<class P, class PR, class PRE=PR> class VectorFunctionModel;

template<class P, class PR, class PRE=PR> struct FunctionModelTraits;

template<class P> using ScalarFunctionModel64Interface = ScalarFunctionModelInterface<P,Precision64>;
template<class P> using VectorFunctionModel64Interface = VectorFunctionModelInterface<P,Precision64>;
template<class P> using ScalarFunctionModel64 = ScalarFunctionModel<P,Precision64>;
template<class P> using VectorFunctionModel64 = VectorFunctionModel<P,Precision64>;

template<class PR> struct FunctionModelTraits<ApproximateTag,PR> {
    typedef FloatApproximation<PR> CoefficientType;
    typedef Void ErrorType;
    typedef FloatApproximation<PR> NumericType;
};

template<class PR, class PRE> struct FunctionModelTraits<ValidatedTag,PR,PRE> {
    typedef FloatValue<PR> CoefficientType;
    typedef FloatError<PRE> ErrorType;
    typedef FloatBounds<PR> NumericType;
};

template<class P, class PR, class PRE=PR> using CanonicalNumericType = typename FunctionModelTraits<P,PR,PRE>::NumericType;
template<class P, class PR> using CanonicalCoefficientType = typename FunctionModelTraits<P,PR>::CoefficientType;
template<class P, class PRE> using CanonicalErrorType = typename FunctionModelTraits<P,PRE,PRE>::ErrorType;

template<class P> using CanonicalNumeric64Type = typename FunctionModelTraits<P,Precision64>::NumericType;
template<class P> using CanonicalCoefficient64Type = typename FunctionModelTraits<P,Precision64>::CoefficientType;
template<class P> using CanonicalError64Type = typename FunctionModelTraits<P,Precision64,Precision64>::ErrorType;

template<class X> using PrecisionType = typename X::PrecisionType;
template<class X> using ErrorPrecisionType = typename X::ErrorPrecisionType;

using ValidatedScalarFunctionModel64Interface = ScalarFunctionModelInterface<ValidatedTag,Precision64>;
using ValidatedVectorFunctionModel64Interface = VectorFunctionModelInterface<ValidatedTag,Precision64>;

using ValidatedScalarFunctionModel64 = ScalarFunctionModel<ValidatedTag,Precision64>;
using ValidatedVectorFunctionModel64 = VectorFunctionModel<ValidatedTag,Precision64>;

template<class P, class PR=Precision64, class PRE=PR> class FunctionModelFactoryInterface;
typedef FunctionModelFactoryInterface<ValidatedTag,Precision64> ValidatedFunctionModel64FactoryInterface;
template<class P, class PR=Precision64, class PRE=PR> class FunctionModelFactory;
typedef FunctionModelFactory<ValidatedTag,Precision64> ValidatedFunctionModel64Factory;
template<class FMF> class FunctionModelCreator;



} // namespace Ariadne

#endif
