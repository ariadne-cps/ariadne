/***************************************************************************
 *            function_concepts.hpp
 *
 *  Copyright  2022  Pieter Collins
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

#ifndef ARIADNE_FUNCTION_CONCEPTS_HPP
#define ARIADNE_FUNCTION_CONCEPTS_HPP


namespace Ariadne {

template<class P, class FLT> class TaylorModel;
template<class FLT> using ApproximateTaylorModel = TaylorModel<ApproximateTag,FLT>;
template<class FLT> using ValidatedTaylorModel = TaylorModel<ValidatedTag,FLT>;

template<class Y> class ElementaryAlgebra;

using ScalarUnivariate = RealScalar(RealScalar);
using VectorUnivariate = RealVector(RealScalar);
using ScalarMultivariate = RealScalar(RealVector);
using VectorMultivariate = RealVector(RealVector);


template<class FLT, template<class>class A, template<class>class R, class T> concept CallableOn = requires(FLT f, A<T> a) {
    { f(a) } -> SameAs<R<T>>;
};

template<class FLT, template<class>class Argument, template<class> class Result> concept IsApproximateCallable = requires(FLT const f) {
    { f(declval<Argument<FloatDPApproximation>>()) } -> SameAs<Result<FloatDPApproximation>>;
    { f(declval<Argument<FloatMPApproximation>>()) } -> SameAs<Result<FloatMPApproximation>>;
    { f(declval<Argument<Differential<FloatDPApproximation>>>()) } -> SameAs<Result<Differential<FloatDPApproximation>>>;
    { f(declval<Argument<Differential<FloatMPApproximation>>>()) } -> SameAs<Result<Differential<FloatMPApproximation>>>;
    { f(declval<Argument<TaylorModel<ApproximateTag,FloatDP>>>()) } -> SameAs<Result<TaylorModel<ApproximateTag,FloatDP>>>;
    { f(declval<Argument<TaylorModel<ApproximateTag,FloatMP>>>()) } -> SameAs<Result<TaylorModel<ApproximateTag,FloatMP>>>;

    { f(declval<Argument<Formula<ApproximateNumber>>>()) } -> SameAs<Result<Formula<ApproximateNumber>>>;
    { f(declval<Argument<ElementaryAlgebra<ApproximateNumber>>>()) } -> SameAs<Result<ElementaryAlgebra<ApproximateNumber>>>;
};

template<class FLT, template<class>class Argument, template<class> class Result>
concept IsValidatedCallable = IsApproximateCallable<FLT,Argument,Result> and requires(FLT const f) {
    { f(declval<Argument<FloatDPBounds>>()) } -> SameAs<Result<FloatDPBounds>>;
    { f(declval<Argument<FloatMPBounds>>()) } -> SameAs<Result<FloatMPBounds>>;
    { f(declval<Argument<Differential<FloatDPBounds>>>()) } -> SameAs<Result<Differential<FloatDPBounds>>>;
    { f(declval<Argument<Differential<FloatMPBounds>>>()) } -> SameAs<Result<Differential<FloatMPBounds>>>;
    { f(declval<Argument<ValidatedTaylorModel<FloatDP>>>()) } -> SameAs<Result<ValidatedTaylorModel<FloatDP>>>;
    { f(declval<Argument<ValidatedTaylorModel<FloatMP>>>()) } -> SameAs<Result<ValidatedTaylorModel<FloatMP>>>;
    { f(declval<Argument<ValidatedTaylorModel<FloatDPUpperInterval>>>()) } -> SameAs<Result<ValidatedTaylorModel<FloatDPUpperInterval>>>;
    { f(declval<Argument<ValidatedTaylorModel<FloatMPUpperInterval>>>()) } -> SameAs<Result<ValidatedTaylorModel<FloatMPUpperInterval>>>;

    { f(declval<Argument<Formula<ValidatedNumber>>>()) } -> SameAs<Result<Formula<ValidatedNumber>>>;
    { f(declval<Argument<ElementaryAlgebra<ValidatedNumber>>>()) } -> SameAs<Result<ElementaryAlgebra<ValidatedNumber>>>;

    { f(declval<Argument<FloatDP>>()) } -> ConvertibleTo<Result<FloatDPBounds>>;
    { f(declval<Argument<FloatMP>>()) } -> ConvertibleTo<Result<FloatMPBounds>>;
};

template<class FLT, template<class>class Argument, template<class> class Result>
concept IsEffectiveCallable = IsValidatedCallable<FLT,Argument,Result> and requires(FLT const f) {
    { f(declval<Argument<Real>>()) } -> SameAs<Result<Real>>;
    { f(declval<Argument<ElementaryAlgebra<Real>>>()) } -> SameAs<Result<ElementaryAlgebra<Real>>>;
    { f(declval<Argument<Formula<Real>>>()) } -> SameAs<Result<Formula<Real>>>;
    { f(declval<Argument<ElementaryAlgebra<EffectiveNumber>>>()) } -> SameAs<Result<ElementaryAlgebra<EffectiveNumber>>>;
    { f(declval<Argument<Formula<EffectiveNumber>>>()) } -> SameAs<Result<Formula<EffectiveNumber>>>;
};

template<class FLT, class P, template<class>class Argument, template<class>class Result>
concept IsCallable
    = ((not Same<P,ApproximateTag>) or IsApproximateCallable<FLT,Argument,Result>)
        and ((not Same<P,ValidatedTag>) or IsValidatedCallable<FLT,Argument,Result>)
        and ((not Same<P,EffectiveTag>) or IsEffectiveCallable<FLT,Argument,Result>);

template<class FLT, class P, class SIG>
concept IsCallableSignature = IsCallable<FLT,P, SignatureTraits<SIG>::template Argument, SignatureTraits<SIG>::template Result>;

template<class FLT, class SIG> concept IsBaseFunction = requires(FLT const& f)
{
    { f.argument_size() } -> SameAs<typename SignatureTraits<SIG>::ArgumentSizeType>;
    { f.result_size() } -> SameAs<typename SignatureTraits<SIG>::ResultSizeType>;
    { f.domain() } -> SameAs<typename SignatureTraits<SIG>::EntireDomainType>;
};

template<class FLT, class SIG> concept IsApproximateFunction = IsBaseFunction<FLT,SIG> and IsCallableSignature<FLT,ApproximateTag,SIG>;
template<class FLT, class SIG> concept IsValidatedFunction = IsApproximateFunction<FLT,SIG> and IsCallableSignature<FLT,ValidatedTag,SIG>;
template<class FLT, class SIG> concept IsEffectiveFunction = IsValidatedFunction<FLT,SIG> and IsCallableSignature<FLT,EffectiveTag,SIG>;

template<class FLT, class P, class SIG> concept AFunction
    = ((not Same<P,ApproximateTag>) or IsApproximateFunction<FLT,SIG>)
        and ((not Same<P,ValidatedTag>) or IsValidatedFunction<FLT,SIG>)
        and ((not Same<P,EffectiveTag>) or IsEffectiveFunction<FLT,SIG>);

// FIXME: Provided to avoid AFunction concept checking conformance before argument classes are defined.
template<class FLT, class SIG> concept IsFunctionClass = Same<FLT,Function<ApproximateTag,SIG>> or Same<FLT,Function<ValidatedTag,SIG>> or Same<FLT,Function<EffectiveTag,SIG>> or Same<FLT,FunctionPatch<ApproximateTag,SIG>> or Same<FLT,FunctionPatch<ValidatedTag,SIG>>;


template<class FLT, class P, class SIG> concept IsEntireFunction = AFunction<FLT,P,SIG> and requires(FLT const& f) {
    { f.domain() } -> ConvertibleTo<typename SignatureTraits<SIG>::EntireDomainType>;
};

template<class FLT, class P, class SIG> concept IsFunctionPatch = AFunction<FLT,P,SIG> and requires(FLT const& f) {
    { f.domain() } -> ConvertibleTo<typename SignatureTraits<SIG>::BoundedDomainType>;
};



} // namespace Ariadne

#endif // ARIADNE_FUNCTION_CONCEPTS_HPP
