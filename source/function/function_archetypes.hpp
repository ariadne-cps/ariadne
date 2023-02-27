/***************************************************************************
 *            function_archetypes.hpp
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

#ifndef ARIADNE_FUNCTION_ARCHETYPES_HPP
#define ARIADNE_FUNCTION_ARCHETYPES_HPP


namespace Ariadne {

template<class P, class SIG> class FunctionArchetype;

template<class SIG> class FunctionArchetype<ApproximateTag,SIG>
{
    using P=ApproximateTag;
  public:
    using DomainType = typename SignatureTraits<SIG>::DomainType;
    using CodomainType = typename SignatureTraits<SIG>::CodomainType;
    using ArgumentSizeType = typename SignatureTraits<SIG>::ArgumentSizeType;
    using ResultSizeType = typename SignatureTraits<SIG>::ResultSizeType;
    using ArgumentIndexType = typename SignatureTraits<SIG>::ArgumentIndexType;
    template<class X> using Argument = typename SignatureTraits<SIG>::template Argument<X>;
    template<class X> using Result = typename SignatureTraits<SIG>::template Result<X>;

    ArgumentSizeType argument_size() const { assert(false); }
    ResultSizeType result_size() const { assert(false); }
    DomainType domain() const { assert(false); }

    Result<FloatDPApproximation> operator() (const Argument<FloatDPApproximation>& x) const { assert(false); }
    Result<FloatMPApproximation> operator() (const Argument<FloatMPApproximation>& x) const { assert(false); }
    Result<Differential<FloatDPApproximation>> operator() (const Argument< Differential<FloatDPApproximation> >& x) const { assert(false); }
    Result<Differential<FloatMPApproximation>> operator() (const Argument< Differential<FloatMPApproximation> >& x) const { assert(false); }
    Result<TaylorModel<ApproximateTag,FloatDP>> operator() (const Argument< TaylorModel<ApproximateTag,FloatDP> >& x) const { assert(false); }
    Result<TaylorModel<ApproximateTag,FloatMP>> operator() (const Argument< TaylorModel<ApproximateTag,FloatMP> >& x) const { assert(false); }
    Result<Formula<ApproximateNumber>> operator() (const Argument< Formula<ApproximateNumber> >& x) const { assert(false); }
    Result<ElementaryAlgebra<ApproximateNumber>> operator() (const Argument< ElementaryAlgebra<ApproximateNumber> >& x) const { assert(false); }

    friend FunctionArchetype<P,SIG> derivative(FunctionArchetype<P,SIG>, SizeType) { assert(false); }

    friend OutputStream& operator<<(OutputStream& os, FunctionArchetype<P,SIG> const& f) { assert(false); }

    template<class X> Result<X> _call(Argument<X> const& x) const { return (*this)(x); }
};



template<class SIG> class FunctionArchetype<ValidatedTag,SIG>
    : public FunctionArchetype<ApproximateTag,SIG>
{
    using P=ValidatedTag;
  public:

    using FunctionArchetype<ApproximateTag,SIG>::operator();
    using FunctionArchetype<ApproximateTag,SIG>::_call;

    Scalar<FloatDPBounds> operator() (const Vector<FloatDPBounds>& x) const { assert(false); }
    Scalar<FloatMPBounds> operator() (const Vector<FloatMPBounds>& x) const { assert(false); }
    Scalar<Differential<FloatDPBounds>> operator() (const Vector< Differential<FloatDPBounds> >& x) const { assert(false); }
    Scalar<Differential<FloatMPBounds>> operator() (const Vector< Differential<FloatMPBounds> >& x) const { assert(false); }
    Scalar<TaylorModel<ValidatedTag,FloatDP>> operator() (const Vector< TaylorModel<ValidatedTag,FloatDP> >& x) const { assert(false); }
    Scalar<TaylorModel<ValidatedTag,FloatMP>> operator() (const Vector< TaylorModel<ValidatedTag,FloatMP> >& x) const { assert(false); }
    Scalar<TaylorModel<ValidatedTag,FloatDPUpperInterval>> operator() (const Vector<TaylorModel<ValidatedTag,FloatDPUpperInterval>>& x) const { assert(false); }
    Scalar<TaylorModel<ValidatedTag,FloatMPUpperInterval>> operator() (const Vector<TaylorModel<ValidatedTag,FloatMPUpperInterval>>& x) const { assert(false); }

    Scalar<Formula<ValidatedNumber>> operator() (const Vector< Formula<ValidatedNumber> >& x) const { assert(false); }
    Scalar<ElementaryAlgebra<ValidatedNumber>> operator() (const Vector< ElementaryAlgebra<ValidatedNumber> >& x) const { assert(false); }

    Scalar<Function<ValidatedTag,SIG>> operator() (const Vector< Function<ValidatedTag,SIG> >& x) const { assert(false); }

    inline Scalar<FloatDPBounds> operator() (const Vector<FloatDP>& x) const { assert(false); }
    inline Scalar<FloatMPBounds> operator() (const Vector<FloatMP>& x) const { assert(false); }

    friend FunctionArchetype<P,SIG> derivative(FunctionArchetype<P,SIG>, SizeType) { assert(false); }
};

template<class SIG> class FunctionArchetype<EffectiveTag,SIG>
    : public FunctionArchetype<ValidatedTag,SIG>
{
    using P=EffectiveTag;
  public:
    using FunctionArchetype<ValidatedTag,SIG>::operator();
    using FunctionArchetype<ValidatedTag,SIG>::_call;

    Scalar<Real> operator() (const Vector<Real>& x) const { assert(false); }
    Scalar<ElementaryAlgebra<Real>> operator() (const Vector<ElementaryAlgebra<Real>>& x) const { assert(false); }
    Scalar<Formula<Real>> operator() (const Vector<Formula<Real>>& x) const { assert(false); }
    Scalar<ElementaryAlgebra<EffectiveNumber>> operator() (const Vector<ElementaryAlgebra<EffectiveNumber>>& x) const { assert(false); }
    Scalar<Formula<EffectiveNumber>> operator() (const Vector<Formula<EffectiveNumber>>& x) const { assert(false); }

    friend FunctionArchetype<P,SIG> derivative(FunctionArchetype<P,SIG>, SizeType) { assert(false); }
};




template<class S> using DimensionOf = std::remove_cv_t<decltype(declval<S>().dimension())>;

template<class P, class SIG, class FLT, class FLTE> class UnitFunctionModelArchetype;

template<class FLT, class FLTE> class UnitFunctionModelArchetype<ValidatedTag, RealScalar(RealVector), FLT, FLTE>
{
    using P=ValidatedTag; using RES=RealScalar; using ARG=RealVector; using SIG=RES(ARG);
    using IND=MultiIndex;
    typedef UnitFunctionModelArchetype<P,SIG,FLT,FLTE> SelfType;
  public:
    typedef P Paradigm;
    typedef Bounds<FLT> NumericType;
    typedef ValidatedNumber GenericNumericType;
    typedef FLT CoefficientType;
    typedef Error<FLTE> ErrorType;
        typedef FLTE ErrorValueType;
    typedef typename FLT::PrecisionType PrecisionType;
    typedef typename FLTE::PrecisionType ErrorPrecisionType;
        typedef FLT RawFloatType;

    template<class X> using Argument = typename SignatureTraits<SIG>::template Argument<X>;
    template<class X> using Result = typename SignatureTraits<SIG>::template Result<X>;

    typedef SizeType ArgumentSizeType;
    typedef SizeType ArgumentIndexType;
    typedef BoxDomainType DomainType;
    typedef IntervalDomainType CodomainType;
    typedef Interval<UpperBound<FLT>> RangeType;
        typedef MultiIndex IndexType;
        typedef Expansion<MultiIndex,FLT> ExpansionType;
        typedef FLT ValueType;
    typedef Positive<UpperBound<FLT>> NormType;

    typedef Sweeper<FLT> PropertiesType;

    typedef Function<P,RES(ARG)> FunctionType;
    typedef Function<P,RealScalar(ARG)> ScalarFunctionType;
    typedef Function<P,RealVector(ARG)> VectorFunctionType;
    static SelfType zero(ArgumentSizeType, PropertiesType);
        //static SelfType constant(ArgumentSizeType, NumericType, PropertiesType);
    static SelfType constant(ArgumentSizeType, GenericNumericType, PropertiesType);
    static SelfType scaling(ArgumentSizeType, ArgumentIndexType, IntervalDomainType, PropertiesType);
    static SelfType unit_ball(ArgumentSizeType, PropertiesType);
    static Vector<SelfType> constants(ArgumentSizeType, Vector<NumericType>, PropertiesType);
    static Vector<SelfType> scalings(BoxDomainType,PropertiesType);
    static Vector<SelfType> unscalings(BoxDomainType,PropertiesType);

        explicit UnitFunctionModelArchetype();
        explicit UnitFunctionModelArchetype(SizeType);
    UnitFunctionModelArchetype(SelfType const& fm);
//    UnitFunctionModelArchetype(UnitFunctionModelArchetype&& fm);
    UnitFunctionModelArchetype& operator=(SelfType const& fm);
//    UnitFunctionModelArchetype& operator=(UnitFunctionModelArchetype&& fm);

        UnitFunctionModelArchetype(Expansion<MultiIndex,ExactDouble>, ExactDouble, PropertiesType);
        UnitFunctionModelArchetype(Expansion<MultiIndex,FLT>, FLTE, PropertiesType);
    UnitFunctionModelArchetype(ExpansionType, ErrorType, PropertiesType);
    UnitFunctionModelArchetype(ArgumentSizeType, PropertiesType);
    UnitFunctionModelArchetype(DomainType, FunctionType, PropertiesType);
//    UnitFunctionModelArchetype(, PropertiesType);
    UnitFunctionModelArchetype& operator=(GenericNumericType);

    SelfType zero_element() const;

    ArgumentSizeType argument_size() const;
    DomainType domain() const;
    CodomainType codomain() const;
    RangeType range() const;
        ExpansionType const& expansion() const;
        ExpansionType& expansion();
        DegreeType degree() const;
        SizeType number_of_nonzeros() const;
    CoefficientType const& operator[] (IndexType) const;
    CoefficientType& operator[] (IndexType);
        CoefficientType const& value() const;
        Void set_value(CoefficientType);
    ErrorType const& error() const;
    ErrorType& error();
    Void set_error(ErrorType);
    PropertiesType properties() const;
    Void set_properties(PropertiesType) const;
    PrecisionType precision() const;
    friend NormType norm(SelfType const&);

    Void simplify();
    Void simplify(PropertiesType);
    Void clobber();

    friend SelfType operator+(SelfType);
    friend SelfType operator-(SelfType);
    friend SelfType operator+(SelfType, SelfType);
    friend SelfType operator-(SelfType, SelfType);
    friend SelfType operator*(SelfType, SelfType);
    friend SelfType operator/(SelfType, SelfType);
    friend SelfType& operator+=(SelfType&, SelfType);
    friend SelfType& operator-=(SelfType&, SelfType);
    friend SelfType& operator*=(SelfType&, SelfType);
    friend SelfType& operator/=(SelfType&, SelfType);
    friend SelfType operator+(SelfType, NumericType);
    friend SelfType operator-(SelfType, NumericType);
    friend SelfType operator*(SelfType, NumericType);
    friend SelfType operator/(SelfType, NumericType);
    friend SelfType operator+(NumericType, SelfType);
    friend SelfType operator-(NumericType, SelfType);
    friend SelfType operator*(NumericType, SelfType);
    friend SelfType operator/(NumericType, SelfType);
    friend SelfType& operator+=(SelfType&, NumericType);
    friend SelfType& operator-=(SelfType&, NumericType);
    friend SelfType& operator*=(SelfType&, NumericType);
    friend SelfType& operator/=(SelfType&, NumericType);
    friend SelfType operator+(SelfType, GenericNumericType);
    friend SelfType operator-(SelfType, GenericNumericType);
    friend SelfType operator*(SelfType, GenericNumericType);
    friend SelfType operator/(SelfType, GenericNumericType);
    friend SelfType operator+(GenericNumericType, SelfType);
    friend SelfType operator-(GenericNumericType, SelfType);
    friend SelfType operator*(GenericNumericType, SelfType);
    friend SelfType operator/(GenericNumericType, SelfType);

    friend SelfType add(SelfType, SelfType);
    friend SelfType sub(SelfType, SelfType);
    friend SelfType mul(SelfType, SelfType);
    friend SelfType div(SelfType, SelfType);
    friend SelfType add(SelfType, NumericType);
    friend SelfType sub(SelfType, NumericType);
    friend SelfType mul(SelfType, NumericType);
    friend SelfType div(SelfType, NumericType);
    friend SelfType max(SelfType, NumericType);
    friend SelfType min(SelfType, NumericType);
    friend SelfType add(NumericType, SelfType);
    friend SelfType sub(NumericType, SelfType);
    friend SelfType mul(NumericType, SelfType);
    friend SelfType div(NumericType, SelfType);
    friend SelfType max(NumericType, SelfType);
    friend SelfType min(NumericType, SelfType);
    friend SelfType add(SelfType, GenericNumericType);
    friend SelfType sub(SelfType, GenericNumericType);
    friend SelfType mul(SelfType, GenericNumericType);
    friend SelfType div(SelfType, GenericNumericType);
    friend SelfType max(SelfType, GenericNumericType);
    friend SelfType min(SelfType, GenericNumericType);
    friend SelfType add(GenericNumericType, SelfType);
    friend SelfType sub(GenericNumericType, SelfType);
    friend SelfType mul(GenericNumericType, SelfType);
    friend SelfType div(GenericNumericType, SelfType);
    friend SelfType max(GenericNumericType, SelfType);
    friend SelfType min(GenericNumericType, SelfType);

    friend SelfType nul(SelfType);
    friend SelfType pos(SelfType);
    friend SelfType neg(SelfType);
    friend SelfType sqr(SelfType);
    friend SelfType hlf(SelfType);
    friend SelfType rec(SelfType);
    friend SelfType pow(SelfType, Int);
    friend SelfType sqrt(SelfType);
    friend SelfType exp(SelfType);
    friend SelfType log(SelfType);
    friend SelfType sin(SelfType);
    friend SelfType cos(SelfType);
    friend SelfType tan(SelfType);
    friend SelfType asin(SelfType);
    friend SelfType acos(SelfType);
    friend SelfType atan(SelfType);

    friend SelfType abs(SelfType);
    friend SelfType max(SelfType, SelfType);
    friend SelfType min(SelfType, SelfType);

#warning Maybe change to explicit types
    template<class X> Result<X> operator() (Argument<X>) const;

    template<class X> Result<X> evaluate(Argument<X>) const;
//    template<class X> friend Result<X> evaluate(SelfType const&, Argument<X>);
    template<class X> friend Result<X> evaluate(SelfType const&, Vector<X> const&);
    template<class X> friend Vector<Result<X>> evaluate(Vector<SelfType> const&, Vector<X> const&);
            template<class X> Result<X> _call(Argument<X>) const;
    friend SelfType evaluate(ScalarUnivariateFunction<P>, SelfType const&);

    friend SelfType partial_evaluate(SelfType, ArgumentIndexType, NumericType);

    friend SelfType derivative(SelfType, ArgumentIndexType);
    friend SelfType antiderivative(SelfType, ArgumentIndexType);
//    friend SelfType antiderivative(SelfType, ArgumentIndexType, NumericType);

    friend Scalar<SelfType> compose(Scalar<SelfType>, Vector<SelfType>);
    friend Vector<SelfType> compose(Vector<SelfType>, Vector<SelfType>);

    friend Scalar<SelfType> embed(SizeType, Scalar<SelfType>, SizeType);
    friend Vector<SelfType> embed(SizeType, Vector<SelfType>, SizeType);

    friend Vector<SelfType> combine(Scalar<SelfType>, Scalar<SelfType>);
    friend Vector<SelfType> combine(Scalar<SelfType>, Vector<SelfType>);
    friend Vector<SelfType> combine(Vector<SelfType>, Scalar<SelfType>);
    friend Vector<SelfType> combine(Vector<SelfType>, Vector<SelfType>);

    friend SelfType refinement(SelfType, SelfType);
    friend Bool refines(SelfType, SelfType);
    friend Bool inconsistent(SelfType, SelfType);

    friend OutputStream& operator<<(OutputStream&, SelfType const&);
};

template<class S> using DimensionOf = std::remove_cv_t<decltype(declval<S>().dimension())>;

auto evaluate(Function<ValidatedTag,RealScalar(RealVector)> const&,
              Vector<UnitFunctionModelArchetype<ValidatedTag,RealScalar(RealVector),FloatDP,FloatDP>>const&)
                  -> Scalar<UnitFunctionModelArchetype<ValidatedTag,RealScalar(RealVector),FloatDP,FloatDP>>;
auto evaluate(Function<ValidatedTag,RealVector(RealVector)> const&,
              Vector<UnitFunctionModelArchetype<ValidatedTag,RealScalar(RealVector),FloatDP,FloatDP>>const&)
                  -> Vector<UnitFunctionModelArchetype<ValidatedTag,RealScalar(RealVector),FloatDP,FloatDP>>;

template<class X> auto
evaluate(Scalar<UnitFunctionModelArchetype<ValidatedTag,RealScalar(RealVector),FloatDP,FloatDP>>const& f, Vector<X> const& x) -> Scalar<X>;
template<class X> auto
evaluate(Vector<UnitFunctionModelArchetype<ValidatedTag,RealScalar(RealVector),FloatDP,FloatDP>>const& f, Vector<X> const& x) -> Vector<X>;

template<class P, class FLT, class FLTE> auto
partial_evaluate(Vector<UnitFunctionModelArchetype<P,RealScalar(RealVector),FLT,FLTE>>, SizeType, Bounds<FLT>)
    -> Vector<UnitFunctionModelArchetype<P,RealScalar(RealVector),FLT,FLTE>>;
template<class P, class FLT, class FLTE> auto
embed(SizeType, Vector<UnitFunctionModelArchetype<P,RealScalar(RealVector),FLT,FLTE>>, SizeType)
    -> Vector<UnitFunctionModelArchetype<P,RealScalar(RealVector),FLT,FLTE>>;

template<class P, class FLT, class FLTE> decltype(auto)
evaluate(Function<P,RealScalar(RealVector)> const& f,
         Vector<UnitFunctionModelArchetype<P,RealScalar(RealVector),FLT,FLTE>>const& x) {
    return x[0];
}
template<class P, class FLT, class FLTE> decltype(auto)
evaluate(Function<P,RealVector(RealVector)> const& f,
         Vector<UnitFunctionModelArchetype<P,RealScalar(RealVector),FLT,FLTE>>const& x) {
    return x;
}

template<class P, class FLT, class FLTE, class X> Vector<X>
evaluate(Vector<UnitFunctionModelArchetype<P,RealScalar(RealVector),FLT,FLTE>>const& f, Vector<X> const& x);
template<class P, class FLT, class FLTE> Vector<UnitFunctionModelArchetype<P,RealScalar(RealVector),FLT,FLTE>>
partial_evaluate(Vector<UnitFunctionModelArchetype<P,RealScalar(RealVector),FLT,FLTE>>, SizeType, Bounds<FLT>);
template<class P, class FLT, class FLTE> Vector<UnitFunctionModelArchetype<P,RealScalar(RealVector),FLT,FLTE>>
embed(UnitBox, Vector<UnitFunctionModelArchetype<P,RealScalar(RealVector),FLT,FLTE>>, UnitBox);

} // namespace Ariadne




#include "function_concepts.hpp"

namespace Ariadne {

using ScalarMultivariate=RealScalar(RealVector);

static_assert(IsApproximateFunction<FunctionArchetype<ApproximateTag,ScalarMultivariate>,ScalarMultivariate>);
static_assert(IsValidatedFunction<FunctionArchetype<ValidatedTag,ScalarMultivariate>,ScalarMultivariate>);
static_assert(IsEffectiveFunction<FunctionArchetype<EffectiveTag,ScalarMultivariate>,ScalarMultivariate>);

static_assert(AFunction<FunctionArchetype<ApproximateTag,ScalarMultivariate>,ApproximateTag,ScalarMultivariate>);
static_assert(AFunction<FunctionArchetype<ValidatedTag,ScalarMultivariate>,ValidatedTag,ScalarMultivariate>);
static_assert(AFunction<FunctionArchetype<EffectiveTag,ScalarMultivariate>,EffectiveTag,ScalarMultivariate>);

} // namespace Ariadne

#endif // ARIADNE_FUNCTION_ARCHETYPES_HPP
