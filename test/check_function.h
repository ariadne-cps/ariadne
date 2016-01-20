/***************************************************************************
 *            check_function.h
 *
 *  Copyright 2009-13  Pieter Collins
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

#include "test.h"
#include "test/utility.h"
#include "test/check_algebra.h"

#include "numeric/operators.h"
#include "function/function.h"

namespace Ariadne {
template<class X> String class_name();
template<> String class_name<ApproximateTag>() { return "ApproximateTag"; }
template<> String class_name<ValidatedTag>() { return "ValidatedTag"; }
template<> String class_name<EffectiveTag>() { return "EffectiveTag"; }
} // namespace Ariadne

using namespace Ariadne;

typedef ExactIntervalType IntervalDomain;
typedef ExactBoxType BoxDomain;
template<class P> using Function = ScalarFunction<P>;
typedef ScalarFunction<EffectiveTag> EffectiveFunction;

template<class F, class R=Return<DontCare>, class = Fallback> struct HasCodomainMethod : False { };
template<class F, class R> struct HasCodomainMethod<F, R, EnableIf<IsDefined<decltype(declval<F>().codomain())>,Fallback>> : True { };

template<class F, class R=Return<DontCare>, class = Fallback> struct HasDomainMethod : False { };
template<class F, class R> struct HasDomainMethod<F, R, EnableIf<IsDefined<decltype(declval<F>().domain())>,Fallback>> : True { };

template<class F, class A, class R=Return<DontCare>, class = Fallback> struct HasCallMethod : False { };
template<class F, class A, class R>
    struct HasCallMethod<F,A,Return<R>,EnableIf<IsConvertible<decltype(declval<F>()(declval<A>())),R>,Fallback>> : True { };

template<class F, class A, class R=Return<DontCare>, class = Fallback> struct HasGradientMethod : False { };
template<class F, class A, class R>
    struct HasGradientMethod<F,A,Return<R>,EnableIf<IsConvertible<decltype(declval<F>().gradient(declval<A>())),R>,Fallback>> : True { };

template<class F, class A, class R=Return<DontCare>, class = Fallback> struct HasDifferentialMethod : False { };
template<class F, class A, class R>
    struct HasDifferentialMethod<F,A,Return<R>,EnableIf<IsConvertible<decltype(declval<F>().differential(declval<A>(),declval<DegreeType>())),R>,Fallback>> : True { };

template<class F, class A, class R=Return<DontCare>, class = Fallback> struct HasEvaluate : False { };
template<class F, class A, class R>
    struct HasEvaluate<F,A,Return<R>,EnableIf<IsConvertible<decltype(evaluate(declval<F>(),declval<A>())),R>,Fallback>> : True { };

template<class F, class G, class R=Return<DontCare>, class = Fallback> struct HasCompose : False { };
template<class F, class G, class R>
    struct HasCompose<F,G,Return<R>,EnableIf<IsConvertible<decltype(compose(declval<F>(),declval<G>())),R>,Fallback>> : True { };

#define ARIADNE_HAS_TYPEDEF(typename_to_check) \
    template<class A, class = Fallback> struct Has##typename_to_check : False { }; \
    template<class A> struct Has##typename_to_check<A, EnableIf<IsDefined<typename A::typename_to_check>,Fallback>> : True { }; \

#define ARIADNE_HAS_METHOD(method_to_check) \
    template<class A, class = Fallback> struct Has_##method_to_check : False { }; \
    template<class A> struct Has_##method_to_check<A, EnableIf<IsDefined<decltype(declval<A>().method_to_check())>,Fallback>> : True { }; \

ARIADNE_HAS_TYPEDEF(Paradigm);
ARIADNE_HAS_TYPEDEF(DomainType);
ARIADNE_HAS_TYPEDEF(CodomainType);

template<class F> class CheckFunctionConcept : public CheckAlgebraConcept<F>
{
    typedef F FunctionType;
    typedef typename F::NumericType NumberType;
  public:
    CheckFunctionConcept() : CheckAlgebraConcept<F>() { }
    void check();
  private:
    void check_create_concept();
    void check_evaluable_concept();
    void check_differentiable_concept();
    void check_integrable_concept();
    void check_composable_concept();
};

template<class F> void CheckFunctionConcept<F>::check()
{
    this->check_division_algebra_concept();
    this->check_transcendental_operators_concept();

    ARIADNE_TEST_CALL(check_create_concept());
    ARIADNE_TEST_CALL(check_evaluable_concept());
    ARIADNE_TEST_CALL(check_differentiable_concept());
    ARIADNE_TEST_CALL(check_integrable_concept());
    ARIADNE_TEST_CALL(check_composable_concept());
}


template<class F> void CheckFunctionConcept<F>::check_create_concept()
{
    typedef typename F::DomainType D;
    typedef typename F::NumericType X;

    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(F::zero(declval<SizeType>())),F>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(F::constant(declval<SizeType>(),declval<X>())),F>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(F::coordinate(declval<SizeType>(),declval<SizeType>())),F>);

    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(declval<F>().create_zero()),F>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(declval<F>().create_constant(-3)),F>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(declval<F>().create_coordinate(2u)),F>);

   // ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(declval<F>().create(declval<Function<P>())),F>);
}

template<class F> void CheckFunctionConcept<F>::check_evaluable_concept()
{
    ARIADNE_TEST_STATIC_ASSERT(HasParadigm<F>);
    ARIADNE_TEST_STATIC_ASSERT(HasDomainType<F>);

    typedef typename F::Paradigm P;
    typedef typename F::DomainType D;
    typedef typename F::CodomainType C;

    ARIADNE_TEST_STATIC_ASSERT(IsConvertible<D,BoxDomain>);
    ARIADNE_TEST_STATIC_ASSERT(IsConvertible<C,IntervalDomain>);

    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(declval<F>().domain()),D>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(declval<F>().codomain()),C>);
    ARIADNE_TEST_STATIC_ASSERT(HasDomainMethod<F>);
    ARIADNE_TEST_STATIC_ASSERT(HasCodomainMethod<F>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(declval<F>().argument_size()),SizeType>);

    if(IsWeaker<ApproximateTag,P>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<ApproximateNumericType>, Return<ApproximateNumericType>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<ApproximateFloat64>, Return<ApproximateFloat64>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<ApproximateNumericType>, Return<ApproximateNumericType>>);
    }

    if(IsWeaker<ValidatedTag,P>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<ValidatedNumericType>, Return<ValidatedNumericType>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<BoundedFloat64>, Return<BoundedFloat64>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<BoundedFloat64>, Return<BoundedFloat64>>);
    }

    if(IsWeaker<EffectiveTag,P>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<EffectiveNumericType>, Return<EffectiveNumericType>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<Real>, Return<Real>>);
    }

}

template<class F> void CheckFunctionConcept<F>::check_differentiable_concept()
{
    typedef typename F::Paradigm P;

    this->check_differential_algebra_concept();

    ARIADNE_TEST_STATIC_ASSERT(HasDerivative<F,SizeType,Return<F>>);

    if(IsWeaker<ApproximateTag,P>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<Differential<ApproximateFloat64>>, Return<Differential<ApproximateFloat64>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasGradientMethod<F,Vector<ApproximateFloat64>, Return<Covector<ApproximateFloat64>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasDifferentialMethod<F,Vector<ApproximateFloat64>, Return<Differential<ApproximateFloat64>>>);
    }

    if(IsWeaker<ValidatedTag,P>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<Differential<BoundedFloat64>>, Return<Differential<BoundedFloat64>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasGradientMethod<F,Vector<BoundedFloat64>, Return<Covector<BoundedFloat64>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasDifferentialMethod<F,Vector<BoundedFloat64>, Return<Differential<BoundedFloat64>>>);
    }

    if(IsWeaker<EffectiveTag,P>::value) {
    }

}

template<class F> void CheckFunctionConcept<F>::check_integrable_concept()
{
    ARIADNE_TEST_STATIC_ASSERT(HasAntiderivative<F,SizeType, Return<F>>);
}

template<class F> void CheckFunctionConcept<F>::check_composable_concept()
{
    typedef typename F::Paradigm P;
    ARIADNE_TEST_STATIC_ASSERT(HasCompose<F,Vector<F>, Return<F>>);

    if(IsWeaker<P,ApproximateTag>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCompose<ScalarFunction<ApproximateTag>,Vector<F>, Return<F>>);
    }

    if(IsWeaker<P,ValidatedTag>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCompose<ScalarFunction<ValidatedTag>,Vector<F>, Return<F>>);
    }

    if(IsWeaker<P,EffectiveTag>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCompose<ScalarFunction<EffectiveTag>,Vector<F>, Return<F>>);
    }
}





template<class VF> class CheckVectorFunctionConcept : CheckVectorConcept<VF>
{
    typedef VF FunctionType;
    void check();
    void test();
    void check_vector_concept();
    void check_evaluable_concept();
};


template<class F> void CheckVectorFunctionConcept<F>::test()
{
    ARIADNE_TEST_CALL(check_vector_concept());
    ARIADNE_TEST_CALL(check_evaluable_concept());
}

template<class F> void CheckVectorFunctionConcept<F>::check_evaluable_concept()
{
    ARIADNE_TEST_STATIC_ASSERT(HasParadigm<F>);
    ARIADNE_TEST_STATIC_ASSERT(HasDomainType<F>);

    typedef typename F::Paradigm P;
    typedef typename F::DomainType D;
    typedef typename F::CodomainType C;

    ARIADNE_TEST_STATIC_ASSERT(IsConvertible<D,BoxDomain>);
    ARIADNE_TEST_STATIC_ASSERT(IsConvertible<C,BoxDomain>);

    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(declval<F>().domain()),D>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(declval<F>().codomain()),C>);
    ARIADNE_TEST_STATIC_ASSERT(HasDomainMethod<F>);
    ARIADNE_TEST_STATIC_ASSERT(HasCodomainMethod<F>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(declval<F>().argument_size()),SizeType>);

    if(IsWeaker<ApproximateTag,P>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<ApproximateNumericType>, Return<Vector<ApproximateNumericType>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<ApproximateFloat64>, Return<Vector<ApproximateFloat64>>>);

        ARIADNE_TEST_STATIC_ASSERT(HasEvaluate<F,Vector<ApproximateNumericType>, Return<Vector<ApproximateNumericType>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasEvaluate<F,Vector<ApproximateFloat64>, Return<Vector<ApproximateFloat64>>>);
    }

    if(IsWeaker<ValidatedTag,P>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<ValidatedNumericType>, Return<Vector<ValidatedNumericType>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<BoundedFloat64>, Return<Vector<BoundedFloat64>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<BoundedFloat64>, Return<Vector<BoundedFloat64>>>);

        ARIADNE_TEST_STATIC_ASSERT(HasEvaluate<F,Vector<ValidatedNumericType>, Return<Vector<ValidatedNumericType>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasEvaluate<F,Vector<BoundedFloat64>, Return<Vector<BoundedFloat64>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasEvaluate<F,Vector<BoundedFloat64>, Return<Vector<BoundedFloat64>>>);

    }

    if(IsWeaker<EffectiveTag,P>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<EffectiveNumericType>, Return<Vector<EffectiveNumericType>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<Real>, Return<Vector<Real>>>);

        ARIADNE_TEST_STATIC_ASSERT(HasEvaluate<F,Vector<EffectiveNumericType>, Return<Vector<EffectiveNumericType>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasEvaluate<F,Vector<Real>, Return<Vector<Real>>>);

    }

}


/*
template<class T1, class T2> inline const bool same_type(const T1& t1, const T2& t2) { return IsSame<T1,T2>::value; }

void CheckFunctionConcept::check_operators_concept()
{
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( SymbolicFunction<ValidatedTag>, SymbolicFunction<ValidatedTag>, SymbolicFunction<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( SymbolicFunction<ValidatedTag>, VectorFunctionElement<ValidatedTag>, Function<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( SymbolicFunction<ValidatedTag>, Function<ValidatedTag>, Function<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( SymbolicFunction<ValidatedTag>, FunctionModel<ValidatedTag>, FunctionModel<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( SymbolicFunction<ValidatedTag>, VectorFunctionModelElement<ValidatedTag>, FunctionModel<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( SymbolicFunction<ValidatedTag>, PolynomialFunctionModel<ValidatedTag>, PolynomialFunctionModel<ValidatedTag> );

    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionElement<ValidatedTag>, SymbolicFunction<ValidatedTag>, Function<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionElement<ValidatedTag>, VectorFunctionElement<ValidatedTag>, Function<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionElement<ValidatedTag>, Function<ValidatedTag>, Function<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionElement<ValidatedTag>, FunctionModel<ValidatedTag>, FunctionModel<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionElement<ValidatedTag>, VectorFunctionModelElement<ValidatedTag>, FunctionModel<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionElement<ValidatedTag>, PolynomialFunctionModel<ValidatedTag>, PolynomialFunctionModel<ValidatedTag> );

    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( Function<ValidatedTag>, SymbolicFunction<ValidatedTag>, Function<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( Function<ValidatedTag>, VectorFunctionElement<ValidatedTag>, Function<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( Function<ValidatedTag>, Function<ValidatedTag>, Function<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( Function<ValidatedTag>, FunctionModel<ValidatedTag>, FunctionModel<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( Function<ValidatedTag>, VectorFunctionModelElement<ValidatedTag>, FunctionModel<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( Function<ValidatedTag>, PolynomialFunctionModel<ValidatedTag>, PolynomialFunctionModel<ValidatedTag> );

    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( FunctionModel<ValidatedTag>, SymbolicFunction<ValidatedTag>, FunctionModel<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( FunctionModel<ValidatedTag>, VectorFunctionElement<ValidatedTag>, FunctionModel<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( FunctionModel<ValidatedTag>, Function<ValidatedTag>, FunctionModel<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( FunctionModel<ValidatedTag>, FunctionModel<ValidatedTag>, FunctionModel<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( FunctionModel<ValidatedTag>, VectorFunctionModelElement<ValidatedTag>, FunctionModel<ValidatedTag> );
    //ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( FunctionModel<ValidatedTag>, PolynomialFunctionModel<ValidatedTag>, PolynomialFunctionModel<ValidatedTag> );

    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionModelElement<ValidatedTag>, SymbolicFunction<ValidatedTag>, FunctionModel<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionModelElement<ValidatedTag>, VectorFunctionElement<ValidatedTag>, FunctionModel<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionModelElement<ValidatedTag>, Function<ValidatedTag>, FunctionModel<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionModelElement<ValidatedTag>, FunctionModel<ValidatedTag>, FunctionModel<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionModelElement<ValidatedTag>, VectorFunctionModelElement<ValidatedTag>, FunctionModel<ValidatedTag> );
    //ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionModelElement<ValidatedTag>, PolynomialFunctionModel<ValidatedTag>, PolynomialFunctionModel<ValidatedTag> );

    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( PolynomialFunctionModel<ValidatedTag>, SymbolicFunction<ValidatedTag>, PolynomialFunctionModel<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( PolynomialFunctionModel<ValidatedTag>, VectorFunctionElement<ValidatedTag>, PolynomialFunctionModel<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( PolynomialFunctionModel<ValidatedTag>, Function<ValidatedTag>, PolynomialFunctionModel<ValidatedTag> );
    //ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( PolynomialFunctionModel<ValidatedTag>, FunctionModel<ValidatedTag>, PolynomialFunctionModel<ValidatedTag> );
    //ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( PolynomialFunctionModel<ValidatedTag>, VectorFunctionModelElement<ValidatedTag>, PolynomialFunctionModel<ValidatedTag> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( PolynomialFunctionModel<ValidatedTag>, PolynomialFunctionModel<ValidatedTag>, PolynomialFunctionModel<ValidatedTag> );
}

void CheckFunctionConcept::check_mixed_operators_concept()
{
    // Function * Function -> Function
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ExactTag>>()+declval<Function<ExactTag>>()),Function<ExactTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ExactTag>>()+declval<Function<ValidatedTag>>()),Function<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ExactTag>>()+declval<Function<ApproximateTag>>()),Function<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ValidatedTag>>()+declval<Function<ValidatedTag>>()),Function<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ValidatedTag>>()+declval<Function<ValidatedTag>>()),Function<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ValidatedTag>>()+declval<Function<ApproximateTag>>()),Function<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ApproximateTag>>()+declval<Function<ExactTag>>()),Function<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ApproximateTag>>()+declval<Function<ValidatedTag>>()),Function<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ApproximateTag>>()+declval<Function<ApproximateTag>>()),Function<ApproximateTag>>::value) )

    // Function * Number -> Function
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ExactTag>>()+declval<Number<ExactTag>>()),Function<ExactTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ExactTag>>()+declval<Number<ValidatedTag>>()),Function<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ValidatedTag>>()+declval<Number<ExactTag>>()),Function<ValidatedTag>>::value) )

    // Number * Function -> Function
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<ExactTag>>()+declval<Function<ExactTag>>()),Function<ExactTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<ExactTag>>()+declval<Function<ValidatedTag>>()),Function<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<ValidatedTag>>()+declval<Function<ExactTag>>()),Function<ValidatedTag>>::value) )

    // SymbolicFunction * SymbolicFunction -> SymbolicFunction
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ExactTag>>()+declval<SymbolicFunction<ExactTag>>()),SymbolicFunction<ExactTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ExactTag>>()+declval<SymbolicFunction<ValidatedTag>>()),SymbolicFunction<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ValidatedTag>>()+declval<SymbolicFunction<ExactTag>>()),SymbolicFunction<ValidatedTag>>::value) )

    // SymbolicFunction * Number -> SymbolicFunction
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ExactTag>>()+declval<Number<ExactTag>>()),SymbolicFunction<ExactTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ExactTag>>()+declval<Number<ValidatedTag>>()),SymbolicFunction<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ValidatedTag>>()+declval<Number<ExactTag>>()),SymbolicFunction<ValidatedTag>>::value) )

    // Number * SymbolicFunction -> SymbolicFunction
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<ExactTag>>()+declval<SymbolicFunction<ExactTag>>()),SymbolicFunction<ExactTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<ExactTag>>()+declval<SymbolicFunction<ValidatedTag>>()),SymbolicFunction<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<ValidatedTag>>()+declval<SymbolicFunction<ExactTag>>()),SymbolicFunction<ValidatedTag>>::value) )

    // Function * SymbolicFunction -> Function
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ExactTag>>()+declval<SymbolicFunction<ExactTag>>()),Function<ExactTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ExactTag>>()+declval<SymbolicFunction<ValidatedTag>>()),Function<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ValidatedTag>>()+declval<SymbolicFunction<ExactTag>>()),Function<ValidatedTag>>::value) )

    // SymbolicFunction * Function -> Function
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ExactTag>>()+declval<Function<ExactTag>>()),Function<ExactTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ExactTag>>()+declval<Function<ValidatedTag>>()),Function<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ValidatedTag>>()+declval<Function<ExactTag>>()),Function<ValidatedTag>>::value) )


    // FunctionModel * FunctionModel -> FunctionModel; no mixed operators
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ValidatedTag>>()+declval<FunctionModel<ValidatedTag>>()),FunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ApproximateTag>>()+declval<FunctionModel<ValidatedTag>>()),FunctionModel<ApproximateTag>>::value) )

    // FunctionModel * Function -> FunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ValidatedTag>>()+declval<Function<ExactTag>>()),FunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ValidatedTag>>()+declval<Function<ValidatedTag>>()),FunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ApproximateTag>>()+declval<Function<ExactTag>>()),FunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ApproximateTag>>()+declval<Function<ValidatedTag>>()),FunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ApproximateTag>>()+declval<Function<ApproximateTag>>()),FunctionModel<ApproximateTag>>::value) )

    // Function * FunctionModel -> FunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ExactTag>>()+declval<FunctionModel<ValidatedTag>>()),FunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ValidatedTag>>()+declval<FunctionModel<ValidatedTag>>()),FunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ExactTag>>()+declval<FunctionModel<ApproximateTag>>()),FunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ValidatedTag>>()+declval<FunctionModel<ApproximateTag>>()),FunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ApproximateTag>>()+declval<FunctionModel<ApproximateTag>>()),FunctionModel<ApproximateTag>>::value) )

    // FunctionModel * SymbolicFunction -> FunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ValidatedTag>>()+declval<SymbolicFunction<ExactTag>>()),FunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ValidatedTag>>()+declval<SymbolicFunction<ValidatedTag>>()),FunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ApproximateTag>>()+declval<SymbolicFunction<ExactTag>>()),FunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ApproximateTag>>()+declval<SymbolicFunction<ValidatedTag>>()),FunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ApproximateTag>>()+declval<SymbolicFunction<ApproximateTag>>()),FunctionModel<ApproximateTag>>::value) )

    // SymbolicFunction * FunctionModel -> FunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ExactTag>>()+declval<FunctionModel<ValidatedTag>>()),FunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ValidatedTag>>()+declval<FunctionModel<ValidatedTag>>()),FunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ExactTag>>()+declval<FunctionModel<ApproximateTag>>()),FunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ValidatedTag>>()+declval<FunctionModel<ApproximateTag>>()),FunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ApproximateTag>>()+declval<FunctionModel<ApproximateTag>>()),FunctionModel<ApproximateTag>>::value) )

    // FunctionModel * Number -> FunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ValidatedTag>>()+declval<Number<ExactTag>>()),FunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ValidatedTag>>()+declval<Number<ValidatedTag>>()),FunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ApproximateTag>>()+declval<Number<ExactTag>>()),FunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ApproximateTag>>()+declval<Number<ValidatedTag>>()),FunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ApproximateTag>>()+declval<Number<ApproximateTag>>()),FunctionModel<ApproximateTag>>::value) )

    // Number * FunctionModel -> FunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<ExactTag>>()+declval<FunctionModel<ValidatedTag>>()),FunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<ValidatedTag>>()+declval<FunctionModel<ValidatedTag>>()),FunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<ExactTag>>()+declval<FunctionModel<ApproximateTag>>()),FunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<ValidatedTag>>()+declval<FunctionModel<ApproximateTag>>()),FunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<ApproximateTag>>()+declval<FunctionModel<ApproximateTag>>()),FunctionModel<ApproximateTag>>::value) )

    // PolynomialFunctionModel * PolynomialFunctionModel -> PolynomialFunctionModel; no mixed operators
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ValidatedTag>>()+declval<PolynomialFunctionModel<ValidatedTag>>()),PolynomialFunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ApproximateTag>>()+declval<PolynomialFunctionModel<ApproximateTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )
    // FunctionModel * PolynomialFunctionModel -> FunctionModel; no mixed operators
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ValidatedTag>>()+declval<PolynomialFunctionModel<ValidatedTag>>()),FunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<ApproximateTag>>()+declval<PolynomialFunctionModel<ApproximateTag>>()),FunctionModel<ApproximateTag>>::value) )

    // PolynomialFunctionModel * FunctionModel -> FunctionModel; no mixed operators
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ValidatedTag>>()+declval<FunctionModel<ValidatedTag>>()),FunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ApproximateTag>>()+declval<FunctionModel<ApproximateTag>>()),FunctionModel<ApproximateTag>>::value) )
    // Function * PolynomialFunctionModel -> PolynomialFunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ExactTag>>()+declval<PolynomialFunctionModel<ValidatedTag>>()),PolynomialFunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ValidatedTag>>()+declval<PolynomialFunctionModel<ValidatedTag>>()),PolynomialFunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ExactTag>>()+declval<PolynomialFunctionModel<ApproximateTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ValidatedTag>>()+declval<PolynomialFunctionModel<ApproximateTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<ApproximateTag>>()+declval<PolynomialFunctionModel<ApproximateTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )

    // PolynomialFunctionModel * Function -> PolynomialFunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ValidatedTag>>()+declval<Function<ExactTag>>()),PolynomialFunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ValidatedTag>>()+declval<Function<ValidatedTag>>()),PolynomialFunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ApproximateTag>>()+declval<Function<ExactTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ApproximateTag>>()+declval<Function<ValidatedTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ApproximateTag>>()+declval<Function<ApproximateTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )

    // SymbolicFunction * PolynomialFunctionModel -> PolynomialFunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ExactTag>>()+declval<PolynomialFunctionModel<ValidatedTag>>()),PolynomialFunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ValidatedTag>>()+declval<PolynomialFunctionModel<ValidatedTag>>()),PolynomialFunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ExactTag>>()+declval<PolynomialFunctionModel<ApproximateTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ValidatedTag>>()+declval<PolynomialFunctionModel<ApproximateTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<ApproximateTag>>()+declval<PolynomialFunctionModel<ApproximateTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )

    // PolynomialFunctionModel * SymbolicFunction -> PolynomialFunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ValidatedTag>>()+declval<SymbolicFunction<ExactTag>>()),PolynomialFunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ValidatedTag>>()+declval<SymbolicFunction<ValidatedTag>>()),PolynomialFunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ApproximateTag>>()+declval<SymbolicFunction<ExactTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ApproximateTag>>()+declval<SymbolicFunction<ValidatedTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ApproximateTag>>()+declval<SymbolicFunction<ApproximateTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )

    // Number * PolynomialFunctionModel -> PolynomialFunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<ExactTag>>()+declval<PolynomialFunctionModel<ValidatedTag>>()),PolynomialFunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<ValidatedTag>>()+declval<PolynomialFunctionModel<ValidatedTag>>()),PolynomialFunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<ExactTag>>()+declval<PolynomialFunctionModel<ApproximateTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<ValidatedTag>>()+declval<PolynomialFunctionModel<ApproximateTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<ApproximateTag>>()+declval<PolynomialFunctionModel<ApproximateTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )

    // PolynomialFunctionModel * Number -> PolynomialFunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ValidatedTag>>()+declval<Number<ExactTag>>()),PolynomialFunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ValidatedTag>>()+declval<Number<ValidatedTag>>()),PolynomialFunctionModel<ValidatedTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ApproximateTag>>()+declval<Number<ExactTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ApproximateTag>>()+declval<Number<ValidatedTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<ApproximateTag>>()+declval<Number<ApproximateTag>>()),PolynomialFunctionModel<ApproximateTag>>::value) )


    // SymbolicFunction -> SymbolicFunction
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(-declval<SymbolicFunction<ExactTag>>()),SymbolicFunction<ExactTag>>::value) )
    // Function -> Function
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(-declval<Function<ExactTag>>()),Function<ExactTag>>::value) )
    // FunctionModel -> FunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(-declval<FunctionModel<ValidatedTag>>()),FunctionModel<ValidatedTag>>::value) )
    // PolynomialFunctionModel -> PolynomialFunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(-declval<PolynomialFunctionModel<ValidatedTag>>()),PolynomialFunctionModel<ValidatedTag>>::value) )

}

OutputStream& operator<<(OutputStream& os, const std::type_info& info) {
    return os << "std::type_info(name=" << info.name() << ")"; }

template<class T> FunctionInterface<Void> const& reference(const T& f) {
    return static_cast<const FunctionInterface<Void>&>(f);
}

void CheckFunctionConcept::check_mixed_operators()
{
    Number<ExactTag> en(1);
    Number<ValidatedTag> vn(2);
    Number<ApproximateTag> an(3);

    Vector<Number<ExactTag>> ev(2);
    Vector<Number<ValidatedTag>> vv(2);
    Vector<Number<ApproximateTag>> av(2);

    Number<ExactTag> ec(1.25);
    Number<ValidatedTag> vc(2.25);
    Number<ApproximateTag> ac(3.25);

    SymbolicFunction<ExactTag> esf=SymbolicFunction<ExactTag>::constant(2,ec);
    SymbolicFunction<ValidatedTag> vsf=SymbolicFunction<ValidatedTag>::constant(2,vc);
    SymbolicFunction<ApproximateTag> asf=SymbolicFunction<ApproximateTag>::constant(2,ac);

    Sweeper swp;
    BoxDomain dom=BoxDomain(2,IntervalDomain(-1,+1));
    PolynomialFunctionModel<ValidatedTag> vtf(dom,swp);
    PolynomialFunctionModel<ApproximateTag> atf(dom,swp);
    Function<ValidatedTag> vvtf(vtf);
    Function<ApproximateTag> avtf(vtf);
    Function<ApproximateTag> aatf(atf);

    FunctionModel<ValidatedTag> vtfm(vtf);
    FunctionModel<ApproximateTag> atfm(atf);

    Function<ExactTag> eesf(esf);
    Function<ValidatedTag> vesf(eesf);
    Function<ApproximateTag> aesf(vesf);
    Function<ValidatedTag> vvsf(vsf);
    Function<ApproximateTag> avsf(vsf);
    Function<ApproximateTag> aasf(asf);

    // Simple tests for checking addition of concrete symbolic functions
    ARIADNE_TEST_EQUALS(typeid(reference(esf+esf)),typeid(esf));
    ARIADNE_TEST_EQUALS(typeid(reference(esf+vsf)),typeid(vsf));
    ARIADNE_TEST_EQUALS(typeid(reference(esf+en)),typeid(esf));
    ARIADNE_TEST_EQUALS(typeid(reference(esf+vn)),typeid(vsf));

    ARIADNE_TEST_EQUALS(typeid(reference(vsf+esf)),typeid(vsf));
    ARIADNE_TEST_EQUALS(typeid(reference(vsf+vsf)),typeid(vsf));
    ARIADNE_TEST_EQUALS(typeid(reference(vsf+en)),typeid(vsf));
    ARIADNE_TEST_EQUALS(typeid(reference(vsf+vn)),typeid(vsf));

    // Hard checks for static vs dynamic typing
    ARIADNE_ASSERT(typeid(reference(aesf+vvsf))==typeid(vsf) || typeid(reference(aesf+vvsf))==typeid(asf));
    if(typeid(reference(aesf+vvsf))==typeid(asf)) {
        ARIADNE_TEST_WARN("Additions of SymbolicFunctions through Function handle gives SymbolicFunction with paradigm determined by static typing.");
    }

    // Checks correct dynamic type for addition of symbolic functions via handle
    // Note that method currently dispatches on *static* type
    ARIADNE_TEST_EQUALS(typeid(reference(eesf+eesf)),typeid(esf));
    ARIADNE_TEST_CHECK_WARN(typeid(reference(eesf+vesf)),typeid(esf));
    ARIADNE_TEST_EQUALS(typeid(reference(eesf+vvsf)),typeid(vsf));
    ARIADNE_TEST_EQUALS(typeid(reference(eesf+esf)),typeid(esf));
    ARIADNE_TEST_EQUALS(typeid(reference(eesf+vsf)),typeid(vsf));
    ARIADNE_TEST_EQUALS(typeid(reference(eesf+en)),typeid(esf));
    ARIADNE_TEST_EQUALS(typeid(reference(eesf+vn)),typeid(vsf));

    ARIADNE_TEST_CHECK_WARN(typeid(reference(vesf+eesf)),typeid(esf));
    ARIADNE_TEST_EQUALS(typeid(reference(vesf+vesf)),typeid(vsf));
    ARIADNE_TEST_CHECK_WARN(typeid(reference(vesf+esf)),typeid(esf));
    ARIADNE_TEST_EQUALS(typeid(reference(vesf+vsf)),typeid(vsf));
    ARIADNE_TEST_CHECK_WARN(typeid(reference(vesf+en)),typeid(esf));
    ARIADNE_TEST_EQUALS(typeid(reference(vesf+vn)),typeid(vsf));

    ARIADNE_TEST_EQUALS(typeid(reference(vvsf+eesf)),typeid(vsf));
    ARIADNE_TEST_EQUALS(typeid(reference(vvsf+vesf)),typeid(vsf));
    ARIADNE_TEST_EQUALS(typeid(reference(vvsf+vvsf)),typeid(vsf));
    ARIADNE_TEST_EQUALS(typeid(reference(vvsf+esf)),typeid(vsf));
    ARIADNE_TEST_EQUALS(typeid(reference(vvsf+vsf)),typeid(vsf));
    ARIADNE_TEST_EQUALS(typeid(reference(vvsf+en)),typeid(vsf));
    ARIADNE_TEST_EQUALS(typeid(reference(vvsf+vn)),typeid(vsf));

    ARIADNE_ASSERT(typeid(reference(eesf+vesf))==typeid(esf) || typeid(reference(eesf+vesf))==typeid(vsf));

    // Check Taylor function model
    ARIADNE_TEST_EQUALS(typeid(reference(vtf+vtf)),typeid(vtf));
//AMBIGUOUS        ARIADNE_TEST_EQUALS(typeid(reference(vtf+vtfm)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vtf+vvtf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vtf+eesf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vtf+esf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vtf+vesf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vtf+vvsf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vtf+vsf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vtf+en)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vtf+vn)),typeid(vtf));

    // Check Taylor function through function model interface
//AMBIGUOUS    ARIADNE_TEST_EQUALS(typeid(reference(vtfm+vtf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vtfm+vvtf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vtfm+eesf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vtfm+vesf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vtfm+vvsf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vtfm+esf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vtfm+vsf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vtfm+en)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vtfm+vn)),typeid(vtf));

    // Check Taylor function model through function interface
    ARIADNE_TEST_EQUALS(typeid(reference(vvtf+vtf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vvtf+vvtf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vvtf+eesf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vvtf+vesf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vvtf+vvsf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vvtf+esf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vvtf+vsf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vvtf+vtf)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vvtf+en)),typeid(vtf));
    ARIADNE_TEST_EQUALS(typeid(reference(vvtf+vn)),typeid(vtf));


}

void CheckFunctionConcept::check_mixed_evaluation()
{
    Number<ExactTag> en(1);
    Number<ValidatedTag> vn(2);
    Number<ApproximateTag> an(3);

    Vector<Number<ExactTag>> ev(2);
    Vector<Number<ValidatedTag>> vv(2);
    Vector<Number<ApproximateTag>> av(2);

    Number<ExactTag> ec(1.25);
    Number<ValidatedTag> vc(2.25);
    Number<ApproximateTag> ac(3.25);

    SymbolicFunction<ExactTag> esf=SymbolicFunction<ExactTag>::constant(2,ec);
    SymbolicFunction<ValidatedTag> vsf=SymbolicFunction<ValidatedTag>::constant(2,vc);
    SymbolicFunction<ApproximateTag> asf=SymbolicFunction<ApproximateTag>::constant(2,ac);

    Sweeper swp;
    BoxDomain dom=BoxDomain(2,IntervalDomain(-1,+1));
    PolynomialFunctionModel<ValidatedTag> vtf(dom,swp);
    PolynomialFunctionModel<ApproximateTag> atf(dom,swp);
    Function<ValidatedTag> vvtf(vtf);
    Function<ApproximateTag> avtf(vtf);
    Function<ApproximateTag> aatf(atf);

    FunctionModel<ValidatedTag> vtfm(vtf);
    FunctionModel<ApproximateTag> atfm(atf);

    Function<ExactTag> eesf(esf);
    Function<ValidatedTag> vesf(eesf);
    Function<ApproximateTag> aesf(vesf);
    Function<ValidatedTag> vvsf(vsf);
    Function<ApproximateTag> avsf(vsf);
    Function<ApproximateTag> aasf(asf);

    ARIADNE_TEST_EXECUTE(evaluate(eesf,ev));
    ARIADNE_TEST_EXECUTE(evaluate(eesf,vv));
    ARIADNE_TEST_EXECUTE(evaluate(eesf,av));
    ARIADNE_TEST_EXECUTE(evaluate(vesf,ev));
    ARIADNE_TEST_EXECUTE(evaluate(vesf,vv));
    ARIADNE_TEST_EXECUTE(evaluate(vesf,av));
    ARIADNE_TEST_EXECUTE(evaluate(vvsf,ev));
    ARIADNE_TEST_EXECUTE(evaluate(vvsf,vv));
    ARIADNE_TEST_EXECUTE(evaluate(vvsf,av));

    ARIADNE_TEST_EXECUTE(evaluate(vvtf,ev));
    ARIADNE_TEST_EXECUTE(evaluate(vvtf,vv));
    ARIADNE_TEST_EXECUTE(evaluate(vvtf,av));
}

void CheckFunctionConcept::check_scalar_function()
{
    // Check constructors of constants and coordinates
    ARIADNE_TEST_NAMED_CONSTRUCT(ExactFunction,o,constant(3,1.0_x));
    ARIADNE_TEST_NAMED_CONSTRUCT(ExactFunction,x,coordinate(3,0));
    ARIADNE_TEST_NAMED_CONSTRUCT(ExactFunction,y,coordinate(3,1));

    // Check constructors of constants and coordinates give elements of type SymbolicFunction
    ARIADNE_TEST_ASSERT(dynamic_cast<const ExactSymbolicFunction*>(o.raw_pointer()));
    ARIADNE_TEST_ASSERT(dynamic_cast<const ExactSymbolicFunction*>(x.raw_pointer()));
    ARIADNE_TEST_ASSERT(dynamic_cast<const ExactSymbolicFunction*>(y.raw_pointer()));

    // Check type hierarchy of function interface
    ARIADNE_TEST_ASSERT(dynamic_cast<const ExactFunctionInterface*>(o.raw_pointer()));
    ARIADNE_TEST_ASSERT(dynamic_cast<const ValidatedFunctionInterface*>(o.raw_pointer()));
    ARIADNE_TEST_ASSERT(dynamic_cast<const ApproximateFunctionInterface*>(o.raw_pointer()));

    // Check creation of arithmetical functions
    ARIADNE_TEST_CONSTRUCT(ExactFunction,f,(o+x*y));
    ARIADNE_TEST_CONSTRUCT(Vector<ApproximateFloat64>,p,({2.0,3.0,5.0}));
    ARIADNE_TEST_EQUALS(f(p),7.0);

    // Check evaluation on vectors
    ARIADNE_TEST_EQUALS(f(Vector<Real>{2,3,5}),7);
    ARIADNE_TEST_EQUALS(f(Vector<BoundedFloat64>{2.0_x,3.0_x,5.0_x}),7.0_x);
    ARIADNE_TEST_EQUALS(f(Vector<ApproximateFloat64>{2.0,3.0,5.0}),7.0);

    ARIADNE_TEST_EQUALS(evaluate(f,Vector<Real>{2,3,5}),7);
    ARIADNE_TEST_EQUALS(evaluate(f,Vector<BoundedFloat64>{2.0_x,3.0_x,5.0_x}),7.0_x);
    ARIADNE_TEST_EQUALS(evaluate(f,Vector<ApproximateFloat64>{2.0,3.0,5.0}),7.0);

    // Check mixed arithmetic
    Vector<ApproximateFloat64> z(3);
    ARIADNE_TEST_NAMED_CONSTRUCT(ExactFunction,ef,constant(3,1));
    ARIADNE_TEST_NAMED_CONSTRUCT(ValidatedFunction,vf,constant(3,2));
    ARIADNE_TEST_NAMED_CONSTRUCT(ApproximateFunction,af,constant(3,4));
    ARIADNE_TEST_STATIC_ASSERT((IsSame<decltype(ef+ef),ExactFunction>::value));
    ARIADNE_TEST_STATIC_ASSERT((IsSame<decltype(ef+vf),ValidatedFunction>::value));
    ARIADNE_TEST_STATIC_ASSERT((IsSame<decltype(ef+af),ApproximateFunction>::value));
    ARIADNE_TEST_STATIC_ASSERT((IsSame<decltype(vf+ef),ValidatedFunction>::value));
    ARIADNE_TEST_STATIC_ASSERT((IsSame<decltype(vf+vf),ValidatedFunction>::value));
    ARIADNE_TEST_STATIC_ASSERT((IsSame<decltype(vf+af),ApproximateFunction>::value));
    ARIADNE_TEST_STATIC_ASSERT((IsSame<decltype(af+ef),ApproximateFunction>::value));
    ARIADNE_TEST_STATIC_ASSERT((IsSame<decltype(af+vf),ApproximateFunction>::value));
    ARIADNE_TEST_STATIC_ASSERT((IsSame<decltype(af+af),ApproximateFunction>::value));
    ARIADNE_TEST_ASSERT(!dynamic_cast<const ExactFunctionInterface*>((ef+vf).raw_pointer()));
    ARIADNE_TEST_ASSERT(dynamic_cast<const ValidatedFunctionInterface*>((ef+vf).raw_pointer()));
    ARIADNE_TEST_ASSERT(dynamic_cast<const ApproximateFunctionInterface*>((ef+vf).raw_pointer()));
    ARIADNE_TEST_ASSERT(!dynamic_cast<const ExactFunctionInterface*>((ef+af).raw_pointer()));
    ARIADNE_TEST_ASSERT(!dynamic_cast<const ValidatedFunctionInterface*>((ef+af).raw_pointer()));
    ARIADNE_TEST_ASSERT(dynamic_cast<const ApproximateFunctionInterface*>((ef+af).raw_pointer()));
    ARIADNE_TEST_ASSERT(!dynamic_cast<const ExactFunctionInterface*>((vf+af).raw_pointer()));
    ARIADNE_TEST_ASSERT(!dynamic_cast<const ValidatedFunctionInterface*>((vf+af).raw_pointer()));
    ARIADNE_TEST_ASSERT(dynamic_cast<const ApproximateFunctionInterface*>((vf+af).raw_pointer()));
    ARIADNE_TEST_EQUALS(evaluate(ef+vf,z),3);
    ARIADNE_TEST_EXECUTE(ef+af);
    ARIADNE_TEST_EXECUTE(vf+ef);
    ARIADNE_TEST_EXECUTE(vf+af);
    ARIADNE_TEST_EXECUTE(af+ef);
    ARIADNE_TEST_EXECUTE(af+vf);

    // Check unary operators
    ARIADNE_TEST_EXECUTE(neg(f));
    ARIADNE_TEST_EXECUTE(rec(f));
    ARIADNE_TEST_EXECUTE(pow(f,5u));
    ARIADNE_TEST_EXECUTE(pow(f,-3));
    ARIADNE_TEST_EXECUTE(sqrt(f));
    ARIADNE_TEST_EXECUTE(exp(f));
    ARIADNE_TEST_EXECUTE(log(f));
    ARIADNE_TEST_EXECUTE(sin(f));
    ARIADNE_TEST_EXECUTE(cos(f));
    ARIADNE_TEST_EXECUTE(tan(f));

    ExactFunction df=f.derivative(1);
    ARIADNE_TEST_PRINT(df);
    ARIADNE_TEST_EQUAL(df(p),2.0);

    f=(1-x*x-y/2);

}

void CheckFunctionConcept::check_vector_function()
{
    ARIADNE_TEST_NAMED_CONSTRUCT(ExactVectorFunction,f,identity(3));

    ARIADNE_TEST_CONSTRUCT(ExactFunction,fa0,(f.get(0)));
    ARIADNE_TEST_CONSTRUCT(ExactFunction,fa1,(f[1]));

    // Regression tests for element
    ARIADNE_TEST_PRINT(f[0]);
    const ExactVectorFunction& fcr=f;
    ARIADNE_TEST_PRINT(fcr[0]);
    ExactVectorFunction& fr=f;
    ARIADNE_TEST_PRINT(fr[0]);

    ARIADNE_TEST_EQUAL(f[0](Vector<ApproximateFloat64>{2.0,3.0,5.0}),2.0);
    ARIADNE_TEST_EXECUTE(f[0]=f[1]);
    ARIADNE_TEST_EQUAL(f[0](Vector<ApproximateFloat64>{2.0,3.0,5.0}),3.0);

    ExactVectorFunction x=ExactVectorFunction::identity(2);
    ExactFunction x0=x[0];
    //ExactFunction fvi=(1-x[0]*x[0]-x[1]/2);

}





void CheckFunctionConcept::check_conversions()
{
    Function<ExactTag> ef(2);
    Function<ValidatedTag> vf(2);
    Function<ApproximateTag> af(2);

    vf=Function<ValidatedTag>(ef);
    af=Function<ApproximateTag>(ef);
    af=Function<ApproximateTag>(vf);

    vf=ef;
    af=ef;
    af=vf;

    std::stringstream cnul;

    const FunctionInterface<ExactTag>& efie=ef;
    const FunctionInterface<ValidatedTag>& vfie=ef;
    const FunctionInterface<ValidatedTag>& vfiv=vf;
    const FunctionInterface<ApproximateTag>& afie=ef;
    const FunctionInterface<ApproximateTag>& afiv=vf;
    const FunctionInterface<ApproximateTag>& afia=af;
    cnul << efie << vfie << vfiv << afie << afiv << afia;

    const ScalarEvaluableInterface<ExactTag>& egie=ef;
    const ScalarEvaluableInterface<ValidatedTag>& vgie=vf;
    const ScalarEvaluableInterface<ValidatedTag>& vgiv=vf;
    const ScalarEvaluableInterface<ApproximateTag>& agie=ef;
    const ScalarEvaluableInterface<ApproximateTag>& agiv=vf;
    const ScalarEvaluableInterface<ApproximateTag>& agia=af;
    cnul << egie << vgie << vgiv << agie << agiv << agia;
}

void CheckFunctionConcept::check_differentiation()
{
    ExactFunction z=ExactFunction::constant(2,0.0_x);
    ExactFunction o=ExactFunction::constant(2,1.0_x);
    ExactFunction x=ExactFunction::coordinate(2,0);
    ExactFunction y=ExactFunction::coordinate(2,1);

    ExactFunction af=3*x-2*y+1;
    ExactFunction daf=af.derivative(1);
    ARIADNE_TEST_EQUAL(evaluate(daf,Vector<ApproximateFloat64>{2.4,1.3}),-2.0);

    ARIADNE_TEST_EQUAL(evaluate(x.derivative(0),Vector<ApproximateFloat64>{2.4,1.3}),1.0);
    ARIADNE_TEST_EQUAL(evaluate(x.derivative(0),Vector<BoundedFloat64>{2.4_x,1.3_x}),1.0_x);
    ARIADNE_TEST_EQUAL(evaluate(x.derivative(1),Vector<ApproximateFloat64>{2.4, 1.3}),0.0);

}
*/

