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

#include "expression/operators.h"
#include "function/function.h"

namespace Ariadne {
template<class X> String class_name();
template<> String class_name<Approximate>() { return "Approximate"; }
template<> String class_name<Validated>() { return "Validated"; }
template<> String class_name<Effective>() { return "Effective"; }
} // namespace Ariadne

using namespace Ariadne;

typedef ExactInterval IntervalDomain;
typedef ExactBox BoxDomain;
template<class P> using Function = ScalarFunction<P>;
typedef ScalarFunction<Effective> EffectiveFunction;

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

    if(IsWeaker<Approximate,P>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<ApproximateNumber>, Return<ApproximateNumber>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<ApproximateFloat64>, Return<ApproximateFloat64>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<ApproximateNumber>, Return<ApproximateNumber>>);
    }

    if(IsWeaker<Validated,P>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<ValidatedNumber>, Return<ValidatedNumber>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<BoundedFloat64>, Return<BoundedFloat64>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<ValidatedFloat64>, Return<ValidatedFloat64>>);
    }

    if(IsWeaker<Effective,P>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<EffectiveNumber>, Return<EffectiveNumber>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<Real>, Return<Real>>);
    }

}

template<class F> void CheckFunctionConcept<F>::check_differentiable_concept()
{
    typedef typename F::Paradigm P;

    this->check_differential_algebra_concept();

    ARIADNE_TEST_STATIC_ASSERT(HasDerivative<F,SizeType,Return<F>>);

    if(IsWeaker<Approximate,P>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<Differential<ApproximateFloat64>>, Return<Differential<ApproximateFloat64>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasGradientMethod<F,Vector<ApproximateFloat64>, Return<Covector<ApproximateFloat64>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasDifferentialMethod<F,Vector<ApproximateFloat64>, Return<Differential<ApproximateFloat64>>>);
    }

    if(IsWeaker<Validated,P>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<Differential<BoundedFloat64>>, Return<Differential<BoundedFloat64>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasGradientMethod<F,Vector<BoundedFloat64>, Return<Covector<BoundedFloat64>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasDifferentialMethod<F,Vector<BoundedFloat64>, Return<Differential<BoundedFloat64>>>);
    }

    if(IsWeaker<Effective,P>::value) {
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

    if(IsWeaker<P,Approximate>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCompose<ScalarFunction<Approximate>,Vector<F>, Return<F>>);
    }

    if(IsWeaker<P,Validated>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCompose<ScalarFunction<Validated>,Vector<F>, Return<F>>);
    }

    if(IsWeaker<P,Effective>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCompose<ScalarFunction<Effective>,Vector<F>, Return<F>>);
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

    if(IsWeaker<Approximate,P>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<ApproximateNumber>, Return<Vector<ApproximateNumber>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<ApproximateFloat64>, Return<Vector<ApproximateFloat64>>>);

        ARIADNE_TEST_STATIC_ASSERT(HasEvaluate<F,Vector<ApproximateNumber>, Return<Vector<ApproximateNumber>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasEvaluate<F,Vector<ApproximateFloat64>, Return<Vector<ApproximateFloat64>>>);
    }

    if(IsWeaker<Validated,P>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<ValidatedNumber>, Return<Vector<ValidatedNumber>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<BoundedFloat64>, Return<Vector<BoundedFloat64>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<ValidatedFloat64>, Return<Vector<ValidatedFloat64>>>);

        ARIADNE_TEST_STATIC_ASSERT(HasEvaluate<F,Vector<ValidatedNumber>, Return<Vector<ValidatedNumber>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasEvaluate<F,Vector<BoundedFloat64>, Return<Vector<BoundedFloat64>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasEvaluate<F,Vector<ValidatedFloat64>, Return<Vector<ValidatedFloat64>>>);

    }

    if(IsWeaker<Effective,P>::value) {
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<EffectiveNumber>, Return<Vector<EffectiveNumber>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasCallMethod<F,Vector<Real>, Return<Vector<Real>>>);

        ARIADNE_TEST_STATIC_ASSERT(HasEvaluate<F,Vector<EffectiveNumber>, Return<Vector<EffectiveNumber>>>);
        ARIADNE_TEST_STATIC_ASSERT(HasEvaluate<F,Vector<Real>, Return<Vector<Real>>>);

    }

}


/*
template<class T1, class T2> inline const bool same_type(const T1& t1, const T2& t2) { return IsSame<T1,T2>::value; }

void CheckFunctionConcept::check_operators_concept()
{
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( SymbolicFunction<Validated>, SymbolicFunction<Validated>, SymbolicFunction<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( SymbolicFunction<Validated>, VectorFunctionElement<Validated>, Function<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( SymbolicFunction<Validated>, Function<Validated>, Function<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( SymbolicFunction<Validated>, FunctionModel<Validated>, FunctionModel<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( SymbolicFunction<Validated>, VectorFunctionModelElement<Validated>, FunctionModel<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( SymbolicFunction<Validated>, PolynomialFunctionModel<Validated>, PolynomialFunctionModel<Validated> );

    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionElement<Validated>, SymbolicFunction<Validated>, Function<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionElement<Validated>, VectorFunctionElement<Validated>, Function<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionElement<Validated>, Function<Validated>, Function<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionElement<Validated>, FunctionModel<Validated>, FunctionModel<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionElement<Validated>, VectorFunctionModelElement<Validated>, FunctionModel<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionElement<Validated>, PolynomialFunctionModel<Validated>, PolynomialFunctionModel<Validated> );

    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( Function<Validated>, SymbolicFunction<Validated>, Function<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( Function<Validated>, VectorFunctionElement<Validated>, Function<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( Function<Validated>, Function<Validated>, Function<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( Function<Validated>, FunctionModel<Validated>, FunctionModel<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( Function<Validated>, VectorFunctionModelElement<Validated>, FunctionModel<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( Function<Validated>, PolynomialFunctionModel<Validated>, PolynomialFunctionModel<Validated> );

    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( FunctionModel<Validated>, SymbolicFunction<Validated>, FunctionModel<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( FunctionModel<Validated>, VectorFunctionElement<Validated>, FunctionModel<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( FunctionModel<Validated>, Function<Validated>, FunctionModel<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( FunctionModel<Validated>, FunctionModel<Validated>, FunctionModel<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( FunctionModel<Validated>, VectorFunctionModelElement<Validated>, FunctionModel<Validated> );
    //ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( FunctionModel<Validated>, PolynomialFunctionModel<Validated>, PolynomialFunctionModel<Validated> );

    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionModelElement<Validated>, SymbolicFunction<Validated>, FunctionModel<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionModelElement<Validated>, VectorFunctionElement<Validated>, FunctionModel<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionModelElement<Validated>, Function<Validated>, FunctionModel<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionModelElement<Validated>, FunctionModel<Validated>, FunctionModel<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionModelElement<Validated>, VectorFunctionModelElement<Validated>, FunctionModel<Validated> );
    //ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( VectorFunctionModelElement<Validated>, PolynomialFunctionModel<Validated>, PolynomialFunctionModel<Validated> );

    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( PolynomialFunctionModel<Validated>, SymbolicFunction<Validated>, PolynomialFunctionModel<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( PolynomialFunctionModel<Validated>, VectorFunctionElement<Validated>, PolynomialFunctionModel<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( PolynomialFunctionModel<Validated>, Function<Validated>, PolynomialFunctionModel<Validated> );
    //ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( PolynomialFunctionModel<Validated>, FunctionModel<Validated>, PolynomialFunctionModel<Validated> );
    //ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( PolynomialFunctionModel<Validated>, VectorFunctionModelElement<Validated>, PolynomialFunctionModel<Validated> );
    ARIADNE_TEST_ARITHMETIC_RESULT_TYPE( PolynomialFunctionModel<Validated>, PolynomialFunctionModel<Validated>, PolynomialFunctionModel<Validated> );
}

void CheckFunctionConcept::check_mixed_operators_concept()
{
    // Function * Function -> Function
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Exact>>()+declval<Function<Exact>>()),Function<Exact>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Exact>>()+declval<Function<Validated>>()),Function<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Exact>>()+declval<Function<Approximate>>()),Function<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Validated>>()+declval<Function<Validated>>()),Function<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Validated>>()+declval<Function<Validated>>()),Function<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Validated>>()+declval<Function<Approximate>>()),Function<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Approximate>>()+declval<Function<Exact>>()),Function<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Approximate>>()+declval<Function<Validated>>()),Function<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Approximate>>()+declval<Function<Approximate>>()),Function<Approximate>>::value) )

    // Function * Number -> Function
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Exact>>()+declval<Number<Exact>>()),Function<Exact>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Exact>>()+declval<Number<Validated>>()),Function<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Validated>>()+declval<Number<Exact>>()),Function<Validated>>::value) )

    // Number * Function -> Function
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<Exact>>()+declval<Function<Exact>>()),Function<Exact>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<Exact>>()+declval<Function<Validated>>()),Function<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<Validated>>()+declval<Function<Exact>>()),Function<Validated>>::value) )

    // SymbolicFunction * SymbolicFunction -> SymbolicFunction
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Exact>>()+declval<SymbolicFunction<Exact>>()),SymbolicFunction<Exact>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Exact>>()+declval<SymbolicFunction<Validated>>()),SymbolicFunction<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Validated>>()+declval<SymbolicFunction<Exact>>()),SymbolicFunction<Validated>>::value) )

    // SymbolicFunction * Number -> SymbolicFunction
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Exact>>()+declval<Number<Exact>>()),SymbolicFunction<Exact>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Exact>>()+declval<Number<Validated>>()),SymbolicFunction<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Validated>>()+declval<Number<Exact>>()),SymbolicFunction<Validated>>::value) )

    // Number * SymbolicFunction -> SymbolicFunction
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<Exact>>()+declval<SymbolicFunction<Exact>>()),SymbolicFunction<Exact>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<Exact>>()+declval<SymbolicFunction<Validated>>()),SymbolicFunction<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<Validated>>()+declval<SymbolicFunction<Exact>>()),SymbolicFunction<Validated>>::value) )

    // Function * SymbolicFunction -> Function
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Exact>>()+declval<SymbolicFunction<Exact>>()),Function<Exact>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Exact>>()+declval<SymbolicFunction<Validated>>()),Function<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Validated>>()+declval<SymbolicFunction<Exact>>()),Function<Validated>>::value) )

    // SymbolicFunction * Function -> Function
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Exact>>()+declval<Function<Exact>>()),Function<Exact>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Exact>>()+declval<Function<Validated>>()),Function<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Validated>>()+declval<Function<Exact>>()),Function<Validated>>::value) )


    // FunctionModel * FunctionModel -> FunctionModel; no mixed operators
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Validated>>()+declval<FunctionModel<Validated>>()),FunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Approximate>>()+declval<FunctionModel<Validated>>()),FunctionModel<Approximate>>::value) )

    // FunctionModel * Function -> FunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Validated>>()+declval<Function<Exact>>()),FunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Validated>>()+declval<Function<Validated>>()),FunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Approximate>>()+declval<Function<Exact>>()),FunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Approximate>>()+declval<Function<Validated>>()),FunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Approximate>>()+declval<Function<Approximate>>()),FunctionModel<Approximate>>::value) )

    // Function * FunctionModel -> FunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Exact>>()+declval<FunctionModel<Validated>>()),FunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Validated>>()+declval<FunctionModel<Validated>>()),FunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Exact>>()+declval<FunctionModel<Approximate>>()),FunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Validated>>()+declval<FunctionModel<Approximate>>()),FunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Approximate>>()+declval<FunctionModel<Approximate>>()),FunctionModel<Approximate>>::value) )

    // FunctionModel * SymbolicFunction -> FunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Validated>>()+declval<SymbolicFunction<Exact>>()),FunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Validated>>()+declval<SymbolicFunction<Validated>>()),FunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Approximate>>()+declval<SymbolicFunction<Exact>>()),FunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Approximate>>()+declval<SymbolicFunction<Validated>>()),FunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Approximate>>()+declval<SymbolicFunction<Approximate>>()),FunctionModel<Approximate>>::value) )

    // SymbolicFunction * FunctionModel -> FunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Exact>>()+declval<FunctionModel<Validated>>()),FunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Validated>>()+declval<FunctionModel<Validated>>()),FunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Exact>>()+declval<FunctionModel<Approximate>>()),FunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Validated>>()+declval<FunctionModel<Approximate>>()),FunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Approximate>>()+declval<FunctionModel<Approximate>>()),FunctionModel<Approximate>>::value) )

    // FunctionModel * Number -> FunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Validated>>()+declval<Number<Exact>>()),FunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Validated>>()+declval<Number<Validated>>()),FunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Approximate>>()+declval<Number<Exact>>()),FunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Approximate>>()+declval<Number<Validated>>()),FunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Approximate>>()+declval<Number<Approximate>>()),FunctionModel<Approximate>>::value) )

    // Number * FunctionModel -> FunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<Exact>>()+declval<FunctionModel<Validated>>()),FunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<Validated>>()+declval<FunctionModel<Validated>>()),FunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<Exact>>()+declval<FunctionModel<Approximate>>()),FunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<Validated>>()+declval<FunctionModel<Approximate>>()),FunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<Approximate>>()+declval<FunctionModel<Approximate>>()),FunctionModel<Approximate>>::value) )

    // PolynomialFunctionModel * PolynomialFunctionModel -> PolynomialFunctionModel; no mixed operators
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Validated>>()+declval<PolynomialFunctionModel<Validated>>()),PolynomialFunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Approximate>>()+declval<PolynomialFunctionModel<Approximate>>()),PolynomialFunctionModel<Approximate>>::value) )
    // FunctionModel * PolynomialFunctionModel -> FunctionModel; no mixed operators
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Validated>>()+declval<PolynomialFunctionModel<Validated>>()),FunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<FunctionModel<Approximate>>()+declval<PolynomialFunctionModel<Approximate>>()),FunctionModel<Approximate>>::value) )

    // PolynomialFunctionModel * FunctionModel -> FunctionModel; no mixed operators
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Validated>>()+declval<FunctionModel<Validated>>()),FunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Approximate>>()+declval<FunctionModel<Approximate>>()),FunctionModel<Approximate>>::value) )
    // Function * PolynomialFunctionModel -> PolynomialFunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Exact>>()+declval<PolynomialFunctionModel<Validated>>()),PolynomialFunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Validated>>()+declval<PolynomialFunctionModel<Validated>>()),PolynomialFunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Exact>>()+declval<PolynomialFunctionModel<Approximate>>()),PolynomialFunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Validated>>()+declval<PolynomialFunctionModel<Approximate>>()),PolynomialFunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Function<Approximate>>()+declval<PolynomialFunctionModel<Approximate>>()),PolynomialFunctionModel<Approximate>>::value) )

    // PolynomialFunctionModel * Function -> PolynomialFunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Validated>>()+declval<Function<Exact>>()),PolynomialFunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Validated>>()+declval<Function<Validated>>()),PolynomialFunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Approximate>>()+declval<Function<Exact>>()),PolynomialFunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Approximate>>()+declval<Function<Validated>>()),PolynomialFunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Approximate>>()+declval<Function<Approximate>>()),PolynomialFunctionModel<Approximate>>::value) )

    // SymbolicFunction * PolynomialFunctionModel -> PolynomialFunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Exact>>()+declval<PolynomialFunctionModel<Validated>>()),PolynomialFunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Validated>>()+declval<PolynomialFunctionModel<Validated>>()),PolynomialFunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Exact>>()+declval<PolynomialFunctionModel<Approximate>>()),PolynomialFunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Validated>>()+declval<PolynomialFunctionModel<Approximate>>()),PolynomialFunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<SymbolicFunction<Approximate>>()+declval<PolynomialFunctionModel<Approximate>>()),PolynomialFunctionModel<Approximate>>::value) )

    // PolynomialFunctionModel * SymbolicFunction -> PolynomialFunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Validated>>()+declval<SymbolicFunction<Exact>>()),PolynomialFunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Validated>>()+declval<SymbolicFunction<Validated>>()),PolynomialFunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Approximate>>()+declval<SymbolicFunction<Exact>>()),PolynomialFunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Approximate>>()+declval<SymbolicFunction<Validated>>()),PolynomialFunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Approximate>>()+declval<SymbolicFunction<Approximate>>()),PolynomialFunctionModel<Approximate>>::value) )

    // Number * PolynomialFunctionModel -> PolynomialFunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<Exact>>()+declval<PolynomialFunctionModel<Validated>>()),PolynomialFunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<Validated>>()+declval<PolynomialFunctionModel<Validated>>()),PolynomialFunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<Exact>>()+declval<PolynomialFunctionModel<Approximate>>()),PolynomialFunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<Validated>>()+declval<PolynomialFunctionModel<Approximate>>()),PolynomialFunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<Number<Approximate>>()+declval<PolynomialFunctionModel<Approximate>>()),PolynomialFunctionModel<Approximate>>::value) )

    // PolynomialFunctionModel * Number -> PolynomialFunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Validated>>()+declval<Number<Exact>>()),PolynomialFunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Validated>>()+declval<Number<Validated>>()),PolynomialFunctionModel<Validated>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Approximate>>()+declval<Number<Exact>>()),PolynomialFunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Approximate>>()+declval<Number<Validated>>()),PolynomialFunctionModel<Approximate>>::value) )
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(declval<PolynomialFunctionModel<Approximate>>()+declval<Number<Approximate>>()),PolynomialFunctionModel<Approximate>>::value) )


    // SymbolicFunction -> SymbolicFunction
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(-declval<SymbolicFunction<Exact>>()),SymbolicFunction<Exact>>::value) )
    // Function -> Function
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(-declval<Function<Exact>>()),Function<Exact>>::value) )
    // FunctionModel -> FunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(-declval<FunctionModel<Validated>>()),FunctionModel<Validated>>::value) )
    // PolynomialFunctionModel -> PolynomialFunctionModel
    ARIADNE_TEST_STATIC_ASSERT( (IsSame<decltype(-declval<PolynomialFunctionModel<Validated>>()),PolynomialFunctionModel<Validated>>::value) )

}

OutputStream& operator<<(OutputStream& os, const std::type_info& info) {
    return os << "std::type_info(name=" << info.name() << ")"; }

template<class T> FunctionInterface<Void> const& reference(const T& f) {
    return static_cast<const FunctionInterface<Void>&>(f);
}

void CheckFunctionConcept::check_mixed_operators()
{
    Number<Exact> en(1);
    Number<Validated> vn(2);
    Number<Approximate> an(3);

    Vector<Number<Exact>> ev(2);
    Vector<Number<Validated>> vv(2);
    Vector<Number<Approximate>> av(2);

    Number<Exact> ec(1.25);
    Number<Validated> vc(2.25);
    Number<Approximate> ac(3.25);

    SymbolicFunction<Exact> esf=SymbolicFunction<Exact>::constant(2,ec);
    SymbolicFunction<Validated> vsf=SymbolicFunction<Validated>::constant(2,vc);
    SymbolicFunction<Approximate> asf=SymbolicFunction<Approximate>::constant(2,ac);

    Sweeper swp;
    BoxDomain dom=BoxDomain(2,IntervalDomain(-1,+1));
    PolynomialFunctionModel<Validated> vtf(dom,swp);
    PolynomialFunctionModel<Approximate> atf(dom,swp);
    Function<Validated> vvtf(vtf);
    Function<Approximate> avtf(vtf);
    Function<Approximate> aatf(atf);

    FunctionModel<Validated> vtfm(vtf);
    FunctionModel<Approximate> atfm(atf);

    Function<Exact> eesf(esf);
    Function<Validated> vesf(eesf);
    Function<Approximate> aesf(vesf);
    Function<Validated> vvsf(vsf);
    Function<Approximate> avsf(vsf);
    Function<Approximate> aasf(asf);

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
    Number<Exact> en(1);
    Number<Validated> vn(2);
    Number<Approximate> an(3);

    Vector<Number<Exact>> ev(2);
    Vector<Number<Validated>> vv(2);
    Vector<Number<Approximate>> av(2);

    Number<Exact> ec(1.25);
    Number<Validated> vc(2.25);
    Number<Approximate> ac(3.25);

    SymbolicFunction<Exact> esf=SymbolicFunction<Exact>::constant(2,ec);
    SymbolicFunction<Validated> vsf=SymbolicFunction<Validated>::constant(2,vc);
    SymbolicFunction<Approximate> asf=SymbolicFunction<Approximate>::constant(2,ac);

    Sweeper swp;
    BoxDomain dom=BoxDomain(2,IntervalDomain(-1,+1));
    PolynomialFunctionModel<Validated> vtf(dom,swp);
    PolynomialFunctionModel<Approximate> atf(dom,swp);
    Function<Validated> vvtf(vtf);
    Function<Approximate> avtf(vtf);
    Function<Approximate> aatf(atf);

    FunctionModel<Validated> vtfm(vtf);
    FunctionModel<Approximate> atfm(atf);

    Function<Exact> eesf(esf);
    Function<Validated> vesf(eesf);
    Function<Approximate> aesf(vesf);
    Function<Validated> vvsf(vsf);
    Function<Approximate> avsf(vsf);
    Function<Approximate> aasf(asf);

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
    ARIADNE_TEST_EQUALS(f(Vector<ValidatedFloat64>{2.0_x,3.0_x,5.0_x}),7.0_x);
    ARIADNE_TEST_EQUALS(f(Vector<ApproximateFloat64>{2.0,3.0,5.0}),7.0);

    ARIADNE_TEST_EQUALS(evaluate(f,Vector<Real>{2,3,5}),7);
    ARIADNE_TEST_EQUALS(evaluate(f,Vector<ValidatedFloat64>{2.0_x,3.0_x,5.0_x}),7.0_x);
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
    Function<Exact> ef(2);
    Function<Validated> vf(2);
    Function<Approximate> af(2);

    vf=Function<Validated>(ef);
    af=Function<Approximate>(ef);
    af=Function<Approximate>(vf);

    vf=ef;
    af=ef;
    af=vf;

    std::stringstream cnul;

    const FunctionInterface<Exact>& efie=ef;
    const FunctionInterface<Validated>& vfie=ef;
    const FunctionInterface<Validated>& vfiv=vf;
    const FunctionInterface<Approximate>& afie=ef;
    const FunctionInterface<Approximate>& afiv=vf;
    const FunctionInterface<Approximate>& afia=af;
    cnul << efie << vfie << vfiv << afie << afiv << afia;

    const ScalarEvaluableInterface<Exact>& egie=ef;
    const ScalarEvaluableInterface<Validated>& vgie=vf;
    const ScalarEvaluableInterface<Validated>& vgiv=vf;
    const ScalarEvaluableInterface<Approximate>& agie=ef;
    const ScalarEvaluableInterface<Approximate>& agiv=vf;
    const ScalarEvaluableInterface<Approximate>& agia=af;
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
    ARIADNE_TEST_EQUAL(evaluate(x.derivative(0),Vector<ValidatedFloat64>{2.4_x,1.3_x}),1.0_x);
    ARIADNE_TEST_EQUAL(evaluate(x.derivative(1),Vector<ApproximateFloat64>{2.4, 1.3}),0.0);

}
*/

