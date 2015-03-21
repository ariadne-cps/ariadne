/***************************************************************************
 *            check_algebra.h
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
#include "numeric/operators.h"

namespace Ariadne {
template<class X> String vector_class_name() { return String("Vector<") + class_name<X>() + ">"; }
}

using namespace Ariadne;


template<class V, class R=Return<DontCare>, class = Fallback> struct HasZeroElementMethod : False { };
template<class V, class R> struct HasZeroElementMethod<V, Return<R>, EnableIf<IsConvertible<decltype(declval<V>().zero_element()),R>,Fallback>> : True { };

template<class V, class R=Return<DontCare>, class = Fallback> struct HasSizeMethod : False { };
template<class V, class R> struct HasSizeMethod<V, Return<R>, EnableIf<IsConvertible<decltype(declval<V>().size()),R>,Fallback>> : True { };

template<class V, class I, class R=DontCare, class = Fallback> struct HasSubscriptingMethod : False { };
template<class V, class I, class R> struct HasSubscriptingMethod<V, I, Return<R>, EnableIf<IsConvertible<decltype(declval<V>()[declval<I>()]),R>,Fallback>> : True { };

template<class A1, class A2, class R=Return<DontCare>, class = Fallback> struct HasJoin : False { };
template<class A1, class A2, class R> struct HasJoin<A1, A2, Return<R>, EnableIf<IsConvertible<decltype(join(declval<A1>(),declval<A2>())),R>,Fallback>> : True { };

template<class V, class R=DontCare, class = Fallback> struct HasNorm : False { };
template<class V, class R> struct HasNorm<V, Return<R>, EnableIf<IsConvertible<decltype(norm(declval<V>())),R>,Fallback>> : True { };


template<class A, class R=DontCare, class = Fallback> struct HasNormMethod : False { };
template<class A, class R> struct HasNormMethod<A, Return<R>, EnableIf<IsConvertible<decltype(declval<A>().norm()),R>,Fallback>> : True { };

template<class A, class K=SizeType, class R=DontCare, class = Fallback> struct HasDerivative : False { };
template<class A, class K, class R> struct HasDerivative<A, K, Return<R>, EnableIf<IsConvertible<decltype(derivative(declval<A>(),declval<K>())),R>,Fallback>> : True { };

template<class A, class K=SizeType, class R=DontCare, class = Fallback> struct HasAntiderivative : False { };
template<class A, class K, class R> struct HasAntiderivative<A, K, Return<R>, EnableIf<IsConvertible<decltype(derivative(declval<A>(),declval<K>())),R>,Fallback>> : True { };

template<class V, class = Fallback> struct HasScalarType : False { }; \
template<class V> struct HasScalarType<V, EnableIf<IsDefined<typename V::ScalarType>,Fallback>> : True { }; \

template<class V, class = Fallback> struct HasNumericType : False { }; \
template<class V> struct HasNumericType<V, EnableIf<IsDefined<typename V::HasNumericType>,Fallback>> : True { }; \


template<class V> class CheckVectorConcept {
  public:
    CheckVectorConcept() { }
  public:
    void check();
    bool has_scalar_type();
    void check_vector_concept();
    void check_norm_concept();
};


template<class V> void CheckVectorConcept<V>::check() {
    ARIADNE_TEST_PRINT(vector_class_name<typename V::ScalarType>());
    if (has_scalar_type()) {
        check_vector_concept();
        check_norm_concept();
    }
}

template<class V> bool CheckVectorConcept<V>::has_scalar_type() {
    return HasScalarType<V>::value;
}

template<class V> void CheckVectorConcept<V>::check_vector_concept() {
    typedef typename V::ScalarType S;

    ARIADNE_TEST_STATIC_ASSERT( HasSizeMethod<const V,SizeType> );
    ARIADNE_TEST_STATIC_ASSERT( HasZeroElementMethod<const V,S> );
    ARIADNE_TEST_STATIC_ASSERT( HasSubscriptingMethod<const V,SizeType,S> );
    ARIADNE_TEST_STATIC_ASSERT( HasSubscriptingMethod<V,SizeType,S> );
    ARIADNE_TEST_STATIC_ASSERT( HasOperator<OperatorPlus,V> );
    ARIADNE_TEST_STATIC_ASSERT( HasOperator<OperatorPlus,V,Return<V>> );
    ARIADNE_TEST_STATIC_ASSERT( HasOperator<OperatorMinus,V,Return<V>> );
    ARIADNE_TEST_STATIC_ASSERT( HasOperator<OperatorPlus,V,V,Return<V>> );
    ARIADNE_TEST_STATIC_ASSERT( HasOperator<OperatorMinus,V,V,Return<V>> );
    ARIADNE_TEST_STATIC_ASSERT( HasOperator<OperatorTimes,S,V,Return<V>> );
    ARIADNE_TEST_STATIC_ASSERT( HasOperator<OperatorTimes,V,S,Return<V>> );
    ARIADNE_TEST_STATIC_ASSERT( HasOperator<OperatorDivides,V,S,Return<V>> );

    ARIADNE_TEST_STATIC_ASSERT( HasJoin<V,V,Return<V>> );
    ARIADNE_TEST_STATIC_ASSERT( HasJoin<V,S,Return<V>> );
    ARIADNE_TEST_STATIC_ASSERT( HasJoin<S,V,Return<V>> );
    ARIADNE_TEST_STATIC_ASSERT( HasJoin<S,S,Return<V>> );
}

template<class V> void CheckVectorConcept<V>::check_norm_concept() {
    typedef decltype(mag(declval<V>()[0])) MagType;
    ARIADNE_TEST_STATIC_ASSERT( HasNorm<V,MagType> );
}




template<class A> class CheckAlgebraConcept
{
  public:
    CheckAlgebraConcept() { }
  public:
    bool has_numeric_type();
    void check_algebra_concept();
    void check_unital_algebra_concept();
    void check_division_algebra_concept();
    void check_banach_algebra_concept();
    void check_graded_algebra_concept();
    void check_transcendental_operators_concept();
    void check_differential_algebra_concept();
};

template<class A> bool CheckAlgebraConcept<A>::has_numeric_type()
{
    ARIADNE_TEST_STATIC_ASSERT(HasNumericType<A>);
    return HasNumericType<A>::value;
}
namespace Ariadne { template<> String class_name<EffectiveScalarFunction>() { return "EffectiveScalarFunction"; } }
template<class A> void CheckAlgebraConcept<A>::check_algebra_concept()
{
    typedef typename A::NumericType X;

    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(declval<A>().create_zero()),A>);

    ARIADNE_TEST_STATIC_ASSERT(IsConvertible<A,A>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorPlus,A>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorPlus,A,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorMinus,A,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorPlus,A,A>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorPlus,A,A,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorMinus,A,A,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorTimes,A,A,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorTimes,A,X,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorTimes,X,A,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorDivides,A,X,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<Sqr,A,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<Pow,A,Nat,Return<A>>);
}
class A { }; class B { public: B& operator=(A a) { return *this; } };
template<class A> void CheckAlgebraConcept<A>::check_unital_algebra_concept()
{
    this->check_algebra_concept();

    typedef typename A::NumericType X;

    ARIADNE_TEST_STATIC_ASSERT(IsAssignable<A,X>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(declval<A>().create_constant(declval<X>())),A>);

    ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorPlus,A,X,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorPlus,X,A,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorMinus,A,X,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorMinus,X,A,Return<A>>);
}

template<class A> void CheckAlgebraConcept<A>::check_division_algebra_concept()
{
    this->check_unital_algebra_concept();

    typedef typename A::NumericType X;

    ARIADNE_TEST_STATIC_ASSERT(HasOperator<Rec,A,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorDivides,A,A,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorDivides,X,A,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<Pow,A,Int,Return<A>>);
}

template<class A> void CheckAlgebraConcept<A>::check_transcendental_operators_concept()
{
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<Sqrt,A,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<Exp,A,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<Log,A,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<Sin,A,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<Cos,A,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<Tan,A,Return<A>>);
    ARIADNE_TEST_STATIC_ASSERT(HasOperator<Atan,A,Return<A>>);
}

template<class A> void CheckAlgebraConcept<A>::check_banach_algebra_concept()
{
    ARIADNE_TEST_STATIC_ASSERT(IsSomething<decltype(declval<A>().norm())>);
    ARIADNE_TEST_STATIC_ASSERT(IsSomething<decltype(declval<A>().unit_coefficient())>);
}

template<class A> void CheckAlgebraConcept<A>::check_graded_algebra_concept()
{
    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(declval<A>().degree()),DegreeType>);
}

template<class A> void CheckAlgebraConcept<A>::check_differential_algebra_concept()
{
//    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(declval<A>().smoothness()),DegreeType>);
    ARIADNE_TEST_STATIC_ASSERT(HasDerivative<A,SizeType,Return<A>>);
}
