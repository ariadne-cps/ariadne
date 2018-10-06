/***************************************************************************
 *            check_algebra.hpp
 *
 *  Copyright 2009-17  Pieter Collins
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

#include "../test.hpp"
#include "../utility.hpp"
#include "utility/metaprogramming.hpp"
#include "numeric/operators.hpp"

namespace Ariadne {
template<class X> String vector_class_name() { return String("Vector<") + class_name<X>() + ">"; }
}

using namespace Ariadne;


template<class V> using DefinedScalarType = typename V::ScalarType;
template<class V> struct HasScalarType : Has<DefinedScalarType,V> { };

template<class A> using DefinedNumericType = typename A::NumericType;
template<class V> struct HasNumericType : Has<DefinedNumericType,V> { };

template<class T> using ZeroElementReturnType = decltype(declval<T>().zero_element());
template<class V, class R=Return<DontCare>> struct HasZeroElementMethod : Is<R,ZeroElementReturnType,V> { };

template<class T> using SizeReturnType = decltype(declval<T>().size());
template<class V, class R=Return<DontCare>> struct HasSizeMethod : Is<R,SizeReturnType,V> { };

template<class V, class I> using SubscriptingReturnType = decltype(declval<const V>().operator[](declval<I>()));
template<class V, class I, class R=Return<DontCare>> struct HasSubscriptingMethod : Is<R,SubscriptingReturnType,V,I> { };

template<class A1, class A2> using JoinReturnType = decltype(join(declval<A1>(),declval<A2>()));
template<class A1, class A2, class R=Return<DontCare>> struct HasJoin : Is<R,JoinReturnType,A1,A2> { };

template<class V> using NormFunctionReturnType = decltype(norm(declval<V>()));
template<class V, class R=Return<DontCare>> struct HasNorm : Is<R,NormFunctionReturnType,V> { };

template<class V> using NormMethodReturnType = decltype(declval<V>().norm());
template<class V, class R=Return<DontCare>> struct HasNormMethod : Is<R,NormMethodReturnType,V> { };

template<class A, class K> using DerivativeReturnType = decltype(derivative(declval<A>(),declval<K>()));
template<class A, class K=SizeType, class R=Return<DontCare>> struct HasDerivative : Is<R,DerivativeReturnType,A,K> { };

template<class A, class K> using AntiderivativeReturnType = decltype(antiderivative(declval<A>(),declval<K>()));
template<class A, class K=SizeType, class R=Return<DontCare>> struct HasAntiderivative : Is<R,AntiderivativeReturnType,A,K> { };



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

namespace Ariadne {
template<class X> String algebra_class_name() { return String("Algebra<") + class_name<X>() + ">"; }
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
    ARIADNE_TEST_STATIC_ASSERT(HasNormMethod<A>);
//    ARIADNE_TEST_STATIC_ASSERT(IsSomething<decltype(declval<A>().unit_coefficient())>);
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
