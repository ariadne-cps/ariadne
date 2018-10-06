/***************************************************************************
 *            metaprogramming.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file metaprogramming.hpp
 *  \brief Classes for template metaprogramming.
 */

#ifndef ARIADNE_METAPROGRAMMING_HPP
#define ARIADNE_METAPROGRAMMING_HPP

#include <type_traits>

namespace Ariadne {

typedef void Void;
typedef bool Bool;


using std::declval;

using True = std::true_type;
using False = std::false_type;

template<bool b, class... PS> struct ValueAnd;
template<class... PS> struct And;
template<class P1, class... PS> struct And<P1,PS...> : ValueAnd<P1::value,PS...> { };
template<> struct And<> : True { };
template<class... PS> struct ValueAnd<false, PS...> : False { };
template<class... PS> struct ValueAnd<true, PS...> : And<PS...> { };

template<class P1, class P2, class P3=void> struct Or { static const bool value = P1::value || P2::value || P3::value; };
template<class P1, class P2> struct Or<P1,P2> { static const bool value = P1::value || P2::value; };
template<class P> struct Not { static const bool value = !P::value; };

template<bool B, class T, class F> struct IfThenElse_Bool;
template<class T, class F> struct IfThenElse_Bool<true,T,F> { typedef T Type; };
template<class T, class F> struct IfThenElse_Bool<false,T,F> { typedef F Type; };
template<class P, class T, class F> using IfThenElse = typename IfThenElse_Bool<P::value,T,F>::Type;

typedef bool Dummy; static const bool dummy=false;
//typedef int Dummy; static const int dummy=1;

template<class P, typename V=Dummy> using EnableIf = typename std::enable_if<P::value,V>::type;
template<class P, typename V=Dummy> using DisableIf = typename std::enable_if<not P::value,V>::type;

template<class T, T V> using IntegralConstant = std::integral_constant<T,V>;
template<class T1, class T2> using IsSame = std::is_same<T1,T2>;
template<class T, class U> using IsConvertible = std::is_convertible<T,U>;
template<class T, class... U> using IsConstructible = std::is_constructible<T,U...>;
template<class T, class U> using IsStrictlyConstructible = And<IsConstructible<T,U>,Not<IsConvertible<U,T>>>;
template<class T> using IsDefaultConstructible = std::is_default_constructible<T>;
template<class T, class U> using IsAssignable = std::is_assignable<T,U>;
template<class T, class U> using IsBaseOf = std::is_base_of<T,U>;
//template<class F, class... AS> using IsInvocable = std::is_invokable<F,AS...>;
//template<class R, class F, class... AS> using IsInvocableReturning = std::is_invokable_r<R,F,AS...>;
template<class T> using IsBuiltinArithmetic = std::is_arithmetic<T>;
template<class T> using IsBuiltinIntegral = std::is_integral<T>;
template<class T> using IsBuiltinFloatingPoint = std::is_floating_point<T>;
template<class T> using IsBuiltinSigned = std::is_signed<T>;
template<class T> using IsBuiltinUnsigned = std::is_unsigned<T>;
template<class T> using IsBuiltinSignedIntegral = std::integral_constant<bool,std::is_integral<T>::value and std::is_signed<T>::value>;
template<class T> using IsBuiltinUnsignedIntegral = std::integral_constant<bool,std::is_integral<T>::value and std::is_unsigned<T>::value>;

template<class SIG> using ResultOf = typename std::result_of<SIG>::type;

template<class T1, class T2, class T3> using AreSame = And<IsSame<T1,T2>,IsSame<T2,T3>>;

struct Fallback { };
struct DontCare { template<class T> DontCare(T); };

template<class R> struct Return { typedef R ReturnType; };
template<class RET> using ReturnType = typename RET::ReturnType;

struct Any { };

template<template<class...>class F, class... T> False _has(...);
template<template<class>class F, class T, class = F<T>> True _has(int);
template<template<class,class>class F, class T1, class T2, class = F<T1,T2>> True _has(int,int);
template<template<class,class,class>class F, class T1, class T2, class T3, class = F<T1,T2,T3>> True _has(int,int,int);
template<template<class...>class F, class... T> struct Has : decltype(_has<F,T...>()) { };
template<template<class>class F, class T> struct Has<F,T> : decltype(_has<F,T>(1)) { };
template<template<class,class>class F, class T1, class T2> struct Has<F,T1,T2> : decltype(_has<F,T1,T2>(1,2)) { };
template<template<class,class,class>class F, class T1, class T2, class T3> struct Has<F,T1,T2,T3> : decltype(_has<F,T1,T2,T3>(1,2,3)) { };

template<class T> struct Self { typedef T Type; };
template<class T> using SelfType = typename Self<T>::Type;

template<class F, class A, class R=DontCare, class=Fallback> struct IsCallable : False { };
template<class F, class A, class R> struct IsCallable<F,A,R, EnableIf<IsConvertible<decltype(declval<F>()(declval<A>())),R>,Fallback>> : True { };

template<class T1, class T2> struct First { typedef T1 Type; };
//template<class T1, class... Tps> struct First { typedef T1 Type; };
template<class T1, class T2, class... Tps> struct Second { typedef T2 Type; };

template<class X> using RemoveConst = typename std::remove_const<X>::type;
template<class X> using RemoveReference = typename std::remove_reference<X>::type;

template<class X> using NegationType = decltype(-declval<X>());
template<class X1, class X2=X1> using SumType = decltype(declval<X1>()+declval<X2>());
template<class X1, class X2=X1> using DifferenceType = decltype(declval<X1>()-declval<X2>());
template<class X1, class X2=X1> using ProductType = decltype(declval<X1>()*declval<X2>());
template<class X1, class X2=X1> using QuotientType = decltype(declval<X1>()/declval<X2>());
template<class X1, class X2=X1> using ArithmeticType = SumType<ProductType<X1,X2>>;
template<class X1, class X2=X1> using EqualsType = decltype(declval<X1>()==declval<X2>());
template<class X1, class X2=X1> using LessType = decltype(declval<X1>()< declval<X2>());
template<class X1, class X2=X1> using ComparisonType = LessType<X1,X2>;

template<class X1, class X2=X1> using EqualityType = decltype(declval<X1>()==declval<X2>());
template<class X1, class X2=X1> using InequalityType = decltype(declval<X1>()!=declval<X2>());

template<class X1, class X2=X1> using InplaceSumType = RemoveReference<decltype(declval<X1&>()+=declval<X2>())>;
template<class X1, class X2=X1> using InplaceDifferenceType = RemoveReference<decltype(declval<X1&>()-=declval<X2>())>;
template<class X1, class X2=X1> using InplaceProductType = RemoveReference<decltype(declval<X1&>()*=declval<X2>())>;
template<class X1, class X2=X1> using InplaceQuotientType = RemoveReference<decltype(declval<X1&>()/=declval<X2>())>;

template<class A1, class A2> struct HasEquality {
    template<class AA1, class AA2, class=decltype(std::declval<AA1>()==std::declval<AA2>())> static std::true_type test(int);
    template<class AA1, class AA2> static std::false_type test(...);
    static const bool value = decltype(test<A1,A2>(1))::value;
};

template<class A1, class A2> struct CanAdd {
    template<class AA1, class AA2, class=decltype(std::declval<AA1>()+std::declval<AA2>())> static std::true_type test(int);
    template<class AA1, class AA2> static std::false_type test(...);
    static const bool value = decltype(test<A1,A2>(1))::value;
};

template<class A1, class A2> struct CanSubtract {
    template<class AA1, class AA2, class=decltype(std::declval<AA1>()-std::declval<AA2>())> static std::true_type test(int);
    template<class AA1, class AA2> static std::false_type test(...);
    static const bool value = decltype(test<A1,A2>(1))::value;
};

template<class A1, class A2> struct CanMultiply {
    template<class AA1, class AA2, class=decltype(std::declval<AA1>()*std::declval<AA2>())> static std::true_type test(int);
    template<class AA1, class AA2> static std::false_type test(...);
    static const bool value = decltype(test<A1,A2>(1))::value;
};

template<class A1, class A2> struct CanDivide {
    template<class AA1, class AA2, class=decltype(std::declval<AA1>()/std::declval<AA2>())> static std::true_type test(int);
    template<class AA1, class AA2> static std::false_type test(...);
    static const bool value = decltype(test<A1,A2>(1))::value;
};

//template<class R, class F, class... AS> using IsInvocableReturning = std::is_invokable_r<R,F,AS...>;
template<class F, class... AS> struct IsInvocable;
template<class R, class F, class... AS> struct IsInvocableReturning;

template<class F, class A> struct IsInvocable<F,A> {
    template<class FF, class AA, class=decltype(std::declval<FF>()(std::declval<AA>()))>
        static std::true_type test(int);
    template<class FF, class AA>
        static std::false_type test(...);
    static const bool value = decltype(test<F,A>(1))::value;
};

template<class R, class F, class A> struct IsInvocableReturning<R,F,A> {
    template<class RR, class FF, class AA, class=decltype(std::declval<RR>()=std::declval<FF>()(std::declval<AA>()))>
        static std::true_type test(int);
    template<class RR, class FF, class AA>
        static std::false_type test(...);
    static const bool value = decltype(test<R,F,A>(1))::value;
};

template<class X> struct HasGenericType {
    template<class XX, class=typename XX::GenericType> static std::true_type test(int);
    template<class XX> static std::false_type test(...);
    static const bool value = decltype(test<X>(1))::value;
};

template<class X> struct HasPrecisionType {
    template<class XX, class=typename XX::PrecisionType> static std::true_type test(int);
    template<class XX> static std::false_type test(...);
    static const bool value = decltype(test<X>(1))::value;
};


// The following class only accepts an exact match as an argument
template<class T> class SuppressConversions {
    T const& _t;
  public:
    template<class TT, EnableIf<IsSame<TT,T>> =dummy> SuppressConversions(TT const& t) : _t(t) { }
    explicit operator T const& () const { return _t; }
};
} // namespace Ariadne

#endif
