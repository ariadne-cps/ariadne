/***************************************************************************
 *            test/utility.h
 *
 *  Copyright 2007-13  Pieter Collins
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

#ifndef ARIADNE_TEST_UTILITY_H
#define ARIADNE_TEST_UTILITY_H

namespace Ariadne {
//namespace Test {

// Compare with double's as if they were ExactFloat64

template<class X> auto operator==(X const& x, double d) -> decltype(x==ExactFloat64(d)) { return x==ExactFloat64(d); }
template<class X> auto operator!=(X const& x, double d) -> decltype(x!=ExactFloat64(d)) { return x!=ExactFloat64(d); }
template<class X> auto operator< (X const& x, double d) -> decltype(x< ExactFloat64(d)) { return x< ExactFloat64(d); }
template<class X> auto operator> (X const& x, double d) -> decltype(x> ExactFloat64(d)) { return x> ExactFloat64(d); }
template<class X> auto operator<=(X const& x, double d) -> decltype(x<=ExactFloat64(d)) { return x<=ExactFloat64(d); }
template<class X> auto operator>=(X const& x, double d) -> decltype(x>=ExactFloat64(d)) { return x>=ExactFloat64(d); }


#define ARIADNE_HAS_TYPEDEF(typename_to_check) \
    template<class A, class = Fallback> struct Has##typename_to_check : False { }; \
    template<class A> struct Has##typename_to_check<A, EnableIf<IsDefined<typename A::typename_to_check>,Fallback>> : True { }; \

#define ARIADNE_HAS_METHOD(method_to_check) \
    template<class A, class = Fallback> struct Has_##method_to_check : False { }; \
    template<class A> struct Has_##method_to_check<A, EnableIf<IsDefined<decltype(declval<A>().method_to_check())>,Fallback>> : True { }; \

struct OperatorPositive { template<class A> auto operator()(A a) const -> decltype(+a) { return +a; } };
struct OperatorNegative { template<class A> auto operator()(A a) const -> decltype(-a) { return -a; } };
struct OperatorNegate { template<class A> auto operator()(A a) const -> decltype(-a) { return -a; } };
struct OperatorUnaryPlus { template<class A> auto operator()(A a) const -> decltype(+a) { return +a; } };
struct OperatorUnaryMinus { template<class A> auto operator()(A a) const -> decltype(-a) { return -a; } };
struct OperatorPlus { template<class A1, class A2> auto operator()(A1 a1, A2 a2) const -> decltype(a1+a2) { return a1+a2; } };
struct OperatorMinus { template<class A1, class A2> auto operator()(A1 a1, A2 a2) const -> decltype(a1-a2) { return a1-a2; } };
struct OperatorTimes { template<class A1, class A2> auto operator()(A1 a1, A2 a2) const -> decltype(a1*a2) { return a1*a2; } };
struct OperatorDivides { template<class A1, class A2> auto operator()(A1 a1, A2 a2) const -> decltype(a1/a2) { return a1/a2; } };

struct OperatorAnd { template<class A1, class A2> auto operator()(A1 a1, A2 a2) const -> decltype(a1 && a2) { return a1 && a2; } };
struct OperatorOr { template<class A1, class A2> auto operator()(A1 a1, A2 a2) const -> decltype(a1 || a2) { return a1 || a2; } };
struct OperatorNot { template<class A> auto operator()(A a) const -> decltype(!a) { return !a; } };

struct OperatorEqual { template<class A1, class A2> auto operator()(A1 a1, A2 a2) const -> decltype(a1==a2) { return a1==a2; } };
struct OperatorLess { template<class A1, class A2> auto operator()(A1 a1, A2 a2) const -> decltype(a1< a2) { return a1< a2; } };

template<class T1, class T2, class = Fallback> struct IsEqualityComparible : False { };
template<class T1, class T2> struct IsEqualityComparible<T1,T2,EnableIf<IsDefined<decltype(declval<T1&>()==declval<T2&>())>,Fallback>> : True { };

template<class T1, class T2, class = Fallback> struct IsLessThanCompartible : False { };
template<class T1, class T2> struct IsLessThanCompartible<T1,T2,EnableIf<IsDefined<decltype(declval<T1>()<declval<T2>())>,Fallback>> : True { };

template<class R, class Op, class A, class = Fallback> struct HasUnaryOperatorReturning : False { };
template<class R, class Op, class A> struct HasUnaryOperatorReturning<R,Op,A,EnableIf<IsConvertible<ResultOf<Op(A)>,R>,Fallback>> : True { };
template<class R, class Op, class A1, class A2, class = Fallback>  struct HasBinaryOperatorReturning : False { };
template<class R, class Op, class A1, class A2> struct HasBinaryOperatorReturning<R,Op,A1,A2,EnableIf<IsConvertible<ResultOf<Op(A1,A2)>,R>,Fallback>> : True { };

template<class R, class Op, class A1, class A2=Void, class = Fallback> struct HasOperatorReturning : False { };
template<class R, class Op, class A> struct HasOperatorReturning<R,Op,A,Void,EnableIf<IsConvertible<ResultOf<Op(A)>,R>,Fallback>> : True { };
template<class R, class Op, class A1, class A2> struct HasOperatorReturning<R,Op,A1,A2,EnableIf<IsConvertible<ResultOf<Op(A1,A2)>,R>,Fallback>> : True { };


//template<Op, class A, class = Fallback> struct HasUnaryOperator : False { };
//template<class Op, class A> struct HasUnaryOperator<Op,A,EnableIf<IsConvertible<ResultOf<Op(A)>,DontCare>,Fallback>> : True { };
//template<class Op, class A1, class A2, class = Fallback>  struct HasBinaryOperator : False { };
//template<class Op, class A1, class A2> struct HasBinaryOperator<R,Op,A1,A2,EnableIf<IsConvertible<ResultOf<Op(A1,A2)>,R>,Fallback>> : True { };

template<class Op, class A> struct HasUnaryOperator : HasUnaryOperatorReturning<DontCare,Op,A> { };
template<class Op, class A1, class A2> struct HasBinaryOperator : HasBinaryOperatorReturning<DontCare,Op,A1,A2> { };
template<class Op, class A1, class A2=Void, class R=Void> struct HasOperator;
template<class Op, class A1, class A2> struct HasOperator<Op,A1,A2,Void> : HasBinaryOperator<Op,A1,A2> { };
template<class Op, class A> struct HasOperator<Op,A,Void> : HasUnaryOperator<Op,A> { };
template<class Op, class A1, class A2, class R> struct HasOperator<Op,A1,A2,Return<R>> : HasBinaryOperatorReturning<R,Op,A1,A2> { };
template<class Op, class A, class R> struct HasOperator<Op,A,Return<R>> : HasUnaryOperatorReturning<R,Op,A> { };

//} // namespace Test
//using namespace Test;
} // namespace Ariadne

#endif // ARIADNE_TEST_UTILITY_H

