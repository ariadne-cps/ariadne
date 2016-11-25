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

// Compare with double's as if they were Float64Value

template<class X> decltype(auto) operator==(X const& x, double d) { return x==Float64Value(d); }
template<class X> decltype(auto) operator!=(X const& x, double d) { return x!=Float64Value(d); }
template<class X> decltype(auto) operator< (X const& x, double d) { return x< Float64Value(d); }
template<class X> decltype(auto) operator> (X const& x, double d) { return x> Float64Value(d); }
template<class X> decltype(auto) operator<=(X const& x, double d) { return x<=Float64Value(d); }
template<class X> decltype(auto) operator>=(X const& x, double d) { return x>=Float64Value(d); }

struct OperatorPlus {
    template<class A> decltype(auto) operator()(A a) const { return +a; }
    template<class A1, class A2> decltype(auto) operator()(A1 a1, A2 a2) const { return a1+a2; }
};
struct OperatorMinus {
    template<class A> decltype(auto) operator()(A a) const { return -a; }
    template<class A1, class A2> decltype(auto) operator()(A1 a1, A2 a2) const { return a1-a2; }
};
struct OperatorTimes {
    template<class A1, class A2> decltype(auto) operator()(A1 a1, A2 a2) const { return a1*a2; }
};
struct OperatorDivides {
    template<class A1, class A2> decltype(auto) operator()(A1 a1, A2 a2) const { return a1/a2; }
};

struct OperatorAnd { template<class A1, class A2> decltype(auto) operator()(A1 a1, A2 a2) const { return a1 && a2; } };
struct OperatorOr { template<class A1, class A2> decltype(auto) operator()(A1 a1, A2 a2) const { return a1 || a2; } };
struct OperatorNot { template<class A> decltype(auto) operator()(A a) const { return !a; } };

struct OperatorEqual { template<class A1, class A2> decltype(auto) operator()(A1 a1, A2 a2) const { return a1==a2; } };
struct OperatorEquals { template<class A1, class A2> decltype(auto) operator()(A1 a1, A2 a2) const { return a1==a2; } };
struct OperatorLess { template<class A1, class A2> decltype(auto) operator()(A1 a1, A2 a2) const { return a1< a2; } };

namespace Detail {

template<class OP, class... AS> struct SafeTypedef {
    template<class O, class = ResultOf<O(AS...)>> static ResultOf<O(AS...)> safe(int);
    template<class O> static Fallback safe(...);
    typedef decltype(safe<OP>(1)) Type;
};

template<class R, template<class>class T, class A1, class = decltype(declval<R>()=declval<T<A1>>())> True is(int);
template<class R, template<class>class T, class A1> False is(...);
template<class R, template<class,class>class T, class A1, class A2, class = decltype(declval<R>()=declval<T<A1,A2>>())> True is(int,int);
template<class R, template<class,class>class T, class A1, class A2> False is(...);
template<class R, template<class,class,class>class T, class A1, class A2, class A3, class = decltype(declval<R>()=declval<T<A1,A2,A3>>())> True is(int,int,int);
template<class R, template<class,class,class>class T, class A1, class A2, class A3> False is(...);


} //namespace Detail

template<class R, template<class...>class T, class... AS> struct Is;
template<class R, template<class>class T, class A> struct Is<Return<R>,T,A> : decltype(Detail::is<R,T,A>(1)) { };
template<class R, template<class,class>class T, class A1, class A2> struct Is<Return<R>,T,A1,A2> : decltype(Detail::is<R,T,A1,A2>(1,2)) { };
template<class R, template<class,class,class>class T, class A1, class A2, class A3> struct Is<Return<R>,T,A1,A2,A3> : decltype(Detail::is<R,T,A1,A2,A3>(1,2,3)) { };


template<class OP, class... AS> using ReturnType = decltype(declval<OP>().operator()(declval<AS...>()));
template<class OP, class A1> using UnaryReturnType = decltype(declval<OP>().operator()(declval<A1>()));
template<class OP, class A1, class A2> using BinaryReturnType = decltype(declval<OP>().operator()(declval<A1>(),declval<A2>()));

template<class OP, class... AS> struct HasOperator;
template<class OP, class A1> struct HasOperator<OP,A1> : Is<Return<DontCare>, UnaryReturnType, OP, A1> { };
template<class OP, class A1, class A2> struct HasOperator<OP,A1,A2> : Is<Return<DontCare>, BinaryReturnType, OP, A1, A2> { };
template<class OP, class A1, class R> struct HasOperator<OP,A1,Return<R>> : Is<Return<R>, UnaryReturnType, OP, A1> { };
template<class OP, class A1, class A2, class R> struct HasOperator<OP,A1,A2,Return<R>> : Is<Return<R>, BinaryReturnType, OP, A1, A2> { };
template<class R, class OP, class... AS> struct HasOperatorReturning : Is<Return<R>, ReturnType, OP, AS...> { };

template<class OP, class... AS> using SafeType = typename Detail::SafeTypedef<OP,AS...>::Type;

template<class T> struct Expect { template<class U, EnableIf<IsSame<T,U>> =dummy> Expect<T>& operator=(U const&) { return *this; } };

//} // namespace Test
//using namespace Test;
} // namespace Ariadne

#endif // ARIADNE_TEST_UTILITY_H

