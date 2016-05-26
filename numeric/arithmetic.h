/***************************************************************************
 *            arithmetic.h
 *
 *  Copyright 2008-16  Pieter Collins
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

/*! \file arithmetic.h
 *  \brief Arithmetic declarations
 */

#ifndef ARIADNE_ARITHMETIC_H
#define ARIADNE_ARITHMETIC_H

#include "utility/metaprogramming.h"
#include "logical.decl.h"
#include "sign.h"

namespace Ariadne {

using Nat = uint;
using Int = int;

class Natural;
class Integer;

template<class T, class NT=T> struct Directed {
    friend NT neg(T const&);
};

template<class T, class AT=T> struct Additive {
    friend AT add(T const&, T const&);
};

template<class T, class MT=T> struct Multiplicative {
    friend MT sqr(T const&);
    friend MT mul(T const&, T const&);
    friend T pow(T const&, Nat);
};

template<class T, class NT=T> struct Abelian {
    friend T nul(T const&);
    friend T pos(T const&);
    friend NT neg(T const&);
    friend T add(T const&, T const&);
    friend T sub(T const&, T const&);
};

template<class T> struct SemiRing {
    friend T add(T const&, T const&);
    friend T mul(T const&, T const&);
    friend T sqr(T const&);
};

template<class T, class NT=T> struct Ring : Abelian<T,NT> {
    friend T sqr(T const&);
    friend T mul(T const&, T const&);
    friend T pow(T const&, Nat);
};

template<class T, class NT=T> struct DiadicRing : Ring<T,NT> {
    friend T hlf(T const&);
};

template<class T, class NT=T, class QT=NT> struct Field : DiadicRing<T,NT> {
    friend QT rec(T const&);
    friend QT div(T const&, T const&);
    friend QT pow(T const&, Int);
};

template<class T, class PT=T> struct Monotone {
    friend PT sqrt(PT const&);
    friend PT exp(T const&);
    friend T log(PT const&);
    friend T asin(T const&);
    friend T atan(T const&);
};

template<class T, class PT=T> struct MonotoneSemiRing : SemiRing<T>, SemiRing<PT>, Monotone<T,PT> { };
template<class T> struct MonotoneSemiRing<T> : SemiRing<T>, Monotone<T,T> { };

template<class T, class NT, class PT> struct DirectedRing : Abelian<T,NT>, SemiRing<PT>, Monotone<T,PT> { };

template<class T, class PT=T> struct Transcendental : Monotone<T,PT> {
    friend T sin(T const&);
    friend T cos(T const&);
    friend T tan(T const&);
};

template<class T, class PT=T> struct TranscendentalField : Field<T>, SemiRing<PT>, Transcendental<T,PT> { };
template<class T> struct TranscendentalField<T> : Field<T>, Transcendental<T> { };

template<class T> struct Lattice {
    friend T max(T const&, T const&);
    friend T min(T const&, T const&);
};


template<class T, class PT=T> struct DirectedLattice : Lattice<T> {
    friend PT abs(T const&);
};


template<class T, class A> struct Apartness {
    friend A neq(T const&, T const&);
};

template<class T, class O, class A=O> struct Ordered : Apartness<T,A> {
    friend O leq(T const&, T const&);
};

template<class T, class PT=T> struct DeclareReal {
    T neg(T const&);
    T add(T const&, T const&);
    T sub(T const&, T const&);
    T mul(T const&, T const&);
    T pow(T const&, Integer const&);
    T rec(T const&);
    PT sqrt(PT const&);
    PT exp(T const&);
    T log(PT const&);
    T sin(T const&);
    T cos(T const&);
    T tan(T const&);
    T atan(T const&);

    T max(T const&, T const&);
    //PT max(T const&, PT const&);
    //PT max(PT const&, T const&);
    T min(T const&, T const&);
    PT abs(T const&);
    PT dist(T const&, T const&);
};

template<class T, class NT, class PT> struct DeclareDirectedReal {
    NT neg(T const&);
    T add(T const&, T const&);
    T sub(T const&, NT const&);
    PT sqrt(PT const&);
    PT exp(T const&);
    T log(PT const&);
    T atan(T const&);

    T max(T const&, T const&);
    PT max(T const&, PT const&);
    PT max(PT const&, T const&);
    T min(T const&, T const&);
};

template<class PT> struct DeclarePositiveReal {
    PT add(PT const&, PT const&);
    PT mul(PT const&, PT const&);
    PT div(PT const&, PT const&);
    PT rec(PT const&);
    PT sqrt(PT const&);
    PT atan(PT const&);

    PT max(PT const&, PT const&);
    PT min(PT const&, PT const&);
};

template<class T, class QT> struct DeclarePositiveDirectedReal {
    T add(T const&, T const&);
    T mul(T const&, T const&);
    T pow(T const&, Natural const&);
    T div(T const&, QT const&);
    QT rec(T const&);
    T sqrt(T const&);
    T atan(T const&);

    T max(T const&, T const&);
    T min(T const&, T const&);
};

template<class T, class NT=T, class QT=NT> struct ArithmeticOperators {
    friend T operator+(T const& t);
    friend NT operator-(T const& t);
    friend T operator+(T const& t1, T const& t2);
    friend T operator-(T const& t1, NT const& t2);
    friend T operator*(T const& t1, T const& t2);
    friend T operator/(T const& t1, NT const& t2);
    friend T& operator+=(T& t1, T const& t2);
    friend T& operator-=(T& t1, NT const& t2);
    friend T& operator*=(T& t1, T const& t2);
    friend T& operator/=(T& t1, NT const& t2);
};

template<class T, class NT=T> struct DefineRingOperators {
    friend T operator+(T const& t) { return pos(t); }
    friend NT operator-(T const& t) { return neg(t); }
    friend T operator+(T const& t1, T const& t2) { return add(t1,t2); }
    friend T operator-(T const& t1, NT const& t2) { return sub(t1,t2); }
    friend T operator*(T const& t1, T const& t2) { return mul(t1,t2); }
    friend T& operator+=(T& t1, T const& t2) { return t1=add(t1,t2); }
    friend T& operator-=(T& t1, NT const& t2) { return t1=sub(t1,t2); }
    friend T& operator*=(T& t1, T const& t2) { return t1=mul(t1,t2); }
};

template<class T, class NT=T, class QT=NT> struct DefineArithmeticOperators {
    friend T operator+(T const& t) { return pos(t); }
    friend NT operator-(T const& t) { return neg(t); }
    friend T operator+(T const& t1, T const& t2) { return add(t1,t2); }
    friend T operator-(T const& t1, NT const& t2) { return sub(t1,t2); }
    friend T operator*(T const& t1, T const& t2) { return mul(t1,t2); }
    friend T operator/(T const& t1, NT const& t2) { return div(t1,t2); }
    friend T& operator+=(T& t1, T const& t2) { return t1=add(t1,t2); }
    friend T& operator-=(T& t1, NT const& t2) { return t1=sub(t1,t2); }
    friend T& operator*=(T& t1, T const& t2) { return t1=mul(t1,t2); }
    friend T& operator/=(T& t1, NT const& t2) { return t1=div(t1,t2); }
};

template<class T, class O, class A=O> struct ComparisonOperators {
    typedef decltype(not declval<A>()) E; typedef decltype(not declval<O>()) S;
    friend E operator==(T const& t1, T const& t2);
    friend A operator!=(T const& t1, T const& t2);
    friend O operator<=(T const& t1, T const& t2);
    friend O operator>=(T const& t1, T const& t2);
    friend S operator< (T const& t1, T const& t2);
    friend S operator> (T const& t1, T const& t2);
};

template<class T, class O, class A=O> struct DefineComparisonOperators {
    typedef decltype(not declval<A>()) E; typedef decltype(not declval<O>()) S;
    friend E operator==(T const& t1, T const& t2) { return not neq(t1,t2); }
    friend A operator!=(T const& t1, T const& t2) { return neq(t1,t2); }
    friend O operator<=(T const& t1, T const& t2) { return leq(t1,t2); }
    friend O operator>=(T const& t1, T const& t2) { return leq(t2,t1); }
    friend S operator< (T const& t1, T const& t2) { return not leq(t2,t1); }
    friend S operator> (T const& t1, T const& t2) { return not leq(t1,t2); }
};

template<class T> struct DefineComparisonOperators<T,Boolean> {
    typedef Boolean B;
    friend Comparison cmp(T const& t1, T const& t2);
    friend B operator==(T const& t1, T const& t2) { return cmp(t1,t2)==Comparison::EQUAL; }
    friend B operator!=(T const& t1, T const& t2) { return cmp(t1,t2)!=Comparison::EQUAL; }
    friend B operator<=(T const& t1, T const& t2) { return cmp(t1,t2)!=Comparison::GREATER; }
    friend B operator>=(T const& t1, T const& t2) { return cmp(t1,t2)!=Comparison::LESS; }
    friend B operator< (T const& t1, T const& t2) { return cmp(t1,t2)==Comparison::LESS; }
    friend B operator> (T const& t1, T const& t2) { return cmp(t1,t2)==Comparison::GREATER; }
};


} // namespace Ariadne

#endif
