/***************************************************************************
 *            algebra_operations.h
 *
 *  Copyright 2010-15  Pieter Collins
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

/*! \file algebra_operations.h
 *  \brief Provide common algebra operations automatically.
 */

#ifndef ARIADNE_ALGEBRA_OPERATIONS_H
#define ARIADNE_ALGEBRA_OPERATIONS_H

#include "numeric/operators.h"
#include "numeric/mixins.h"

namespace Ariadne {



template<class A, class X=typename A::NumericType> struct AlgebraOperations {
  public:
    template<class OP> static A _create(OP, A const&);
    template<class OP> static A _create(OP, A const&, A const&);
    template<class OP> static A _create(OP, X const&, A const&);
    template<class OP> static A _create(OP, A const&, X const&);

    static A _nul(A const& a);
    static A _pos(A const& a);
    static A _neg(A const& a);
    static A _half(A const& a);
    static A _sqr(A const& a);
    static A _rec(A const& a);
    static A _add(A const& a1, A const& a2);
    static A _sub(A const& a1, A const& a2);
    static A _mul(A const& a1, A const& a2);
    static A _div(A const& a1, A const& a2);
    static A _add(A const& a1, X const& x2);
    static A _sub(A const& a1, X const& x2);
    static A _mul(A const& a1, X const& x2);
    static A _div(A const& a1, X const& x2);
    static A _add(X const& x1, A const& a2);
    static A _sub(X const& x1, A const& a2);
    static A _mul(X const& x1, A const& a2);
    static A _div(X const& x1, A const& a2);
    static A _pow(A const& a, Nat m);
    static A _pow(A const& a, Int n);
    static A _sqrt(A const& a);
    static A _exp(A const& a);
    static A _log(A const& a);
    static A _sin(A const& a);
    static A _cos(A const& a);
    static A _tan(A const& a);
    static A _asin(A const& a);
    static A _acos(A const& a);
    static A _atan(A const& a);
};

template<class A> struct GradedAlgebraOperations {
    typedef typename A::NumericType X;
    static A _sqr(A const& a) {
        return a*a; }
    static A _div(A const& a1, A const& a2) {
        return mul(a1,rec(a2)); }
    static A _pow(A const& a, Nat m) {
        A s=a; A r=s.create_constant(1); while(m) { if(m%2) { r*=s; } s=sqr(s); m/=2; } return r; }
    static A _pow(A const& a, Int n) {
        if(n>=0) { return _pow(a,Nat(n)); } else { return rec(_pow(a,Nat(-n))); } }
    static A _compose(const Series<X>& f, const A& a) { return A::_compose(f,a); }
    static A _rec(const A& a) { return _compose(Series<X>::rec(a.value()),a); }
    static A _sqrt(const A& a) { return _compose(Series<X>::sqrt(a.value()),a); }
    static A _exp(const A& a) { return _compose(Series<X>::exp(a.value()),a); }
    static A _log(const A& a) { return _compose(Series<X>::log(a.value()),a); }
    static A _sin(const A& a) { return _compose(Series<X>::sin(a.value()),a); }
    static A _cos(const A& a) { return _compose(Series<X>::cos(a.value()),a); }
    static A _tan(const A& a) { return _compose(Series<X>::tan(a.value()),a); }
    static A _atan(const A& a) { return _compose(Series<X>::atan(a.value()),a); }
};

template<class A> struct NormedAlgebraOperations {
    typedef typename A::NumericType X;
    static A _sqr(A const& a) {
        return a*a; }
    static A _pow(A const& a, Nat m) {
        A s=a; A r=s.create_constant(1); while(m) { if(m%2) { r*=s; } s=sqr(s); m/=2; } return r; }
    static A _pow(A const& a, Int n) {
        if(n>=0) { return _pow(a,Nat(n)); } else { return rec(_pow(a,Nat(-n))); } }
    static A _div(A const& a1, A const& a2) {
        return mul(a1,rec(a2)); }
    static A _rec(const A& a);
    static A _sqrt(const A& a);
    static A _exp(const A& a);
    static A _log(const A& a);
    static A _sin(const A& a);
    static A _cos(const A& a);
    static A _tan(const A& a);
    static A _asin(const A& a);
    static A _acos(const A& a);
    static A _atan(const A& a);
};

template<class A, class X> struct DeclareAlgebraOperators
{
    friend A operator+(A const& a);
    friend A operator-(A const& a);
    friend A operator+(A const& a1, A const& a2);
    friend A operator-(A const& a1, A const& a2);
    friend A operator*(A const& a1, A const& a2);
    friend A& operator+=(A& a1, A const& a2);
    friend A& operator-=(A& a1, A const& a2);
    friend A& operator*=(A& a1, A const& a2);
    friend A operator+(A const& a1, X const& x2);
    friend A operator-(A const& a1, X const& x2);
    friend A operator*(A const& a1, X const& x2);
    friend A operator/(A const& a1, X const& x2);
    friend A& operator+=(A& a1, X const& x2);
    friend A& operator-=(A& a1, X const& x2);
    friend A& operator*=(A& a1, X const& x2);
    friend A& operator/=(A& a1, X const& x2);
    friend A operator+(X const& x1, A const& a2);
    friend A operator-(X const& x1, A const& a2);
    friend A operator*(X const& x1, A const& a2);
};

template<class A, class X> struct DeclareAlgebraOperations
    : DeclareAlgebraOperators<A,X>
{
    friend A nul(A const& a);
    friend A pos(A const& a);
    friend A neg(A const& a);
    friend A half(A const& a);
    friend A sqr(A const& a);
    friend A add(A const& a1, A const& a2);
    friend A sub(A const& a1, A const& a2);
    friend A mul(A const& a1, A const& a2);

    friend A add(A const& a1, X const& x2);
    friend A sub(A const& a1, X const& x2);
    friend A mul(A const& a1, X const& x2);
    friend A div(A const& a1, X const& x2);
    friend A add(X const& x1, A const& a2);
    friend A sub(X const& x1, A const& a2);
    friend A mul(X const& x1, A const& a2);
    friend A pow(A const& a, Nat m);
    friend A pow(A const& a, Int n);
};

template<class A, class X> struct DeclareDivisionAlgebraOperations : DeclareAlgebraOperations<A,X> {
    friend A operator/(A const& a1, A const& a2);
    friend A& operator/=(A& a1, A const& a2);
    friend A operator/(X const& x1, A const& a2);

    friend A rec(A const& a);
    friend A div(A const& a1, A const& a2);
    friend A div(X const& x1, A const& a2);
    friend A pow(A const& a, Int n);
};

template<class A, class X> struct DeclareTranscendentalAlgebraOperations : DeclareDivisionAlgebraOperations<A,X> {
    friend A sqrt(A const& a);
    friend A exp(A const& a);
    friend A log(A const& a);
    friend A sin(A const& a);
    friend A cos(A const& a);
    friend A tan(A const& a);
    friend A asin(A const& a);
    friend A acos(A const& a);
    friend A atan(A const& a);
};

template<class A, class Y> struct DispatchMixedAlgebraNumberOperations;

template<class A> struct DispatchMixedAlgebraNumberOperations<A,int> {
    template<class N, class T> using IfInt = typename EnableIf<IsIntegral<N>,T>::Type;
    template<class N> friend IfInt<N,A> operator+(A a1, N x2);
    friend A& operator+=(A& a1, int x2) { return a1+=ValidatedNumericType(x2); }
    friend A& operator-=(A& a1, int x2) { return a1-=ValidatedNumericType(x2); }
    friend A& operator*=(A& a1, int x2) { return a1*=ValidatedNumericType(x2); }
    friend A& operator/=(A& a1, int x2) { return a1/=ValidatedNumericType(x2); }
    friend A operator+(A a1, int x2) { a1+=x2; return std::move(a1); }
    friend A operator-(A a1, int x2) { a1-=x2; return std::move(a1); }
    friend A operator*(A a1, int x2) { a1*=x2; return std::move(a1); }
    friend A operator/(A a1, int x2) { a1/=x2; return std::move(a1); }
    friend A operator+(int x1, A a2) { return a2+x1; }
    friend A operator-(int x1, A a2) { return (-a2)+x1; }
};


template<class A, class X> struct DispatchAlgebraOperators
{
    typedef AlgebraOperations<A,X> OperationsType;
   public:
    friend A operator+(A const& a) { return pos(a); }
    friend A operator-(A const& a) { return neg(a); }
    friend A operator+(A const& a1, A const& a2) { return add(a1,a2); }
    friend A operator-(A const& a1, A const& a2) { return sub(a1,a2); }
    friend A operator*(A const& a1, A const& a2) { return mul(a1,a2); }
    friend A& operator+=(A& a1, A const& a2) { return a1=add(a1,a2); }
    friend A& operator-=(A& a1, A const& a2) { return a1=sub(a1,a2); }
    friend A& operator*=(A& a1, A const& a2) { return a1=mul(a1,a2); }
    friend A operator+(A const& a1, X const& x2) { return add(a1,x2); }
    friend A operator-(A const& a1, X const& x2) { return sub(a1,x2); }
    friend A operator*(A const& a1, X const& x2) { return mul(a1,x2); }
    friend A operator/(A const& a1, X const& x2) { return div(a1,x2); }
    friend A& operator+=(A& a1, X const& x2) { return a1=add(a1,x2); }
    friend A& operator-=(A& a1, X const& x2) { return a1=sub(a1,x2); }
    friend A& operator*=(A& a1, X const& x2) { return a1=mul(a1,x2); }
    friend A& operator/=(A& a1, X const& x2) { return a1=div(a1,x2); }
    friend A operator+(X const& x1, A const& a2) { return add(x1,a2); }
    friend A operator-(X const& x1, A const& a2) { return sub(x1,a2); }
    friend A operator*(X const& x1, A const& a2) { return mul(x1,a2); }
};

template<class A, class X> struct DispatchAlgebraOperations
    : DispatchAlgebraOperators<A,X>
{
    typedef AlgebraOperations<A,X> OperationsType;
   public:
    friend A nul(A const& a) { return OperationsType()._nul(a); }
    friend A pos(A const& a) { return OperationsType()._pos(a); }
    friend A neg(A const& a) { return OperationsType()._neg(a); }
    friend A half(A const& a) { return OperationsType()._half(a); }
    friend A sqr(A const& a) { return OperationsType()._sqr(a); }
    friend A add(A const& a1, A const& a2) { return OperationsType()._add(a1, a2); }
    friend A sub(A const& a1, A const& a2) { return OperationsType()._sub(a1, a2); }
    friend A mul(A const& a1, A const& a2) { return OperationsType()._mul(a1, a2); }

    friend A add(A const& a1, X const& x2) { return OperationsType()._add(a1, x2); }
    friend A sub(A const& a1, X const& x2) { return OperationsType()._add(a1, neg(x2)); }
    friend A mul(A const& a1, X const& x2) { return OperationsType()._mul(a1, x2); }
    friend A div(A const& a1, X const& x2) { return OperationsType()._mul(a1, rec(x2)); }
    friend A add(X const& x1, A const& a2) { return OperationsType()._add(a2, x1); }
    friend A sub(X const& x1, A const& a2) { return OperationsType()._add(neg(a2), x1); }
    friend A mul(X const& x1, A const& a2) { return OperationsType()._mul(a2, x1); }
    friend A pow(A const& a, Nat m) { return OperationsType()._pow(a, m); }
};

template<class A, class X> struct DispatchSymbolicAlgebraOperations
    : DispatchAlgebraOperators<A,X>
{
    typedef A OperationsType;
    //typedef AlgebraOperations<A,X> OperationsType;
  public:
    friend A operator/(A const& a1, A const& a2) { return div(a1,a2); }
    friend A& operator/=(A& a1, A const& a2) { return a1=div(a1,a2); }
    friend A operator/(X const& x1, A const& a2) { return div(x1,a2); }
  public:
    friend A pos(A const& a) { return OperationsType::_create(Pos(),a); }
    friend A neg(A const& a) { return OperationsType::_create(Neg(),a); }
    friend A add(A const& a1, A const& a2) { return OperationsType::_create(Add(),a1,a2); }
    friend A sub(A const& a1, A const& a2) { return OperationsType::_create(Sub(),a1,a2); }
    friend A mul(A const& a1, A const& a2) { return OperationsType::_create(Mul(),a1,a2); }
    friend A div(A const& a1, A const& a2) { return OperationsType::_create(Div(),a1,a2); }

    friend A add(A const& a1, X const& x2) { return OperationsType::_create(Add(),a1,x2); }
    friend A sub(A const& a1, X const& x2) { return OperationsType::_create(Sub(),a1,x2); }
    friend A mul(A const& a1, X const& x2) { return OperationsType::_create(Mul(),a1,x2); }
    friend A div(A const& a1, X const& x2) { return OperationsType::_create(Div(),a1,x2); }
    friend A add(X const& x1, A const& a2) { return OperationsType::_create(Add(),x1,a2); }
    friend A sub(X const& x1, A const& a2) { return OperationsType::_create(Sub(),x1,a2); }
    friend A mul(X const& x1, A const& a2) { return OperationsType::_create(Mul(),x1,a2); }
    friend A div(X const& x1, A const& a2) { return OperationsType::_create(Div(),x1,a2); }

    friend A rec(A const& a) { return OperationsType::_create(Rec(),a); }
    friend A sqrt(A const& a) { return OperationsType::_create(Sqrt(),a); }
    friend A exp(A const& a) { return OperationsType::_create(Exp(),a); }
    friend A log(A const& a) { return OperationsType::_create(Log(),a); }
    friend A sin(A const& a) { return OperationsType::_create(Sin(),a); }
    friend A cos(A const& a) { return OperationsType::_create(Cos(),a); }
    friend A tan(A const& a) { return OperationsType::_create(Tan(),a); }
    friend A asin(A const& a) { return OperationsType::_create(Asin(),a); }
    friend A acos(A const& a) { return OperationsType::_create(Acos(),a); }
    friend A atan(A const& a) { return OperationsType::_create(Atan(),a); }
};

template<class A, class X> struct DispatchTranscendentalAlgebraOperations : DispatchAlgebraOperations<A,X> {
    typedef AlgebraOperations<A,X> OperationsType;
  public:
    friend A operator/(A const& a1, X const& x2);
    friend A div(A const& a1, X const& x2);
    friend A pow(A const& a, Nat m);

    friend A operator/(A const& a1, A const& a2) { return div(a1,a2); }
    friend A& operator/=(A& a1, A const& a2) { return a1=div(a1,a2); }
    friend A operator/(X const& x1, A const& a2) { return div(x1,a2); }
    friend A div(A const& a1, A const& a2) { return OperationsType::_div(a1, a2); }

    friend A rec(A const& a) { return OperationsType()._rec(a); }
    friend A div(X const& x1, A const& a2) { return OperationsType()._mul(rec(a2), x1); }
    friend A pow(A const& a, Int n) { return OperationsType()._pow(a, n); }

    friend A sqrt(A const& a) { return OperationsType()._sqrt(a); }
    friend A exp(A const& a) { return OperationsType()._exp(a); }
    friend A log(A const& a) { return OperationsType()._log(a); }
    friend A sin(A const& a) { return OperationsType()._sin(a); }
    friend A cos(A const& a) { return OperationsType()._cos(a); }
    friend A tan(A const& a) { return OperationsType()._tan(a); }
    friend A asin(A const& a) { return OperationsType()._asin(a); }
    friend A acos(A const& a) { return OperationsType()._acos(a); }
    friend A atan(A const& a) { return OperationsType()._atan(a); }
};

template<class A, class X> struct DispatchOrderedAlgebraOperations : DispatchAlgebraOperations<A,X> {
    typedef AlgebraOperations<A,X> OperationsType;
    friend A max(A const& a1, A const& a2) { return OperationsType()._max(a1, a2); }
    friend A min(A const& a1, A const& a2) { return OperationsType()._min(a1, a2); }
    friend A abs(A const& a) { return OperationsType()._abs(a); }
};


template<class X> class Algebra;
template<class X> class NormedAlgebra;
template<class X> class GradedAlgebra;
template<class X> class SymbolicAlgebra;

template<class T> using NumericType = typename T::NumericType;

template<class A> struct IsAlgebra { static const Bool value = false; };
//template<class X> struct IsAlgebra< Algebra<X> > { static const Bool value = true; };
template<class X> struct IsAlgebra< GradedAlgebra<X> > { static const Bool value = true; };
template<class X> struct IsAlgebra< NormedAlgebra<X> > { static const Bool value = true; };
template<class X> struct IsAlgebra< SymbolicAlgebra<X> > { static const Bool value = true; };
template<class A> using EnableIfAlgebra = EnableIf<IsAlgebra<A>,A>;

template<class A> struct IsNormedAlgebra { static const Bool value = false; };
template<class X> struct IsNormedAlgebra< NormedAlgebra<X> > { static const Bool value = true; };
template<class A> using EnableIfNormedAlgebra = EnableIf<IsNormedAlgebra<A>,A>;

template<class A> struct IsGradedAlgebra { static const Bool value = false; };
template<class X> struct IsGradedAlgebra< GradedAlgebra<X> > { static const Bool value = true; };
template<class A> using EnableIfGradedAlgebra = EnableIf<IsGradedAlgebra<A>,A>;

template<class A> struct IsSymbolicAlgebra { static const Bool value = false; };
template<class X> struct IsSymbolicAlgebra< SymbolicAlgebra<X> > { static const Bool value = true; };
template<class A> using EnableIfSymbolicAlgebra = EnableIf<IsSymbolicAlgebra<A>,A>;


// The following operators should be defined in order to define other algebra operations
template<class A> EnableIfAlgebra<A> operator+(const A& a1, const A& a2);
template<class A> EnableIfAlgebra<A> operator*(const A& a1, const A& a2);
template<class A> EnableIfAlgebra<A>& operator+=(A& a, const NumericType<A>& c);
template<class A> EnableIfAlgebra<A>& operator*=(A& a, const NumericType<A>& c);

template<class A> inline EnableIfAlgebra<A> operator+(A a) { return std::move(a); }
template<class A> inline EnableIfAlgebra<A> operator-(A a) { a.imul(-1); return std::move(a); }
template<class A> inline EnableIfAlgebra<A> operator+(const A& a1, const A& a2) { A r=a1; r.isma(+1,a2); return std::move(r); }
template<class A> inline EnableIfAlgebra<A> operator-(const A& a1, const A& a2) { return a1+(-a2); }
template<class A> inline EnableIfAlgebra<A> operator*(const A& a1, const A& a2) { A r=a1.create(); r.ifma(a1,a2); return std::move(r); }

template<class A> inline EnableIfAlgebra<A>& operator+=(A& a1, const A& a2) { return a1=a1+a2; }
template<class A> inline EnableIfAlgebra<A>& operator-=(A& a1, const A& a2) { return a1=a1-a2; }
template<class A> inline EnableIfAlgebra<A>& operator*=(A& a1, const A& a2) { return a1=a1*a2; }

template<class A> inline EnableIfAlgebra<A>& operator+=(A& a, const NumericType<A>& c) { a.iadd(c); return a; }
template<class A> inline EnableIfAlgebra<A>& operator-=(A& a, const NumericType<A>& c) { return a+=neg(c); }
template<class A> inline EnableIfAlgebra<A>& operator*=(A& a, const NumericType<A>& c) { a.imul(c); return a; }
template<class A> inline EnableIfAlgebra<A>& operator/=(A& a, const NumericType<A>& c) { return a*=rec(c); }

template<class A> inline EnableIfAlgebra<A> operator+(A a, const NumericType<A>& c) { a+=c; return std::move(a); }
template<class A> inline EnableIfAlgebra<A> operator-(A a, const NumericType<A>& c) { a-=c; return std::move(a); }
template<class A> inline EnableIfAlgebra<A> operator*(A a, const NumericType<A>& c) { a*=c; return std::move(a); }
template<class A> inline EnableIfAlgebra<A> operator/(A a, const NumericType<A>& c) { a/=c; return std::move(a); }
template<class A> inline EnableIfAlgebra<A> operator+(const NumericType<A>& c, A a) { a+=c; return std::move(a); }
template<class A> inline EnableIfAlgebra<A> operator-(const NumericType<A>& c, A a) { a-=c; return neg(std::move(a)); }
template<class A> inline EnableIfAlgebra<A> operator*(const NumericType<A>& c, A a) { a*=c; return std::move(a); }

template<class A> inline EnableIfAlgebra<A> add(const A& a1, const A& a2) { return a1+a2; }
template<class A> inline EnableIfAlgebra<A> sub(const A& a1, const A& a2) { return a1-a2; }
template<class A> inline EnableIfAlgebra<A> mul(const A& a1, const A& a2) { return a1*a2; }
template<class A> inline EnableIfAlgebra<A> pos(const A& a) { return +a; }
template<class A> inline EnableIfAlgebra<A> neg(const A& a) { return -a; }
template<class A> inline EnableIfAlgebra<A> sqr(const A& a) { return a*a; }

template<class A>  inline EnableIfAlgebra<A> pow(const A& x, Nat m) {
    A s=x; A r=s.create_constant(1); while(m) { if(m%2) { r*=s; } s=sqr(s); m/=2; } return r; }


// Operations requiring reciprocal
template<class A> inline EnableIfAlgebra<A> pow(const A& a, Int n) { return n<0 ? rec(pow(a,Nat(-n))) : pow(a,Nat(n)); }
template<class A> inline EnableIfAlgebra<A> operator/(const A& a1, const A& a2) { return a1*rec(a2); }
template<class A> inline EnableIfAlgebra<A> operator/(const NumericType<A>& c, const A& a) { return c*rec(a); }

template<class A> EnableIfNormedAlgebra<A> rec(const A& a);
template<class A> EnableIfNormedAlgebra<A> sqrt(const A& a);
template<class A> EnableIfNormedAlgebra<A> exp(const A& a);
template<class A> EnableIfNormedAlgebra<A> log(const A& a);
template<class A> EnableIfNormedAlgebra<A> sin(const A& a);
template<class A> EnableIfNormedAlgebra<A> cos(const A& a);
template<class A> EnableIfNormedAlgebra<A> tan(const A& a);
template<class A> EnableIfNormedAlgebra<A> atan(const A& a);

template<class A> EnableIfGradedAlgebra<A> compose(const Series<NumericType<A>>& x, const A& y);
template<class A> EnableIfGradedAlgebra<A> rec(const A& a) { return compose(Series<NumericType<A>>::rec(a.value()),a); }
template<class A> EnableIfGradedAlgebra<A> sqrt(const A& a) { return compose(Series<NumericType<A>>::sqrt(a.value()),a); }
template<class A> EnableIfGradedAlgebra<A> exp(const A& a) { return compose(Series<NumericType<A>>::exp(a.value()),a); }
template<class A> EnableIfGradedAlgebra<A> log(const A& a) { return compose(Series<NumericType<A>>::log(a.value()),a); }
template<class A> EnableIfGradedAlgebra<A> sin(const A& a) { return compose(Series<NumericType<A>>::sin(a.value()),a); }
template<class A> EnableIfGradedAlgebra<A> cos(const A& a) { return compose(Series<NumericType<A>>::cos(a.value()),a); }
template<class A> EnableIfGradedAlgebra<A> tan(const A& a) { return compose(Series<NumericType<A>>::tan(a.value()),a); }
template<class A> EnableIfGradedAlgebra<A> atan(const A& a) { return compose(Series<NumericType<A>>::atan(a.value()),a); }

template<class A> EnableIfSymbolicAlgebra<A> rec(const A& a) { return apply(Rec(),a); }
template<class A> EnableIfSymbolicAlgebra<A> sqrt(const A& a) { return apply(Sqrt(),a); }
template<class A> EnableIfSymbolicAlgebra<A> exp(const A& a) { return apply(Exp(),a); }
template<class A> EnableIfSymbolicAlgebra<A> log(const A& a) { return apply(Log(),a); }
template<class A> EnableIfSymbolicAlgebra<A> sin(const A& a) { return apply(Sin(),a); }
template<class A> EnableIfSymbolicAlgebra<A> cos(const A& a) { return apply(Cos(),a); }
template<class A> EnableIfSymbolicAlgebra<A> tan(const A& a) { return apply(Tan(),a); }
template<class A> EnableIfSymbolicAlgebra<A> atan(const A& a) { return apply(Atan(),a); }


} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_OPERATIONS_H */
