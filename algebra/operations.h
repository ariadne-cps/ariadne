/***************************************************************************
 *            algebra/operations.h
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

/*! \file operations.h
 *  \brief Provide common algebra operations automatically.
 */

#ifndef ARIADNE_ALGEBRA_OPERATIONS_H
#define ARIADNE_ALGEBRA_OPERATIONS_H

#include "numeric/operators.h"
#include "numeric/arithmetic.h"

namespace Ariadne {


template<class A, class X=typename A::NumericType> class AlgebraOperations {
    template<class OP> static A _create(OP, A const&);
    template<class OP> static A _create(OP, A const&, A const&);
    template<class OP> static A _create(OP, X const&, A const&);
    template<class OP> static A _create(OP, A const&, X const&);
  public:
    static A _nul(A const& a);
    static A _pos(A const& a);
    static A _neg(A const& a);
    static A _hlf(A const& a);
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

template<class A> class GradedAlgebraOperations {
    typedef typename A::NumericType X;
  public:
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

template<class A> class NormedAlgebraOperations {
    typedef typename A::NumericType X;
  public:
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

template<class V, class X> struct DeclareVectorOperators
{
    friend V operator+(V const& v);
    friend V operator-(V const& v);
    friend V operator+(V const& v1, V const& v2);
    friend V operator-(V const& v1, V const& v2);
    friend V operator*(V const& v1, X const& x2);
    friend V operator/(V const& v1, X const& x2);
    friend V operator*(X const& x1, V const& v2);
    friend V& operator+=(V& v1, V const& v2);
    friend V& operator-=(V& v1, V const& v2);
    friend V& operator*=(V& v1, X const& x2);
    friend V& operator/=(V& v1, X const& x2);
};

template<class V, class A, class X> struct DeclareVectorAlgebraOperators
    : DeclareVectorOperators<V,A>
{
    friend V operator*(V const& v1, X const& x2);
    friend V operator/(V const& v1, X const& x2);
    friend V operator*(X const& x1, V const& v2);
    friend V& operator*=(V& v1, X const& x2);
    friend V& operator/=(V& v1, X const& x2);
};

template<class A, class X> struct DeclareAlgebraOperators
    : DeclareVectorOperators<A,X>
{
    friend A operator*(A const& a1, A const& a2);
    friend A& operator*=(A& a1, A const& a2);
    friend A operator+(A const& a1, X const& x2);
    friend A operator-(A const& a1, X const& x2);
    friend A operator+(X const& x1, A const& a2);
    friend A operator-(X const& x1, A const& a2);
    friend A& operator+=(A& a1, X const& x2);
    friend A& operator-=(A& a1, X const& x2);
};

template<class A, class X> struct DeclareAlgebraOperations
    : DeclareAlgebraOperators<A,X>
{
    friend A nul(A const& a);
    friend A pos(A const& a);
    friend A neg(A const& a);
    friend A hlf(A const& a);
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

template<class A, class X, class Y> struct DispatchConcreteGenericAlgebraNumberOperations {
    friend X create(Y const& y, A const& a) { return X(y,a.precision()); }
    friend A& operator+=(A& a1, Y const& y2) { X x2=create(y2,a1); return a1+=x2; }
    friend A& operator-=(A& a1, Y const& y2) { X x2=create(y2,a1); return a1-=x2; }
    friend A& operator*=(A& a1, Y const& y2) { X x2=create(y2,a1); return a1*=x2; }
    friend A& operator/=(A& a1, Y const& y2) { X x2=create(y2,a1); return a1/=x2; }
    friend A operator+(A a1, Y const& y2) { a1+=y2; return std::move(a1); }
    friend A operator-(A a1, Y const& y2) { a1-=y2; return std::move(a1); }
    friend A operator*(A a1, Y const& y2) { a1*=y2; return std::move(a1); }
    friend A operator/(A a1, Y const& y2) { a1/=y2; return std::move(a1); }
    friend A operator+(Y const& y1, A a2) { return a2+y1; }
    friend A operator-(Y const& y1, A a2) { return (-a2)+y1; }
    friend A operator*(Y const& y1, A a2) { return a2*y1; }
    friend A operator/(Y const& y1, A a2) { return rec(a2)*y1; }
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
    friend A hlf(A const& a) { return OperationsType()._hlf(a); }
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
    friend A sqr(A const& a) { return OperationsType::_create(Sqr(),a); }
    friend A add(A const& a1, A const& a2) { return OperationsType::_create(Add(),a1,a2); }
    friend A sub(A const& a1, A const& a2) { return OperationsType::_create(Sub(),a1,a2); }
    friend A mul(A const& a1, A const& a2) { return OperationsType::_create(Mul(),a1,a2); }
    friend A div(A const& a1, A const& a2) { return OperationsType::_create(Div(),a1,a2); }
    friend A pow(A const& a, Int n) { return OperationsType::_create(Pow(),a,n); }

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
    // BUG: Cannot declare these due to bug in Clang
    // friend A operator/(A const& a1, X const& x2);
    // friend A div(A const& a1, X const& x2);
    // friend A pow(A const& a, Nat m);

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


} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_OPERATIONS_H */
