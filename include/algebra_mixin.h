/***************************************************************************
 *            algebra_mixin.h
 *
 *  Copyright 2010  Pieter Collins
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

/*! \file algebra_mixin.h
 *  \brief Mixin class providing operations for (Banach) algebras.
 */

#ifndef ARIADNE_ALGEBRA_MIXIN_H
#define ARIADNE_ALGEBRA_MIXIN_H

#include "algebra_interface.h"

namespace Ariadne {

template<class A, class X> class AlgebraOperators { };
template<class A, class X> class NormedAlgebraOperators { };
template<class A, class X> class GradedAlgebraOperators { };
template<class A, class X> class SymbolicAlgebraOperators { };

template<class A, class X> class AlgebraMixin
    : public virtual AlgebraInterface<X>
    , public AlgebraOperators<A,X>
{
    typedef X NumericType;
  public:
    virtual AlgebraMixin<A,X>* _create() const { return new A(static_cast<const A&>(*this)*static_cast<X>(0)); }
    virtual AlgebraMixin<A,X>* _clone() const { return new A(static_cast<const A&>(*this)); }
    virtual void _iadd(const X& c) { static_cast<A*>(this)->A::iadd(c); }
    virtual void _imul(const X& c) { static_cast<A*>(this)->A::imul(c); }
    virtual AlgebraInterface<X>* _add(const AlgebraInterface<X>& x) const {
        const A& a1=static_cast<const A&>(*this); const A& a2=dynamic_cast<const A&>(x); return new A(a1+a2); }
    virtual AlgebraInterface<X>* _mul(const AlgebraInterface<X>& x) const {
        const A& a1=static_cast<const A&>(*this); const A& a2=dynamic_cast<const A&>(x); return new A(a1*a2); }
    virtual void _isma(const X& c, const AlgebraInterface<X>& x) {
        static_cast<A*>(this)->A::isma(c,dynamic_cast<const A&>(x)); }
    virtual void _ifma(const AlgebraInterface<X>& x1, const AlgebraInterface<X>& x2)  {
        static_cast<A*>(this)->A::ifma(dynamic_cast<const A&>(x1),dynamic_cast<const A&>(x2)); }
    virtual std::ostream& write(std::ostream& os) const { os << static_cast<const A&>(*this); return os; }
};

template<class A, class X> class GradedAlgebraMixin
    : public virtual GradedAlgebraInterface<X>
    , public AlgebraMixin<A,X>
{
    virtual GradedAlgebraMixin<A,X>* _apply(const Series<X>& f) const { return new A(compose(f,static_cast<const A&>(*this))); }
};

template<class A, class X> class SymbolicAlgebraMixin
    : public virtual SymbolicAlgebraInterface<X>
    , public AlgebraMixin<A,X>
{
    // virtual SymbolicAlgebraMixin<A,X>* _apply(Operator op) { return new A(op,static_cast<const A&>(*this)); }
};

template<class A, class X> inline A operator+(const AlgebraOperators<A,X>& a) {
    A r=static_cast<const A&>(a); return r; }
template<class A, class X> inline A operator-(const AlgebraOperators<A,X>& a) {
    A r=static_cast<const A&>(a); r.A::imul(-1); return r; }

template<class A, class X> inline A operator+(const AlgebraOperators<A,X>& a1, const AlgebraOperators<A,X>& a2) {
    A r=static_cast<const A&>(a1); r.isma(+1,static_cast<const A&>(a2)); return r; }
template<class A, class X> inline A operator-(const AlgebraOperators<A,X>& a1, const AlgebraOperators<A,X>& a2) {
    A r=static_cast<const A&>(a1); r.isma(-1,static_cast<const A&>(a2)); return r; }
template<class A, class X> inline A operator*(const AlgebraOperators<A,X>& a1, const AlgebraOperators<A,X>& a2) {
    A r(a1*0); r.ifma(static_cast<const A&>(a1),static_cast<const A&>(a2)); return r; }

template<class A, class X> inline A& operator+=(AlgebraOperators<A,X>& a1, const AlgebraOperators<A,X>& a2) {
    A& r=static_cast<A&>(a1); r.A::isma(+1,static_cast<const A&>(a2)); return r; }
template<class A, class X> inline A& operator-=(AlgebraOperators<A,X>& a1, const AlgebraOperators<A,X>& a2) {
    A& r=static_cast<A&>(a1); r.A::isma(-1,static_cast<const A&>(a2)); return r; }

template<class A, class X, class C> inline typename EnableIfNumeric<C,A>::Type
operator+(const AlgebraOperators<A,X>& a1, const C& c2) {
    A r=static_cast<const A&>(a1); r.A::iadd(static_cast<X>(c2)); return r; }
template<class A, class X, class C> inline typename EnableIfNumeric<C,A>::Type
operator-(const AlgebraOperators<A,X>& a1, const C& c2) {
    A r=static_cast<const A&>(a1); r.A::iadd(neg(static_cast<X>(c2))); return r; }
template<class A, class X, class C> inline typename EnableIfNumeric<C,A>::Type
operator*(const AlgebraOperators<A,X>& a1, const C& c2) {
    A r=static_cast<const A&>(a1); r.A::imul(static_cast<X>(c2)); return r; }
template<class A, class X, class C> inline typename EnableIfNumeric<C,A>::Type
operator/(const AlgebraOperators<A,X>& a1, const C& c2) {
    A r=static_cast<const A&>(a1); r.A::imul(rec(static_cast<X>(c2))); return r; }
template<class A, class X, class C> inline typename EnableIfNumeric<C,A>::Type
operator+(const C& c1, const AlgebraOperators<A,X>& a2) {
    A r=static_cast<const A&>(a2); r.A::iadd(static_cast<X>(c1)); return r; }
template<class A, class X, class C> inline typename EnableIfNumeric<C,A>::Type
operator-(const C& c1, const AlgebraOperators<A,X>& a2) {
    A r=neg(static_cast<const A&>(a2)); r.A::iadd(static_cast<X>(c1)); return r; }
template<class A, class X, class C> inline typename EnableIfNumeric<C,A>::Type
operator*(const C& c1, const AlgebraOperators<A,X>& a2) {
    A r=static_cast<const A&>(a2); r.A::imul(c1); return r; }

template<class A, class X, class C> inline typename EnableIfNumeric<C,A&>::Type
operator+=(AlgebraOperators<A,X>& a, const C& c) { A& r=static_cast<A&>(a); r.A::iadd(static_cast<X>(c)); return r; }
template<class A, class X, class C> inline typename EnableIfNumeric<C,A&>::Type
operator-=(AlgebraOperators<A,X>& a, const C& c) { A& r=static_cast<A&>(a); r.A::iadd(neg(static_cast<X>(c))); return r; }
template<class A, class X, class C> inline typename EnableIfNumeric<C,A&>::Type
operator*=(AlgebraOperators<A,X>& a, const C& c) { A& r=static_cast<A&>(a); r.A::imul(static_cast<X>(c)); return r; }
template<class A, class X, class C> inline typename EnableIfNumeric<C,A&>::Type
operator/=(AlgebraOperators<A,X>& a, const C& c) { A& r=static_cast<A&>(a); r.A::imul(rec(static_cast<X>(c))); return r; }



template<class A, class X> inline A neg(const AlgebraOperators<A,X>& x) {
    return -static_cast<const A&>(x); }
template<class A, class X> A inline sqr(const AlgebraOperators<A,X>& x) {
    return static_cast<const A&>(x)*static_cast<const A&>(x); }
template<class A, class X> A pow(const AlgebraOperators<A,X>& x, uint m) {
    A s=static_cast<const A&>(x); A r=s.create(); r+=static_cast<X>(1); while(m) { if(m%2) { r=r*s; } s=sqr(s); m/=2; } return r; }


template<class A, class X> inline A operator/(const NormedAlgebraOperators<A,X>& a1, const NormedAlgebraOperators<A,X>& a2) {
    return static_cast<const A&>(a1)*rec(static_cast<const A&>(a2)); }
template<class A, class X> inline A operator/(const X& c1, const NormedAlgebraOperators<A,X>& a2) {
    return c1*rec(static_cast<const A&>(a2)); }


template<class A, class X> A rec(const NormedAlgebraOperators<A,X>& a);
template<class A, class X> A sqrt(const NormedAlgebraOperators<A,X>& a);
template<class A, class X> A exp(const NormedAlgebraOperators<A,X>& a);
template<class A, class X> A log(const NormedAlgebraOperators<A,X>& a);
template<class A, class X> A sin(const NormedAlgebraOperators<A,X>& a);
template<class A, class X> A cos(const NormedAlgebraOperators<A,X>& a);
template<class A, class X> A tan(const NormedAlgebraOperators<A,X>& a);


} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_MIXIN_H */