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

namespace Ariadne {

template<class A, class X> class AlgebraMixin { };
template<class A, class X> class NormedAlgebraMixin : public AlgebraMixin<A,X> { };

template<class A, class X> inline A operator+(const AlgebraMixin<A,X>& a) {
    return static_cast<const A&>(a); }
template<class A, class X> inline A operator-(const AlgebraMixin<A,X>& a) {
    A r=static_cast<const A&>(a); r.imul(-1); return r; }

template<class A, class X> inline A operator+(const AlgebraMixin<A,X>& a1, const AlgebraMixin<A,X>& a2) {
    return static_cast<const A&>(a1).sma(+1,static_cast<const A&>(a2)); }
template<class A, class X> inline A operator-(const AlgebraMixin<A,X>& a1, const AlgebraMixin<A,X>& a2) {
    return static_cast<const A&>(a1).sma(-1,static_cast<const A&>(a2)); }
template<class A, class X> inline A operator*(const AlgebraMixin<A,X>& a1, const AlgebraMixin<A,X>& a2) {
    A r=static_cast<const A&>(a1).null(); r.ifma(static_cast<const A&>(a1),static_cast<const A&>(a2)); return r; }

template<class A, class X> inline A& operator+=(AlgebraMixin<A,X>& a1, const AlgebraMixin<A,X>& a2) {
    A& r=static_cast<A&>(a1); r.isma(+1,static_cast<const A&>(a2)); return r; }
template<class A, class X> inline A& operator-=(AlgebraMixin<A,X>& a1, const AlgebraMixin<A,X>& a2) {
    A& r=static_cast<A&>(a1); r.isma(-1,static_cast<const A&>(a2)); return r; }

template<class A, class X, class C> inline
typename enable_if< is_numeric<C>::value, A>::type
operator+(const AlgebraMixin<A,X>& a1, const C& c2) {
    A r=static_cast<const A&>(a1); r.iadd(static_cast<X>(c2)); return r; }
template<class A, class X> inline A operator-(const AlgebraMixin<A,X>& a1, const X& c2) {
    A r=static_cast<const A&>(a1); r.iadd(neg(c2)); return r; }
template<class A, class X> inline A operator*(const AlgebraMixin<A,X>& a1, const X& c2) {
    A r=static_cast<const A&>(a1); r.imul(c2); return r; }
template<class A, class X> inline A operator/(const AlgebraMixin<A,X>& a1, const X& c2) {
    A r=static_cast<const A&>(a1); r.imul(rec(c2)); return r; }
template<class A, class X> inline A operator+(const X& c1, const AlgebraMixin<A,X>& a2) {
    A r=static_cast<const A&>(a2); r.iadd(c1); return r; }
template<class A, class X> inline A operator-(const X& c1, const AlgebraMixin<A,X>& a2) {
    A r=neg(static_cast<const A&>(a2)); r.iadd(c1); return r; }
template<class A, class X> inline A operator*(const X& c1, const AlgebraMixin<A,X>& a2) {
    A r=static_cast<const A&>(a2); r.imul(c1); return r; }

template<class X, class C> inline enable_if< is_numeric<C>, Algebra<X>& >::type
operator+=(Algebra<X>& x, const C& c) { x.iadd(static_cast<X>(c)); return x; }
template<class X, class C> inline enable_if< is_numeric<C>, Algebra<X>& >::type
operator*=(Algebra<X>& x, const X& c) { x.imul(static_cast<X>(c)); return x; }

A& operator+=(AlgebraMixin<A,X>& a, const X& c) {
    A& r=static_cast<A&>(a); r.iadd(c); return r; }
template<class A, class X> inline A& operator-=(AlgebraMixin<A,X>& a, const X& c) {
    A& r=static_cast<A&>(a); r.iadd(neg(c)); return r; }
template<class A, class X> inline A& operator*=(AlgebraMixin<A,X>& a, const X& c) {
    A& r=static_cast<A&>(a); r.imul(c); return r; }
template<class A, class X> inline A& operator/=(AlgebraMixin<A,X>& a, const X& c) {
    A& r=static_cast<A&>(a); r.imul(rec(c)); return r; }


template<class A, class X> inline A neg(const AlgebraMixin<A,X>& x) {
    return -static_cast<const A&>(x); }
template<class A, class X> A inline sqr(const AlgebraMixin<A,X>& x) {
    return static_cast<const A&>(x)*static_cast<const A&>(x); }
template<class A, class X> A pow(const AlgebraMixin<A,X>& x, uint m) {
    A s=static_cast<const A&>(x); A r=s.null(); r+=static_cast<X>(1); while(m) { if(m%2) { r=r*s; } s=sqr(s); m/=2; } return r; }



template<class A, class X> inline A operator/(const NormedAlgebraMixin<A,X>& a1, const NormedAlgebraMixin<A,X>& a2) {
    return static_cast<const A&>(a1)*rec(static_cast<const A&>(a2)); }
template<class A, class X> inline A operator/(const X& c1, const NormedAlgebraMixin<A,X>& a2) {
    return c1*rec(static_cast<const A&>(a2)); }


template<class A, class X> A rec(const NormedAlgebraMixin<A,X>& a);
template<class A, class X> A sqrt(const NormedAlgebraMixin<A,X>& a);
template<class A, class X> A exp(const NormedAlgebraMixin<A,X>& a);
template<class A, class X> A log(const NormedAlgebraMixin<A,X>& a);
template<class A, class X> A sin(const NormedAlgebraMixin<A,X>& a);
template<class A, class X> A cos(const NormedAlgebraMixin<A,X>& a);
template<class A, class X> A tan(const NormedAlgebraMixin<A,X>& a);

} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_MIXIN_H */