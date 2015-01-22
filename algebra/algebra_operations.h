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

#include "expression/operators.h"

namespace Ariadne {

template<class X> class Algebra;
template<class X> class NormedAlgebra;
template<class X> class GradedAlgebra;
template<class X> class SymbolicAlgebra;

template<class T> using NumericType = typename T::NumericType;

template<class A> struct IsAlgebra { static const Bool value = false; };
template<class X> struct IsAlgebra< Algebra<X> > { static const Bool value = true; };
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
    A s=x; A r=s.create_zero(); r+=1; while(m) { if(m%2) { r*=s; } s=sqr(s); m/=2; } return r; }

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
