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

#include "algebra/algebra_interface.h"
#include "expression/operators.h"

namespace Ariadne {

template<class X> class Differential;
template<class X> class TaylorModel;

template<class X> class Algebra;
template<class X> class NormedAlgebra;
template<class X> class GradedAlgebra;
template<class X> class SymbolicAlgebra;

template<class A> struct IsAlgebra { static const bool value = false; };
template<class X> struct IsAlgebra< Algebra<X> > { static const bool value = true; };
template<class X> struct IsAlgebra< GradedAlgebra<X> > { static const bool value = true; };
template<class X> struct IsAlgebra< NormedAlgebra<X> > { static const bool value = true; };
template<class X> struct IsAlgebra< SymbolicAlgebra<X> > { static const bool value = true; };
//template<class X> struct IsAlgebra< Differential<X> > { static const bool value = true; };
template<class A> using EnableIfAlgebra = EnableIf<IsAlgebra<A>,A>;
template<class A, class X> struct IsAlgebraOverX : public And< IsAlgebra<A>, IsNumeric<X> > { };
template<class A, class X> using EnableIfAlgebraOverX = EnableIf<IsAlgebraOverX<A,X>,A>;

template<class A, class X> class AlgebraMixin
    : public virtual AlgebraInterface<X>
{
    typedef X NumericType;
  public:
    virtual AlgebraInterface<X>* _create() const { return new A(static_cast<const A&>(*this).A::create()); }
    virtual AlgebraInterface<X>* _clone() const { return new A(static_cast<const A&>(*this)); }
    virtual void _iadd(const X& c) { static_cast<A*>(this)->A::iadd(c); }
    virtual void _imul(const X& c) { static_cast<A*>(this)->A::imul(c); }
    virtual void _isma(const X& c, const AlgebraInterface<X>& x) {
        static_cast<A*>(this)->A::isma(c,dynamic_cast<const A&>(x)); }
    virtual void _ifma(const AlgebraInterface<X>& x1, const AlgebraInterface<X>& x2)  {
        static_cast<A*>(this)->A::ifma(dynamic_cast<const A&>(x1),dynamic_cast<const A&>(x2)); }
    virtual std::ostream& write(std::ostream& os) const { os << static_cast<const A&>(*this); return os; }
};

template<class A, class X> class NormedAlgebraMixin
    : public virtual NormedAlgebraInterface<X>
    , public AlgebraMixin<A,X>
{
    virtual NormedAlgebraInterface<X>* _create_ball(ErrorType r) const { return new A(static_cast<const A&>(*this).A::create_ball(r)); }
    virtual NormedAlgebraInterface<X>* _create() const { return new A(static_cast<const A&>(*this).A::create()); }
    virtual NormedAlgebraInterface<X>* _clone() const { return new A(static_cast<const A&>(*this)); }
};

template<class A, class X> class GradedAlgebraMixin
    : public virtual GradedAlgebraInterface<X>
    , public AlgebraMixin<A,X>
{
    virtual GradedAlgebraMixin<A,X>* _create() const { return new A(static_cast<const A&>(*this).A::create()); }
    virtual GradedAlgebraMixin<A,X>* _clone() const { return new A(static_cast<const A&>(*this)); }
    virtual GradedAlgebraMixin<A,X>* _apply(const Series<X>& f) const { return new A(compose(f,static_cast<const A&>(*this))); }
};

template<class A, class X> class SymbolicAlgebraMixin
    : public virtual SymbolicAlgebraInterface<X>
    , public AlgebraMixin<A,X>
{
    virtual SymbolicAlgebraMixin<A,X>* _create() const { return new A(static_cast<const A&>(*this).A::create()); }
    virtual SymbolicAlgebraMixin<A,X>* _clone() const { return new A(static_cast<const A&>(*this)); }
    virtual SymbolicAlgebraMixin<A,X>* _apply(OperatorCode op) { return new A(op,static_cast<const A&>(*this)); }
};

template<class A> inline EnableIfAlgebra<A> operator+(const A& a) { A r=a; return r; }
template<class A> inline EnableIfAlgebra<A> operator-(const A& a) { A r=a; r.imul(-1); return r; }
template<class A> inline EnableIfAlgebra<A> operator+(const A& a1, const A& a2) { A r=a1; r.isma(+1,a2); return r; }
template<class A> inline EnableIfAlgebra<A> operator-(const A& a1, const A& a2) { A r=a1; r.isma(-1,a2); return r; }
template<class A> inline EnableIfAlgebra<A> operator*(const A& a1, const A& a2) { A r=a1.create(); r.ifma(a1,a2); return r; }

template<class A> inline EnableIfAlgebra<A>& operator+=(A& a1, const A& a2) { a1.isma(+1,a2); return a1; }
template<class A> inline EnableIfAlgebra<A>& operator-=(A& a1, const A& a2) { a1.isma(-1,a2); return a1; }
template<class A> inline EnableIfAlgebra<A>& operator*=(A& a1, const A& a2) { A r = a1.create(); r.ifma(a1,a2); a1=r; return a1; }

template<class A, class X> inline EnableIfAlgebraOverX<A,X>&
operator+=(A& a, const X& c) { a.iadd(numeric_cast<typename A::NumericType>(c)); return a; }
template<class A, class X> EnableIfAlgebraOverX<A,X>&
operator-=(A& a, const X& c) { a.iadd(neg(numeric_cast<typename A::NumericType>(c))); return a; }
template<class A, class X> EnableIfAlgebraOverX<A,X>&
operator*=(A& a, const X& c) { a.imul(numeric_cast<typename A::NumericType>(c)); return a; }
template<class A, class X> EnableIfAlgebraOverX<A,X>&
operator/=(A& a, const X& c) { a.imul(rec(numeric_cast<typename A::NumericType>(c))); return a; }

template<class A, class X> inline EnableIfAlgebraOverX<A,X> operator+(const A& a, const X& c) {
    A r=a; r+=c; return r; }
template<class A, class X> inline EnableIfAlgebraOverX<A,X> operator-(const A& a, const X& c) {
    A r=a; r+=neg(numeric_cast<typename A::NumericType>(c)); return r; }
template<class A, class X> inline EnableIfAlgebraOverX<A,X> operator*(const A& a, const X& c) {
    A r=a; r*=c; return r; }
template<class A, class X> inline EnableIfAlgebraOverX<A,X> operator/(const A& a, const X& c) {
    A r=a; r*=rec(numeric_cast<typename A::NumericType>(c)); return r; }
template<class A, class X> inline EnableIfAlgebraOverX<A,X> operator+(const X& c, const A& a) {
    A r=a; r+=c; return r; }
template<class A, class X> inline EnableIfAlgebraOverX<A,X> operator-(const X& c, const A& a) {
    A r=neg(a); r+=c; return r; }
template<class A, class X> inline EnableIfAlgebraOverX<A,X> operator*(const X& c, const A& a) {
    A r=a; r*=c; return r; }

template<class A> inline EnableIfAlgebra<A> add(const A& a1, const A& a2) { return a1+a2; }
template<class A> inline EnableIfAlgebra<A> mul(const A& a1, const A& a2) { return a1*a2; }
template<class A> inline EnableIfAlgebra<A> neg(const A& a) { return -a; }
template<class A> inline EnableIfAlgebra<A> sqr(const A& a) { return a*a; }

template<class A>  inline EnableIfAlgebra<A> pow(const A& x, uint m) {
    A s=x; A r=s.create(); r+=1; while(m) { if(m%2) { r*=s; } s=sqr(s); m/=2; } return r; }

template<class A> inline EnableIfAlgebra<A> operator/(const A& a1, const A& a2) { return a1*rec(a2); }
template<class A, class X> inline EnableIfAlgebraOverX<A,X> operator/(const X& c, const A& a) { return c*rec(a); }

template<class A> struct IsNormedAlgebra { static const bool value = false; };
template<class X> struct IsNormedAlgebra< NormedAlgebra<X> > { static const bool value = true; };
template<class A> using EnableIfNormedAlgebra = EnableIf<IsNormedAlgebra<A>,A>;

template<class A> EnableIfNormedAlgebra<A> rec(const A& a);
template<class A> EnableIfNormedAlgebra<A> pow(const A& a, int n) { return n<0 ? rec(pow(a,uint(-n))) : pow(a,uint(n)); }
template<class A> EnableIfNormedAlgebra<A> sqrt(const A& a);
template<class A> EnableIfNormedAlgebra<A> exp(const A& a);
template<class A> EnableIfNormedAlgebra<A> log(const A& a);
template<class A> EnableIfNormedAlgebra<A> sin(const A& a);
template<class A> EnableIfNormedAlgebra<A> cos(const A& a);
template<class A> EnableIfNormedAlgebra<A> tan(const A& a);

template<class A> EnableIfNormedAlgebra<A> asin(const A& a);
template<class A> EnableIfNormedAlgebra<A> acos(const A& a);
template<class A> EnableIfNormedAlgebra<A> atan(const A& a);

template<class A> struct IsGradedAlgebra { static const bool value = false; };
template<class X> struct IsGradedAlgebra< GradedAlgebra<X> > { static const bool value = true; };
//template<class X> struct IsGradedAlgebra< Differential<X> > { static const bool value = true; };
template<class A> using EnableIfGradedAlgebra = EnableIf<IsGradedAlgebra<A>,A>;

template<class A> EnableIfGradedAlgebra<A> compose(const Series<typename A::NumericType>& x, const A& y);

template<class A> EnableIfGradedAlgebra<A> sqrt(const A& a) { return compose(Series<typename A::NumericType>::sqrt(a.dgree(),a.value()),a); }
template<class A> EnableIfGradedAlgebra<A> exp(const A& a) { return compose(Series<typename A::NumericType>::exp(a.dgree(),a.value()),a); }
template<class A> EnableIfGradedAlgebra<A> log(const A& a) { return compose(Series<typename A::NumericType>::log(a.dgree(),a.value()),a); }
template<class A> EnableIfGradedAlgebra<A> sin(const A& a) { return compose(Series<typename A::NumericType>::sin(a.dgree(),a.value()),a); }
template<class A> EnableIfGradedAlgebra<A> cos(const A& a) { return compose(Series<typename A::NumericType>::cos(a.dgree(),a.value()),a); }
template<class A> EnableIfGradedAlgebra<A> tan(const A& a) { return compose(Series<typename A::NumericType>::tan(a.dgree(),a.value()),a); }

template<class A> struct IsSymbolicAlgebra { static const bool value = false; };
template<class X> struct IsSymbolicAlgebra< SymbolicAlgebra<X> > { static const bool value = true; };
template<class A> using EnableIfSymbolicAlgebra = EnableIf<IsSymbolicAlgebra<A>,A>;

template<class A> EnableIfSymbolicAlgebra<A> sqrt(const A& a) { return A(SQRT,a); }
template<class A> EnableIfSymbolicAlgebra<A> exp(const A& a) { return A(EXP,a); }
template<class A> EnableIfSymbolicAlgebra<A> log(const A& a) { return A(LOG,a); }
template<class A> EnableIfSymbolicAlgebra<A> sin(const A& a) { return A(SIN,a); }
template<class A> EnableIfSymbolicAlgebra<A> cos(const A& a) { return A(COS,a); }
template<class A> EnableIfSymbolicAlgebra<A> tan(const A& a) { return A(TAN,a); }


} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_MIXIN_H */
