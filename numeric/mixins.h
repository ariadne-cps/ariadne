/***************************************************************************
 *            numeric/mixins.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file numeric/mixins.h
 *  \brief
 */



#ifndef ARIADNE_NUMERIC_MIXINS_H
#define ARIADNE_NUMERIC_MIXINS_H

namespace Ariadne {

typedef Nat Nat;
typedef Int Int;

template<class X, class NX=X, class R=X> class Operations {
  public:
    static X _nul(X const& x);
    static X _pos(X const& x);
    static X _neg(NX const& x);
    static X _half(X const& x);
    static R _sqr(X const& x);
    static R _rec(NX const& x);
    static R _add(X const& x1, X const& x2);
    static R _sub(X const& x1, NX const& x2);
    static R _mul(X const& x1, X const& x2);
    static R _div(X const& x1, NX const& x2);
    static R _pow(X const& x, Nat m);
    static R _pow(X const& x, Int n);
    static R _sqrt(X const& x);
    static R _exp(X const& x);
    static R _log(X const& x);
    static R _sin(X const& x);
    static R _cos(X const& x);
    static R _tan(X const& x);
    static R _asin(X const& x);
    static R _acos(X const& x);
    static R _atan(X const& x);
    static X _max(X const& x1, X const& x2);
    static X _min(X const& x1, X const& x2);
    static X _abs(X const& x);
};


//! \ingroup NumericAlgebraSubModule
//! \brief Declare operators corresponding to \a X being a (mathematical) ring.
//! \brief \sa DeclareOrderedOperations, DeclareComparisons
//! \details Inheriting from this class declares operators which are found by argument-dependent lookup.
template<class X> struct DeclareRingOperations {
    //@{
    //! \name Named operations.

    //! \brief Positive \a +x. Usually constructs a direct copy.
    friend X pos(X const& x);
    //! \brief Negative \a -x. Usually constructs a direct copy with reversed sign.
    friend X neg(X const& x);
    //! \brief Square \a x<sup>2</sup>. May dispatch to <code>mul(x,x)</code>, but,
    //! especially in the case of interval arithmetic, may be provided separately for more accuracy.
    friend X sqr(X const& x);

    //! \brief Sum \a x1+x2.
    friend X add(X const& x1, X const& x2);
    //! \brief Difference \a x1-x2. May be implemented in terms of add() and neg() as <code>add(x1,neg(x2))</code>.
    friend X sub(X const& x1, X const& x2);
    //! \brief Product \a x1*x2.
    friend X mul(X const& x1, X const& x2);

    //! \brief Power \a x<sup>m</sup>. May be implemented in terms of mul() and/or sqr().
    friend X pow(X const& x, Nat m);

    //@}


    //@{
    //! \name  Standard overloadable operators. */

    //! \brief Positive \a +x. Usually dispatches to <code>pos(x)</code>
    friend X operator+(X const& x);
    //! \brief Negative \a -x. Usually dispatches to <code>neg(x)</code>
    friend X operator-(X const& x);
    //! \brief Sum \a x1+x2. Usually dispatches to <code>sum(x1,x2)</code>
    friend X operator+(X const& x1, X const& x2);
    //! \brief Difference \a x1-x2. Usually dispatches to <code>sub(x1,x2)</code>
    friend X operator-(X const& x1, X const& x2);
    //! \brief Product \a x1*x2. Usually dispatches to <code>mul(x1,x2)</code>
    friend X operator*(X const& x1, X const& x2);

    //@}
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare operators corresponding to \a X being a (mathematical) field.
//! Includes all operations supported in DeclareRingOperations<X>, and those listed below.
template<class X> struct DeclareFieldOperations
    : DeclareRingOperations<X>
{
    //! \brief Reciprocal \a 1/x.
    friend X rec(X const& x);
    //! \brief Quotient \a x1/x2. May be implemented in terms of mul() and rec() as <code>mul(x1,rec(x2))</code>.
    friend X div(X const& x1, X const& x2);
    //! \brief Power \a x<sup>n</sup>. May be implemented using rec() and pow(X,Nat) as  <code> n>=0 ? pow(x,Nat(n)) : rec(pow(x,Nat(-n)))</code>
    friend X pow(X const& x, Int n);

    //! \brief Quotient \a x1*x2. Usually dispatches to <code>div(x1,x2)</code>
    friend X operator/(X const&, X const&);
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare operators corresponding to \a X being a mutable field.
template<class X> struct DeclareInplaceFieldOperations
    : DeclareFieldOperations<X>
{
    //! \brief Inplace addition. May be implemented as \c x1=x1+x2.
    friend X& operator+=(X& x1, X const& x2);
    //! \brief Inplace subtraction. May be implemented as \c x1=x1-x2.
    friend X& operator-=(X& x1, X const& x2);
    //! \brief Inplace multiplication. May be implemented as \c x1=x1*x2.
    friend X& operator*=(X& x1, X const& x2);
    //! \brief Inplace division. May be implemented as \c x1=x1/x2.
    friend X& operator/=(X& x1, X const& x2);
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare operators corresponding to \a X being an ordered mutable field.
//! \sa DeclareRingOperations, DeclareFieldOperations, DeclareComparisons
template<class X> class DeclareOrderedOperations
{
    //! \brief The absolute value of \a x.
    friend X abs(X const& x);
    //! \brief The mimimum of \a x1 and \a x2.
    friend X min(X const& x1, X const& x2);
    //! \brief The maximum of \a x1 and \a x2.
    friend X max(X const& x1, X const& x2);
    //friend U mag(X const& x);
    //friend L mig(X const& x);
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare elementary algebraic and transcendental operations.
template<class X> class DeclareTranscendentalOperations
{
    //! \brief The square root of \a x, √\a x. Requires \c x>=0.
    friend X sqrt(X const& x);
    //! \brief The natural exponent of \a x, \em e<sup>x</sup>.
    friend X exp(X const& x);
    //! \brief The natural logarithm of \a x. Requires \c x>=0.
    friend X log(X const& x);
    //! \brief The sine of \a x.
    friend X sin(X const& x);
    //! \brief The cosine of \a x.
    friend X cos(X const& x);
    //! \brief The tangent of \a x, sin(\a x)/cos(\a x) \f$.
    friend X tan(X const& x);
    //! \brief The arc-tangent of \a x.
    friend X atan(X const& x);
};

template<class X> class DeclareAnalyticOperations
    : DeclareFieldOperations<X>, DeclareTranscendentalOperations<X>
{
};





//! \ingroup NumericAlgebraSubModule
//! \brief Declare comparison operations.
//! \sa DeclareRingOperations, DeclareComparisons
template<class X1, class X2, class R> struct DeclareFieldComparisons {
    //! \brief Tests if \a x1 is equal to \a x2.
    friend R eq(X1 x1, X2 x2);
    //! \brief Tests if \a x1 is less than \a x2.
    friend R lt(X1 x1, X2 x2);
    //! \brief Tests if \a x1 is equal to \a x2. If equality is undecidable, may not return \a true even if values are equal.
    friend R operator==(X1 x1, X2 x2);
    //! \brief Tests if \a x1 is strictly less than \a x2.
    friend R operator< (X1 x1, X2 x2);
    //! \brief Tests if \a x1 is not equal to \a x2.
    friend R operator!=(X1 x1, X2 x2);
    //! \brief Tests if \a x1 is strictly greater than \a x2.
    friend R operator> (X1 x1, X2 x2);
    //! \brief Tests if \a x1 is less than or equal to \a x2.
    friend R operator<=(X1 x1, X2 x2);
    //! \brief Tests if \a x1 is greater than or equal to \a x2.
    friend R operator>=(X1 x1, X2 x2);
};

template<class X1, class X2, class R> struct ProvideConvertedOperations;

template<class X1, class X2, class R> struct ProvideConvertedOperations<X1,X2,Return<R>> {
    friend R operator+(X1 x1, X2 x2) { return R(x1)+R(x2); }
    friend R operator-(X1 x1, X2 x2) { return R(x1)-R(x2); }
    friend R operator*(X1 x1, X2 x2) { return R(x1)*R(x2); }
    friend R operator/(X1 x1, X2 x2) { return R(x1)/R(x2); }
    friend R operator+(X2 x1, X1 x2) { return R(x1)+R(x2); }
    friend R operator-(X2 x1, X1 x2) { return R(x1)-R(x2); }
    friend R operator*(X2 x1, X1 x2) { return R(x1)*R(x2); }
    friend R operator/(X2 x1, X1 x2) { return R(x1)/R(x2); }

    typedef decltype(declval<R>()< declval<R>()) L;
    typedef decltype(declval<R>()==declval<R>()) E;
};

template<class X1, class X2> struct ProvideConvertedOperations<X1,X2,Return<X1>> {
    typedef X1 R;
    friend R operator+(X1 x1, X2 x2) { return x1+R(x2); }
    friend R operator-(X1 x1, X2 x2) { return x1-R(x2); }
    friend R operator*(X1 x1, X2 x2) { return x1*R(x2); }
    friend R operator/(X1 x1, X2 x2) { return x1/R(x2); }
    friend R operator+(X2 x1, X1 x2) { return R(x1)+x2; }
    friend R operator-(X2 x1, X1 x2) { return R(x1)-x2; }
    friend R operator*(X2 x1, X1 x2) { return R(x1)*x2; }
    friend R operator/(X2 x1, X1 x2) { return R(x1)/x2; }

    friend R add(X1 x1, X2 x2) { return x1+R(x2); }
    friend R sub(X1 x1, X2 x2) { return x1-R(x2); }
    friend R mul(X1 x1, X2 x2) { return x1*R(x2); }
    friend R div(X1 x1, X2 x2) { return x1/R(x2); }
    friend R add(X2 x1, X1 x2) { return R(x1)+x2; }
    friend R sub(X2 x1, X1 x2) { return R(x1)-x2; }
    friend R mul(X2 x1, X1 x2) { return R(x1)*x2; }
    friend R div(X2 x1, X1 x2) { return R(x1)/x2; }
};

template<class X, class R> struct ProvideConvertedOperations<X,X,Return<R>> {
    friend R operator+(X const& x1, X const& x2) { return R(x1)+R(x2); }
    friend R operator-(X const& x1, X const& x2) { return R(x1)-R(x2); }
    friend R operator*(X const& x1, X const& x2) { return R(x1)*R(x2); }
    friend R operator/(X const& x1, X const& x2) { return R(x1)/R(x2); }
};

template<class X> struct ProvideFieldOperators : DeclareFieldOperations<X> {
    friend X operator+(X const& x) { return pos(x); }
    friend X operator-(X const& x) { return neg(x); }
    friend X operator+(X const& x1, X const& x2) { return add(x1,x2); }
    friend X operator-(X const& x1, X const& x2) { return sub(x1,x2); }
    friend X operator*(X const& x1, X const& x2) { return mul(x1,x2); }
    friend X operator/(X const& x1, X const& x2) { return div(x1,x2); }
    friend X& operator+=(X& x1, X const& x2) { return x1=add(x1,x2); }
    friend X& operator-=(X& x1, X const& x2) { return x1=sub(x1,x2); }
    friend X& operator*=(X& x1, X const& x2) { return x1=mul(x1,x2); }
    friend X& operator/=(X& x1, X const& x2) { return x1=div(x1,x2); }
};



template<class X, class R=X> struct DeclareInplaceOperators {
};
template<class X> struct DeclareInplaceOperators<X,X> {
    friend X& operator+=(X& x1, X const& x2);
    friend X& operator-=(X& x1, X const& x2);
    friend X& operator*=(X& x1, X const& x2);
    friend X& operator/=(X& x1, X const& x2);
};

template<class X, class R=X> struct DeclareArithmeticOperators
    : DeclareInplaceOperators<X,R>
{
    friend X operator+(X const& x);
    friend X operator-(X const& x);
    friend R operator+(X const& x1, X const& x2);
    friend R operator-(X const& x1, X const& x2);
    friend R operator*(X const& x1, X const& x2);
    friend R operator/(X const& x1, X const& x2);
};

template<class X, class NX, class R=X> struct DeclareDirectedInplaceOperators {
};
template<class X, class NX> struct DeclareDirectedInplaceOperators<X,NX,X> {
    friend X operator+=(X const& x1, X const& x2);
    friend X operator-=(X const& x1, NX const& x2);
};

template<class X, class NX, class R=X, class NR=NX> struct DeclareDirectedArithmeticOperators
    : DeclareDirectedInplaceOperators<X,NX,R>, DeclareDirectedInplaceOperators<NX,X,NR>
{
    friend X operator+(X const& x);
    friend X operator-(NX const& x);
    friend R operator+(X const& x1, X const& x2);
    friend R operator-(X const& x1, NX const& x2);
    friend NX operator+(NX const& x);
    friend NX operator-(X const& x);
    friend NR operator+(NX const& x1, NX const& x2);
    friend NR operator-(NX const& x1, X const& x2);
};

template<class X, class NX, class R=X, class NR=NX> struct DeclareDirectedNumericOperations
    : DeclareDirectedArithmeticOperators<X,NX,R>
{
    friend X nul(X const& x);
    friend X pos(X const& x);
    friend X neg(NX const& x);
    friend X hlf(X const& x);
    friend R add(X const& x1, X const& x2);
    friend R sub(X const& x1, NX const& x2);
    friend R sqrt(X const& x);
    friend R exp(X const& x);
    friend R log(X const& x);
    friend R atan(X const& x);
    friend X max(X const& x1, X const& x2);
    friend X min(X const& x1, X const& x2);
    friend X abs(X const& x);
};

template<class X, class R=X> struct DeclareNumericOperations
    : DeclareArithmeticOperators<X,R>
{
    friend X nul(X const& x);
    friend X pos(X const& x);
    friend X neg(X const& x);
    friend X half(X const& x);
    friend R sqr(X const& x);
    friend R rec(X const& x);
    friend R add(X const& x1, X const& x2);
    friend R sub(X const& x1, X const& x2);
    friend R mul(X const& x1, X const& x2);
    friend R div(X const& x1, X const& x2);
    friend R pow(X const& x, Nat m);
    friend R pow(X const& x, Int n);
    friend R sqrt(X const& x);
    friend R exp(X const& x);
    friend R log(X const& x);
    friend R sin(X const& x);
    friend R cos(X const& x);
    friend R tan(X const& x);
    friend R asin(X const& x);
    friend R acos(X const& x);
    friend R atan(X const& x);
    friend X max(X const& x1, X const& x2);
    friend X min(X const& x1, X const& x2);
};


template<class X, class R=X> struct ProvideInplaceOperators {
};
template<class X> struct ProvideInplaceOperators<X,X> {
    friend X& operator+=(X& x1, X const& x2) { return x1=add(x1,x2); }
    friend X& operator-=(X& x1, X const& x2) { return x1=sub(x1,x2); }
    friend X& operator*=(X& x1, X const& x2) { return x1=mul(x1,x2); }
    friend X& operator/=(X& x1, X const& x2) { return x1=div(x1,x2); }
};

template<class X, class R=X> struct ProvideArithmeticOperators
    : ProvideInplaceOperators<X,R>
{
    friend X operator+(X const& x) { return pos(x); }
    friend X operator-(X const& x) { return neg(x); }
    friend R operator+(X const& x1, X const& x2) { return add(x1,x2); }
    friend R operator-(X const& x1, X const& x2) { return sub(x1,x2); }
    friend R operator*(X const& x1, X const& x2) { return mul(x1,x2); }
    friend R operator/(X const& x1, X const& x2) { return div(x1,x2); }
};

template<class X, class NX> class ProvideDirectedArithmeticOperators
{
    friend X neg(NX const&);
    friend NX neg(X const&);
    friend NX operator-(X const&);
    friend X operator+(X const& x) { return pos(x); }
    friend X operator-(NX const& x) { return neg(x); }
    friend X operator+(X const& x1, X const& x2) { return add(x1,x2); }
    friend X operator-(X const& x1, NX const& x2) { return sub(x1,x2); }
    friend X& operator+=(X& x1, X const& x2) { return x1=add(x1,x2); }
    friend X& operator-=(X& x1, NX const& x2) { return x1=sub(x1,x2); }
};


template<class X, class Y> struct ProvideMixedConcreteGenericOperators {
    friend X operator+(const X& x1, const Y& y2) { return x1+x1.create(y2); }
    friend X operator-(const X& x1, const Y& y2) { return x1-x1.create(y2); }
    friend X operator*(const X& x1, const Y& y2) { return x1*x1.create(y2); }
    friend X operator/(const X& x1, const Y& y2) { return x1/x1.create(y2); }
    friend X operator+(const Y& y1, const X& x2) { return x2.create(y1)+x2; }
    friend X operator-(const Y& y1, const X& x2) { return x2.create(y1)-x2; }
    friend X operator*(const Y& y1, const X& x2) { return x2.create(y1)*x2; }
    friend X operator/(const Y& y1, const X& x2) { return x2.create(y1)/x2; }
};

template<class X, class QX=X, class R=X, class QR=QX> class ProvidePositiveDirectedArithmeticOperators
{
    friend R rec(QX const&);
    friend R operator+(X const& x1, X const& x2) { return add(x1,x2); }
    friend R operator*(X const& x1, X const& x2) { return mul(x1,x2); }
    friend R operator/(X const& x1, QX const& x2) { return div(x1,x2); }
    friend X& operator+=(X& x1, X const& x2) { return x1=add(x1,x2); }
    friend X& operator*=(X& x1, X const& x2) { return x1=mul(x1,x2); }
    friend X& operator/=(X& x1, QX const& x2) { return x1=div(x1,x2); }

    friend QR rec(X const&);
    friend QR operator/(QX const&, X const&);
};

template<class X, class R=X> struct DispatchNumericOperations
    : ProvideArithmeticOperators<X,R>
{
    friend X nul(X const& x) { return Operations<X,X,R>::_nul(x); }
    friend X pos(X const& x) { return Operations<X,X,R>::_pos(x); }
    friend X neg(X const& x) { return Operations<X,X,R>::_neg(x); }
    friend X half(X const& x) { return Operations<X,X,R>::_half(x); }
    friend R sqr(X const& x) { return Operations<X,X,R>::_sqr(x); }
    friend R rec(X const& x) { return Operations<X,X,R>::_rec(x); }
    friend R add(X const& x1, X const& x2) { return Operations<X,X,R>::_add(x1,x2); }
    friend R sub(X const& x1, X const& x2) { return Operations<X,X,R>::_sub(x1,x2); }
    friend R mul(X const& x1, X const& x2) { return Operations<X,X,R>::_mul(x1,x2); }
    friend R div(X const& x1, X const& x2) { return Operations<X,X,R>::_div(x1,x2); }
    friend R pow(X const& x, Nat m) { return Operations<X,X,R>::_pow(x,m); }
    friend R pow(X const& x, Int n) { return Operations<X,X,R>::_pow(x,n); }
    friend R sqrt(X const& x) { return Operations<X,X,R>::_sqrt(x); }
    friend R exp(X const& x) { return Operations<X,X,R>::_exp(x); }
    friend R log(X const& x) { return Operations<X,X,R>::_log(x); }
    friend R sin(X const& x) { return Operations<X,X,R>::_sin(x); }
    friend R cos(X const& x) { return Operations<X,X,R>::_cos(x); }
    friend R tan(X const& x) { return Operations<X,X,R>::_tan(x); }
    friend R asin(X const& x) { return Operations<X,X,R>::_asin(x); }
    friend R acos(X const& x) { return Operations<X,X,R>::_acos(x); }
    friend R atan(X const& x) { return Operations<X,X,R>::_atan(x); }
    friend X max(X const& x1, X const& x2) { return Operations<X,X,R>::_max(x1,x2); }
    friend X min(X const& x1, X const& x2) { return Operations<X,X,R>::_min(x1,x2); }
    friend X abs(X const& x) { return Operations<X,X,R>::_abs(x); }
};

template<class X, class NX, class R=X, class NR=NX> struct DispatchDirectedNumericOperations
    : ProvideDirectedArithmeticOperators<X,NX>
{
    friend X nul(X const& x) { return Operations<X,NX>::_nul(x); }
    friend X pos(X const& x) { return Operations<X,NX>::_pos(x); }
    friend X neg(NX const& x) { return Operations<X,NX>::_neg(x); }
    friend X hlf(X const& x) { return Operations<X,NX>::_hlf(x); }
    friend R add(X const& x1, X const& x2) { return Operations<X,NX>::_add(x1,x2); }
    friend R sub(X const& x1, NX const& x2) { return Operations<X,NX>::_sub(x1,x2); }
    friend R sqrt(X const& x) { return Operations<X,NX>::_sqrt(x); }
    friend R exp(X const& x) { return Operations<X,NX>::_exp(x); }
    friend R log(X const& x) { return Operations<X,NX>::_log(x); }
    friend R atan(X const& x) { return Operations<X,NX>::_atan(x); }
    friend X max(X const& x1, X const& x2) { return Operations<X,NX>::_max(x1,x2); }
    friend X min(X const& x1, X const& x2) { return Operations<X,NX>::_min(x1,x2); }
    friend X abs(X const& x) { return Operations<X,NX>::_abs(x); }

};

template<class X, class QX, class R=X, class QR=QX> struct DispatchPositiveDirectedNumericOperations
    : ProvidePositiveDirectedArithmeticOperators<X,QX,R,QR>
{
    friend X nul(X const& x) { return Operations<X,QX>::_nul(x); }
    friend X hlf(X const& x) { return Operations<X,QX>::_hlf(x); }
    friend R sqr(X const& x) { return Operations<X,QX>::_sqr(x); }
    friend R rec(QX const& x) { return Operations<X,QX>::_rec(x); }
    friend R add(X const& x1, X const& x2) { return Operations<X,QX>::_add(x1,x2); }
    friend R mul(X const& x1, X const& x2) { return Operations<X,QX>::_mul(x1,x2); }
    friend R div(X const& x1, QX const& x2) { return Operations<X,QX>::_div(x1,x2); }
    friend R pow(X const& x, Nat m) { return Operations<X,QX>::_pow(x,m); }
    friend R sqrt(X const& x) { return Operations<X,QX>::_sqrt(x); }
    friend R atan(X const& x) { return Operations<X,QX>::_atan(x); }
    friend X max(X const& x1, X const& x2) { return Operations<X,QX>::_max(x1,x2); }
    friend X min(X const& x1, X const& x2) { return Operations<X,QX>::_min(x1,x2); }
    friend X abs(X const& x) { return Operations<X,QX>::_abs(x); }
};

template<class X> struct DispatchNumericComparisons
{
    //typedef decltype(Operators<X>::lt(declval<X>(),declval<X>())) R;
    //friend R eq(X const& x1, X const& x2) { return eq(x1,x2); }
    //friend R lt(X const& x1, X const& x2) { return eq(x1,x2); }
    typedef decltype(lt(declval<X>(),declval<X>())) R;
    friend R operator==(X const& x1, X const& x2) { return eq(x1,x2); }
    friend R operator< (X const& x1, X const& x2) { return lt(x1,x2); }
    friend R operator!=(X const& x1, X const& x2) { return !(x1==x2); }
    friend R operator> (X const& x1, X const& x2) { return  (x2< x1); }
    friend R operator<=(X const& x1, X const& x2) { return !(x2< x1); }
    friend R operator>=(X const& x1, X const& x2) { return !(x1< x2); }
};

} // namespace Ariadne

#endif
