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

template<class A, class X=Void> class Operators;

template<class X> class Operators<X> {
  public:
    static X add(X const& x1, X const& x2);
    static X sub(X const& x1, X const& x2);
    static X mul(X const& x1, X const& x2);
    static X div(X const& x1, X const& x2);
    static X pow(X const& x, Int n);
    static X pos(X const& x);
    static X sqr(X const& x);
    static X neg(X const& x);
    static X rec(X const& x);
    static X sqrt(X const& x);
    static X exp(X const& x);
    static X log(X const& x);
    static X sin(X const& x);
    static X cos(X const& x);
    static X tan(X const& x);
    static X atan(X const& x);
    static X max(X const& x1, X const& x2);
    static X min(X const& x1, X const& x2);
    static X abs(X const& x);
};


template<class X> struct DeclareOperations { };


//! \ingroup NumericAlgebraSubModule
//! \brief Declare operators corresponding to \a X being a (mathematical) ring.
//! \brief \sa DeclareOrderedOperations, DeclareComparisons
//! \details Inheriting from this class declares operators which are found by argument-dependent lookup.
template<class X> struct DeclareRingOperations {
    //@{
    //! \name Named operations.

    //! \brief Positive \a +x. Usually constructs a direct copy.
    friend X pos(X);
    //! \brief Negative \a -x. Usually constructs a direct copy with reversed sign.
    friend X neg(X);
    //! \brief Square \a x<sup>2</sup>. May dispatch to <code>mul(x,x)</code>, but,
    //! especially in the case of interval arithmetic, may be provided separately for more accuracy.
    friend X sqr(X);

    //! \brief Sum \a x1+x2.
    friend X add(X x1, X x2);
    //! \brief Difference \a x1-x2. May be implemented in terms of add() and neg() as <code>add(x1,neg(x2))</code>.
    friend X sub(X x1, X x2);
    //! \brief Product \a x1*x2.
    friend X mul(X x1, X x2);

    //! \brief Power \a x<sup>m</sup>. May be implemented in terms of mul() and/or sqr().
    friend X pow(X x, Nat m);

    //@}


    //@{
    //! \name  Standard overloadable operators. */

    //! \brief Positive \a +x. Usually dispatches to <code>pos(x)</code>
    friend X operator+(X x);
    //! \brief Negative \a -x. Usually dispatches to <code>neg(x)</code>
    friend X operator-(X x);
    //! \brief Sum \a x1+x2. Usually dispatches to <code>sum(x1,x2)</code>
    friend X operator+(X x1, X x2);
    //! \brief Difference \a x1-x2. Usually dispatches to <code>sub(x1,x2)</code>
    friend X operator-(X x1, X x2);
    //! \brief Product \a x1*x2. Usually dispatches to <code>mul(x1,x2)</code>
    friend X operator*(X x1, X x2);

    //@}
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare operators corresponding to \a X being a (mathematical) field.
//! Includes all operations supported in DeclareRingOperations<X>, and those listed below.
template<class X> struct DeclareFieldOperations
    : DeclareRingOperations<X>
{
    //! \brief Reciprocal \a 1/x.
    friend X rec(X x);
    //! \brief Quotient \a x1/x2. May be implemented in terms of mul() and rec() as <code>mul(x1,rec(x2))</code>.
    friend X div(X x1, X x2);
    //! \brief Power \a x<sup>n</sup>. May be implemented using rec() and pow(X,Nat) as  <code> n>=0 ? pow(x,Nat(n)) : rec(pow(x,Nat(-n)))</code>
    friend X pow(X x, Int n);

    //! \brief Quotient \a x1*x2. Usually dispatches to <code>div(x1,x2)</code>
    friend X operator/(X, X);
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare operators corresponding to \a X being a mutable field.
template<class X> struct DeclareInplaceFieldOperations
    : DeclareFieldOperations<X>
{
    //! \brief Inplace addition. May be implemented as \c x1=x1+x2.
    friend X& operator+=(X& x1, X x2);
    //! \brief Inplace subtraction. May be implemented as \c x1=x1-x2.
    friend X& operator-=(X& x1, X x2);
    //! \brief Inplace multiplication. May be implemented as \c x1=x1*x2.
    friend X& operator*=(X& x1, X x2);
    //! \brief Inplace division. May be implemented as \c x1=x1/x2.
    friend X& operator/=(X& x1, X x2);
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare operators corresponding to \a X being an ordered mutable field.
//! \sa DeclareRingOperations, DeclareFieldOperations, DeclareComparisons
template<class X> class DeclareOrderedOperations
{
    //! \brief The absolute value of \a x.
    friend X abs(X x);
    //! \brief The mimimum of \a x1 and \a x2.
    friend X min(X x1, X x2);
    //! \brief The maximum of \a x1 and \a x2.
    friend X max(X x1, X x2);
    //friend U mag(X);
    //friend L mig(X);
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare elementary algebraic and transcendental operations.
template<class X> class DeclareTranscendentalOperations
{
    //! \brief The square root of \a x, √\a x. Requires \c x>=0.
    friend X sqrt(X);
    //! \brief The natural exponent of \a x, \em e<sup>x</sup>.
    friend X exp(X);
    //! \brief The natural logarithm of \a x. Requires \c x>=0.
    friend X log(X);
    //! \brief The sine of \a x.
    friend X sin(X);
    //! \brief The cosine of \a x.
    friend X cos(X);
    //! \brief The tangent of \a x, sin(\a x)/cos(\a x) \f$.
    friend X tan(X);
    //! \brief The arc-tangent of \a x.
    friend X atan(X);
};

template<class X> class DeclareAnalyticOperations
    : DeclareTranscendentalOperations<X>, DeclareFieldOperations<X>
{
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare comparison operations.
//! \sa DeclareRingOperations, DeclareComparisons
template<class X1, class X2, class R> struct DeclareComparisons {
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
    friend R operator+(X x1, X x2) { return R(x1)+R(x2); }
    friend R operator-(X x1, X x2) { return R(x1)-R(x2); }
    friend R operator*(X x1, X x2) { return R(x1)*R(x2); }
    friend R operator/(X x1, X x2) { return R(x1)/R(x2); }
};

template<class X1, class X2, class R> struct ProvideDerivedComparisons : DeclareComparisons<X1,X2,R> {
    friend R operator==(X1 x1, X2 x2) { return eq(x1,x2); }
    friend R operator< (X1 x1, X2 x2) { return lt(x1,x2); }
    friend R operator!=(X1 x1, X2 x2) { return !(x1==x2); }
    friend R operator> (X1 x1, X2 x2) { return  (x2< x1); }
    friend R operator<=(X1 x1, X2 x2) { return !(x2< x1); }
    friend R operator>=(X1 x1, X2 x2) { return !(x1< x2); }
};

template<class X> struct ProvideFieldOperators : DeclareFieldOperations<X> {
    friend X operator+(X x) { return pos(x); }
    friend X operator-(X x) { return neg(x); }
    friend X operator+(X x1, X x2) { return add(x1,x2); }
    friend X operator-(X x1, X x2) { return sub(x1,x2); }
    friend X operator*(X x1, X x2) { return mul(x1,x2); }
    friend X operator/(X x1, X x2) { return div(x1,x2); }
    friend X& operator+=(X& x1, X x2) { return x1=add(x1,x2); }
    friend X& operator-=(X& x1, X x2) { return x1=sub(x1,x2); }
    friend X& operator*=(X& x1, X x2) { return x1=mul(x1,x2); }
    friend X& operator/=(X& x1, X x2) { return x1=div(x1,x2); }
};

template<class X, class R=LogicalType<Paradigm<X>>> struct DeclareNumericOperations
    : public DeclareAnalyticOperations<X>, DeclareOrderedOperations<X>, DeclareComparisons<X,X,R>, DeclareInplaceFieldOperations<X>
{
};

template<class X> struct DispatchNumericOperations : ProvideFieldOperators<X>
{
    friend X add(X x1, X x2) { return Operators<X>::add(x1,x2); }
    friend X sub(X x1, X x2) { return Operators<X>::sub(x1,x2); }
    friend X mul(X x1, X x2) { return Operators<X>::mul(x1,x2); }
    friend X div(X x1, X x2) { return Operators<X>::div(x1,x2); }
    friend X pow(X x, Nat m) { return Operators<X>::pow(x,m); }
    friend X pow(X x, Int n) { return Operators<X>::pow(x,n); }
    friend X pos(X x) { return Operators<X>::pos(x); }
    friend X neg(X x) { return Operators<X>::neg(x); }
    friend X sqr(X x) { return Operators<X>::sqr(x); }
    friend X rec(X x) { return Operators<X>::rec(x); }
    friend X sqrt(X x) { return Operators<X>::sqrt(x); }
    friend X exp(X x) { return Operators<X>::exp(x); }
    friend X log(X x) { return Operators<X>::log(x); }
    friend X sin(X x) { return Operators<X>::sin(x); }
    friend X cos(X x) { return Operators<X>::cos(x); }
    friend X tan(X x) { return Operators<X>::tan(x); }
    friend X atan(X x) { return Operators<X>::atan(x); }
    friend X max(X x1, X x2) { return Operators<X>::max(x1,x2); }
    friend X min(X x1, X x2) { return Operators<X>::min(x1,x2); }
    friend X abs(X x) { return Operators<X>::abs(x); }
};

template<class X> struct DispatchNumericComparisons
{
    //typedef decltype(Operators<X>::lt(declval<X>(),declval<X>())) R;
    //friend R eq(X x1, X x2) { return eq(x1,x2); }
    //friend R lt(X x1, X x2) { return eq(x1,x2); }
    typedef decltype(lt(declval<X>(),declval<X>())) R;
    friend R operator==(X x1, X x2) { return eq(x1,x2); }
    friend R operator< (X x1, X x2) { return lt(x1,x2); }
    friend R operator!=(X x1, X x2) { return !(x1==x2); }
    friend R operator> (X x1, X x2) { return  (x2< x1); }
    friend R operator<=(X x1, X x2) { return !(x2< x1); }
    friend R operator>=(X x1, X x2) { return !(x1< x2); }
};

} // namespace Ariadne

#endif
