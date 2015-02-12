/***************************************************************************
 *            numeric/number.h
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

/*! \file numeric/number.h
 *  \brief Generic numbers
 */



#ifndef ARIADNE_NUMBER_H
#define ARIADNE_NUMBER_H

// Use friend declarations to provide number operators
#define ARIADNE_FRIEND_NUMBER
// Use template functions to provide number operators
// #define ARIADNE_TEMPLATE_NUMBER

#include "utility/handle.h"
#include "numeric/paradigm.h"
#include "utility/prototype.h"

#include "number.decl.h"
#include "float.decl.h"
#include "logical.decl.h"

#include "number_interface.h"

#include "mixins.h"
#include "integer.h"
#include "rational.h"
#include "real.h"

#include "number_interface.h"

namespace Ariadne {

/************ Number *********************************************************/

template<class X> struct IsNumber;

class NumberInterface;
template<class X> class NumberWrapper;

template<class P> class Number;
template<class P> struct IsNumber<Number<P>> : True { };

template<class N1, class N2, EnableIf<And<IsNumber<N1>,IsNumber<N2>>> =dummy> inline N1& operator+=(N1& n1, const N2& n2) { n1=n1+n2; return n1; }
template<class N1, class N2, EnableIf<And<IsNumber<N1>,IsNumber<N2>>> =dummy> inline N1& operator-=(N1& n1, const N2& n2) { n1=n1-n2; return n1; }
template<class N1, class N2, EnableIf<And<IsNumber<N1>,IsNumber<N2>>> =dummy> inline N1& operator*=(N1& n1, const N2& n2) { n1=n1*n2; return n1; }
template<class N1, class N2, EnableIf<And<IsNumber<N1>,IsNumber<N2>>> =dummy> inline N1& operator/=(N1& n1, const N2& n2) { n1=n1/n2; return n1; }


/*
struct DefineBuiltinFloatOperators {
    using ApN = Number<Approximate>;

    template<class X, class D, EnableIf<IsFloatingPoint<D>> =dummy> friend auto operator+(X x, D d) -> decltype(add(x,ApN(d))) { return add(x,ApN(d)); }
    template<class X, class D, EnableIf<IsFloatingPoint<D>> =dummy> friend auto operator+(D d, X x) -> decltype(add(ApN(d),x)) { return add(ApN(d),x); }

    template<class X, class D, EnableIf<IsFloatingPoint<D>> =dummy> friend auto operator-(X x, D d) -> decltype(add(x,ApN(d))) { return sub(x,ApN(d)); }
    template<class X, class D, EnableIf<IsFloatingPoint<D>> =dummy> friend auto operator-(D d, X x) -> decltype(add(ApN(d),x)) { return sub(ApN(d),x); }

    template<class X, class D, EnableIf<IsFloatingPoint<D>> =dummy> friend auto operator*(X x, D d) -> decltype(add(x,ApN(d))) { return mul(x,ApN(d)); }
    template<class X, class D, EnableIf<IsFloatingPoint<D>> =dummy> friend auto operator*(D d, X x) -> decltype(add(ApN(d),x)) { return mul(ApN(d),x); }

    template<class X, class D, EnableIf<IsFloatingPoint<D>> =dummy> friend auto operator/(X x, D d) -> decltype(add(x,ApN(d))) { return div(x,ApN(d)); }
    template<class X, class D, EnableIf<IsFloatingPoint<D>> =dummy> friend auto operator/(D d, X x) -> decltype(add(ApN(d),x)) { return div(ApN(d),x); }
};
*/

#ifdef ARIADNE_FRIEND_NUMBER

//template<class P1, class P2> struct DisableIfWeaker { typedef typename std::enable_if<not std::is_convertible<P2,P1>::value,Dummy>::type Type; };
template<class P1, class P2> using DisableIfWeaker = DisableIf<IsWeaker<P1,P2>>;

//! \ingroup NumericModule
//! \brief Generic numbers with computational paradigm \a P, which may be %Effective, %Validated, %Upper, %Lower or %Approximate.
// Number
template<class P> class Number
{
    static_assert(IsParadigm<P>::value,"P must be a paradigm");
    template<class PP> friend class Number;
    template<class X> using IsGettableAs = And<IsNumber<X>,IsWeaker<typename X::Paradigm,P>,Not<IsSame<typename X::Paradigm,Exact>>>;
  private:
    typedef Opposite<P> NP;
    typedef Weaker<P,NP> SP;
    typedef Widen<P> WP;
    //friend class DispatchGenericField<Number<P>>;
  private:
    Handle<NumberInterface> _handle;
    NumberInterface const* pointer() const { return _handle.pointer(); }
    NumberInterface const& ref() const { return _handle.ref(); }
  public:
    explicit Number(NumberInterface* p) : _handle(p) { }
  private:
    explicit Number(Handle<NumberInterface> h) : _handle(h) { }
    Handle<NumberInterface> handle() const { return this->_handle; }
  public:
    typedef P Paradigm;

    Number() : Number(Integer(0)) { }

    // Construct from a Number of a weaker paradigm
    template<class PP, EnableIf<IsWeaker<P,PP>> = dummy> Number(const Number<PP>& y) : Number<P>(y.handle()) { }

    // Disable construction from a Number of a non-weaker paradigm
   // template<class PP, DisableIfWeaker<P,PP> =dummy> Number(const Number<PP>& y) = delete;

    // Construct from a builtin integer
    template<class N, EnableIf<IsIntegral<N>> =dummy> Number(const N& n) : Number<P>(Integer(n)) { }
    // Construct from a builtin floating-point number
    template<class X, EnableIf<And<IsSame<P,Approximate>,IsFloatingPoint<X>>> =dummy> Number(const X& x) : Number<P>(Float64Template<P>(x)) { }

    // Construct from a type which is convertible to another Number type.
    template<class X, EnableIf<IsWeaker<P,ParadigmTag<X>>> =dummy, EnableIf<IsConvertible<X,Number<ParadigmTag<X>>>> = dummy>
        Number<P>(X const & x) : Number<P>(x.operator Number<ParadigmTag<X>>()) { }

    //! \brief Get the value of the number as a double-precision floating-point type
    template<class WP, EnableIf<IsWeaker<WP,P>> =dummy>
    Float<WP,Precision64> get() const { return pointer()->_get(WP()); }
    //! \brief Get the value of the number as a double-precision floating-point type
    template<class WP, EnableIf<IsWeaker<WP,P>> =dummy>
    Float<WP,Precision64> get(WP par, Precision64 const& prec) const { return pointer()->_get(WP()); }
    //! \brief Get the value of the number as a multiple-precision floating-point type
    template<class WP, EnableIf<IsWeaker<WP,P>> =dummy>
    Float<WP,PrecisionMP> get(WP par, PrecisionMP const& prec) const { return pointer()->_get(WP(),prec); }


    friend Number<P> operator+(Number<P> y) { return pos(y); }
    friend Number<NP> operator-(Number<P> y) { return neg(y); }
    friend Number<P> operator+(Number<P> y1, Number<P> y2) { return add(y1,y2); }
    friend Number<P> operator-(Number<P> y1, Number<NP> y2) { return sub(y1,y2); }
    friend Number<P> operator*(Number<P> y1, Number<P> y2) { return mul(y1,y2); }
    friend Number<P> operator/(Number<P> y1, Number<NP> y2) { return div(y1,y2); }
    friend Number<P>& operator+=(Number<P>& y1, Number<P> y2) { return y1=y1+y2; }
    friend Number<P>& operator-=(Number<P>& y1, Number<NP> y2) { return y1=y1-y2; }
    friend Number<P>& operator*=(Number<P>& y1, Number<P> y2) { return y1=y1*y2; }
    friend Number<P>& operator/=(Number<P>& y1, Number<NP> y2) { return y1=y1/y2; }

    friend Number<NP> operator-(Number<NP> y1, Number<P> y2);
    friend Number<NP> operator/(Number<NP> y1, Number<P> y2);

    friend Number<NP> pos(Number<P> y) { return Number<NP>(y.ref()._pos()); }
    friend Number<NP> neg(Number<P> y) { return Number<NP>(y.ref()._neg()); }
    friend Number<NP> sqr(Number<P> y) { return Number<NP>(y.ref()._sqr()); }
    friend Number<NP> rec(Number<P> y) { return Number<NP>(y.ref()._rec()); }
    friend Number<P> add(Number<P> y1, Number<P> y2) { return Number<P>(y1.ref()._add(y2.ref())); }
    friend Number<P> sub(Number<P> y1, Number<NP> y2) { return Number<P>(y1.ref()._sub(y2.ref())); }
    friend Number<P> mul(Number<P> y1, Number<P> y2) { return Number<P>(y1.ref()._mul(y2.ref())); }
    friend Number<P> div(Number<P> y1, Number<NP> y2) { return Number<P>(y1.ref()._div(y2.ref())); }

    friend Number<P> sqrt(Number<P> y) { return Number<P>(y.ref()._sqrt()); }
    friend Number<P> exp(Number<P> y) { return Number<P>(y.ref()._exp()); }
    friend Number<P> log(Number<P> y) { return Number<P>(y.ref()._log()); }
    friend Number<SP> sin(Number<P> y) { return Number<SP>(y.ref()._sin()); }
    friend Number<SP> cos(Number<P> y) { return Number<SP>(y.ref()._cos()); }
    friend Number<P> tan(Number<P> y) { return Number<P>(y.ref()._tan()); }
    friend Number<P> atan(Number<P> y) { return Number<P>(y.ref()._atan()); }

    friend Number<P> pow(Number<P> y, Nat m) { return Number<P>(y.ref()._pow(m)); }
    friend Number<SP> pow(Number<P> y, Int n) { return Number<SP>(y.ref()._pow(n)); }

    friend Number<P> abs(Number<P> y) { return Number<P>(y.ref()._abs()); }
    friend Number<P> max(Number<P> y1, Number<P> y2) { return Number<P>(y1.ref()._min(y2.ref())); }
    friend Number<P> min(Number<P> y1, Number<P> y2) { return Number<P>(y1.ref()._max(y2.ref())); }

    friend Logical<Equality<P>> operator==(Number<P> y1, Number<NP> y2) {
        return Logical<Equality<P>>(y1.ref()._equals(y2.ref())); }
    friend Logical<LessThan<P>> operator< (Number<P> y1, Number<Negated<P>> y2) {
        return Logical<LessThan<P>>(y1.ref()._less(y2.ref())); }
    friend Logical<LessThan<Negated<P>>> operator> (Number<P> y1, Number<Negated<P>> y2) { return (y2<y1); }
    friend Logical<Negated<Equality<P>>> operator!=(Number<P> y1, Number<Negated<P>> y2) { return !(y1==y2); }
    friend Logical<LessThan<P>> operator<=(Number<P> y1, Number<Negated<P>> y2) { return !(y1>y2); }
    friend Logical<LessThan<Negated<P>>> operator>=(Number<P> y1, Number<Negated<P>> y2) { return !(y2<y1); }

    friend OutputStream& operator<<(OutputStream& os, Number<P> y) { return os << y.ref(); }

    friend Number<SP> operator+(Number<SP>, Number<SP>);
    friend Number<SP> operator-(Number<SP>, Number<SP>);
    friend Number<SP> operator*(Number<SP>, Number<SP>);
    friend Number<SP> operator/(Number<SP>, Number<SP>);
};

template<class R> struct IsConcreteNumber : IsConvertible<R,Real> { };

template<class R, class P, EnableIf<IsConcreteNumber<R>> =dummy> auto operator+(R const& r1, Number<P> const& y2) -> decltype(Number<Paradigm<R>>(r1)+y2) {
    return Number<Paradigm<R>>(r1)+y2; }
template<class R, class P, EnableIf<IsConcreteNumber<R>> =dummy> auto operator+(Number<P> const& y1, R const& r2) -> decltype(y1+Number<Paradigm<R>>(r2)) {
    return y1+Number<Paradigm<R>>(r2); }

//! \brief Get the value of the number represented by \a X with the same precision paramters as \a p.
template<class X, class P, EnableIf<IsWeaker<ParadigmTag<X>,P>> =dummy>
inline X extract(Number<P> y);

#endif

#ifdef ARIADNE_TEMPLATE_NUMBER

//! \ingroup NumericModule
//! \brief Generic numbers with computational paradigm \a P, which may be %Effective, %Validated, %Upper, %Lower or %Approximate.
// Number
template<class P> class Number
    : public Handle<NumberInterface>
{
    static_assert(IsConvertible<P,Approximate>::value,"P must be a paradigm");
    template<class PP> friend class Number;
    template<class X> using IsGettableAs = And<IsNumber<X>,IsWeaker<typename X::Paradigm,P>,Not<IsSame<typename X::Paradigm,Exact>>>;
  public:
    explicit Number(Handle<NumberInterface> h) : Handle<NumberInterface>(h) { }
    Handle<NumberInterface> handle() const { return *this; }
  public:
    typedef P Paradigm;

    Number();

    template<class PP, EnableIf<IsWeaker<P,PP>> = dummy> Number(const Number<PP>& n) : Number<P>(n.handle()) { }
    template<class N, EnableIf<IsIntegral<N>> =dummy> Number(const N& n) : Number<P>(Integer(n)) { }

    // Construct from raw double
    template<class X, EnableIf<And<IsSame<P,Approximate>,IsFloatingPoint<X>>> =dummy> explicit Number(const X& n);
};

// Unary and binary operators
template<class P> Number<P> inline operator+(Number<P> y) {
    return pos(y); }
template<class P> Number<P> inline operator-(Number<P> y) {
    return neg(y); }

template<class P1, class P2> inline Number<Weaker<P1,P2>> operator+(Number<P1> y1, Number<P2> y2) {
    typedef Weaker<P1,P2> P0; return add(Number<P0>(y1),Number<P0>(y2)); }
template<class P1, class P2> inline Number<Weaker<P1,P2>> operator-(Number<P1> y1, Number<P2> y2) {
    typedef Weaker<P1,P2> P0; return sub(Number<P0>(y1),Number<P0>(y2)); }
template<class P1, class P2> inline Number<Weaker<P1,P2>> operator*(Number<P1> y1, Number<P2> y2) {
    typedef Weaker<P1,P2> P0; return mul(Number<P0>(y1),Number<P0>(y2)); }
template<class P1, class P2> inline Number<Weaker<P1,P2>> operator/(Number<P1> y1, Number<P2> y2) {
    typedef Weaker<P1,P2> P0; return div(Number<P0>(y1),Number<P0>(y2)); }

// Overloads
template<class P1, class P2, EnableIf<And<IsSame<P1,Upper>,IsSame<P2,Lower>>> =dummy>
    inline Number<Approximate> operator+(Number<P1> y1, Number<P2> y2) { return Number<Approximate>(y1)+Number<Approximate>(y2); }
template<class P1, class P2, EnableIf<And<IsSame<P1,Lower>,IsSame<P2,Upper>>> =dummy>
    inline Number<Approximate> operator+(Number<P1> y1, Number<P2> y2) { return Number<Approximate>(y1)+Number<Approximate>(y2); }
template<class P1, class P2, EnableIf<And<IsSame<P1,Upper>,IsSame<P2,Upper>>> =dummy>
    inline Number<Approximate> operator-(Number<P1> y1, Number<P2> y2) { return Number<Approximate>(y1)-Number<Approximate>(y2); }
template<class P1, class P2, EnableIf<And<IsSame<P1,Lower>,IsSame<P2,Lower>>> =dummy>
    inline Number<Approximate> operator-(Number<P1> y1, Number<P2> y2) { return Number<Approximate>(y1)-Number<Approximate>(y2); }

// Comparison operators
template<class P> inline Logical<Equality<P>> operator==(Number<P> y1, Number<Negated<P>> y2) {
    return Logical<Equality<P>>(y1.ref()._eq(y2.ref())); }
template<class P> inline Logical<LessThan<P>> operator< (Number<P> y1, Number<Negated<P>> y2) {
    return Logical<LessThan<P>>(y1.ref()._lt(y2.ref())); }
template<class P> inline Logical<LessThan<Negated<P>>> operator> (Number<P> y1, Number<Negated<P>> y2) { return (y2<y1); }
template<class P> inline Logical<Negated<Equality<P>>> operator!=(Number<P> y1, Number<Negated<P>> y2) { return !(y1==y2); }
template<class P> inline Logical<LessThan<P>> operator<=(Number<P> y1, Number<Negated<P>> y2) { return !(y1>y2); }
template<class P> inline Logical<LessThan<Negated<P>>> operator>=(Number<P> y1, Number<Negated<P>> y2) { return !(y2<y1); }

template<class P> inline OutputStream& operator<<(OutputStream& os, Number<P> y) {
    return os << y.handle(); }

template<class P, class X, EnableIf<And<IsFloat<X>,IsWeaker<Paradigm<X>,P>>> =dummy>
    inline X extract(Number<P> y, Prototype<X> const& p) { return y.handle().get(p); }

template<class P> inline Number<P> abs(Number<P> y);
template<class P1, class P2> inline Number<Weaker<P1,P2>> max(Number<P1> y1, Number<P2> y2);
template<class P1, class P2> inline Number<Weaker<P1,P2>> min(Number<P1> y1, Number<P2> y2);

template<class P> inline Number<P> add(Number<P> y1, Number<P> y2) { return Number<P>(y1.ref()._add(y2.ref())); }
template<class P> inline Number<P> sub(Number<P> y1, Number<Negated<P>> y2) { return Number<P>(y1.ref()._sub(y2.ref())); }
template<class P> inline Number<P> mul(Number<P> y1, Number<P> y2) { return Number<P>(y1.ref()._mul(y2.ref())); }
template<class P> inline Number<P> div(Number<P> y1, Number<Negated<P>> y2) { return Number<P>(y1.ref()._div(y2.ref())); }
template<class P> inline Number<P> pos(Number<P> y) { return Number<P>(y.ref()._pos()); }
template<class P> inline Number<P> sqr(Number<P> y) { return Number<P>(y.ref()._sqr()); }
template<class P> inline Number<P> neg(Number<P> y) { return Number<P>(y.ref()._neg()); }
template<class P> inline Number<P> rec(Number<P> y) { return Number<P>(y.ref()._rec()); }
template<class P> inline Number<P> sqrt(Number<P> y) { return Number<P>(y.ref()._sqrt()); }
template<class P> inline Number<P> exp(Number<P> y) { return Number<P>(y.ref()._exp()); }
template<class P> inline Number<P> log(Number<P> y) { return Number<P>(y.ref()._log()); }
template<class P> inline Number<P> sin(Number<P> y) { return Number<P>(y.ref()._sin()); }
template<class P> inline Number<P> cos(Number<P> y) { return Number<P>(y.ref()._cos()); }
template<class P> inline Number<P> tan(Number<P> y) { return Number<P>(y.ref()._tan()); }
template<class P> inline Number<P> atan(Number<P> y) { return Number<P>(y.ref()._atan()); }

template<class P> inline Number<P> pow(Number<P> y, Integer n) { return Number<P>(y.ref()._pow(n)); }

template<class P, class N, EnableIf<IsIntegral<N>> =dummy> inline Number<P> pow(Number<P> y1, N n2) {
    return Number<P>(y1.ref()._pow(n2)); }

template<class P> inline Logical<Equality<P>> eq(Number<P> y1, Number<Negated<P>> y2) {
    return Logical<P>(y1.ref()._eq(y2.ref())); }

template<class P> OutputStream& operator<<(OutputStream& os, Number<P> y) { return y._ptr->_write(os); }

#endif /* ARIADNE_FRIEND_NUMBER */

} // namespace Ariadne

#endif
