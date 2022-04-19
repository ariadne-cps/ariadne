/***************************************************************************
 *            numeric/float_value.hpp
 *
 *  Copyright  2008-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file numeric/float_value.hpp
 *  \brief Exact Floating-point representations of real numbers.
 */

#ifndef ARIADNE_FLOAT_VALUE_HPP
#define ARIADNE_FLOAT_VALUE_HPP

#include "logical.decl.hpp"
#include "number.decl.hpp"
#include "float.decl.hpp"

#include "float_traits.hpp"
#include "float_operations.hpp"

#include "integer.hpp"
#include "dyadic.hpp"

#include "floatdp.hpp"
#include "floatmp.hpp"
#include "float_bounds.hpp"

namespace Ariadne {

struct DefaultTag;

static_assert(not GenericNumber<FloatValue<DoublePrecision>>);
static_assert(not GenericNumber<FloatValue<MultiplePrecision>>);


template<class F> struct IsARawFloat { static const bool value = false; };
template<> struct IsARawFloat<FloatDP> { static const bool value = true; };
template<> struct IsARawFloat<FloatMP> { static const bool value = true; };

template<class F> concept ARawFloat = IsARawFloat<F>::value;

template<ARawFloat F> Bounds<F> operator+(F const& x1, F const& x2);
template<ARawFloat F> Bounds<F> operator-(F const& x1, F const& x2);
template<ARawFloat F> Bounds<F> operator*(F const& x1, F const& x2);
template<ARawFloat F> Bounds<F> operator/(F const& x1, F const& x2);

template<ARawFloat F> Integer cast_integer(Value<F> const& x) {
    Dyadic w(x); Integer z=round(w); ARIADNE_ASSERT_MSG(z==w,"Cannot cast non-integral value "<<z<<" to an Integer"); return z;
}

template<ARawFloat F> inline FloatFactory<PrecisionType<F>> factory(Value<F> const& flt) { return FloatFactory<PrecisionType<F>>(flt.precision()); }
template<class PR> inline FloatValue<PR> FloatFactory<PR>::create(Dyadic const& y, ExactTag) { return FloatValue<PR>(y,_pr); }
template<class PR> inline FloatValue<PR> FloatFactory<PR>::create(Integer const& y, ExactTag) { return FloatValue<PR>(y,_pr); }
template<class PR> inline FloatValue<PR> FloatFactory<PR>::create(ExactDouble const& y) { return FloatValue<PR>(y,_pr); }

template<class PR> template<BuiltinSignedIntegral N> inline
FloatValue<PR> FloatFactory<PR>::create(N const& y) { return FloatValue<PR>(y,_pr); }

template<class PR> template<BuiltinUnsignedIntegral M> inline
PositiveFloatValue<PR> FloatFactory<PR>::create(M const& y) { return PositiveFloatValue<PR>(y,_pr); }

template<ARawFloat F> class Positive<Value<F>> : public F {
    using PR = typename Value<F>::PrecisionType;
  public:
    Positive() : Value<F>() { }
    explicit Positive(PR const& pr) : Value<F>(pr) { }
    template<BuiltinUnsignedIntegral M> Positive(M m, PR pr) : Value<F>(m,pr) { }
    Positive(TwoExp const& ex, PR pr) : Value<F>(ex,pr) { }
    explicit Positive(Dyadic const& w, PR pr) : Value<F>(w,pr) { }
    explicit Positive(Value<F> const& x) : Value<F>(x) { }
  public:
    friend Positive<Bounds<F>> operator+(PositiveValue<F> const& v1, PositiveValue<F> const& v2) {
        return cast_positive(static_cast<Value<F>const&>(v1)+static_cast<Value<F>const&>(v2)); }
    friend Positive<Bounds<F>> operator*(PositiveValue<F> const& v1, PositiveValue<F> const& v2) {
        return cast_positive(static_cast<Value<F>const&>(v1)*static_cast<Value<F>const&>(v2)); }

    friend Positive<Value<F>> nul(Positive<Value<F>> const& x) { return PositiveValue<F>(nul(static_cast<F const&>(x))); }
    friend Positive<Value<F>> pos(Positive<Value<F>> const& x) { return PositiveValue<F>(pos(static_cast<F const&>(x))); }
    friend Positive<Value<F>> hlf(Positive<Value<F>> const& x) { return PositiveValue<F>(hlf(static_cast<F const&>(x))); }
    friend Positive<Bounds<F>> pow(Positive<Value<F>> const& x, Nat m) {
        return pow(Positive<Bounds<F>>(x),m); }
    friend Positive<Bounds<F>> pow(Positive<Value<F>> const& x, Int n) {
        return pow(Positive<Bounds<F>>(x),n); }
};

template<class F> inline PositiveValue<F> cast_positive(Value<F> const& x) {
    return PositiveValue<F>(x);
}
template<class F> inline Positive<F> cast_positive(Positive<F> const& x) {
    return x;
}
template<class F> inline F const& cast_unsigned(Positive<F> const& x) {
    return static_cast<F const&>(x);
}


template<ARawFloat F> Positive<Bounds<F>> operator+(Positive<F> const& x1, Positive<F> const& x2) {
    return cast_positive(cast_unsigned(x1)+cast_unsigned(x2)); }
template<ARawFloat F> Positive<Bounds<F>> operator*(Positive<F> const& x1, Positive<F> const& x2) {
    return cast_positive(cast_unsigned(x1)*cast_unsigned(x2)); }
template<ARawFloat F> Positive<Bounds<F>> operator/(Positive<F> const& x1, Positive<F> const& x2) {
    return cast_positive(cast_unsigned(x1)/cast_unsigned(x2)); }
template<ARawFloat F> Positive<F> max(Positive<F> const& x1, Positive<F> const& x2) {
    return cast_positive(max(cast_unsigned(x1),cast_unsigned(x2))); }
template<ARawFloat F> Positive<F> max(Positive<F> const& x1, F const& x2) {
    return cast_positive(max(cast_unsigned(x1),x2)); }
template<ARawFloat F> Positive<F> max(F const& x1, Positive<F> const& x2) {
    return cast_positive(max(x1,cast_unsigned(x2))); }
template<ARawFloat F> Positive<F> min(Positive<F> const& x1, Positive<F> const& x2) {
    return cast_positive(min(cast_unsigned(x1),cast_unsigned(x2))); }

inline Bounds<FloatDP> add(Value<FloatDP> const& x1, Value<FloatDP> const& x2) { return add(Bounds<FloatDP>(x1),Bounds<FloatDP>(x2)); }
inline Bounds<FloatDP> sub(Value<FloatDP> const& x1, Value<FloatDP> const& x2) { return sub(Bounds<FloatDP>(x1),Bounds<FloatDP>(x2)); }
inline Bounds<FloatDP> mul(Value<FloatDP> const& x1, Value<FloatDP> const& x2) { return mul(Bounds<FloatDP>(x1),Bounds<FloatDP>(x2)); }
inline Bounds<FloatDP> div(Value<FloatDP> const& x1, Value<FloatDP> const& x2) { return div(Bounds<FloatDP>(x1),Bounds<FloatDP>(x2)); }
inline Bounds<FloatDP> pow(Value<FloatDP> const& x, Int n) { return pow(Bounds<FloatDP>(x),n); }
inline Bounds<FloatDP> sqr(Value<FloatDP> const& x) { return sqr(Bounds<FloatDP>(x)); }
inline Bounds<FloatDP> rec(Value<FloatDP> const& x) { return rec(Bounds<FloatDP>(x)); }
inline Bounds<FloatDP> sqrt(Value<FloatDP> const& x) { return sqrt(Bounds<FloatDP>(x)); }
inline Bounds<FloatDP> exp(Value<FloatDP> const& x) { return exp(Bounds<FloatDP>(x)); }
inline Bounds<FloatDP> log(Value<FloatDP> const& x) { return log(Bounds<FloatDP>(x)); }
inline Bounds<FloatDP> sin(Value<FloatDP> const& x) { return sin(Bounds<FloatDP>(x)); }
inline Bounds<FloatDP> cos(Value<FloatDP> const& x) { return cos(Bounds<FloatDP>(x)); }
inline Bounds<FloatDP> tan(Value<FloatDP> const& x) { return tan(Bounds<FloatDP>(x)); }
inline Bounds<FloatDP> asin(Value<FloatDP> const& x) { return asin(Bounds<FloatDP>(x)); }
inline Bounds<FloatDP> acos(Value<FloatDP> const& x) { return acos(Bounds<FloatDP>(x)); }
inline Bounds<FloatDP> atan(Value<FloatDP> const& x) { return atan(Bounds<FloatDP>(x)); }

inline Boolean eq(Value<FloatDP> const& x1, Value<FloatDP> const& x2) { return x1.dbl == x2.dbl; }
inline Boolean lt(Value<FloatDP> const& x1, Value<FloatDP> const& x2) { return x1.dbl <  x2.dbl; }

inline Bounds<FloatDP> operator+(Value<FloatDP> const& x1, Value<FloatDP> const& x2) { return add(x1,x2); }
inline Bounds<FloatDP> operator-(Value<FloatDP> const& x1, Value<FloatDP> const& x2) { return sub(x1,x2); }
inline Bounds<FloatDP> operator*(Value<FloatDP> const& x1, Value<FloatDP> const& x2) { return mul(x1,x2); }
inline Bounds<FloatDP> operator/(Value<FloatDP> const& x1, Value<FloatDP> const& x2) { return div(x1,x2); }

Bounds<FloatMP> add(Value<FloatMP> const& x1, Value<FloatMP> const& x2);
Bounds<FloatMP> sub(Value<FloatMP> const& x1, Value<FloatMP> const& x2);
Bounds<FloatMP> mul(Value<FloatMP> const& x1, Value<FloatMP> const& x2);
Bounds<FloatMP> div(Value<FloatMP> const& x1, Value<FloatMP> const& x2);
inline Bounds<FloatMP> pow(Value<FloatMP> const& x, Int n) { return pow(Bounds<FloatMP>(x),n); }
inline Bounds<FloatMP> sqr(Value<FloatMP> const& x) { return sqr(Bounds<FloatMP>(x)); }
inline Bounds<FloatMP> rec(Value<FloatMP> const& x) { return rec(Bounds<FloatMP>(x)); }
inline Bounds<FloatMP> sqrt(Value<FloatMP> const& x) { return sqrt(Bounds<FloatMP>(x)); }
inline Bounds<FloatMP> exp(Value<FloatMP> const& x) { return exp(Bounds<FloatMP>(x)); }
inline Bounds<FloatMP> log(Value<FloatMP> const& x) { return log(Bounds<FloatMP>(x)); }
inline Bounds<FloatMP> sin(Value<FloatMP> const& x) { return sin(Bounds<FloatMP>(x)); }
inline Bounds<FloatMP> cos(Value<FloatMP> const& x) { return cos(Bounds<FloatMP>(x)); }
inline Bounds<FloatMP> tan(Value<FloatMP> const& x) { return tan(Bounds<FloatMP>(x)); }
inline Bounds<FloatMP> asin(Value<FloatMP> const& x) { return asin(Bounds<FloatMP>(x)); }
inline Bounds<FloatMP> acos(Value<FloatMP> const& x) { return acos(Bounds<FloatMP>(x)); }
inline Bounds<FloatMP> atan(Value<FloatMP> const& x) { return atan(Bounds<FloatMP>(x)); }

inline Boolean eq(Value<FloatMP> const& x1, Value<FloatMP> const& x2) { return x1 == x2; }
inline Boolean lt(Value<FloatMP> const& x1, Value<FloatMP> const& x2) { return x1 <  x2; }

inline Bounds<FloatMP> operator+(Value<FloatMP> const& x1, Value<FloatMP> const& x2) { return add(x1,x2); }
inline Bounds<FloatMP> operator-(Value<FloatMP> const& x1, Value<FloatMP> const& x2) { return sub(x1,x2); }
inline Bounds<FloatMP> operator*(Value<FloatMP> const& x1, Value<FloatMP> const& x2) { return mul(x1,x2); }
inline Bounds<FloatMP> operator/(Value<FloatMP> const& x1, Value<FloatMP> const& x2) { return div(x1,x2); }



inline Bounds<FloatDP> add(Int const& n1, Value<FloatDP> const& x2) { return add(Value<FloatDP>(n1,x2.precision()),x2); }
inline Bounds<FloatDP> sub(Int const& n1, Value<FloatDP> const& x2) { return sub(Value<FloatDP>(n1,x2.precision()),x2); }
inline Bounds<FloatDP> mul(Int const& n1, Value<FloatDP> const& x2) { return mul(Value<FloatDP>(n1,x2.precision()),x2); }
inline Bounds<FloatDP> div(Int const& n1, Value<FloatDP> const& x2) { return div(Value<FloatDP>(n1,x2.precision()),x2); }
inline Bounds<FloatDP> add(Value<FloatDP> const& x1, Int const& n2) { return add(x1,Value<FloatDP>(n2,x1.precision())); }
inline Bounds<FloatDP> sub(Value<FloatDP> const& x1, Int const& n2) { return sub(x1,Value<FloatDP>(n2,x1.precision())); }
inline Bounds<FloatDP> mul(Value<FloatDP> const& x1, Int const& n2) { return mul(x1,Value<FloatDP>(n2,x1.precision())); }
inline Bounds<FloatDP> div(Value<FloatDP> const& x1, Int const& n2) { return div(x1,Value<FloatDP>(n2,x1.precision())); }

inline Value<FloatDP> max(Int const& n1, Value<FloatDP> const& x2) { return max(Value<FloatDP>(n1,x2.precision()),x2); }
inline Value<FloatDP> min(Int const& n1, Value<FloatDP> const& x2) { return min(Value<FloatDP>(n1,x2.precision()),x2); }
inline Value<FloatDP> max(Value<FloatDP> const& x1, Int const& n2) { return max(x1,Value<FloatDP>(n2,x1.precision())); }
inline Value<FloatDP> min(Value<FloatDP> const& x1, Int const& n2) { return min(x1,Value<FloatDP>(n2,x1.precision())); }

inline Bounds<FloatDP> operator+(Value<FloatDP> x1, Int n2) { return add(x1,n2); }
inline Bounds<FloatDP> operator-(Value<FloatDP> x1, Int n2) { return sub(x1,n2); }
inline Bounds<FloatDP> operator*(Value<FloatDP> x1, Int n2) { return mul(x1,n2); }
inline Bounds<FloatDP> operator/(Value<FloatDP> x1, Int n2) { return div(x1,n2); }
inline Bounds<FloatDP> operator+(Int n1, Value<FloatDP> x2) { return add(n1,x2); }
inline Bounds<FloatDP> operator-(Int n1, Value<FloatDP> x2) { return sub(n1,x2); }
inline Bounds<FloatDP> operator*(Int n1, Value<FloatDP> x2) { return mul(n1,x2); }
inline Bounds<FloatDP> operator/(Int n1, Value<FloatDP> x2) { return div(n1,x2); }

inline Boolean operator==(Value<FloatDP> x1, Int n2) { return x1==Value<FloatDP>(n2,x1.precision()); }
inline Boolean operator!=(Value<FloatDP> x1, Int n2) { return x1!=Value<FloatDP>(n2,x1.precision()); }
inline Boolean operator> (Value<FloatDP> x1, Int n2) { return x1> Value<FloatDP>(n2,x1.precision()); }
inline Boolean operator< (Value<FloatDP> x1, Int n2) { return x1< Value<FloatDP>(n2,x1.precision()); }
inline Boolean operator>=(Value<FloatDP> x1, Int n2) { return x1>=Value<FloatDP>(n2,x1.precision()); }
inline Boolean operator<=(Value<FloatDP> x1, Int n2) { return x1<=Value<FloatDP>(n2,x1.precision()); }
inline Boolean operator==(Int n1, Value<FloatDP> x2) { return Value<FloatDP>(n1,x2.precision())==x2; }
inline Boolean operator!=(Int n1, Value<FloatDP> x2) { return Value<FloatDP>(n1,x2.precision())!=x2; }
inline Boolean operator> (Int n1, Value<FloatDP> x2) { return Value<FloatDP>(n1,x2.precision())> x2; }
inline Boolean operator< (Int n1, Value<FloatDP> x2) { return Value<FloatDP>(n1,x2.precision())< x2; }
inline Boolean operator>=(Int n1, Value<FloatDP> x2) { return Value<FloatDP>(n1,x2.precision())>=x2; }
inline Boolean operator<=(Int n1, Value<FloatDP> x2) { return Value<FloatDP>(n1,x2.precision())<=x2; }

inline Boolean operator==(Value<FloatDP> x1, Nat m2) { return x1==Value<FloatDP>(m2,x1.precision()); }
inline Boolean operator!=(Value<FloatDP> x1, Nat m2) { return x1!=Value<FloatDP>(m2,x1.precision()); }
inline Boolean operator> (Value<FloatDP> x1, Nat m2) { return x1> Value<FloatDP>(m2,x1.precision()); }
inline Boolean operator< (Value<FloatDP> x1, Nat m2) { return x1< Value<FloatDP>(m2,x1.precision()); }
inline Boolean operator>=(Value<FloatDP> x1, Nat m2) { return x1>=Value<FloatDP>(m2,x1.precision()); }
inline Boolean operator<=(Value<FloatDP> x1, Nat m2) { return x1<=Value<FloatDP>(m2,x1.precision()); }
inline Boolean operator==(Nat m1, Value<FloatDP> x2) { return Value<FloatDP>(m1,x2.precision())==x2; }
inline Boolean operator!=(Nat m1, Value<FloatDP> x2) { return Value<FloatDP>(m1,x2.precision())!=x2; }
inline Boolean operator> (Nat m1, Value<FloatDP> x2) { return Value<FloatDP>(m1,x2.precision())> x2; }
inline Boolean operator< (Nat m1, Value<FloatDP> x2) { return Value<FloatDP>(m1,x2.precision())< x2; }
inline Boolean operator>=(Nat m1, Value<FloatDP> x2) { return Value<FloatDP>(m1,x2.precision())>=x2; }
inline Boolean operator<=(Nat m1, Value<FloatDP> x2) { return Value<FloatDP>(m1,x2.precision())<=x2; }

/* Comparisons with Dbl already present
inline Boolean operator==(Value<FloatDP> x1, Int n2) { return x1==Value<FloatDP>(n2,x1.precision()); }
inline Boolean operator!=(Value<FloatDP> x1, Int n2) { return x1!=Value<FloatDP>(n2,x1.precision()); }
inline Boolean operator< (Value<FloatDP> x1, Int n2) { return x1< Value<FloatDP>(n2,x1.precision()); }
inline Boolean operator> (Value<FloatDP> x1, Int n2) { return x1> Value<FloatDP>(n2,x1.precision()); }
inline Boolean operator<=(Value<FloatDP> x1, Int n2) { return x1<=Value<FloatDP>(n2,x1.precision()); }
inline Boolean operator>=(Value<FloatDP> x1, Int n2) { return x1>=Value<FloatDP>(n2,x1.precision()); }
inline Boolean operator==(Int n1, Value<FloatDP> x2) { return Value<FloatDP>(n1,x2.precision())==x2; }
inline Boolean operator!=(Int n1, Value<FloatDP> x2) { return Value<FloatDP>(n1,x2.precision())!=x2; }
inline Boolean operator< (Int n1, Value<FloatDP> x2) { return Value<FloatDP>(n1,x2.precision())< x2; }
inline Boolean operator> (Int n1, Value<FloatDP> x2) { return Value<FloatDP>(n1,x2.precision())> x2; }
inline Boolean operator<=(Int n1, Value<FloatDP> x2) { return Value<FloatDP>(n1,x2.precision())<=x2; }
inline Boolean operator>=(Int n1, Value<FloatDP> x2) { return Value<FloatDP>(n1,x2.precision())>=x2; }
*/

inline Bounds<FloatDP> operator+(Value<FloatDP> x1, ExactDouble d2) { return add(x1,Value<FloatDP>(d2,x1.precision())); }
inline Bounds<FloatDP> operator-(Value<FloatDP> x1, ExactDouble d2) { return sub(x1,Value<FloatDP>(d2,x1.precision())); }
inline Bounds<FloatDP> operator*(Value<FloatDP> x1, ExactDouble d2) { return mul(x1,Value<FloatDP>(d2,x1.precision())); }
inline Bounds<FloatDP> operator/(Value<FloatDP> x1, ExactDouble d2) { return div(x1,Value<FloatDP>(d2,x1.precision())); }

inline Bounds<FloatMP> add(Int const& n1, Value<FloatMP> const& x2) { return add(Value<FloatMP>(n1,x2.precision()),x2); }
inline Bounds<FloatMP> sub(Int const& n1, Value<FloatMP> const& x2) { return sub(Value<FloatMP>(n1,x2.precision()),x2); }
inline Bounds<FloatMP> mul(Int const& n1, Value<FloatMP> const& x2) { return mul(Value<FloatMP>(n1,x2.precision()),x2); }
inline Bounds<FloatMP> div(Int const& n1, Value<FloatMP> const& x2) { return div(Value<FloatMP>(n1,x2.precision()),x2); }
inline Bounds<FloatMP> add(Value<FloatMP> const& x1, Int const& n2) { return add(x1,Value<FloatMP>(n2,x1.precision())); }
inline Bounds<FloatMP> sub(Value<FloatMP> const& x1, Int const& n2) { return sub(x1,Value<FloatMP>(n2,x1.precision())); }
inline Bounds<FloatMP> mul(Value<FloatMP> const& x1, Int const& n2) { return mul(x1,Value<FloatMP>(n2,x1.precision())); }
inline Bounds<FloatMP> div(Value<FloatMP> const& x1, Int const& n2) { return div(x1,Value<FloatMP>(n2,x1.precision())); }

inline Value<FloatMP> max(Int const& n1, Value<FloatMP> const& x2) { return max(Value<FloatMP>(n1,x2.precision()),x2); }
inline Value<FloatMP> min(Int const& n1, Value<FloatMP> const& x2) { return min(Value<FloatMP>(n1,x2.precision()),x2); }
inline Value<FloatMP> max(Value<FloatMP> const& x1, Int const& n2) { return max(x1,Value<FloatMP>(n2,x1.precision())); }
inline Value<FloatMP> min(Value<FloatMP> const& x1, Int const& n2) { return min(x1,Value<FloatMP>(n2,x1.precision())); }

inline Bounds<FloatMP> operator+(Value<FloatMP> x1, Int n2) { return add(x1,n2); }
inline Bounds<FloatMP> operator-(Value<FloatMP> x1, Int n2) { return sub(x1,n2); }
inline Bounds<FloatMP> operator*(Value<FloatMP> x1, Int n2) { return mul(x1,n2); }
inline Bounds<FloatMP> operator/(Value<FloatMP> x1, Int n2) { return div(x1,n2); }
inline Bounds<FloatMP> operator+(Int n1, Value<FloatMP> x2) { return add(n1,x2); }
inline Bounds<FloatMP> operator-(Int n1, Value<FloatMP> x2) { return sub(n1,x2); }
inline Bounds<FloatMP> operator*(Int n1, Value<FloatMP> x2) { return mul(n1,x2); }
inline Bounds<FloatMP> operator/(Int n1, Value<FloatMP> x2) { return div(n1,x2); }

/* Comparisons with Dbl already present
inline Boolean operator==(Value<FloatMP> x1, Int n2) { return x1==Value<FloatMP>(n2,x1.precision()); }
inline Boolean operator!=(Value<FloatMP> x1, Int n2) { return x1!=Value<FloatMP>(n2,x1.precision()); }
inline Boolean operator< (Value<FloatMP> x1, Int n2) { return x1< Value<FloatMP>(n2,x1.precision()); }
inline Boolean operator> (Value<FloatMP> x1, Int n2) { return x1> Value<FloatMP>(n2,x1.precision()); }
inline Boolean operator<=(Value<FloatMP> x1, Int n2) { return x1<=Value<FloatMP>(n2,x1.precision()); }
inline Boolean operator>=(Value<FloatMP> x1, Int n2) { return x1>=Value<FloatMP>(n2,x1.precision()); }
inline Boolean operator==(Int n1, Value<FloatMP> x2) { return Value<FloatMP>(n1,x2.precision())==x2; }
inline Boolean operator!=(Int n1, Value<FloatMP> x2) { return Value<FloatMP>(n1,x2.precision())!=x2; }
inline Boolean operator< (Int n1, Value<FloatMP> x2) { return Value<FloatMP>(n1,x2.precision())< x2; }
inline Boolean operator> (Int n1, Value<FloatMP> x2) { return Value<FloatMP>(n1,x2.precision())> x2; }
inline Boolean operator<=(Int n1, Value<FloatMP> x2) { return Value<FloatMP>(n1,x2.precision())<=x2; }
inline Boolean operator>=(Int n1, Value<FloatMP> x2) { return Value<FloatMP>(n1,x2.precision())>=x2; }
*/

inline Bounds<FloatMP> operator+(Value<FloatMP> x1, ExactDouble d2) { return add(x1,Value<FloatMP>(d2,x1.precision())); }
inline Bounds<FloatMP> operator-(Value<FloatMP> x1, ExactDouble d2) { return sub(x1,Value<FloatMP>(d2,x1.precision())); }
inline Bounds<FloatMP> operator*(Value<FloatMP> x1, ExactDouble d2) { return mul(x1,Value<FloatMP>(d2,x1.precision())); }
inline Bounds<FloatMP> operator/(Value<FloatMP> x1, ExactDouble d2) { return div(x1,Value<FloatMP>(d2,x1.precision())); }


/*
template<ARawFloat F> Bounds<F> operator+(Int const& n1, Value<F> const& x2);
template<ARawFloat F> Bounds<F> operator-(Int const& n1, Value<F> const& x2);
template<ARawFloat F> Bounds<F> operator*(Int const& n1, Value<F> const& x2);
template<ARawFloat F> Bounds<F> operator/(Int const& n1, Value<F> const& x2);
template<ARawFloat F> Bounds<F> operator+(Value<F> const& x1, Int const& n2);
template<ARawFloat F> Bounds<F> operator-(Value<F> const& x1, Int const& n2);
template<ARawFloat F> Bounds<F> operator*(Value<F> const& x1, Int const& n2);
template<ARawFloat F> Bounds<F> operator/(Value<F> const& x1, Int const& n2);

template<ARawFloat F> Bounds<F> add(Int const& n1, Value<F> const& x2);
template<ARawFloat F> Bounds<F> sub(Int const& n1, Value<F> const& x2);
template<ARawFloat F> Bounds<F> mul(Int const& n1, Value<F> const& x2);
template<ARawFloat F> Bounds<F> div(Int const& n1, Value<F> const& x2);
template<ARawFloat F> Bounds<F> max(Int const& n1, Value<F> const& x2);
template<ARawFloat F> Bounds<F> min(Int const& n1, Value<F> const& x2);
template<ARawFloat F> Bounds<F> add(Value<F> const& x1, Int const& n2);
template<ARawFloat F> Bounds<F> sub(Value<F> const& x1, Int const& n2);
template<ARawFloat F> Bounds<F> mul(Value<F> const& x1, Int const& n2);
template<ARawFloat F> Bounds<F> div(Value<F> const& x1, Int const& n2);
template<ARawFloat F> Bounds<F> max(Value<F> const& x1, Int const& n2);
template<ARawFloat F> Bounds<F> min(Value<F> const& x1, Int const& n2);

template<ARawFloat F> Boolean operator==(Int const& n1, Value<F> const& x2);
template<ARawFloat F> Boolean operator!=(Int const& n1, Value<F> const& x2);
template<ARawFloat F> Boolean operator< (Int const& n1, Value<F> const& x2);
template<ARawFloat F> Boolean operator> (Int const& n1, Value<F> const& x2);
template<ARawFloat F> Boolean operator<=(Int const& n1, Value<F> const& x2);
template<ARawFloat F> Boolean operator>=(Int const& n1, Value<F> const& x2);
template<ARawFloat F> Boolean operator==(Value<F> const& x1, Int const& n2);
template<ARawFloat F> Boolean operator!=(Value<F> const& x1, Int const& n2);
template<ARawFloat F> Boolean operator< (Value<F> const& x1, Int const& n2);
template<ARawFloat F> Boolean operator> (Value<F> const& x1, Int const& n2);
template<ARawFloat F> Boolean operator<=(Value<F> const& x1, Int const& n2);
template<ARawFloat F> Boolean operator>=(Value<F> const& x1, Int const& n2);
*/

class Integer;
template<ARawFloat F> inline Bounds<F> add(Integer const& z1, Value<F> const& x2) { return add(Value<F>(z1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> sub(Integer const& z1, Value<F> const& x2) { return sub(Value<F>(z1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> mul(Integer const& z1, Value<F> const& x2) { return mul(Value<F>(z1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> div(Integer const& z1, Value<F> const& x2) { return div(Value<F>(z1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> add(Value<F> const& x1, Integer const& z2) { return add(x1,Value<F>(z2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> sub(Value<F> const& x1, Integer const& z2) { return sub(x1,Value<F>(z2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> mul(Value<F> const& x1, Integer const& z2) { return mul(x1,Value<F>(z2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> div(Value<F> const& x1, Integer const& z2) { return div(x1,Value<F>(z2,x1.precision())); }

template<ARawFloat F> inline Value<F> max(Integer const& z1, Value<F> const& x2) { return max(Value<F>(z1,x2.precision()),x2); }
template<ARawFloat F> inline Value<F> min(Integer const& z1, Value<F> const& x2) { return min(Value<F>(z1,x2.precision()),x2); }
template<ARawFloat F> inline Value<F> max(Value<F> const& x1, Integer const& z2) { return max(x1,Value<F>(z2,x1.precision())); }
template<ARawFloat F> inline Value<F> min(Value<F> const& x1, Integer const& z2) { return min(x1,Value<F>(z2,x1.precision())); }

template<ARawFloat F> inline Bounds<F> operator+(Value<F> x1, Integer z2) { return add(x1,z2); }
template<ARawFloat F> inline Bounds<F> operator-(Value<F> x1, Integer z2) { return sub(x1,z2); }
template<ARawFloat F> inline Bounds<F> operator*(Value<F> x1, Integer z2) { return mul(x1,z2); }
template<ARawFloat F> inline Bounds<F> operator/(Value<F> x1, Integer z2) { return div(x1,z2); }
template<ARawFloat F> inline Bounds<F> operator+(Integer z1, Value<F> x2) { return add(z1,x2); }
template<ARawFloat F> inline Bounds<F> operator-(Integer z1, Value<F> x2) { return sub(z1,x2); }
template<ARawFloat F> inline Bounds<F> operator*(Integer z1, Value<F> x2) { return mul(z1,x2); }
template<ARawFloat F> inline Bounds<F> operator/(Integer z1, Value<F> x2) { return div(z1,x2); }

template<ARawFloat F> inline Boolean operator==(Value<F> x1, Integer z2) { return Dyadic(x1)==z2; }
template<ARawFloat F> inline Boolean operator!=(Value<F> x1, Integer z2) { return Dyadic(x1)!=z2; }
template<ARawFloat F> inline Boolean operator< (Value<F> x1, Integer z2) { return Dyadic(x1)< z2; }
template<ARawFloat F> inline Boolean operator> (Value<F> x1, Integer z2) { return Dyadic(x1)> z2; }
template<ARawFloat F> inline Boolean operator<=(Value<F> x1, Integer z2) { return Dyadic(x1)<=z2; }
template<ARawFloat F> inline Boolean operator>=(Value<F> x1, Integer z2) { return Dyadic(x1)>=z2; }
template<ARawFloat F> inline Boolean operator==(Integer z1, Value<F> x2) { return z1==Dyadic(x2); }
template<ARawFloat F> inline Boolean operator!=(Integer z1, Value<F> x2) { return z1!=Dyadic(x2); }
template<ARawFloat F> inline Boolean operator< (Integer z1, Value<F> x2) { return z1< Dyadic(x2); }
template<ARawFloat F> inline Boolean operator> (Integer z1, Value<F> x2) { return z1> Dyadic(x2); }
template<ARawFloat F> inline Boolean operator<=(Integer z1, Value<F> x2) { return z1<=Dyadic(x2); }
template<ARawFloat F> inline Boolean operator>=(Integer z1, Value<F> x2) { return z1>=Dyadic(x2); }

class Dyadic;
template<ARawFloat F> inline Bounds<F> add(Dyadic const& w1, Value<F> const& x2) { return add(Value<F>(w1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> sub(Dyadic const& w1, Value<F> const& x2) { return sub(Value<F>(w1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> mul(Dyadic const& w1, Value<F> const& x2) { return mul(Value<F>(w1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> div(Dyadic const& w1, Value<F> const& x2) { return div(Value<F>(w1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> add(Value<F> const& x1, Dyadic const& w2) { return add(x1,Value<F>(w2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> sub(Value<F> const& x1, Dyadic const& w2) { return sub(x1,Value<F>(w2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> mul(Value<F> const& x1, Dyadic const& w2) { return mul(x1,Value<F>(w2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> div(Value<F> const& x1, Dyadic const& w2) { return div(x1,Value<F>(w2,x1.precision())); }

template<ARawFloat F> inline Value<F> max(Dyadic const& w1, Value<F> const& x2) { return max(Value<F>(w1,x2.precision()),x2); }
template<ARawFloat F> inline Value<F> min(Dyadic const& w1, Value<F> const& x2) { return min(Value<F>(w1,x2.precision()),x2); }
template<ARawFloat F> inline Value<F> max(Value<F> const& x1, Dyadic const& w2) { return max(x1,Value<F>(w2,x1.precision())); }
template<ARawFloat F> inline Value<F> min(Value<F> const& x1, Dyadic const& w2) { return min(x1,Value<F>(w2,x1.precision())); }

template<ARawFloat F> inline Bounds<F> operator+(Value<F> x1, Dyadic w2) { return add(x1,w2); }
template<ARawFloat F> inline Bounds<F> operator-(Value<F> x1, Dyadic w2) { return sub(x1,w2); }
template<ARawFloat F> inline Bounds<F> operator*(Value<F> x1, Dyadic w2) { return mul(x1,w2); }
template<ARawFloat F> inline Bounds<F> operator/(Value<F> x1, Dyadic w2) { return div(x1,w2); }
template<ARawFloat F> inline Bounds<F> operator+(Dyadic w1, Value<F> x2) { return add(w1,x2); }
template<ARawFloat F> inline Bounds<F> operator-(Dyadic w1, Value<F> x2) { return sub(w1,x2); }
template<ARawFloat F> inline Bounds<F> operator*(Dyadic w1, Value<F> x2) { return mul(w1,x2); }
template<ARawFloat F> inline Bounds<F> operator/(Dyadic w1, Value<F> x2) { return div(w1,x2); }

template<ARawFloat F> inline Boolean operator==(Value<F> x1, Dyadic w2) { return Dyadic(x1)==w2; }
template<ARawFloat F> inline Boolean operator!=(Value<F> x1, Dyadic w2) { return Dyadic(x1)!=w2; }
template<ARawFloat F> inline Boolean operator< (Value<F> x1, Dyadic w2) { return Dyadic(x1)< w2; }
template<ARawFloat F> inline Boolean operator> (Value<F> x1, Dyadic w2) { return Dyadic(x1)> w2; }
template<ARawFloat F> inline Boolean operator<=(Value<F> x1, Dyadic w2) { return Dyadic(x1)<=w2; }
template<ARawFloat F> inline Boolean operator>=(Value<F> x1, Dyadic w2) { return Dyadic(x1)>=w2; }
template<ARawFloat F> inline Boolean operator==(Dyadic w1, Value<F> x2) { return w1==Dyadic(x2); }
template<ARawFloat F> inline Boolean operator!=(Dyadic w1, Value<F> x2) { return w1!=Dyadic(x2); }
template<ARawFloat F> inline Boolean operator< (Dyadic w1, Value<F> x2) { return w1< Dyadic(x2); }
template<ARawFloat F> inline Boolean operator> (Dyadic w1, Value<F> x2) { return w1> Dyadic(x2); }
template<ARawFloat F> inline Boolean operator<=(Dyadic w1, Value<F> x2) { return w1<=Dyadic(x2); }
template<ARawFloat F> inline Boolean operator>=(Dyadic w1, Value<F> x2) { return w1>=Dyadic(x2); }

class Rational;
template<ARawFloat F> inline Bounds<F> add(Rational const& q1, Value<F> const& x2) { return add(Bounds<F>(q1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> sub(Rational const& q1, Value<F> const& x2) { return sub(Bounds<F>(q1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> mul(Rational const& q1, Value<F> const& x2) { return mul(Bounds<F>(q1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> div(Rational const& q1, Value<F> const& x2) { return div(Bounds<F>(q1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> add(Value<F> const& x1, Rational const& q2) { return add(x1,Bounds<F>(q2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> sub(Value<F> const& x1, Rational const& q2) { return sub(x1,Bounds<F>(q2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> mul(Value<F> const& x1, Rational const& q2) { return mul(x1,Bounds<F>(q2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> div(Value<F> const& x1, Rational const& q2) { return div(x1,Bounds<F>(q2,x1.precision())); }

template<ARawFloat F> inline Bounds<F> max(Rational const& q1, Value<F> const& x2) { return max(Bounds<F>(q1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> min(Rational const& q1, Value<F> const& x2) { return min(Bounds<F>(q1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> max(Value<F> const& x1, Rational const& q2) { return max(x1,Bounds<F>(q2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> min(Value<F> const& x1, Rational const& q2) { return min(x1,Bounds<F>(q2,x1.precision())); }

template<ARawFloat F> inline Bounds<F> operator+(Value<F> x1, Rational q2) { return add(x1,q2); }
template<ARawFloat F> inline Bounds<F> operator-(Value<F> x1, Rational q2) { return sub(x1,q2); }
template<ARawFloat F> inline Bounds<F> operator*(Value<F> x1, Rational q2) { return mul(x1,q2); }
template<ARawFloat F> inline Bounds<F> operator/(Value<F> x1, Rational q2) { return div(x1,q2); }
template<ARawFloat F> inline Bounds<F> operator+(Rational q1, Value<F> x2) { return add(q1,x2); }
template<ARawFloat F> inline Bounds<F> operator-(Rational q1, Value<F> x2) { return sub(q1,x2); }
template<ARawFloat F> inline Bounds<F> operator*(Rational q1, Value<F> x2) { return mul(q1,x2); }
template<ARawFloat F> inline Bounds<F> operator/(Rational q1, Value<F> x2) { return div(q1,x2); }

template<ARawFloat F> ValidatedKleenean operator==(Rational const& q1, Value<F> const& x2) { return Bounds<F>(q1,x2.precision())==x2; }
template<ARawFloat F> ValidatedKleenean operator!=(Rational const& q1, Value<F> const& x2) { return Bounds<F>(q1,x2.precision())!=x2; }
template<ARawFloat F> ValidatedKleenean operator< (Rational const& q1, Value<F> const& x2) { return Bounds<F>(q1,x2.precision())< x2; }
template<ARawFloat F> ValidatedKleenean operator> (Rational const& q1, Value<F> const& x2) { return Bounds<F>(q1,x2.precision())> x2; }
template<ARawFloat F> ValidatedKleenean operator<=(Rational const& q1, Value<F> const& x2) { return Bounds<F>(q1,x2.precision())<=x2; }
template<ARawFloat F> ValidatedKleenean operator>=(Rational const& q1, Value<F> const& x2) { return Bounds<F>(q1,x2.precision())>=x2; }
template<ARawFloat F> ValidatedKleenean operator==(Value<F> const& x1, Rational const& q2) { return x1==Bounds<F>(q2,x1.precision()); }
template<ARawFloat F> ValidatedKleenean operator!=(Value<F> const& x1, Rational const& q2) { return x1!=Bounds<F>(q2,x1.precision()); }
template<ARawFloat F> ValidatedKleenean operator< (Value<F> const& x1, Rational const& q2) { return x1< Bounds<F>(q2,x1.precision()); }
template<ARawFloat F> ValidatedKleenean operator> (Value<F> const& x1, Rational const& q2) { return x1> Bounds<F>(q2,x1.precision()); }
template<ARawFloat F> ValidatedKleenean operator<=(Value<F> const& x1, Rational const& q2) { return x1<=Bounds<F>(q2,x1.precision()); }
template<ARawFloat F> ValidatedKleenean operator>=(Value<F> const& x1, Rational const& q2) { return x1>=Bounds<F>(q2,x1.precision()); }

class Real;
template<ARawFloat F> inline Bounds<F> add(Real const& r1, Value<F> const& x2) { return add(Bounds<F>(r1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> sub(Real const& r1, Value<F> const& x2) { return sub(Bounds<F>(r1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> mul(Real const& r1, Value<F> const& x2) { return mul(Bounds<F>(r1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> div(Real const& r1, Value<F> const& x2) { return div(Bounds<F>(r1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> add(Value<F> const& x1, Real const& r2) { return add(x1,Bounds<F>(r2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> sub(Value<F> const& x1, Real const& r2) { return sub(x1,Bounds<F>(r2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> mul(Value<F> const& x1, Real const& r2) { return mul(x1,Bounds<F>(r2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> div(Value<F> const& x1, Real const& r2) { return div(x1,Bounds<F>(r2,x1.precision())); }

template<ARawFloat F> inline Bounds<F> max(Real const& r1, Value<F> const& x2) { return max(Bounds<F>(r1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> min(Real const& r1, Value<F> const& x2) { return min(Bounds<F>(r1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> max(Value<F> const& x1, Real const& r2) { return max(x1,Bounds<F>(r2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> min(Value<F> const& x1, Real const& r2) { return min(x1,Bounds<F>(r2,x1.precision())); }

template<ARawFloat F> inline Bounds<F> operator+(Value<F> x1, Real r2) { return add(x1,r2); }
template<ARawFloat F> inline Bounds<F> operator-(Value<F> x1, Real r2) { return sub(x1,r2); }
template<ARawFloat F> inline Bounds<F> operator*(Value<F> x1, Real r2) { return mul(x1,r2); }
template<ARawFloat F> inline Bounds<F> operator/(Value<F> x1, Real r2) { return div(x1,r2); }
template<ARawFloat F> inline Bounds<F> operator+(Real r1, Value<F> x2) { return add(r1,x2); }
template<ARawFloat F> inline Bounds<F> operator-(Real r1, Value<F> x2) { return sub(r1,x2); }
template<ARawFloat F> inline Bounds<F> operator*(Real r1, Value<F> x2) { return mul(r1,x2); }
template<ARawFloat F> inline Bounds<F> operator/(Real r1, Value<F> x2) { return div(r1,x2); }

template<ARawFloat F> ValidatedKleenean operator==(Real const& r1, Value<F> const& x2) { return Bounds<F>(r1,x2.precision())==x2; }
template<ARawFloat F> ValidatedKleenean operator!=(Real const& r1, Value<F> const& x2) { return Bounds<F>(r1,x2.precision())!=x2; }
template<ARawFloat F> ValidatedKleenean operator< (Real const& r1, Value<F> const& x2) { return Bounds<F>(r1,x2.precision())< x2; }
template<ARawFloat F> ValidatedKleenean operator> (Real const& r1, Value<F> const& x2) { return Bounds<F>(r1,x2.precision())> x2; }
template<ARawFloat F> ValidatedKleenean operator<=(Real const& r1, Value<F> const& x2) { return Bounds<F>(r1,x2.precision())<=x2; }
template<ARawFloat F> ValidatedKleenean operator>=(Real const& r1, Value<F> const& x2) { return Bounds<F>(r1,x2.precision())>=x2; }
template<ARawFloat F> ValidatedKleenean operator==(Value<F> const& x1, Real const& r2) { return x1==Bounds<F>(r2,x1.precision()); }
template<ARawFloat F> ValidatedKleenean operator!=(Value<F> const& x1, Real const& r2) { return x1!=Bounds<F>(r2,x1.precision()); }
template<ARawFloat F> ValidatedKleenean operator< (Value<F> const& x1, Real const& r2) { return x1< Bounds<F>(r2,x1.precision()); }
template<ARawFloat F> ValidatedKleenean operator> (Value<F> const& x1, Real const& r2) { return x1> Bounds<F>(r2,x1.precision()); }
template<ARawFloat F> ValidatedKleenean operator<=(Value<F> const& x1, Real const& r2) { return x1<=Bounds<F>(r2,x1.precision()); }
template<ARawFloat F> ValidatedKleenean operator>=(Value<F> const& x1, Real const& r2) { return x1>=Bounds<F>(r2,x1.precision()); }

using ValidatedNumber = Number<ValidatedTag>;
template<ARawFloat F> inline Bounds<F> add(ValidatedNumber const& y1, Value<F> const& x2) { return add(Bounds<F>(y1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> sub(ValidatedNumber const& y1, Value<F> const& x2) { return sub(Bounds<F>(y1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> mul(ValidatedNumber const& y1, Value<F> const& x2) { return mul(Bounds<F>(y1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> div(ValidatedNumber const& y1, Value<F> const& x2) { return div(Bounds<F>(y1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> add(Value<F> const& x1, ValidatedNumber const& y2) { return add(x1,Bounds<F>(y2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> sub(Value<F> const& x1, ValidatedNumber const& y2) { return sub(x1,Bounds<F>(y2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> mul(Value<F> const& x1, ValidatedNumber const& y2) { return mul(x1,Bounds<F>(y2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> div(Value<F> const& x1, ValidatedNumber const& y2) { return div(x1,Bounds<F>(y2,x1.precision())); }

template<ARawFloat F> inline Bounds<F> max(ValidatedNumber const& y1, Value<F> const& x2) { return max(Bounds<F>(y1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> min(ValidatedNumber const& y1, Value<F> const& x2) { return min(Bounds<F>(y1,x2.precision()),x2); }
template<ARawFloat F> inline Bounds<F> max(Value<F> const& x1, ValidatedNumber const& y2) { return max(x1,Bounds<F>(y2,x1.precision())); }
template<ARawFloat F> inline Bounds<F> min(Value<F> const& x1, ValidatedNumber const& y2) { return min(x1,Bounds<F>(y2,x1.precision())); }

template<ARawFloat F> inline Bounds<F> operator+(Value<F> x1, ValidatedNumber y2) { return add(x1,y2); }
template<ARawFloat F> inline Bounds<F> operator-(Value<F> x1, ValidatedNumber y2) { return sub(x1,y2); }
template<ARawFloat F> inline Bounds<F> operator*(Value<F> x1, ValidatedNumber y2) { return mul(x1,y2); }
template<ARawFloat F> inline Bounds<F> operator/(Value<F> x1, ValidatedNumber y2) { return div(x1,y2); }
template<ARawFloat F> inline Bounds<F> operator+(ValidatedNumber y1, Value<F> x2) { return add(y1,x2); }
template<ARawFloat F> inline Bounds<F> operator-(ValidatedNumber y1, Value<F> x2) { return sub(y1,x2); }
template<ARawFloat F> inline Bounds<F> operator*(ValidatedNumber y1, Value<F> x2) { return mul(y1,x2); }
template<ARawFloat F> inline Bounds<F> operator/(ValidatedNumber y1, Value<F> x2) { return div(y1,x2); }

template<ARawFloat F> ValidatedKleenean operator==(ValidatedNumber const& y1, Value<F> const& x2) { return Bounds<F>(y1,x2.precision())==x2; }
template<ARawFloat F> ValidatedKleenean operator!=(ValidatedNumber const& y1, Value<F> const& x2) { return Bounds<F>(y1,x2.precision())!=x2; }
template<ARawFloat F> ValidatedKleenean operator< (ValidatedNumber const& y1, Value<F> const& x2) { return Bounds<F>(y1,x2.precision())< x2; }
template<ARawFloat F> ValidatedKleenean operator> (ValidatedNumber const& y1, Value<F> const& x2) { return Bounds<F>(y1,x2.precision())> x2; }
template<ARawFloat F> ValidatedKleenean operator<=(ValidatedNumber const& y1, Value<F> const& x2) { return Bounds<F>(y1,x2.precision())<=x2; }
template<ARawFloat F> ValidatedKleenean operator>=(ValidatedNumber const& y1, Value<F> const& x2) { return Bounds<F>(y1,x2.precision())>=x2; }
template<ARawFloat F> ValidatedKleenean operator==(Value<F> const& x1, ValidatedNumber const& y2) { return x1==Bounds<F>(y2,x1.precision()); }
template<ARawFloat F> ValidatedKleenean operator!=(Value<F> const& x1, ValidatedNumber const& y2) { return x1!=Bounds<F>(y2,x1.precision()); }
template<ARawFloat F> ValidatedKleenean operator< (Value<F> const& x1, ValidatedNumber const& y2) { return x1< Bounds<F>(y2,x1.precision()); }
template<ARawFloat F> ValidatedKleenean operator> (Value<F> const& x1, ValidatedNumber const& y2) { return x1> Bounds<F>(y2,x1.precision()); }
template<ARawFloat F> ValidatedKleenean operator<=(Value<F> const& x1, ValidatedNumber const& y2) { return x1<=Bounds<F>(y2,x1.precision()); }
template<ARawFloat F> ValidatedKleenean operator>=(Value<F> const& x1, ValidatedNumber const& y2) { return x1>=Bounds<F>(y2,x1.precision()); }



template<ARawFloat F> Bounds<F> add(Bounds<F> const& x1, Value<F> const& x2) {
    return Bounds<F>(add(down,x1.lower_raw(),x2),add(up,x1.upper_raw(),x2)); }
template<ARawFloat F> Bounds<F> sub(Bounds<F> const& x1, Value<F> const& x2) {
    return Bounds<F>(sub(down,x1.lower_raw(),x2),sub(up,x1.upper_raw(),x2)); }
template<ARawFloat F> Bounds<F> mul(Bounds<F> const& x1, Value<F> const& x2) {
    return mul(x1,Bounds<F>(x2)); }
template<ARawFloat F> Bounds<F> div(Bounds<F> const& x1, Value<F> const& x2) {
    return div(x1,Bounds<F>(x2)); }
template<ARawFloat F> Bounds<F> add(Value<F> const& x1, Bounds<F> const& x2) {
    return Bounds<F>(add(down,x1,x2.lower_raw()),add(up,x1,x2.upper_raw())); }
template<ARawFloat F> Bounds<F> sub(Value<F> const& x1, Bounds<F> const& x2) {
    return Bounds<F>(sub(down,x1,x2.upper_raw()),sub(up,x1,x2.lower_raw())); }
template<ARawFloat F> Bounds<F> mul(Value<F> const& x1, Bounds<F> const& x2) {
    return mul(Bounds<F>(x1),x2); }
template<ARawFloat F> Bounds<F> div(Value<F> const& x1, Bounds<F> const& x2) {
    return div(Bounds<F>(x1),x2); }

template<ARawFloat F> Bounds<F> max(Bounds<F> const& x1, Value<F> const& x2) {
    return Bounds<F>(max(down,x1.lower_raw(),x2),max(up,x1.upper_raw(),x2)); }
template<ARawFloat F> Bounds<F> min(Bounds<F> const& x1, Value<F> const& x2) {
    return Bounds<F>(min(down,x1.lower_raw(),x2),min(up,x1.upper_raw(),x2)); }
template<ARawFloat F> Bounds<F> max(Value<F> const& x1, Bounds<F> const& x2) {
    return Bounds<F>(max(down,x1,x2.lower_raw()),max(up,x1,x2.upper_raw())); }
template<ARawFloat F> Bounds<F> min(Value<F> const& x1, Bounds<F> const& x2) {
    return Bounds<F>(min(down,x1,x2.lower_raw()),min(up,x1,x2.upper_raw())); }

template<ARawFloat F> Bounds<F> operator+(Bounds<F> const& x1, Value<F> const& x2) { return add(x1,x2); }
template<ARawFloat F> Bounds<F> operator-(Bounds<F> const& x1, Value<F> const& x2) { return sub(x1,x2); }
template<ARawFloat F> Bounds<F> operator*(Bounds<F> const& x1, Value<F> const& x2) { return mul(x1,x2); }
template<ARawFloat F> Bounds<F> operator/(Bounds<F> const& x1, Value<F> const& x2) { return div(x1,x2); }
template<ARawFloat F> Bounds<F> operator+(Value<F> const& x1, Bounds<F> const& x2) { return add(x1,x2); }
template<ARawFloat F> Bounds<F> operator-(Value<F> const& x1, Bounds<F> const& x2) { return sub(x1,x2); }
template<ARawFloat F> Bounds<F> operator*(Value<F> const& x1, Bounds<F> const& x2) { return mul(x1,x2); }
template<ARawFloat F> Bounds<F> operator/(Value<F> const& x1, Bounds<F> const& x2) { return div(x1,x2); }

template<ARawFloat F> ValidatedKleenean operator==(Bounds<F> const& x1, Value<F> const& x2) { return x1==Bounds<F>(x2); }
template<ARawFloat F> ValidatedKleenean operator!=(Bounds<F> const& x1, Value<F> const& x2) { return x1!=Bounds<F>(x2); }
template<ARawFloat F> ValidatedKleenean operator< (Bounds<F> const& x1, Value<F> const& x2) { return x1< Bounds<F>(x2); }
template<ARawFloat F> ValidatedKleenean operator> (Bounds<F> const& x1, Value<F> const& x2) { return x1> Bounds<F>(x2); }
template<ARawFloat F> ValidatedKleenean operator<=(Bounds<F> const& x1, Value<F> const& x2) { return x1<=Bounds<F>(x2); }
template<ARawFloat F> ValidatedKleenean operator>=(Bounds<F> const& x1, Value<F> const& x2) { return x1>=Bounds<F>(x2); }
template<ARawFloat F> ValidatedKleenean operator==(Value<F> const& x1, Bounds<F> const& x2) { return Bounds<F>(x1)==x2; }
template<ARawFloat F> ValidatedKleenean operator!=(Value<F> const& x1, Bounds<F> const& x2) { return Bounds<F>(x1)!=x2; }
template<ARawFloat F> ValidatedKleenean operator< (Value<F> const& x1, Bounds<F> const& x2) { return Bounds<F>(x1)< x2; }
template<ARawFloat F> ValidatedKleenean operator> (Value<F> const& x1, Bounds<F> const& x2) { return Bounds<F>(x1)> x2; }
template<ARawFloat F> ValidatedKleenean operator<=(Value<F> const& x1, Bounds<F> const& x2) { return Bounds<F>(x1)<=x2; }
template<ARawFloat F> ValidatedKleenean operator>=(Value<F> const& x1, Bounds<F> const& x2) { return Bounds<F>(x1)>=x2; }


template<ARawFloat F, ARawFloat FE> Ball<F,FE> add(Ball<F,FE> const& x1, Value<F> const& x2) {
    return add(x1,Ball<F,FE>(x2,x1.error_precision())); }
template<ARawFloat F, ARawFloat FE> Ball<F,FE> sub(Ball<F,FE> const& x1, Value<F> const& x2) {
    return sub(x1,Ball<F,FE>(x2,x1.error_precision())); }
template<ARawFloat F, ARawFloat FE> Ball<F,FE> mul(Ball<F,FE> const& x1, Value<F> const& x2) {
    return mul(x1,Ball<F,FE>(x2,x1.error_precision())); }
template<ARawFloat F, ARawFloat FE> Ball<F,FE> div(Ball<F,FE> const& x1, Value<F> const& x2) {
    return div(x1,Ball<F,FE>(x2,x1.error_precision())); }
template<ARawFloat F, ARawFloat FE> Ball<F,FE> add(Value<F> const& x1, Ball<F,FE> const& x2) {
    return add(Ball<F,FE>(x1,x2.error_precision()),x2); }
template<ARawFloat F, ARawFloat FE> Ball<F,FE> sub(Value<F> const& x1, Ball<F,FE> const& x2) {
    return sub(Ball<F,FE>(x1,x2.error_precision()),x2); }
template<ARawFloat F, ARawFloat FE> Ball<F,FE> mul(Value<F> const& x1, Ball<F,FE> const& x2) {
    return mul(Ball<F,FE>(x1,x2.error_precision()),x2); }
template<ARawFloat F, ARawFloat FE> Ball<F,FE> div(Value<F> const& x1, Ball<F,FE> const& x2) {
    return div(Ball<F,FE>(x1,x2.error_precision()),x2); }

template<ARawFloat F, ARawFloat FE> Ball<F,FE> max(Ball<F,FE> const& x1, Value<F> const& x2) {
    return max(x1,Ball<F,FE>(x2,x1.error_precision())); }
template<ARawFloat F, ARawFloat FE> Ball<F,FE> min(Ball<F,FE> const& x1, Value<F> const& x2) {
    return min(x1,Ball<F,FE>(x2,x1.error_precision())); }
template<ARawFloat F, ARawFloat FE> Ball<F,FE> max(Value<F> const& x1, Ball<F,FE> const& x2) {
    return max(Ball<F,FE>(x1,x2.error_precision()),x2); }
template<ARawFloat F, ARawFloat FE> Ball<F,FE> min(Value<F> const& x1, Ball<F,FE> const& x2) {
    return min(Ball<F,FE>(x1,x2.error_precision()),x2); }

template<ARawFloat F, ARawFloat FE> Ball<F,FE> operator+(Ball<F,FE> const& x1, Value<F> const& x2) { return add(x1,x2); }
template<ARawFloat F, ARawFloat FE> Ball<F,FE> operator-(Ball<F,FE> const& x1, Value<F> const& x2) { return sub(x1,x2); }
template<ARawFloat F, ARawFloat FE> Ball<F,FE> operator*(Ball<F,FE> const& x1, Value<F> const& x2) { return mul(x1,x2); }
template<ARawFloat F, ARawFloat FE> Ball<F,FE> operator/(Ball<F,FE> const& x1, Value<F> const& x2) { return div(x1,x2); }
template<ARawFloat F, ARawFloat FE> Ball<F,FE> operator+(Value<F> const& x1, Ball<F,FE> const& x2) { return add(x1,x2); }
template<ARawFloat F, ARawFloat FE> Ball<F,FE> operator-(Value<F> const& x1, Ball<F,FE> const& x2) { return sub(x1,x2); }
template<ARawFloat F, ARawFloat FE> Ball<F,FE> operator*(Value<F> const& x1, Ball<F,FE> const& x2) { return mul(x1,x2); }
template<ARawFloat F, ARawFloat FE> Ball<F,FE> operator/(Value<F> const& x1, Ball<F,FE> const& x2) { return div(x1,x2); }

template<ARawFloat F, ARawFloat FE> ValidatedKleenean operator==(Ball<F,FE> const& x1, Value<F> const& x2) {
    return x1==Ball<F,FE>(x2,x1.error_precision()); }
template<ARawFloat F, ARawFloat FE> ValidatedKleenean operator!=(Ball<F,FE> const& x1, Value<F> const& x2) {
    return x1!=Ball<F,FE>(x2,x1.error_precision()); }
template<ARawFloat F, ARawFloat FE> ValidatedKleenean operator< (Ball<F,FE> const& x1, Value<F> const& x2) {
    return x1< Ball<F,FE>(x2,x1.error_precision()); }
template<ARawFloat F, ARawFloat FE> ValidatedKleenean operator> (Ball<F,FE> const& x1, Value<F> const& x2) {
    return x1> Ball<F,FE>(x2,x1.error_precision()); }
template<ARawFloat F, ARawFloat FE> ValidatedKleenean operator<=(Ball<F,FE> const& x1, Value<F> const& x2) {
    return x1<=Ball<F,FE>(x2,x1.error_precision()); }
template<ARawFloat F, ARawFloat FE> ValidatedKleenean operator>=(Ball<F,FE> const& x1, Value<F> const& x2) {
    return x1>=Ball<F,FE>(x2,x1.error_precision()); }
template<ARawFloat F, ARawFloat FE> ValidatedKleenean operator==(Value<F> const& x1, Ball<F,FE> const& x2) {
    return Ball<F,FE>(x1,x2.error_precision())==x2; }
template<ARawFloat F, ARawFloat FE> ValidatedKleenean operator!=(Value<F> const& x1, Ball<F,FE> const& x2) {
    return Ball<F,FE>(x1,x2.error_precision())!=x2; }
template<ARawFloat F, ARawFloat FE> ValidatedKleenean operator< (Value<F> const& x1, Ball<F,FE> const& x2) {
    return Ball<F,FE>(x1,x2.error_precision())< x2; }
template<ARawFloat F, ARawFloat FE> ValidatedKleenean operator> (Value<F> const& x1, Ball<F,FE> const& x2) {
    return Ball<F,FE>(x1,x2.error_precision())> x2; }
template<ARawFloat F, ARawFloat FE> ValidatedKleenean operator<=(Value<F> const& x1, Ball<F,FE> const& x2) {
    return Ball<F,FE>(x1,x2.error_precision())<=x2; }
template<ARawFloat F, ARawFloat FE> ValidatedKleenean operator>=(Value<F> const& x1, Ball<F,FE> const& x2) {
    return Ball<F,FE>(x1,x2.error_precision())>=x2; }


template<ARawFloat F, class PRE> Ball<F,RawFloatType<PRE>> add(Value<F> const& x1, Value<F> const& x2, PRE pre) {
    typedef RawFloat<PRE> FE; return Operations<Ball<F,FE>>::_add(x1,x2,pre); }
template<ARawFloat F, class PRE> Ball<F,RawFloatType<PRE>> sub(Value<F> const& x1, Value<F> const& x2, PRE pre) {
    typedef RawFloat<PRE> FE; return Operations<Ball<F,FE>>::_sub(x1,x2,pre); }
template<ARawFloat F, class PRE> Ball<F,RawFloatType<PRE>> mul(Value<F> const& x1, Value<F> const& x2, PRE pre) {
    typedef RawFloat<PRE> FE; return Operations<Ball<F,FE>>::_mul(x1,x2,pre); }
template<ARawFloat F, class PRE> Ball<F,RawFloatType<PRE>> div(Value<F> const& x1, Value<F> const& x2, PRE pre) {
    typedef RawFloat<PRE> FE; return Operations<Ball<F,FE>>::_div(x1,x2,pre); }

template<ARawFloat F> Value<F> operator/(Value<F> const& v1, TwoExp const& e2) { return shft(v1,-e2.exponent()); }
template<ARawFloat F> Value<F>& operator/=(Value<F>& v1, TwoExp const& e2) { return v1=v1/e2; }


template<ARawFloat F> Error<F> operator+(Positive<Value<F>> const& v1, Error<F> const& e2);

template<ARawFloat F> Bool same(Value<F> const& v1, Value<F> const& v2);
template<ARawFloat F> Positive<F> cast_exact(Positive<UpperBound<F>> const&);
template<ARawFloat F> Positive<F> cast_exact(Positive<Approximation<F>> const&);

}


#endif
