/***************************************************************************
 *            float_value.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file float_value.hpp
 *  \brief Exact Floating-point representations of real numbers.
 */

#ifndef ARIADNE_FLOAT_VALUE_HPP
#define ARIADNE_FLOAT_VALUE_HPP

#include "../utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"

#include "float_operations.hpp"

#include "logical.hpp"
#include "builtin.hpp"
#include "twoexp.hpp"

namespace Ariadne {

template<class F> struct NumericTraits<Value<F>> {
    typedef ExactNumber GenericType;
    typedef Value<F> OppositeType;
    typedef PositiveValue<F> PositiveType;
    typedef Boolean LessType;
    typedef Boolean EqualsType;
};

static_assert(not IsGenericNumericType<FloatValue<DoublePrecision>>::value,"");
static_assert(not IsGenericNumericType<FloatValue<MultiplePrecision>>::value,"");

//! \ingroup NumericModule
//! \brief A floating-point number, which is taken to represent the \em exact value of a real quantity.
//! \sa FloatDP , FloatMP, FloatBall, FloatBounds, FloatApproximation.
template<class F> class Value
    : DispatchNumericOperations<Value<F>,Bounds<F>>
    , DispatchComparisonOperations<Value<F>,Boolean>
    , DefineMixedComparisonOperators<Value<F>,ExactNumber,Boolean>
    , DefineMixedComparisonOperators<Value<F>,Rational,Boolean>
    , DefineMixedComparisonOperators<Value<F>,Dyadic,Boolean>
//    , DefineMixedComparisonOperators<Value<F>,Integer,Boolean>
//    , DefineMixedComparisonOperators<Value<F>,Int,Boolean>
//        , public DispatchFloatOperations<Ball<F>>
        , public DispatchFloatOperations<Bounds<F>>
    , DefineConcreteGenericArithmeticOperators<Value<F>>
    , DefineConcreteGenericComparisonOperators<Value<F>>
{
  protected:
    typedef ExactTag P; typedef typename F::PrecisionType PR;
  public:
    typedef ExactTag Paradigm;
    typedef Value<F> NumericType;
    typedef ExactNumber GenericType;
    typedef F RawType;
    typedef PR PrecisionType;
    typedef PR PropertiesType;
  public:
    Value<F>() : _v(0.0) { }
    explicit Value<F>(PrecisionType pr) : _v(0.0,pr) { }
    explicit Value<F>(RawType const& v) : _v(v) { }

    template<class N, EnableIf<IsBuiltinIntegral<N>> = dummy> Value<F>(N n, PR pr) : Value<F>(ExactDouble(n),pr) { }
    Value<F>(ExactDouble d, PR pr);
    Value<F>(const Integer& z, PR pr);
    Value<F>(const TwoExp& t, PR pr);
    Value<F>(const Dyadic& w, PR pr);
    Value<F>(const Value<F>& x, PR pr);

    template<class N, EnableIf<IsBuiltinIntegral<N>> = dummy> Value<F>& operator=(N n) { _v=n; return *this; }
    Value<F>& operator=(const Integer& z);
    Value<F>& operator=(const TwoExp& t);
    Value<F>& operator=(const Dyadic& w);

    operator ExactNumber () const;
    explicit operator Dyadic () const;
    explicit operator Rational () const;

    Ball<F> create(ValidatedNumber const&) const;
//    explicit operator RawType () const { return _v; }

    PrecisionType precision() const { return _v.precision(); }
    PropertiesType properties() const { return _v.precision(); }
    GenericType generic() const { return this->operator GenericType(); }
    RawType const& raw() const { return _v; }
    RawType& raw() { return _v; }
    double get_d() const { return _v.get_d(); }

    template<class FE> FloatBall<PR,Ariadne::PrecisionType<FE>> pm(Error<FE> e) const;
  public:
    friend Value<F> operator*(Value<F> const&, TwoExp const&);
    friend Value<F> operator/(Value<F> const&, TwoExp const&);
    friend Value<F>& operator*=(Value<F>&, TwoExp const&);
    friend Value<F>& operator/=(Value<F>&, TwoExp const&);
    friend Error<F> mag(Value<F> const&);
    friend LowerBound<F> mig(Value<F> const&);
    friend Bool same(Value<F> const&, Value<F> const&);
    friend OutputStream& operator<<(OutputStream&, Value<F> const&);
  public:
    friend Comparison cmp(Value<F> const& x1, Rational const& q2) { return cmp(x1.raw(),q2); }
    friend Comparison cmp(Value<F> const& x1, Dyadic const& w2) { return cmp(x1.raw(),w2); }
    friend Comparison cmp(Value<F> const& x1, Integer const& z2) { return cmp(x1.raw(),z2); }
  public:
    static Nat output_places;
    static Void set_output_places(Nat p) { output_places=p; }
  private: public:
    RawType _v;
  private:
    friend Value<F> shft(Value<F> const& x, Int n) {
        return Value<F>(shft(x.raw(),n)); }
    friend Value<F> operator*(Value<F> const& x, TwoExp const& y) {
        return Value<F>(mul(near,x.raw(),F(y,x.precision()))); }
    friend Value<F> operator/(Value<F> const& x, TwoExp const& y) {
        return Value<F>(div(near,x.raw(),F(y,x.precision()))); }
    friend Value<F>& operator*=(Value<F>& x, TwoExp const& y) { return x=x*y; }
    friend Value<F>& operator/=(Value<F>& x, TwoExp const& y) { return x=x/y; }
    friend OutputStream& operator<<(OutputStream& os, Value<F> const& x) {
        return Operations<Value<F>>::_write(os,x); }
};

/*
template<class PR> class FloatValue : public Value<RawFloatType<PR>> {
    typedef RawFloatType<PR> F;
    static Nat output_places;
    using Value<F>::Value;
};
*/

template<class PR> Value(Dyadic, PR) -> Value<RawFloatType<PR>>;
template<class F> Value(F) -> Value<F>;

template<class F> inline FloatFactory<PrecisionType<F>> factory(Value<F> const& flt) { return FloatFactory<PrecisionType<F>>(flt.precision()); }
template<class PR> inline FloatValue<PR> FloatFactory<PR>::create(Dyadic const& y, ExactTag) { return FloatValue<PR>(y,_pr); }
template<class PR> inline FloatValue<PR> FloatFactory<PR>::create(Integer const& y, ExactTag) { return FloatValue<PR>(y,_pr); }
template<class PR> inline FloatValue<PR> FloatFactory<PR>::create(ExactDouble const& y) { return FloatValue<PR>(y,_pr); }

template<class PR> template<class N, EnableIf<IsBuiltinSignedIntegral<N>>> inline
FloatValue<PR> FloatFactory<PR>::create(N const& y) { return FloatValue<PR>(y,_pr); }

template<class PR> template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>>> inline
PositiveFloatValue<PR> FloatFactory<PR>::create(M const& y) { return PositiveFloatValue<PR>(y,_pr); }

template<class F> class Positive<Value<F>> : public Value<F> {
    using typename Value<F>::PR;
  public:
    Positive<Value<F>>() : Value<F>() { }
    explicit Positive<Value<F>>(PR const& pr) : Value<F>(pr) { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy> Positive<Value<F>>(M m, PR pr) : Value<F>(m,pr) { }
    Positive<Value<F>>(TwoExp const& ex, PR pr) : Value<F>(ex,pr) { }
    explicit Positive<Value<F>>(Dyadic const& w, PR pr) : Value<F>(w,pr) { }
    explicit Positive<Value<F>>(F const& x) : Value<F>(x) { }
    explicit Positive<Value<F>>(Value<F> const& x) : Value<F>(x) { }
  public:
    friend PositiveBounds<F> operator+(PositiveValue<F> const& v1, PositiveValue<F> const& v2) {
        return cast_positive(static_cast<Value<F>const&>(v1)+static_cast<Value<F>const&>(v2)); }
    friend PositiveBounds<F> operator*(PositiveValue<F> const& v1, PositiveValue<F> const& v2) {
        return cast_positive(static_cast<Value<F>const&>(v1)*static_cast<Value<F>const&>(v2)); }

    friend Positive<Value<F>> hlf(Positive<Value<F>> const&);
    friend Positive<Bounds<F>> pow(Positive<Value<F>> const& x, Nat m) {
        return pow(Positive<Bounds<F>>(x),m); }
    friend Positive<Bounds<F>> pow(Positive<Bounds<F>> const& x, Int n) {
        return pow(Positive<Bounds<F>>(x),n); }
};

template<class F> inline PositiveValue<F> cast_positive(Value<F> const& x) {
    return PositiveValue<F>(x); }

static_assert(IsSame<decltype(declval<FloatDPValue>() < declval<Rational>()),Boolean>::value,"");

extern template Ariadne::Nat Ariadne::Value<Ariadne::FloatDP>::output_places;
extern template Ariadne::Nat Ariadne::Value<Ariadne::FloatMP>::output_places;


}

#endif
