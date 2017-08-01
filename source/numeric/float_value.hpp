/***************************************************************************
 *            float_value.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file float_value.hpp
 *  \brief Exact Floating-point representations of real numbers.
 */

#ifndef ARIADNE_FLOAT_VALUE_HPP
#define ARIADNE_FLOAT_VALUE_HPP

#include "utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"

#include "float_operations.hpp"

#include "logical.hpp"
#include "builtin.hpp"
#include "twoexp.hpp"

namespace Ariadne {

template<class PR> struct NumericTraits<FloatValue<PR>> {
    typedef ExactNumber GenericType;
    typedef PositiveFloatValue<PR> PositiveType;
    typedef Boolean LessType;
    typedef Boolean EqualsType;
};

static_assert(not IsGenericNumericType<FloatValue<Precision64>>::value,"");
static_assert(not IsGenericNumericType<FloatValue<PrecisionMP>>::value,"");

//! \ingroup NumericModule
//! \brief A floating-point number, which is taken to represent the \em exact value of a real quantity.
//! \sa Float64 , FloatMP, FloatBall, FloatBounds, FloatApproximation.
template<class PR> class FloatValue
    : DispatchNumericOperations<FloatValue<PR>,FloatBounds<PR>>
    , DispatchComparisonOperations<FloatValue<PR>,Boolean>
    , DefineMixedComparisonOperators<FloatValue<PR>,ExactNumber,Boolean>
    , DefineMixedComparisonOperators<FloatValue<PR>,Rational,Boolean>
//    , DefineMixedComparisonOperators<FloatValue<PR>,Dyadic,Boolean>
//    , DefineMixedComparisonOperators<FloatValue<PR>,Integer,Boolean>
//    , DefineMixedComparisonOperators<FloatValue<PR>,Int,Boolean>
//        , public DispatchFloatOperations<FloatBall<PR>>
        , public DispatchFloatOperations<FloatBounds<PR>>
    , DefineConcreteGenericArithmeticOperators<FloatValue<PR>>
    , DefineConcreteGenericComparisonOperators<FloatValue<PR>>
{
    typedef ExactTag P; typedef RawFloat<PR> FLT;
  public:
    typedef ExactTag Paradigm;
    typedef FloatValue<PR> NumericType;
    typedef ExactNumber GenericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
    typedef PR PropertiesType;
  public:
    FloatValue<PR>() : _v(0.0) { }
    explicit FloatValue<PR>(PrecisionType pr) : _v(0.0,pr) { }
    explicit FloatValue<PR>(RawFloatType const& v) : _v(v) { }

    template<class N, EnableIf<IsBuiltinIntegral<N>> = dummy> FloatValue<PR>(N n, PR pr) : FloatValue<PR>(ExactDouble(n),pr) { }
    FloatValue<PR>(ExactDouble d, PR pr);
    FloatValue<PR>(const Integer& z, PR pr);
    FloatValue<PR>(const TwoExp& t, PR pr);
    FloatValue<PR>(const Dyadic& w, PR pr);
    FloatValue<PR>(const FloatValue<PR>& x, PR pr);

    template<class N, EnableIf<IsBuiltinIntegral<N>> = dummy> FloatValue<PR>& operator=(N n) { _v=n; return *this; }
    FloatValue<PR>& operator=(const Integer& z);
    FloatValue<PR>& operator=(const TwoExp& t);
    FloatValue<PR>& operator=(const Dyadic& w);

    operator ExactNumber () const;
    explicit operator Dyadic () const;
    explicit operator Rational () const;

    FloatBall<PR> create(ValidatedNumber const&) const;
//    explicit operator RawFloatType () const { return _v; }

    PrecisionType precision() const { return _v.precision(); }
    PropertiesType properties() const { return _v.precision(); }
    GenericType generic() const { return this->operator GenericType(); }
    RawFloatType const& raw() const { return _v; }
    RawFloatType& raw() { return _v; }
    double get_d() const { return _v.get_d(); }

    template<class PRE> FloatBall<PR,PRE> pm(FloatError<PRE> e) const;
  public:
    friend FloatValue<PR> operator*(FloatValue<PR> const&, TwoExp const&);
    friend FloatValue<PR> operator/(FloatValue<PR> const&, TwoExp const&);
    friend FloatValue<PR>& operator*=(FloatValue<PR>&, TwoExp const&);
    friend FloatValue<PR>& operator/=(FloatValue<PR>&, TwoExp const&);
    friend FloatError<PR> mag(FloatValue<PR> const&);
    friend FloatLowerBound<PR> mig(FloatValue<PR> const&);
    friend Bool same(FloatValue<PR> const&, FloatValue<PR> const&);
    friend OutputStream& operator<<(OutputStream&, FloatValue<PR> const&);
  public:
    friend Comparison cmp(FloatValue<PR> const& x1, Rational const& q2) { return cmp(x1.raw(),q2); }
    friend Comparison cmp(FloatValue<PR> const& x1, Dyadic const& w2) { return cmp(x1.raw(),w2); }
    friend Comparison cmp(FloatValue<PR> const& x1, Integer const& z2) { return cmp(x1.raw(),z2); }
  public:
    static Nat output_places;
    static Void set_output_places(Nat p) { output_places=p; }
  private: public:
    RawFloatType _v;
  private:
    friend FloatValue<PR> shft(FloatValue<PR> const& x, Int n) {
        return FloatValue<PR>(shft(x.raw(),n)); }
    friend FloatValue<PR> operator*(FloatValue<PR> const& x, TwoExp const& y) {
        return FloatValue<PR>(mul(near,x.raw(),RawFloat<PR>(y,x.precision()))); }
    friend FloatValue<PR> operator/(FloatValue<PR> const& x, TwoExp const& y) {
        return FloatValue<PR>(div(near,x.raw(),RawFloat<PR>(y,x.precision()))); }
    friend FloatValue<PR>& operator*=(FloatValue<PR>& x, TwoExp const& y) { return x=x*y; }
    friend FloatValue<PR>& operator/=(FloatValue<PR>& x, TwoExp const& y) { return x=x/y; }
    friend OutputStream& operator<<(OutputStream& os, FloatValue<PR> const& x) {
        return Operations<FloatValue<PR>>::_write(os,x); }
};

template<class PR> inline FloatFactory<PR> factory(FloatValue<PR> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class PR> inline FloatValue<PR> FloatFactory<PR>::create(Dyadic const& y, ExactTag) { return FloatValue<PR>(y,_pr); }
template<class PR> inline FloatValue<PR> FloatFactory<PR>::create(Integer const& y, ExactTag) { return FloatValue<PR>(y,_pr); }
template<class PR> inline FloatValue<PR> FloatFactory<PR>::create(ExactDouble const& y) { return FloatValue<PR>(y,_pr); }

template<class PR> template<class N, EnableIf<IsBuiltinSignedIntegral<N>>> inline
FloatValue<PR> FloatFactory<PR>::create(N const& y) { return FloatValue<PR>(y,_pr); }

template<class PR> template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>>> inline
PositiveFloatValue<PR> FloatFactory<PR>::create(M const& y) { return PositiveFloatValue<PR>(y,_pr); }

template<class PR> class Positive<FloatValue<PR>> : public FloatValue<PR> {
  public:
    Positive<FloatValue<PR>>() : FloatValue<PR>() { }
    explicit Positive<FloatValue<PR>>(PR const& pr) : FloatValue<PR>(pr) { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy> Positive<FloatValue<PR>>(M m, PR pr) : FloatValue<PR>(m,pr) { }
    Positive<FloatValue<PR>>(TwoExp const& ex, PR pr) : FloatValue<PR>(ex,pr) { }
    explicit Positive<FloatValue<PR>>(Dyadic const& w, PR pr) : FloatValue<PR>(w,pr) { }
    explicit Positive<FloatValue<PR>>(RawFloat<PR> const& x) : FloatValue<PR>(x) { }
    explicit Positive<FloatValue<PR>>(FloatValue<PR> const& x) : FloatValue<PR>(x) { }
  public:
    friend Positive<FloatValue<PR>> hlf(Positive<FloatValue<PR>> const&);
    friend Positive<FloatBounds<PR>> pow(Positive<FloatValue<PR>> const& x, Nat m) {
        return pow(Positive<FloatBounds<PR>>(x),m); }
    friend Positive<FloatBounds<PR>> pow(Positive<FloatBounds<PR>> const& x, Int n) {
        return pow(Positive<FloatBounds<PR>>(x),n); }
};

template<class PR> inline PositiveFloatValue<PR> cast_positive(FloatValue<PR> const& x) {
    return PositiveFloatValue<PR>(x); }

static_assert(IsSame<decltype(declval<Float64Value>() < declval<Rational>()),Boolean>::value,"");

}

#endif
