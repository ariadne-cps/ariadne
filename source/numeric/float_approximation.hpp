/***************************************************************************
 *            numeric/float_approximation.hpp
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

/*! \file numeric/float_approximation.hpp
 *  \brief Floating-point approximations to real numbers.
 */

#ifndef ARIADNE_FLOAT_APPROXIMATION_HPP
#define ARIADNE_FLOAT_APPROXIMATION_HPP

#include "logical.decl.hpp"
#include "number.decl.hpp"
#include "float.decl.hpp"

#include "float_traits.hpp"
#include "float_operations.hpp"

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Floating point number approximations to real numbers supporting approxiamate arithmetic.
//! \details
//! The \c Approximation<F> class represents approximate floating-point numbers.
//! Operations are performed approximately, with no guarantees on the output.
//! \sa Real, FloatDP , FloatMP, FloatValue, FloatBall, FloatBounds.
template<class F> class Approximation
    : public DefineFieldOperators<Approximation<F>>
    , public DefineConcreteGenericOperators<Approximation<F>>
    , public DefineComparisonOperators<Approximation<F>,LessTrait<Approximation<F>>,EqualsTrait<Approximation<F>>>
{
  protected:
    typedef ApproximateTag P; typedef typename F::RoundingModeType RND; typedef typename F::PrecisionType PR;
  public:
    typedef ApproximateTag Paradigm;
    typedef Approximation<F> NumericType;
    typedef ApproximateNumber GenericType;
    typedef F RawType;
    typedef PR PrecisionType;
    typedef PR PropertiesType;
  public:
    Approximation<F>() : _a(0.0) { }
    explicit Approximation<F>(PrecisionType pr) : _a(0.0,pr) { }
    explicit Approximation<F>(RawType const& a) : _a(a) { }

        Approximation<F>(double d, PR pr) : _a(d,near,pr) { }
        Approximation<F>(ExactDouble const& d, PR pr) : _a(d,pr) { }
        Approximation<F>(TwoExp const& t, PR pr) :_a(t,pr) { }
        Approximation<F>(const Integer& z, PR pr) : _a(z,near,pr) { }
        Approximation<F>(const Dyadic& w, PR pr) : _a(w,near,pr) { }
        Approximation<F>(const Decimal& d, PR pr) : _a(d,near,pr) { }
        Approximation<F>(const Rational& q, PR pr) : _a(q,near,pr) { }
        Approximation<F>(const Real& r, PR pr); // : _a(r.get(pr)) { }
        Approximation<F>(const Approximation<F>& x, PR pr) : _a(x.raw(),near,pr) { }
    Approximation<F>(const ApproximateNumber& y, PR pr); // : _a(y.get(ApproximateTag(),pr)) { }

    template<class FF, EnableIf<IsConstructible<F,FF,RND,PR>> =dummy>
        Approximation<F>(const Approximation<FF>& x, PR pr) : _a(x.raw(),near,pr) { }

    Approximation<F>(Error<F> const& x); // FIXME: Remove
    Approximation<F>(Value<F> const& x);
    template<class FE> Approximation<F>(Ball<F,FE> const& x);
    Approximation<F>(Bounds<F> const& x);
    Approximation<F>(UpperBound<F> const& x);
    Approximation<F>(LowerBound<F> const& x);

    template<class N, EnableIf<IsBuiltinIntegral<N>> =dummy> Approximation<F>& operator=(N n) { this->_a=n; return *this; }
    template<class D, EnableIf<IsBuiltinFloatingPoint<D>> =dummy> Approximation<F>& operator=(D x) { this->_a=x; return *this; }
        Approximation<F>& operator=(const LowerBound<F>& x) { return *this=Approximation<F>(x); }
        Approximation<F>& operator=(const UpperBound<F>& x) { return *this=Approximation<F>(x); }
        Approximation<F>& operator=(const Bounds<F>& x) { return *this=Approximation<F>(x); }
        Approximation<F>& operator=(const Value<F>& x) { return *this=Approximation<F>(x); }
    Approximation<F>& operator=(const ApproximateNumber& y) { return *this=Approximation<F>(y,this->precision()); }
    Approximation<F> create(const ApproximateNumber& y) const { return Approximation<F>(y,this->precision()); }

    operator ApproximateNumber () const;

    friend Approximation<F> round(Approximation<F> const& x);

    PrecisionType precision() const { return _a.precision(); }
    PropertiesType properties() const { return _a.precision(); }
    GenericType generic() const { return this->operator GenericType(); }
    explicit operator RawType () const { return this->_a; }
    RawType const& raw() const { return this->_a; }
    RawType& raw() { return this->_a; }
    double get_d() const { return this->_a.get_d(); }
  public:
    friend Approximation<F> floor(Approximation<F> const& x) {
        return Approximation<F>(floor(x._a)); }
    friend Approximation<F> ceil(Approximation<F> const& x) {
        return Approximation<F>(ceil(x._a)); }
    friend Approximation<F> round(Approximation<F> const& x) {
        return Approximation<F>(round(x._a)); }

    friend Approximation<F> abs(Approximation<F> const& x) {
        return Approximation<F>(abs(x._a)); }
    friend Approximation<F> max(Approximation<F> const& x, Approximation<F> const& y) {
        return Approximation<F>(max(x._a,y._a)); }
    friend Approximation<F> min(Approximation<F> const& x, Approximation<F> const& y) {
        return Approximation<F>(min(x._a,y._a)); }
    friend PositiveApproximation<F> mag(Approximation<F> const& x) {
        return PositiveApproximation<F>(abs(x._a)); }
    friend PositiveApproximation<F> mig(Approximation<F> const& x) {
        return PositiveApproximation<F>(abs(x._a)); }

    friend Approximation<F> nul(Approximation<F> const& x) {
        return Approximation<F>(nul(x._a)); }
    friend Approximation<F> pos(Approximation<F> const& x) {
        return Approximation<F>(pos(x._a)); }
    friend Approximation<F> neg(Approximation<F> const& x) {
        return Approximation<F>(neg(x._a)); }
    friend Approximation<F> hlf(Approximation<F> const& x) {
        return Approximation<F>(hlf(x._a)); }
    friend Approximation<F> sqr(Approximation<F> const& x) {
        return Approximation<F>(mul(near,x._a,x._a)); }
    friend Approximation<F> rec(Approximation<F> const& x) {
        return Approximation<F>(div(near,1.0,x._a)); }

    friend Approximation<F> add(Approximation<F> const& x1, Approximation<F> const& x2) {
        return Approximation<F>(add(near,x1._a,x2._a)); }
    friend Approximation<F> sub(Approximation<F> const& x1, Approximation<F> const& x2) {
        return Approximation<F>(sub(near,x1._a,x2._a)); }
    friend Approximation<F> mul(Approximation<F> const& x1, Approximation<F> const& x2) {
        return Approximation<F>(mul(near,x1._a,x2._a)); }
    friend Approximation<F> div(Approximation<F> const& x1, Approximation<F> const& x2) {
        return Approximation<F>(div(near,x1._a,x2._a)); }
    friend Approximation<F> fma(Approximation<F> const& x1, Approximation<F> const& x2, Approximation<F> const& x3) {
        return Approximation<F>(fma(near,x1._a,x2._a,x3._a)); }

    friend Approximation<F> pow(Approximation<F> const& x, Nat m) {
        return Approximation<F>(pow(approx,x._a,static_cast<Int>(m))); }
    friend Approximation<F> pow(Approximation<F> const& x, Int n) {
        return Approximation<F>(pow(approx,x._a,n)); }

    friend Approximation<F> sqrt(Approximation<F> const& x) {
        return Approximation<F>(sqrt(approx,x._a)); }
    friend Approximation<F> exp(Approximation<F> const& x) {
        return Approximation<F>(exp(approx,x._a)); }
    friend Approximation<F> log(Approximation<F> const& x) {
        return Approximation<F>(log(approx,x._a)); }
    friend Approximation<F> sin(Approximation<F> const& x) {
        return Approximation<F>(sin(approx,x._a)); }
    friend Approximation<F> cos(Approximation<F> const& x) {
        return Approximation<F>(cos(approx,x._a)); }
    friend Approximation<F> tan(Approximation<F> const& x) {
        return Approximation<F>(tan(approx,x._a)); }
    friend Approximation<F> asin(Approximation<F> const& x) {
        return Approximation<F>(asin(approx,x._a)); }
    friend Approximation<F> acos(Approximation<F> const& x) {
        return Approximation<F>(acos(approx,x._a)); }
    friend Approximation<F> atan(Approximation<F> const& x) {
        return Approximation<F>(atan(approx,x._a)); }

    friend ApproximateKleenean eq(Approximation<F> const& x1, Approximation<F> const& x2) {
        return x1._a==x2._a; }
    friend ApproximateKleenean lt(Approximation<F> const& x1, Approximation<F> const& x2) {
        return x1._a< x2._a; }

    friend Bool same(Approximation<F> const& x1, Approximation<F> const& x2) {
        return x1._a==x2._a; }

    friend Integer cast_integer(Approximation<F> const& x) {
        return round(static_cast<Dyadic>(x._a)); }

    friend OutputStream& operator<<(OutputStream& os, Approximation<F> const& x) {
        return write(os,x.raw(),DecimalPrecision{Approximation<F>::output_places},to_nearest); }
    friend InputStream& operator>>(InputStream& is, Approximation<F>& x) {
        is >> x._a; return is; }
  public:
    static Nat output_places;
    static Void set_output_places(Nat p) { output_places=p; }
    Approximation<F> pm(Approximation<F> _e) { return *this; }
  public:
    RawType _a;
};

template<class F> template<class FE> Approximation<F>::Approximation(Ball<F,FE> const& x) : Approximation<F>(x.value_raw()) { }

template<class PR> Approximation(ApproximateNumber, PR) -> Approximation<RawFloatType<PR>>;
template<class F> Approximation(F) -> Approximation<F>;

template<class F> inline FloatFactory<PrecisionType<F>> factory(Approximation<F> const& flt) { return FloatFactory<PrecisionType<F>>(flt.precision()); }
template<class PR> inline FloatApproximation<PR> FloatFactory<PR>::create(ApproximateNumber const& y) { return FloatApproximation<PR>(y,_pr); }
template<class PR> inline PositiveFloatApproximation<PR> FloatFactory<PR>::create(PositiveApproximateNumber const& y) { return PositiveFloatApproximation<PR>(y,_pr); }
template<class PR> template<class D, EnableIf<IsBuiltinFloatingPoint<D>>> inline
    FloatApproximation<PR> FloatFactory<PR>::create(D const& y) { return FloatApproximation<PR>(y,_pr); }

template<class F> class Positive<Approximation<F>> : public Approximation<F>
    , public DeclarePositiveFloatOperations<PositiveApproximation<F>>
    , public DefineSemiFieldOperators<PositiveApproximation<F>>
    , public DefineConcreteGenericOperators<PositiveApproximation<F>>
{
    using typename Approximation<F>::PR;
  public:
    Positive<Approximation<F>>() : Approximation<F>() { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy>
        Positive<Approximation<F>>(M m) : Approximation<F>(m) { }
    explicit Positive<Approximation<F>>(PR const& pr) : Approximation<F>(pr) { }
    explicit Positive<Approximation<F>>(F const& x) : Approximation<F>(x) { }
    explicit Positive<Approximation<F>>(Approximation<F> const& x) : Approximation<F>(x) { }
    Positive<Approximation<F>>(PositiveApproximateNumber const& y, PR pr) : Approximation<F>(y,pr) { }
    Positive<Approximation<F>>(PositiveLowerBound<F> const& x) : Approximation<F>(x) { }
    Positive<Approximation<F>>(PositiveUpperBound<F> const& x) : Approximation<F>(x) { }
    Positive<Approximation<F>>(PositiveValue<F> const& x) : Approximation<F>(x) { }
    Positive<Approximation<F>>(Error<F> const& x) : Approximation<F>(x) { }
  public:
    friend PositiveApproximation<F> nul(PositiveApproximation<F> const& x) { return PositiveApproximation<F>(nul(x.raw())); }
    friend PositiveApproximation<F> hlf(PositiveApproximation<F> const& x) { return PositiveApproximation<F>(hlf(x.raw())); }
    friend PositiveApproximation<F> sqr(PositiveApproximation<F> const& x) { return PositiveApproximation<F>(sqr(near,x.raw())); }
    friend PositiveApproximation<F> rec(PositiveApproximation<F> const& x) { return PositiveApproximation<F>(rec(near,x.raw())); }
    friend PositiveApproximation<F> add(PositiveApproximation<F> const& x1, PositiveApproximation<F> const& x2) {
        return PositiveApproximation<F>(add(near,x1.raw(),x2.raw())); }
    friend PositiveApproximation<F> mul(PositiveApproximation<F> const& x1, PositiveApproximation<F> const& x2) {
        return PositiveApproximation<F>(mul(near,x1.raw(),x2.raw())); }
    friend PositiveApproximation<F> div(PositiveApproximation<F> const& x1, PositiveApproximation<F> const& x2) {
        return PositiveApproximation<F>(div(near,x1.raw(),x2.raw())); }
    friend PositiveApproximation<F> pow(PositiveApproximation<F> const& x, Nat m) { return PositiveApproximation<F>(pow(near,x.raw(),static_cast<Int>(m))); }
    friend PositiveApproximation<F> pow(PositiveApproximation<F> const& x, Int n) { return PositiveApproximation<F>(pow(near,x.raw(),n)); }
    friend PositiveApproximation<F> sqrt(PositiveApproximation<F> const& x) { return PositiveApproximation<F>(sqrt(near,x.raw())); }
    friend PositiveApproximation<F> exp(PositiveApproximation<F> const& x) { return PositiveApproximation<F>(exp(near,x.raw())); }
    friend Approximation<F> log(PositiveApproximation<F> const& x) { return Approximation<F>(log(near,x.raw())); }
    friend PositiveApproximation<F> atan(PositiveApproximation<F> const& x) { return PositiveApproximation<F>(atan(near,x.raw())); }
    friend PositiveApproximation<F> max(PositiveApproximation<F> const& x1, PositiveApproximation<F> const& x2) {
        return PositiveApproximation<F>(max(near,x1.raw(),x2.raw())); }
    friend PositiveApproximation<F> max(PositiveApproximation<F> const& x1, Approximation<F> const& x2) {
        return PositiveApproximation<F>(max(near,x1.raw(),x2.raw())); }
    friend PositiveApproximation<F> max(Approximation<F> const& x1, PositiveApproximation<F> const& x2) {
        return PositiveApproximation<F>(max(near,x1.raw(),x2.raw())); }
    friend PositiveApproximation<F> min(PositiveApproximation<F> const& x1, PositiveApproximation<F> const& x2) {
        return PositiveApproximation<F>(min(near,x1.raw(),x2.raw())); }
    friend PositiveApproximation<F> abs(PositiveApproximation<F> const& x) { return PositiveApproximation<F>(abs(x)); }
};

template<class F> inline PositiveApproximation<F> cast_positive(Approximation<F> const& x) {
    return PositiveApproximation<F>(x); }

extern template Ariadne::Nat Ariadne::Approximation<Ariadne::FloatDP>::output_places;
extern template Ariadne::Nat Ariadne::Approximation<Ariadne::FloatMP>::output_places;

}

#endif
