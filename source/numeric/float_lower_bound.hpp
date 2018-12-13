/***************************************************************************
 *            float_lower_bound.hpp
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

/*! \file float_lower_bound.hpp
 *  \brief Floating-point lower bounds for real numbers.
 */

#ifndef FLOAT_LOWER_BOUND_H
#define FLOAT_LOWER_BOUND_H

#include "../utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"

namespace Ariadne {

template<class F> struct NumericTraits<LowerBound<F>> {
    typedef ValidatedLowerNumber GenericType;
    typedef UpperBound<F> OppositeType;
    typedef PositiveLowerBound<F> PositiveType;
    typedef ValidatedUpperKleenean LessType;
    typedef ValidatedNegatedSierpinskian EqualsType;
};

//! \ingroup NumericModule
//! \brief Floating-point lower bounds for real numbers.
//! \sa UpperReal, FloatDP, FloatMP, FloatBounds, FloatUpperBound.
template<class F> class LowerBound
    : public DispatchDirectedFloatOperations<LowerBound<F>>
    , public DispatchFloatOperations<Approximation<F>>
{
  protected:
    typedef LowerTag P; typedef typename F::PrecisionType PR;
  public:
    typedef LowerTag Paradigm;
    typedef LowerBound<F> NumericType;
    typedef ValidatedLowerNumber GenericType;
    typedef F RawType;
    typedef PR PrecisionType;
    typedef PR PropertiesType;
  public:
    LowerBound<F>() : _l(0.0) { }
    explicit LowerBound<F>(PrecisionType pr) : _l(0.0,pr) { }
    explicit LowerBound<F>(RawType const& l) : _l(l) { }

    template<class N, EnableIf<IsBuiltinIntegral<N>> = dummy> LowerBound<F>(N n, PR pr) : LowerBound<F>(ExactDouble(n),pr) { }
    LowerBound<F>(ExactDouble d, PR pr);
        LowerBound<F>(TwoExp t, PR pr);
        LowerBound<F>(const Integer& z, PR pr);
        LowerBound<F>(const Dyadic& w, PR pr);
        LowerBound<F>(const Decimal& d, PR pr);
        LowerBound<F>(const Rational& q, PR pr);
        LowerBound<F>(const Real& r, PR pr);
    LowerBound<F>(const LowerBound<F>& x, PR pr);
    LowerBound<F>(const ValidatedLowerNumber& y, PR pr);

    LowerBound<F>(Bounds<F> const& x);
    LowerBound<F>(Ball<F> const& x);
    LowerBound<F>(Value<F> const& x);

        LowerBound<F>& operator=(const Value<F>& x) { return *this=LowerBound<F>(x); }
    LowerBound<F>& operator=(const ValidatedLowerNumber&);
    LowerBound<F> create(const ValidatedLowerNumber& y) const;
    UpperBound<F> create(const ValidatedUpperNumber& y) const;

    operator ValidatedLowerNumber () const;

    PrecisionType precision() const { return _l.precision(); }
    PropertiesType properties() const { return _l.precision(); }
    GenericType generic() const { return this->operator GenericType(); }
    RawType const& raw() const { return _l; }
    RawType& raw() { return _l; }
    double get_d() const { return _l.get_d(); }
  public: // To be removed
    friend Bool same(LowerBound<F> const&, LowerBound<F> const&);
    friend Bool refines(LowerBound<F> const&, LowerBound<F> const&);
    friend LowerBound<F> refinement(LowerBound<F> const&, LowerBound<F> const&);
  public:
    friend LowerBound<F> operator*(PositiveBounds<F> const& x1, LowerBound<F> const& x2) {
        return LowerBound<F>(mul(down,x2.raw()>=0?x1.lower().raw():x1.upper().raw(),x2.raw())); }
    friend LowerBound<F> operator*(LowerBound<F> const& x1, PositiveBounds<F> const& x2) {
        return LowerBound<F>(mul(down,x1.raw(),x1.raw()>=0?x2.lower().raw():x2.upper().raw())); }
    friend LowerBound<F> operator/(LowerBound<F> const& x1, PositiveBounds<F> const& x2) {
        return LowerBound<F>(div(down,x1.raw(),x1.raw()>=0?x2.upper().raw():x2.lower().raw())); }
    // Needed to prevent ambiguity; useful as implementation is easier than PositiveBounds version.
    friend LowerBound<F> operator*(PositiveValue<F> const& x1, LowerBound<F> const& x2) {
        return LowerBound<F>(mul(down,x1.raw(),x2.raw())); }
    friend LowerBound<F> operator*(LowerBound<F> const& x1, PositiveValue<F> const& x2) {
        return LowerBound<F>(mul(down,x1.raw(),x2.raw())); }
    friend LowerBound<F> operator/(LowerBound<F> const& x1, PositiveValue<F> const& x2) {
        return LowerBound<F>(div(down,x1.raw(),x2.raw())); }
  private: public:
    static Nat output_places;
    RawType _l;
};

template<class PR> LowerBound(ValidatedLowerNumber, PR) -> LowerBound<RawFloatType<PR>>;
template<class F> LowerBound(F) -> LowerBound<F>;

template<class F> inline FloatFactory<PrecisionType<F>> factory(LowerBound<F> const& flt) { return FloatFactory<PrecisionType<F>>(flt.precision()); }
template<class PR> inline FloatLowerBound<PR> FloatFactory<PR>::create(Number<LowerTag> const& y) { return FloatLowerBound<PR>(y,_pr); }

template<class F> class Positive<LowerBound<F>> : public LowerBound<F>
    , public DispatchPositiveDirectedFloatOperations<PositiveLowerBound<F>,PositiveUpperBound<F>>
{
    using typename LowerBound<F>::PR;
  public:
    Positive<LowerBound<F>>() : LowerBound<F>() { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy>
        Positive<LowerBound<F>>(M m) : LowerBound<F>(m) { }
    explicit Positive<LowerBound<F>>(F const& x) : LowerBound<F>(x) { }
    explicit Positive<LowerBound<F>>(LowerBound<F> const& x) : LowerBound<F>(x) { }
    explicit Positive<LowerBound<F>>(ValidatedLowerNumber const& y, PR pr) : LowerBound<F>(y,pr) { }
    Positive<LowerBound<F>>(PositiveValue<F> const& x) : LowerBound<F>(x) { }
    Positive<LowerBound<F>>(PositiveBounds<F> const& x) : LowerBound<F>(x) { }
  public:
    // Needed to prevent ambiguity
    friend PositiveLowerBound<F> operator*(PositiveBounds<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveLowerBound<F>(mul(down,x1.lower().raw(),x2.raw())); }
    friend PositiveLowerBound<F> operator*(PositiveLowerBound<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveLowerBound<F>(mul(down,x1.raw(),x2.lower().raw())); }
    friend PositiveLowerBound<F> operator/(PositiveLowerBound<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveLowerBound<F>(div(down,x1.raw(),x2.upper().raw())); }
    friend PositiveLowerBound<F> operator*(PositiveValue<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveLowerBound<F>(mul(down,x1.raw(),x2.raw())); }
    friend PositiveLowerBound<F> operator*(PositiveLowerBound<F> const& x1, PositiveValue<F> const& x2) {
        return PositiveLowerBound<F>(mul(down,x1.raw(),x2.raw())); }
    friend PositiveLowerBound<F> operator/(PositiveLowerBound<F> const& x1, PositiveValue<F> const& x2) {
        return PositiveLowerBound<F>(div(down,x1.raw(),x2.raw())); }
};

template<class F> inline PositiveLowerBound<F> cast_positive(LowerBound<F> const& x) {
    return PositiveLowerBound<F>(x); }

}

#endif
