/***************************************************************************
 *            float_upper_bound.hpp
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

/*! \file float_upper_bound.hpp
 *  \brief Floating-point upper bounds for real numbers.
 */

#ifndef ARIADNE_FLOAT_UPPER_BOUND_HPP
#define ARIADNE_FLOAT_UPPER_BOUND_HPP

#include "../utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"

namespace Ariadne {

template<class F> struct NumericTraits<UpperBound<F>> {
    typedef ValidatedUpperNumber GenericType;
    typedef LowerBound<F> OppositeType;
    typedef PositiveUpperBound<F> PositiveType;
    typedef ValidatedLowerKleenean LessType;
    typedef ValidatedNegatedSierpinskian EqualsType;
};

//! \ingroup NumericModule
//! \brief Floating-point upper bounds for real numbers.
//! \sa UpperReal, FloatDP, FloatMP, FloatBounds, FloatLowerBound.
template<class F> class UpperBound
    : public DispatchDirectedFloatOperations<UpperBound<F>>
    , public DispatchFloatOperations<Approximation<F>>
{
  protected:
    typedef UpperTag P; typedef typename F::PrecisionType PR;
  public:
    typedef UpperTag Paradigm;
    typedef UpperBound<F> NumericType;
    typedef ValidatedUpperNumber GenericType;
    typedef F RawType;
    typedef PR PrecisionType;
    typedef PR PropertiesType;
  public:
    UpperBound<F>() : _u(0.0) { }
    explicit UpperBound<F>(PrecisionType pr) : _u(0.0,pr) { }
    explicit UpperBound<F>(RawType const& u) : _u(u) { }

    template<class N, EnableIf<IsBuiltinIntegral<N>> = dummy> UpperBound<F>(N n, PR pr) : UpperBound<F>(ExactDouble(n),pr) { }
    UpperBound<F>(ExactDouble d, PR pr);
        UpperBound<F>(TwoExp t, PR pr);
        UpperBound<F>(const Integer& z, PR pr);
        UpperBound<F>(const Dyadic& w, PR pr);
        UpperBound<F>(const Decimal& d, PR pr);
        UpperBound<F>(const Rational& q, PR pr);
        UpperBound<F>(const Real& r, PR pr);
    UpperBound<F>(const UpperBound<F>& x, PR pr);
    UpperBound<F>(const ValidatedUpperNumber& y, PR pr);

    UpperBound<F>(Bounds<F> const& x);
    UpperBound<F>(Ball<F> const& x);
    UpperBound<F>(Value<F> const& x);
    UpperBound<F>(Error<F> const& x); // FIXME: Remove

        UpperBound<F>& operator=(const Value<F>& x) { return *this=UpperBound<F>(x); }
    UpperBound<F>& operator=(const ValidatedUpperNumber& y);
    UpperBound<F> create(const ValidatedUpperNumber& y) const;
    LowerBound<F> create(const ValidatedLowerNumber& y) const;

    operator ValidatedUpperNumber () const;

    PrecisionType precision() const { return _u.precision(); }
    PropertiesType properties() const { return _u.precision(); }
    GenericType generic() const { return this->operator GenericType(); }
    RawType const& raw() const { return _u; }
    RawType& raw() { return _u; }
    double get_d() const { return _u.get_d(); }
  public: // To be removed
    friend Bool same(UpperBound<F> const&, UpperBound<F> const&);
    friend Bool refines(UpperBound<F> const&, UpperBound<F> const&);
    friend UpperBound<F> refinement(UpperBound<F> const&, UpperBound<F> const&);
  public:
    friend UpperBound<F> operator*(PositiveBounds<F> const& x1, UpperBound<F> const& x2) {
        return UpperBound<F>(mul(up,x2.raw()>=0?x1.upper().raw():x1.lower().raw(),x2.raw())); }
    friend UpperBound<F> operator*(UpperBound<F> const& x1, PositiveBounds<F> const& x2) {
        return UpperBound<F>(mul(up,x1.raw(),x1.raw()>=0?x2.upper().raw():x2.lower().raw())); }
    friend UpperBound<F> operator/(UpperBound<F> const& x1, PositiveBounds<F> const& x2) {
        return UpperBound<F>(div(up,x1.raw(),x1.raw()>=0?x2.lower().raw():x2.upper().raw())); }
    // Needed to prevent ambiguity; useful as implementation is easier than PositiveBounds version.
    friend UpperBound<F> operator*(PositiveValue<F> const& x1, UpperBound<F> const& x2) {
        return UpperBound<F>(mul(up,x1.raw(),x2.raw())); }
    friend UpperBound<F> operator*(UpperBound<F> const& x1, PositiveValue<F> const& x2) {
        return UpperBound<F>(mul(up,x1.raw(),x2.raw())); }
    friend UpperBound<F> operator/(UpperBound<F> const& x1, PositiveValue<F> const& x2) {
        return UpperBound<F>(div(up,x1.raw(),x2.raw())); }
  private: public:
    static Nat output_places;
    RawType _u;
};

template<class PR> UpperBound(ValidatedUpperNumber, PR) -> UpperBound<RawFloatType<PR>>;
template<class F> UpperBound(F) -> UpperBound<F>;

template<class F> inline FloatFactory<PrecisionType<F>> factory(UpperBound<F> const& flt) { return FloatFactory<PrecisionType<F>>(flt.precision()); }
template<class PR> inline FloatUpperBound<PR> FloatFactory<PR>::create(Number<UpperTag> const& y) { return FloatUpperBound<PR>(y,_pr); }

template<class F> class Positive<UpperBound<F>> : public UpperBound<F>
    , public DispatchPositiveDirectedFloatOperations<PositiveUpperBound<F>,PositiveLowerBound<F>>
{
    using typename UpperBound<F>::PR;
  public:
    Positive<UpperBound<F>>() : UpperBound<F>() { }
    explicit Positive<UpperBound<F>>(PR const& pr) : UpperBound<F>(pr) { }
    explicit Positive<UpperBound<F>>(F const& x) : UpperBound<F>(x) { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy> Positive<UpperBound<F>>(M m, PR pr) : UpperBound<F>(m,pr) { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy> PositiveValue<F> create(M m) const { return PositiveValue<F>(m,this->precision()); }
    explicit Positive<UpperBound<F>>(UpperBound<F> const& x) : UpperBound<F>(x) { ARIADNE_PRECONDITION_MSG(!(this->_u<0),"x="<<x); }
    explicit Positive<UpperBound<F>>(ValidatedUpperNumber const& y, PR pr) : UpperBound<F>(y,pr) { ARIADNE_PRECONDITION_MSG(!(this->_u<0),"y="<<y); }
    Positive<UpperBound<F>>(PositiveValue<F> const& x) : UpperBound<F>(x) { }
    Positive<UpperBound<F>>(PositiveBounds<F> const& x) : UpperBound<F>(x) { }
  public:
    // Needed to prevent ambiguity
    friend PositiveUpperBound<F> operator*(PositiveBounds<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveUpperBound<F>(mul(up,x1.upper().raw(),x2.raw())); }
    friend PositiveUpperBound<F> operator*(PositiveUpperBound<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveUpperBound<F>(mul(up,x1.raw(),x2.upper().raw())); }
    friend PositiveUpperBound<F> operator/(PositiveUpperBound<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveUpperBound<F>(div(up,x1.raw(),x2.lower().raw())); }
    friend PositiveUpperBound<F> operator*(PositiveValue<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveUpperBound<F>(mul(up,x1.raw(),x2.raw())); }
    friend PositiveUpperBound<F> operator*(PositiveUpperBound<F> const& x1, PositiveValue<F> const& x2) {
        return PositiveUpperBound<F>(mul(up,x1.raw(),x2.raw())); }
    friend PositiveUpperBound<F> operator/(PositiveUpperBound<F> const& x1, PositiveValue<F> const& x2) {
        return PositiveUpperBound<F>(div(up,x1.raw(),x2.raw())); }
};

template<class F> inline PositiveUpperBound<F> cast_positive(UpperBound<F> const& x) {
    return PositiveUpperBound<F>(x); }

}

#endif
