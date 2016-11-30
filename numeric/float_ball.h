/***************************************************************************
 *            float_ball.h
 *
 *  Copyright 2008-16  Pieter Collins
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

/*! \file float_ball.h
 *  \brief Floating-point approximations to real numbers.
 */

#ifndef ARIADNE_FLOAT_BALL_H
#define ARIADNE_FLOAT_BALL_H

#include "utility/macros.h"

#include "number.decl.h"
#include "float.decl.h"

namespace Ariadne {

template<class PR> struct NumericTraits<FloatBall<PR>> {
    typedef ValidatedNumber GenericType;
    typedef PositiveFloatBall<PR> PositiveType;
    typedef ValidatedKleenean LessType;
    typedef ValidatedKleenean EqualsType;
};

//! \ingroup NumericModule
//! \brief Floating point approximations to a real number with guaranteed error bounds.
template<class PR> class FloatBall
    : public DispatchFloatOperations<FloatBall<PR>>
    , public ProvideConvertedFieldOperations<FloatBounds<PR>,FloatBall<PR>>
    , public ProvideConvertedFieldOperations<FloatBall<PR>,FloatValue<PR>>
{
    typedef MetricTag P; typedef RawFloat<PR> FLT;
  public:
    typedef MetricTag Paradigm;
    typedef FloatBall<PR> NumericType;
    typedef ValidatedNumber GenericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    FloatBall<PR>() : _v(0.0), _e(0.0) { }
    explicit FloatBall<PR>(PrecisionType pr) : _v(0.0,pr), _e(0.0,pr) { }
    explicit FloatBall<PR>(RawFloatType const& v) : _v(v), _e(0.0) { }
    explicit FloatBall<PR>(RawFloatType const& v, RawFloatType const& e) : _v(v), _e(e) { }
    FloatBall<PR>(FloatValue<PR> const& value, FloatError<PR> const& error);
    FloatBall<PR>(FloatLowerBound<PR> const& lower, FloatUpperBound<PR> const& upper) =  delete;

    FloatBall<PR>(ExactDouble d, PR pr);
        FloatBall<PR>(const Integer& z, PR pr);
        FloatBall<PR>(const Dyadic& w, PR pr);
        FloatBall<PR>(const Rational& q, PR pr);
        FloatBall<PR>(const Real& r, PR pr);
        FloatBall<PR>(const FloatBall<PR>& x, PR pr);
    FloatBall<PR>(const ValidatedNumber& y, PR pr);

    explicit FloatBall<PR>(FloatBounds<PR> const& x);
    FloatBall<PR>(FloatValue<PR> const& x);

    FloatBall<PR>& operator=(const ValidatedNumber& y);

    operator ValidatedNumber () const;

    FloatBall<PR> create(const ValidatedNumber& y) const;

    FloatLowerBound<PR> const lower() const;
    FloatUpperBound<PR> const upper() const;
    FloatValue<PR> const value() const;
    FloatError<PR> const error() const;

    RawFloatType const lower_raw() const { return sub_down(_v,_e); }
    RawFloatType const upper_raw() const { return add_up(_v,_e); }
    RawFloatType const& value_raw() const { return _v; }
    RawFloatType const& error_raw() const { return _e; }
    double get_d() const { return _v.get_d(); }

    PrecisionType precision() const { return _v.precision(); }
    FloatBall<PR> pm(FloatError<PR> e) const;
    friend FloatApproximation<PR> round(FloatApproximation<PR> const& x);
  public:
    friend PositiveFloatUpperBound<PR> mag(FloatBall<PR> const&);
    friend PositiveFloatLowerBound<PR> mig(FloatBall<PR> const&);
    friend Bool same(FloatBall<PR> const&, FloatBall<PR> const&);
    friend Bool same(FloatBall<PR> const&, FloatBall<PR> const&);
    friend Bool models(FloatBall<PR> const&, FloatValue<PR> const&);
    friend Bool consistent(FloatBall<PR> const&, FloatBall<PR> const&);
    friend Bool refines(FloatBall<PR> const&, FloatBall<PR> const&);
    friend FloatBall<PR> refinement(FloatBall<PR> const&, FloatBall<PR> const&);
  private: public:
    RawFloatType _v, _e;
};

template<class PR> class Positive<FloatBall<PR>> : public FloatBall<PR> {
  public:
    Positive<FloatBall<PR>>() : FloatBounds<PR>() { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy>
        Positive<FloatBall<PR>>(M m, PR pr) : FloatBall<PR>(m,pr) { }
    explicit Positive<FloatBall<PR>>(FloatBall<PR> const& x) : FloatBall<PR>(x) { }
};

template<class PR> inline PositiveFloatBall<PR> cast_positive(FloatBall<PR> const& x) {
    return PositiveFloatBall<PR>(x); }

}

#endif
