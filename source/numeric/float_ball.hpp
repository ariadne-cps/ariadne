/***************************************************************************
 *            float_ball.hpp
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

/*! \file float_ball.hpp
 *  \brief Floating-point approximations to real numbers.
 */

#ifndef ARIADNE_FLOAT_BALL_HPP
#define ARIADNE_FLOAT_BALL_HPP

#include "utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"

#include "float_operations.hpp"

namespace Ariadne {

template<class PR, class PRE> struct NumericTraits<FloatBall<PR,PRE>> {
    typedef ValidatedNumber GenericType;
    typedef PositiveFloatBall<PR,PRE> PositiveType;
    typedef ValidatedKleenean LessType;
    typedef ValidatedKleenean EqualsType;
};

template<class PRE, class PR, EnableIf<IsDefaultConstructible<PRE>> = dummy> PRE _error_precision(PR const&) { return PRE(); }
template<class PRE, class PR, DisableIf<IsDefaultConstructible<PRE>> = dummy, EnableIf<IsConstructible<PRE,PR>> = dummy> PRE _error_precision(PR const& pr) { return PRE(pr); }

//! \ingroup NumericModule
//! \brief Floating point approximations to a real number with guaranteed error bounds.
template<class PR, class PRE> class FloatBall
    : public DispatchFloatOperations<FloatBall<PR,PRE>>
    , public ProvideConvertedFieldOperations<FloatBounds<PR>,FloatBall<PR,PRE>>
    , public ProvideConvertedFieldOperations<FloatBall<PR,PRE>,FloatValue<PR>>
{
    static_assert(IsConstructible<PR,PRE>::value or IsDefaultConstructible<PRE>::value,"");

    typedef ValidatedTag P; typedef RawFloat<PR> FLT; typedef RawFloat<PRE> EFLT;
  public:
    typedef P Paradigm;
    typedef FloatBall<PR,PRE> NumericType;
    typedef Number<P> GenericType;
    typedef FLT RawFloatType;
    typedef EFLT RawErrorFloatType;
    typedef PR PrecisionType;
    typedef PRE ErrorPrecisionType;
    typedef PR PropertiesType;
  public:
    FloatBall<PR,PRE>() : _v(0.0), _e(0.0) { }
    explicit FloatBall<PR,PRE>(PrecisionType pr) : _v(0.0,pr), _e(0.0,_error_precision<PRE>(pr)) { }
    explicit FloatBall<PR,PRE>(PrecisionType pr, ErrorPrecisionType pre) : _v(0.0,pr), _e(0.0,pre) { }
    explicit FloatBall<PR,PRE>(RawFloat<PR> const& v) : _v(v), _e(0.0,_error_precision<PRE>(v.precision())) { }
    explicit FloatBall<PR,PRE>(RawFloat<PR> const& v, PRE pre) : _v(v), _e(0.0,pre) { }
    explicit FloatBall<PR,PRE>(RawFloat<PR> const& v, RawFloat<PRE> const& e) : _v(v), _e(e) { }
    FloatBall<PR,PRE>(FloatValue<PR> const& value, FloatError<PRE> const& error);
    FloatBall<PR,PRE>(FloatLowerBound<PR> const& lower, FloatUpperBound<PR> const& upper) = delete;

    FloatBall<PR,PRE>(ExactDouble d, PR pr);
        FloatBall<PR,PRE>(TwoExp t, PR pr);
        FloatBall<PR,PRE>(const Integer& z, PR pr);
        FloatBall<PR,PRE>(const Dyadic& w, PR pr);
        FloatBall<PR,PRE>(const Decimal& d, PR pr);
        FloatBall<PR,PRE>(const Rational& q, PR pr);
        FloatBall<PR,PRE>(const Real& r, PR pr);
        FloatBall<PR,PRE>(const FloatBall<PR,PRE>& x, PR pr);
    FloatBall<PR,PRE>(const ValidatedNumber& y, PR pr);

    explicit FloatBall<PR,PRE>(FloatBounds<PR> const& x);
    FloatBall<PR,PRE>(FloatValue<PR> const& x);

    FloatBall<PR,PRE>& operator=(const ValidatedNumber& y);

    operator ValidatedNumber () const;

    FloatBall<PR,PRE> create(const ValidatedNumber& y) const;

    FloatLowerBound<PR> const lower() const;
    FloatUpperBound<PR> const upper() const;
    FloatValue<PR> const value() const;
    FloatError<PRE> const error() const;

    RawFloatType const lower_raw() const { return sub(down,_v,_e); }
    RawFloatType const upper_raw() const { return add(up,_v,_e); }
    RawFloatType const& value_raw() const { return _v; }
    RawErrorFloatType const& error_raw() const { return _e; }
    double get_d() const { return _v.get_d(); }

    PrecisionType precision() const { return _v.precision(); }
    ErrorPrecisionType error_precision() const { return _e.precision(); }
    PropertiesType properties() const { return _v.precision(); }
    GenericType generic() const { return this->operator GenericType(); }
    FloatBall<PR,PRE> pm(FloatError<PRE> e) const;
    friend FloatApproximation<PR> round(FloatApproximation<PR> const& x);
  public:
    friend PositiveFloatUpperBound<PR> mag(FloatBall<PR,PRE> const&);
    friend PositiveFloatLowerBound<PR> mig(FloatBall<PR,PRE> const&);
    friend Bool same(FloatBall<PR,PRE> const&, FloatBall<PR,PRE> const&);
    friend Bool models(FloatBall<PR,PRE> const&, FloatValue<PR> const&);
    friend Bool consistent(FloatBall<PR,PRE> const&, FloatBall<PR,PRE> const&);
    friend Bool refines(FloatBall<PR,PRE> const&, FloatBall<PR,PRE> const&);
    friend FloatBall<PR,PRE> refinement(FloatBall<PR,PRE> const&, FloatBall<PR,PRE> const&);
  private: public:
    RawFloat<PR> _v; RawFloat<PRE> _e;
};

template<class PR, class PRE> inline FloatFactory<PR> factory(FloatBall<PR,PRE> const& flt) {
    return FloatFactory<PR>(flt.precision());
}

template<class PR, class PRE> class Positive<FloatBall<PR,PRE>> : public FloatBall<PR,PRE> {
  public:
    Positive<FloatBall<PR,PRE>>() : FloatBounds<PR>() { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy>
        Positive<FloatBall<PR,PRE>>(M m, PR pr) : FloatBall<PR,PRE>(m,pr) { }
    explicit Positive<FloatBall<PR,PRE>>(FloatBall<PR,PRE> const& x) : FloatBall<PR,PRE>(x) { }
};

template<class PR, class PRE> inline PositiveFloatBall<PR,PRE> cast_positive(FloatBall<PR,PRE> const& x) {
    return PositiveFloatBall<PR,PRE>(x); }


}

#endif
