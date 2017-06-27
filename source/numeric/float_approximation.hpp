/***************************************************************************
 *            float_approximation.hpp
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

/*! \file float_approximation.hpp
 *  \brief Floating-point approximations to real numbers.
 */

#ifndef ARIADNE_FLOAT_APPROXIMATION_HPP
#define ARIADNE_FLOAT_APPROXIMATION_HPP

#include "utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"

namespace Ariadne {

template<class PR> struct NumericTraits<FloatApproximation<PR>> {
    typedef ApproximateNumber GenericType;
    typedef PositiveFloatApproximation<PR> PositiveType;
    typedef Fuzzy LessType;
    typedef Fuzzy EqualsType;
};

//! \ingroup NumericModule
//! \brief Floating point number approximations to real numbers supporting approxiamate arithmetic.
//! \details
//! The \c %FloatApproximation<PR> class represents approximate floating-point numbers.
//! Operations are performed approximately, with no guarantees on the output.
//! \sa Real, Float64 , FloatMP, FloatValue, FloatBall, FloatBounds.
template<class PR> class FloatApproximation
    : public DispatchFloatOperations<FloatApproximation<PR>>
{
    typedef ApproximateTag P; typedef RawFloat<PR> FLT;
  public:
    typedef ApproximateTag Paradigm;
    typedef FloatApproximation<PR> NumericType;
    typedef ApproximateNumber GenericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
    typedef PR PropertiesType;
  public:
    FloatApproximation<PR>() : _a(0.0) { }
    explicit FloatApproximation<PR>(PrecisionType pr) : _a(0.0,pr) { }
    explicit FloatApproximation<PR>(RawFloatType const& a) : _a(a) { }

        FloatApproximation<PR>(double d, PR pr);
        FloatApproximation<PR>(ExactDouble d, PR pr);
        FloatApproximation<PR>(const Integer& z, PR pr);
        FloatApproximation<PR>(const Dyadic& w, PR pr);
        FloatApproximation<PR>(const Decimal& d, PR pr);
        FloatApproximation<PR>(const Rational& q, PR pr);
        FloatApproximation<PR>(const Real& r, PR pr);
        FloatApproximation<PR>(const FloatApproximation<PR>& r, PR pr);
    FloatApproximation<PR>(const ApproximateNumber& y, PR pr);

    FloatApproximation<PR>(FloatError<PR> const& x); // FIXME: Remove
    FloatApproximation<PR>(FloatValue<PR> const& x);
    FloatApproximation<PR>(FloatBall<PR> const& x);
    FloatApproximation<PR>(FloatBounds<PR> const& x);
    FloatApproximation<PR>(FloatUpperBound<PR> const& x);
    FloatApproximation<PR>(FloatLowerBound<PR> const& x);

    template<class N, EnableIf<IsBuiltinIntegral<N>> =dummy> FloatApproximation<PR>& operator=(N n) { this->_a=n; return *this; }
    template<class D, EnableIf<IsBuiltinFloatingPoint<D>> =dummy> FloatApproximation<PR>& operator=(D x) { this->_a=x; return *this; }
        FloatApproximation<PR>& operator=(const FloatLowerBound<PR>& x) { return *this=FloatApproximation<PR>(x); }
        FloatApproximation<PR>& operator=(const FloatUpperBound<PR>& x) { return *this=FloatApproximation<PR>(x); }
        FloatApproximation<PR>& operator=(const FloatBounds<PR>& x) { return *this=FloatApproximation<PR>(x); }
        FloatApproximation<PR>& operator=(const FloatValue<PR>& x) { return *this=FloatApproximation<PR>(x); }
    FloatApproximation<PR>& operator=(const ApproximateNumber& y);
    FloatApproximation<PR> create(const ApproximateNumber& y) const;

    operator ApproximateNumber () const;

    friend FloatApproximation<PR> round(FloatApproximation<PR> const& x);

    PrecisionType precision() const { return _a.precision(); }
    PropertiesType properties() const { return _a.precision(); }
    GenericType generic() const { return this->operator GenericType(); }
    explicit operator RawFloatType () const { return this->_a; }
    RawFloatType const& raw() const { return this->_a; }
    RawFloatType& raw() { return this->_a; }
    double get_d() const { return this->_a.get_d(); }
  public:
    friend Bool same(FloatApproximation<PR> const&, FloatApproximation<PR> const&);
    friend PositiveFloatApproximation<PR> mag(FloatApproximation<PR> const&);
  public:
    static Void set_output_places(Nat p) { output_places=p; }
    FloatApproximation<PR> pm(FloatApproximation<PR> _e) { return *this; }
  private: public:
    static Nat output_places;
    RawFloatType _a;
};

template<class PR> class Positive<FloatApproximation<PR>> : public FloatApproximation<PR>
    , public DispatchPositiveFloatOperations<PositiveFloatApproximation<PR>>
{
  public:
    Positive<FloatApproximation<PR>>() : FloatApproximation<PR>() { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy>
        Positive<FloatApproximation<PR>>(M m) : FloatApproximation<PR>(m) { }
    explicit Positive<FloatApproximation<PR>>(RawFloat<PR> const& x) : FloatApproximation<PR>(x) { }
    explicit Positive<FloatApproximation<PR>>(FloatApproximation<PR> const& x) : FloatApproximation<PR>(x) { }
    explicit Positive<FloatApproximation<PR>>(ApproximateNumber const& y, PR pr) : FloatApproximation<PR>(y,pr) { }
    Positive<FloatApproximation<PR>>(PositiveFloatLowerBound<PR> const& x) : FloatApproximation<PR>(x) { }
    Positive<FloatApproximation<PR>>(PositiveFloatUpperBound<PR> const& x) : FloatApproximation<PR>(x) { }
    Positive<FloatApproximation<PR>>(PositiveFloatValue<PR> const& x) : FloatApproximation<PR>(x) { }
    Positive<FloatApproximation<PR>>(FloatError<PR> const& x) : FloatApproximation<PR>(x) { }
  public:
};

template<class PR> inline PositiveFloatApproximation<PR> cast_positive(FloatApproximation<PR> const& x) {
    return PositiveFloatApproximation<PR>(x); }


}

#endif
