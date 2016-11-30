/***************************************************************************
 *            float_bounds.h
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

/*! \file float_bounds.h
 *  \brief Floating-point bounds for real numbers.
 */

#ifndef ARIADNE_FLOAT_BOUNDS_H
#define ARIADNE_FLOAT_BOUNDS_H

#include "utility/macros.h"

#include "number.decl.h"
#include "float.decl.h"

namespace Ariadne {

template<class PR> struct NumericTraits<FloatBounds<PR>> {
    typedef ValidatedNumber GenericType;
    typedef PositiveFloatBounds<PR> PositiveType;
    typedef ValidatedKleenean LessType;
    typedef ValidatedKleenean EqualsType;
};

//! \ingroup NumericModule
//! \brief Validated bounds on a number with floating-point endpoints supporting outwardly-rounded arithmetic.
//! \details
//! Note that direct construction from a floating-point number is prohibited, since <c>%Float64Bounds(3.3)</c> would the singleton interval \f$[3.2999999999999998224,3.2999999999999998224]\f$ (the constant is first interpreted by the C++ compiler to give a C++ \c double, whereas <c>%Float64Bounds(3.3_decimal)</c> yields the interval \f$[3.2999999999999998224,3.3000000000000002665]\f$ enclosing \f$3.3\f$.
//!
//! Comparison tests on \c FloatBounds use the idea that an interval represents a single number with an unknown value.
//! Hence the result is of type \c ValidatedKleenean, which can take values { \c True, \c False, \c Indeterminate }.
//! Hence a test \f$[\underline{x},\overline{x}]\leq [\underline{y},\overline{y}]\f$ returns \c True if \f$\overline{x}\leq \underline{y}\f$, since in this case \f$x\leq x\f$ whenever \f$x_1\in[\underline{x},\overline{x}]\f$ and \f$y\in[\underline{y},\overline{y}]\f$, \c False if \f$\underline{x}>\overline{y}\f$, since in this case we know \f$x>y\f$, and \c Indeterminate otherwise, since in this case we can find \f$x,y\f$ making the result either true or false.
//! In the case of equality, the comparison \f$[\underline{x},\overline{x}]\f$==\f$[\underline{y},\overline{y}]\f$ only returns \c True if both intervals are singletons, since otherwise we can find values making the result either true of false.
//! To test equality of representation, use \c same(x,y)
//!
//! To obtain the lower and upper bounds of the possible values, use \c x.lower() and \c x.upper().
//! To obtain a best estimate of the value, use \c x.value(), which has an error at most \a x.error().
//! If \f$v\f$ and \f$e\f$ are the returned value and error for the bounds \f$[l,u]\f$, then it is guaranteed that \f$v-e\leq l\f$ and \f$v+e\geq u\f$ in exact arithmetic.
//!
//! To test if the bounds contain a number , use \c models(FloatBounds,FloatValue), and to test if bounds are inconsistent use \c inconsistent(x,y), and to test if \c x provides a better approximation, use \c refines(x,y).
//! \sa Float64, FloatMP
//!
//! \par Python interface
//!
//! In the Python interface, %Ariadne validated bounds can be constructed from Python literals of the form \c {a:b} or (deprecated) \c [a,b] .
//! The former is preferred, as it cannot be confused with literals for other classes such as Vector and Array types.
//! Automatic conversion is used to convert FloatBounds literals of the form \c {a,b} to an FloatBounds in functions.
//!
//! Care must be taken when defining intervals using floating-point coefficients, since values are first converted to the nearest
//! representable value by the Python interpreter. <br><br>
//! \code
//!   Float64Bounds({1.1:2.3}) # Create the interval [1.1000000000000001, 2.2999999999999998]
//!   Float64Bounds({2.5:4.25}) # Create the interval [2.5, 4.25], which can be represented exactly
//!   Float64Bounds([2.5,4.25]) # Alternative syntax for creating the interval [2.5, 4.25]
//! \endcode
template<class PR> class FloatBounds
    : public DispatchFloatOperations<FloatBounds<PR>>
//    , public ProvideConvertedFieldOperations<FloatBounds<PR>,FloatValue<PR>>
{
    typedef BoundedTag P; typedef RawFloat<PR> FLT;
  public:
    typedef BoundedTag Paradigm;
    typedef FloatBounds<PR> NumericType;
    typedef ValidatedNumber GenericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
  public:
    FloatBounds<PR>() : _l(0.0), _u(0.0) { }
    explicit FloatBounds<PR>(PrecisionType pr) : _l(0.0,pr), _u(0.0,pr) { }
    explicit FloatBounds<PR>(RawFloatType const& v) : _l(v), _u(v) { }
    explicit FloatBounds<PR>(RawFloatType const& l, RawFloatType const& u) : _l(l), _u(u) { }
    FloatBounds<PR>(FloatLowerBound<PR> const& lower, FloatUpperBound<PR> const& upper);
    FloatBounds<PR>(FloatLowerBound<PR> const& lower, ValidatedUpperNumber const& upper);
    FloatBounds<PR>(ValidatedLowerNumber const& lower, FloatUpperBound<PR> const& upper);
    template<class N1, class N2, EnableIf<And<IsIntegral<N1>,IsIntegral<N2>>> = dummy> FloatBounds<PR>(N1 n1, N2 n2, PR pr) : _l(n1,pr), _u(n2,pr) { }
    FloatBounds<PR>(ExactDouble const& dl, ExactDouble const& du, PrecisionType pr);
    FloatBounds<PR>(Dyadic const& wl, Dyadic const& wu, PrecisionType pr);
    FloatBounds<PR>(Rational const& ql, Rational const& qu, PrecisionType pr);
    FloatBounds<PR>(Pair<ExactDouble,ExactDouble> const& dlu, PrecisionType pr);

    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatBounds<PR>(N n, PR pr) : FloatBounds<PR>(ExactDouble(n),pr) { }
    FloatBounds<PR>(ExactDouble d, PR pr);
        FloatBounds<PR>(const Integer& z, PR pr);
        FloatBounds<PR>(const Dyadic& w, PR pr);
        FloatBounds<PR>(const Rational& q, PR pr);
        FloatBounds<PR>(const Real& x, PR pr);
        FloatBounds<PR>(const FloatBounds<PR>& x, PR pr);
    FloatBounds<PR>(const ValidatedNumber& y, PR pr);

    FloatBounds<PR>(FloatBall<PR> const& x);
    FloatBounds<PR>(FloatValue<PR> const& x);

        FloatBounds<PR>& operator=(const FloatValue<PR>& x) { return *this=FloatBounds<PR>(x); }
    FloatBounds<PR>& operator=(const ValidatedNumber& y);
    template<class N, EnableIf<IsIntegral<N>> = dummy> FloatBounds<PR>& operator=(N n) { return *this=ValidatedNumber(n); }

    operator ValidatedNumber () const;

    FloatBounds<PR> create(const ValidatedNumber& y) const;

    FloatLowerBound<PR> const lower() const;
    FloatUpperBound<PR> const upper() const;
    FloatValue<PR> const value() const;
    FloatError<PR> const error() const;

    RawFloatType const& lower_raw() const { return _l; }
    RawFloatType const& upper_raw() const { return _u; }
    RawFloatType const value_raw() const { return hlf(add_near(_l,_u)); }
    RawFloatType const error_raw() const { RawFloatType v=value_raw(); return max(sub_up(_u,v),sub_up(v,_l)); }
    double get_d() const { return value_raw().get_d(); }

    PrecisionType precision() const { ARIADNE_DEBUG_ASSERT(_l.precision()==_u.precision()); return _u.precision(); }
    GenericType generic() const;

    FloatBounds<PR> pm(FloatError<PR> e) const;

    // DEPRECATED
    explicit operator RawFloatType () const { return value_raw(); }
    friend FloatApproximation<PR> round(FloatApproximation<PR> const& x);
    friend FloatValue<PR> midpoint(FloatBounds<PR> const& x);
  public:
    friend PositiveFloatUpperBound<PR> mag(FloatBounds<PR> const&);
    friend PositiveFloatLowerBound<PR> mig(FloatBounds<PR> const&);
    friend Bool same(FloatBounds<PR> const&, FloatBounds<PR> const&);
    friend Bool models(FloatBounds<PR> const&, FloatValue<PR> const&);
    friend Bool consistent(FloatBounds<PR> const&, FloatBounds<PR> const&);
    friend Bool inconsistent(FloatBounds<PR> const&, FloatBounds<PR> const&);
    friend Bool refines(FloatBounds<PR> const&, FloatBounds<PR> const&);
    friend FloatBounds<PR> refinement(FloatBounds<PR> const&, FloatBounds<PR> const&);
  public:
    static Nat output_places;
    static Void set_output_places(Nat p) { output_places=p; }
  private: public:
    RawFloatType _l, _u;
};

template<class PR> class Positive<FloatBounds<PR>> : public FloatBounds<PR>
    , public DispatchPositiveFloatOperations<PositiveFloatBounds<PR>>
{
  public:
    Positive<FloatBounds<PR>>() : FloatBounds<PR>() { }
    explicit Positive<FloatBounds<PR>>(PR const& pr) : FloatBounds<PR>(pr) { }
    template<class M, EnableIf<IsUnsignedIntegral<M>> =dummy> Positive<FloatBounds<PR>>(M m, PR pr) : FloatBounds<PR>(m,pr) { }
    explicit Positive<FloatBounds<PR>>(RawFloat<PR> const& x) : FloatBounds<PR>(x) { }
    explicit Positive<FloatBounds<PR>>(RawFloat<PR> const& l, RawFloat<PR> const& u) : FloatBounds<PR>(l,u) { }
    explicit Positive<FloatBounds<PR>>(FloatBounds<PR> const& x) : FloatBounds<PR>(x) { }
  public:
};

template<class PR> inline PositiveFloatBounds<PR> cast_positive(FloatBounds<PR> const& x) {
    return PositiveFloatBounds<PR>(x); }

}

#endif
