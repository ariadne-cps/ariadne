/***************************************************************************
 *            float_bounds.hpp
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

/*! \file float_bounds.hpp
 *  \brief Floating-point bounds for real numbers.
 */

#ifndef ARIADNE_FLOAT_BOUNDS_HPP
#define ARIADNE_FLOAT_BOUNDS_HPP

#include "../utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"

#include "float_operations.hpp"
#include "float_factory.hpp"

namespace Ariadne {

template<class F> struct NumericTraits<Bounds<F>> {
    typedef ValidatedNumber GenericType;
    typedef Bounds<F> OppositeType;
    typedef PositiveBounds<F> PositiveType;
    typedef ValidatedKleenean LessType;
    typedef ValidatedKleenean EqualsType;
};

//! \ingroup NumericModule
//! \brief Validated bounds on a number with floating-point endpoints supporting outwardly-rounded arithmetic.
//! \details
//! Note that direct construction from a floating-point number is prohibited, since <c>%FloatDPBounds(3.3)</c> would the singleton interval \f$[3.2999999999999998224,3.2999999999999998224]\f$ (the constant is first interpreted by the C++ compiler to give a C++ \c double, whereas <c>%FloatDPBounds(3.3_decimal)</c> yields the interval \f$[3.2999999999999998224,3.3000000000000002665]\f$ enclosing \f$3.3\f$.
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
//! \sa Real, FloatDP, FloatMP, FloatValue, FloatBall, FloatUpperBound, FloatLowerBound, FloatApproximation.
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
//!   FloatDPBounds({1.1:2.3}) # Create the interval [1.1000000000000001, 2.2999999999999998]
//!   FloatDPBounds({2.5:4.25}) # Create the interval [2.5, 4.25], which can be represented exactly
//!   FloatDPBounds([2.5,4.25]) # Alternative syntax for creating the interval [2.5, 4.25]
//! \endcode
template<class F> class Bounds
    : public DispatchFloatOperations<Bounds<F>>
//    , public ProvideConvertedFieldOperations<Bounds<F>,Value<F>>
{
  protected:
    typedef ValidatedTag P; typedef typename F::RoundingModeType RND; typedef typename F::PrecisionType PR;
  public:
    typedef P Paradigm;
    typedef Bounds<F> NumericType;
    typedef Number<P> GenericType;
    typedef F RawType;
    typedef PR PrecisionType;
    typedef PR PropertiesType;
  public:
    Bounds<F>() : _l(0.0), _u(0.0) { }
    explicit Bounds<F>(PrecisionType pr) : _l(0.0,pr), _u(0.0,pr) { }
    explicit Bounds<F>(RawType const& v) : _l(v), _u(v) { }
    explicit Bounds<F>(RawType const& l, RawType const& u) : _l(l), _u(u) { }
    Bounds<F>(LowerBound<F> const& lower, UpperBound<F> const& upper);
    Bounds<F>(LowerBound<F> const& lower, ValidatedUpperNumber const& upper);
    Bounds<F>(ValidatedLowerNumber const& lower, UpperBound<F> const& upper);
    Bounds<F>(ValidatedLowerNumber const& lower, ValidatedUpperNumber const& upper, PR pr);
    template<class N1, class N2, EnableIf<And<IsBuiltinIntegral<N1>,IsBuiltinIntegral<N2>>> = dummy> Bounds<F>(N1 n1, N2 n2, PR pr) : _l(n1,pr), _u(n2,pr) { }
    Bounds<F>(ExactDouble const& dl, ExactDouble const& du, PrecisionType pr);
    Bounds<F>(Dyadic const& wl, Dyadic const& wu, PrecisionType pr);
    Bounds<F>(Rational const& ql, Rational const& qu, PrecisionType pr);

    template<class FF, EnableIf<IsConstructible<F,FF,RND,PR>> =dummy>
        Bounds<F>(Bounds<FF> const& x, PR pr) : Bounds<F>(F(x.lower_raw(),downward,pr),F(x.upper_raw(),upward,pr)) { }

    template<class N, EnableIf<IsBuiltinIntegral<N>> = dummy> Bounds<F>(N n, PR pr) : Bounds<F>(ExactDouble(n),pr) { }
    Bounds<F>(ExactDouble d, PR pr);
        Bounds<F>(TwoExp t, PR pr);
        Bounds<F>(const Integer& z, PR pr);
        Bounds<F>(const Dyadic& w, PR pr);
        Bounds<F>(const Decimal& d, PR pr);
        Bounds<F>(const Rational& q, PR pr);
        Bounds<F>(const Real& x, PR pr);
        Bounds<F>(const Bounds<F>& x, PR pr);
    Bounds<F>(const ValidatedNumber& y, PR pr);

    template<class FE> Bounds<F>(Ball<F,FE> const& x);
    Bounds<F>(Value<F> const& x);

        Bounds<F>& operator=(const Value<F>& x) { return *this=Bounds<F>(x); }
    Bounds<F>& operator=(const ValidatedNumber& y);
    template<class N, EnableIf<IsBuiltinIntegral<N>> = dummy> Bounds<F>& operator=(N n) { return *this=ValidatedNumber(n); }

    operator ValidatedNumber () const;

    Bounds<F> create(const ValidatedNumber& y) const;

    LowerBound<F> const lower() const;
    UpperBound<F> const upper() const;
    Value<F> const value() const;
    Error<F> const error() const;

    RawType const& lower_raw() const { return _l; }
    RawType const& upper_raw() const { return _u; }
    RawType const value_raw() const { return hlf(add(near,_l,_u)); }
    RawType const error_raw() const { RawType v=value_raw(); return max(sub(up,_u,v),sub(up,v,_l)); }
    double get_d() const { return value_raw().get_d(); }

    PrecisionType precision() const { ARIADNE_DEBUG_ASSERT(_l.precision()==_u.precision()); return _u.precision(); }
    PropertiesType properties() const { ARIADNE_DEBUG_ASSERT(_l.precision()==_u.precision()); return _u.precision(); }
    GenericType generic() const { return this->operator GenericType(); }

    Bounds<F> pm(Error<F> e) const;

    // DEPRECATED
    explicit operator RawType () const { return value_raw(); }
    friend Approximation<F> round(Approximation<F> const& x);
    friend Value<F> midpoint(Bounds<F> const& x);
  public:
    friend PositiveUpperBound<F> mag(Bounds<F> const&);
    friend PositiveLowerBound<F> mig(Bounds<F> const&);
    friend Bool same(Bounds<F> const&, Bounds<F> const&);
    friend Bool models(Bounds<F> const&, Value<F> const&);
    friend Bool consistent(Bounds<F> const&, Bounds<F> const&);
    friend Bool inconsistent(Bounds<F> const&, Bounds<F> const&);
    friend Bool refines(Bounds<F> const&, Bounds<F> const&);
    friend Bounds<F> refinement(Bounds<F> const&, Bounds<F> const&);
    friend Bounds<F> round(Bounds<F> const&);
  public:
    static Nat output_places;
    static Void set_output_places(Nat p) { output_places=p; }
  private: public:
    RawType _l, _u;
};

template<class PR> Bounds(ValidatedNumber, PR) -> Bounds<RawFloatType<PR>>;
template<class PR> Bounds(ValidatedLowerNumber, ValidatedUpperNumber, PR) -> Bounds<RawFloatType<PR>>;
template<class F> Bounds(F,F) -> Bounds<F>;


template<class F> inline FloatFactory<PrecisionType<F>> factory(Bounds<F> const& flt) { return FloatFactory<PrecisionType<F>>(flt.precision()); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::create(Number<ValidatedTag> const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::create(Number<EffectiveTag> const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::create(Number<ExactTag> const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::create(Real const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::create(Rational const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::create(Dyadic const& y) { return FloatBounds<PR>(y,_pr); }
template<class PR> inline FloatBounds<PR> FloatFactory<PR>::create(Integer const& y) { return FloatBounds<PR>(y,_pr); }


template<class F> class Positive<Bounds<F>> : public Bounds<F>
    , public DispatchPositiveFloatOperations<PositiveBounds<F>>
{
    using typename Bounds<F>::PR;
  public:
    Positive<Bounds<F>>() : Bounds<F>() { }
    explicit Positive<Bounds<F>>(PR const& pr) : Bounds<F>(pr) { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy> Positive<Bounds<F>>(M m, PR pr) : Bounds<F>(m,pr) { }
    explicit Positive<Bounds<F>>(F const& x) : Bounds<F>(x) { }
    explicit Positive<Bounds<F>>(F const& l, F const& u) : Bounds<F>(l,u) { }
    explicit Positive<Bounds<F>>(Bounds<F> const& x) : Bounds<F>(x) { }
    Positive<Bounds<F>>(Positive<LowerBound<F>> const& xl, Positive<UpperBound<F>> const& xu) : Bounds<F>(xl,xu) { }
  public:
    Positive<Value<F>> value() const { return cast_positive(this->Bounds<F>::value()); }
    Positive<LowerBound<F>> lower() const { return cast_positive(this->Bounds<F>::lower()); }
    Positive<UpperBound<F>> upper() const { return cast_positive(this->Bounds<F>::upper()); }
};

template<class F> inline PositiveBounds<F> cast_positive(Bounds<F> const& x) {
    return PositiveBounds<F>(x); }

extern template Ariadne::Nat Ariadne::Bounds<Ariadne::FloatDP>::output_places;
extern template Ariadne::Nat Ariadne::Bounds<Ariadne::FloatMP>::output_places;

}


#endif
