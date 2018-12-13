/***************************************************************************
 *            float_approximation.hpp
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

/*! \file float_approximation.hpp
 *  \brief Floating-point approximations to real numbers.
 */

#ifndef ARIADNE_FLOAT_APPROXIMATION_HPP
#define ARIADNE_FLOAT_APPROXIMATION_HPP

#include "../utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"

#include "float_operations.hpp"

namespace Ariadne {

template<class F> struct NumericTraits<Approximation<F>> {
    typedef ApproximateNumber GenericType;
    typedef PositiveApproximation<F> PositiveType;
    typedef Approximation<F> OppositeType;
    typedef Fuzzy LessType;
    typedef Fuzzy EqualsType;
};

//! \ingroup NumericModule
//! \brief Floating point number approximations to real numbers supporting approxiamate arithmetic.
//! \details
//! The \c %Approximation<F> class represents approximate floating-point numbers.
//! Operations are performed approximately, with no guarantees on the output.
//! \sa Real, FloatDP , FloatMP, FloatValue, FloatBall, FloatBounds.
template<class F> class Approximation
    : public DispatchFloatOperations<Approximation<F>>
{
  protected:
    typedef ApproximateTag P; typedef typename F::PrecisionType PR;
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

        Approximation<F>(double d, PR pr);
        Approximation<F>(ExactDouble d, PR pr);
        Approximation<F>(TwoExp t, PR pr);
        Approximation<F>(const Integer& z, PR pr);
        Approximation<F>(const Dyadic& w, PR pr);
        Approximation<F>(const Decimal& d, PR pr);
        Approximation<F>(const Rational& q, PR pr);
        Approximation<F>(const Real& r, PR pr);
        Approximation<F>(const Approximation<F>& r, PR pr);
    Approximation<F>(const ApproximateNumber& y, PR pr);

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
    Approximation<F>& operator=(const ApproximateNumber& y);
    Approximation<F> create(const ApproximateNumber& y) const;

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
    friend Bool same(Approximation<F> const&, Approximation<F> const&);
    friend PositiveApproximation<F> mag(Approximation<F> const&);
    friend PositiveApproximation<F> mig(Approximation<F> const&);
  public:
    static Nat output_places;
    static Void set_output_places(Nat p) { output_places=p; }
    Approximation<F> pm(Approximation<F> _e) { return *this; }
  public:
    RawType _a;
};

template<class PR> Approximation(ApproximateNumber, PR) -> Approximation<RawFloatType<PR>>;
template<class F> Approximation(F) -> Approximation<F>;

template<class F> inline FloatFactory<PrecisionType<F>> factory(Approximation<F> const& flt) { return FloatFactory<PrecisionType<F>>(flt.precision()); }
template<class PR> inline FloatApproximation<PR> FloatFactory<PR>::create(Number<ApproximateTag> const& y) { return FloatApproximation<PR>(y,_pr); }
template<class PR> template<class D, EnableIf<IsBuiltinFloatingPoint<D>>> inline
    FloatApproximation<PR> FloatFactory<PR>::create(D const& y) { return FloatApproximation<PR>(y,_pr); }

template<class F> class Positive<Approximation<F>> : public Approximation<F>
    , public DispatchPositiveFloatOperations<PositiveApproximation<F>>
{
    using typename Approximation<F>::PR;
  public:
    Positive<Approximation<F>>() : Approximation<F>() { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy>
        Positive<Approximation<F>>(M m) : Approximation<F>(m) { }
    explicit Positive<Approximation<F>>(F const& x) : Approximation<F>(x) { }
    explicit Positive<Approximation<F>>(Approximation<F> const& x) : Approximation<F>(x) { }
    explicit Positive<Approximation<F>>(ApproximateNumber const& y, PR pr) : Approximation<F>(y,pr) { }
    Positive<Approximation<F>>(PositiveLowerBound<F> const& x) : Approximation<F>(x) { }
    Positive<Approximation<F>>(PositiveUpperBound<F> const& x) : Approximation<F>(x) { }
    Positive<Approximation<F>>(PositiveValue<F> const& x) : Approximation<F>(x) { }
    Positive<Approximation<F>>(Error<F> const& x) : Approximation<F>(x) { }
  public:
};

template<class F> inline PositiveApproximation<F> cast_positive(Approximation<F> const& x) {
    return PositiveApproximation<F>(x); }

extern template Ariadne::Nat Ariadne::Approximation<Ariadne::FloatDP>::output_places;
extern template Ariadne::Nat Ariadne::Approximation<Ariadne::FloatMP>::output_places;

}

#endif
