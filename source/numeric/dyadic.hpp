/***************************************************************************
 *            dyadic.hpp
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

/*! \file dyadic.hpp
 *  \brief Dyadic numbers.
 */

#ifndef ARIADNE_DYADIC_HPP
#define ARIADNE_DYADIC_HPP

#include "../external/gmp.hpp"

#include <iostream> // For std::floor std::ceil etc
#include <iomanip> // For std::setprecision
#include <cmath> // For std::floor std::ceil etc
#include <algorithm> // For std::max, std::min
#include <limits> // For std::numeric_limits<double>

#include "../numeric/logical.hpp"
#include "../numeric/integer.hpp"
#include "../numeric/arithmetic.hpp"

namespace Ariadne {

class ExactDouble;
class Dyadic;

class InfinityException : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

//! \ingroup NumericModule
//! \related FloatDP, ExactIntervalType
//! \brief A floating-point number, which is taken to represent the \em exact value of a real quantity.
class Dyadic
    : DeclareRingOperations<Dyadic>
    , DeclareLatticeOperations<Dyadic,Dyadic>
    , DeclareComparisonOperations<Dyadic,Boolean,Boolean>
    , DefineRingOperators<Dyadic>
    , DefineComparisonOperators<Dyadic,Boolean,Boolean>
{
  public:
    mpf_t _mpf;
  public:
    typedef ExactTag Paradigm;

    //! \brief Construct a Dyadic number from a GNU mpf object.
    explicit Dyadic (mpf_t mpf);
    //! \brief Construct the Dyadic number \a p/2<sup>q</sup>.
    Dyadic (Integer const& p, Nat q);
    Dyadic (Integer const& p, Int q) = delete; // Can only construct a dyadic p/2^q for a positive (unsigned) value q, such at uint or Nat.
    Dyadic (Integer const& p, Natural q);
    //! \brief Destructor.
    ~Dyadic();
    //! \brief Default constructor creates the number 0 (zero).
    Dyadic();
    //! \brief Copy constructor.
    Dyadic(Dyadic const& n);
    Dyadic(Dyadic&& n);
    //! \brief Assignment constructor.
    Dyadic& operator=(Dyadic const& n);
    Dyadic& operator=(Dyadic&& n);
    //! \brief Convert from a built-in positive integer.
    template<class M, EnableIf<And<IsBuiltinIntegral<M>,IsBuiltinUnsigned<M>>> = dummy> Dyadic(M m);
    //! \brief Convert from a built-in integer.
    template<class N, EnableIf<And<IsBuiltinIntegral<N>,IsBuiltinSigned<N>>> = dummy> Dyadic(N n);
    //! \brief Convert from an exact double-precision number.
    Dyadic(const ExactDouble& d);
    //! \brief Convert from an integer.
    Dyadic(const Integer& z);
    //! \brief Convert from a power of two.
    Dyadic(const TwoExp& t);
    //! \brief Explicit construction from a built-in double-precision value.
    //! \details Tests to ensure that the number is not 'accidentally' created from a rounded version of a string literal,
    //! by comparing the input with it's single-precision approximation.
    explicit Dyadic(double x);
    //! \brief Convert to a generic number.
    operator Number<ExactTag> () const;
    //! \brief A representation of ±∞ or NaN.
    static Dyadic inf(Sign sgn);
    //! \brief A representation of +∞.
    static Dyadic inf();
    //! \brief A representation of NaN (not-a-number).
    static Dyadic nan();

    //! \brief The smallest integer \a p such that \a x=p/2<sup>q</sup>
    Integer mantissa() const;
    //! \brief The (negative) integer \a -q such that \a x=p/2<sup>q</sup>
    Int exponent() const;
    //! \brief A double-precision approximateion.
    double get_d() const;
    mpf_t const& get_mpf() const;
    //! \brief Convert a floating-point literal to Dyadic i.e. long binary format.
    friend Dyadic operator"" _bin(long double x);
    //! \brief Halve the number.
    friend Dyadic hlf(Dyadic const&);
    //| \brief Power of a number (m always positive).
    friend Dyadic pow(Dyadic const& x, Int m);
    friend Dyadic pow(Dyadic const& x, Nat m);

    friend Rational rec(Rational const&);
    friend Rational div(Rational const&, Rational const&);
    friend Rational operator/(Rational const&, Rational const&);

    friend Real sqrt(Real const&);
    friend Real exp(Real const&);
    friend Real log(Real const&);
    friend Real sin(Real const&);
    friend Real cos(Real const&);
    friend Real tan(Real const&);
    friend Real atan(Real const&);

    //! \brief The sign of the number.
    friend Sign sgn(Dyadic const&);
    //! \brief Round down to the nearest lower integer.
    friend Integer floor(Dyadic const&);
    //! \brief Round to the nearest integer. Rounding of halves is implementation-dependent.
    friend Integer round(Dyadic const&);
    //! \brief Round up to the nearest higher integer.
    friend Integer ceil(Dyadic const&);
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, Dyadic const& x);

    //! \brief Tests whether the value is NaN (not-a-number).
    friend Bool is_nan(Dyadic const& w);
    //! \brief Tests whether the value is ±∞.
    friend Bool is_inf(Dyadic const& w);
    //! \brief Tests whether the value is finite.
    friend Bool is_finite(Dyadic const& w);
    //! \brief Tests whether the value is zero.
    friend Bool is_zero(Dyadic const& w);
};

template<class M, EnableIf<And<IsBuiltinIntegral<M>,IsBuiltinUnsigned<M>>>> inline Dyadic::Dyadic(M m) : Dyadic(Integer(m)) { }
template<class N, EnableIf<And<IsBuiltinIntegral<N>,IsBuiltinSigned<N>>>> inline Dyadic::Dyadic(N n) : Dyadic(Integer(n)) { }

template<> class Bounds<Dyadic> {
    Dyadic _l, _u;
  public:
    Bounds<Dyadic>(Dyadic w) : _l(w), _u(w) { }
    Bounds<Dyadic>(Dyadic l, Dyadic u) : _l(l), _u(u) { }
    template<class X, EnableIf<IsConstructible<Dyadic,X>> =dummy> Bounds<Dyadic>(Bounds<X> const& x)
        : DyadicBounds(Dyadic(x.lower_raw()),Dyadic(x.upper_raw())) { }
    Bounds<Dyadic> pm(Dyadic e) { return DyadicBounds(_l-e,_u+e); }
    Dyadic lower_raw() const { return _l; }
    Dyadic upper_raw() const { return _u; }
    friend Bounds<Dyadic> add(DyadicBounds const& w1, DyadicBounds& w2) { return DyadicBounds(w1._l+w2._l,w1._u+w2._u); }
    friend Bounds<Dyadic> sub(DyadicBounds const& w1, DyadicBounds& w2) { return DyadicBounds(w1._l-w2._u,w1._u-w2._l); }
    friend DyadicBounds abs(DyadicBounds const& w) { return DyadicBounds(max(min(w._l,-w._u),0),max(-w._l,w._u)); }
    friend ValidatedKleenean operator<(DyadicBounds const& w1, DyadicBounds const& w2) {
        if (w1._u<w2._l) { return true; } else if (w1._l >= w2._u) { return false; } else { return indeterminate; } }
    friend DyadicBounds refinement(DyadicBounds const& w1, DyadicBounds const& w2) { return DyadicBounds(max(w1._l,w2._l),min(w1._u,w2._u)); }
    friend OutputStream& operator<<(OutputStream& os, DyadicBounds y) { return os << "[" << y._l << ":" << y._u << "]"; }
};

using DyadicBounds = Bounds<Dyadic>;


inline Dyadic operator"" _dyadic(long double x) { return Dyadic(static_cast<double>(x)); }
inline Dyadic operator"" _dy(long double x) { return operator"" _dyadic(x); }
inline Dyadic operator"" _q2(long double x) { return operator"" _dyadic(x); }
inline Dyadic operator"" _bin(long double x) { return operator"" _dyadic(x); }

Comparison cmp(Dyadic const& x1, Dyadic const& x2);
Dyadic make_dyadic(unsigned long long int n);

} // namespace Ariadne

#endif
