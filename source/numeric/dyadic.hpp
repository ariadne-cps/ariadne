/***************************************************************************
 *            dyadic.hpp
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

/*! \file dyadic.hpp
 *  \brief Dyadic numbers.
 */

#ifndef ARIADNE_DYADIC_HPP
#define ARIADNE_DYADIC_HPP

#include "external/gmp.hpp"

#include <iostream> // For std::floor std::ceil etc
#include <iomanip> // For std::setprecision
#include <cmath> // For std::floor std::ceil etc
#include <algorithm> // For std::max, std::min
#include <limits> // For std::numeric_limits<double>

#include "numeric/logical.hpp"
#include "numeric/integer.hpp"
#include "numeric/arithmetic.hpp"

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
    Dyadic (Integer const& p, Int q) = delete;
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
    //! \brief The smallest integer \a p such that \a x=p/2<sup>q</sup>
    Integer mantissa() const;
    //! \brief The (negative) integer \a -q such that \a x=p/2<sup>q</sup>
    Int exponent() const;
    //! \brief A double-precision approximateion.
    double get_d() const;
    mpf_t const& get_mpf() const;
    //! \brief Convert a floating-point literal to Dyadic i.e. long binary format.
    friend Dyadic operator"" _bin(long double x);
    //! \brief Halve a number.
    friend Dyadic hlf(Dyadic const&);
    friend Rational rec(Rational const&);
    friend Rational div(Rational const&, Rational const&);
    friend Integer round(Dyadic const&);
    //! \brief Convert a floating-point literal to Dyadic i.e. long binary format.
    friend OutputStream& operator<<(OutputStream& os, Dyadic const& x);
};

template<class M, EnableIf<And<IsBuiltinIntegral<M>,IsBuiltinUnsigned<M>>>> inline Dyadic::Dyadic(M m) : Dyadic(Integer(m)) { }
template<class N, EnableIf<And<IsBuiltinIntegral<N>,IsBuiltinSigned<N>>>> inline Dyadic::Dyadic(N n) : Dyadic(Integer(n)) { }


inline Dyadic operator"" _dyadic(long double x) { return Dyadic(static_cast<double>(x)); }
inline Dyadic operator"" _dy(long double x) { return operator"" _dyadic(x); }
inline Dyadic operator"" _q2(long double x) { return operator"" _dyadic(x); }
inline Dyadic operator"" _bin(long double x) { return operator"" _dyadic(x); }



} // namespace Ariadne

#endif
