/***************************************************************************
 *            dyadic.h
 *
 *  Copyright 2008-14  Pieter Collins
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

/*! \file dyadic.h
 *  \brief Dyadic numbers.
 */

#ifndef ARIADNE_DYADIC_H
#define ARIADNE_DYADIC_H

#include "external/gmp.h"

#include <iostream> // For std::floor std::ceil etc
#include <iomanip> // For std::setprecision
#include <cmath> // For std::floor std::ceil etc
#include <algorithm> // For std::max, std::min
#include <limits> // For std::numeric_limits<double>

#include "numeric/logical.h"
#include "numeric/integer.h"
#include "numeric/arithmetic.h"

namespace Ariadne {

class Dyadic;

//! \ingroup NumericModule
//! \related Float64, ExactIntervalType
//! \brief A floating-point number, which is taken to represent the \em exact value of a real quantity.
class Dyadic
    : Ring<Dyadic>
    , DirectedLattice<Dyadic>
    , Ordered<Dyadic,Boolean>
    , DefineRingOperators<Dyadic>
    , DefineComparisonOperators<Dyadic,Boolean>
{
  public:
    mpf_t _mpf;
  public:
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
    template<class M, EnableIf<And<IsIntegral<M>,IsUnsigned<M>>> = dummy> Dyadic(M m);
    //! \brief Convert from a built-in integer.
    template<class N, EnableIf<And<IsIntegral<N>,IsSigned<N>>> = dummy> Dyadic(N n);
    //! \brief Convert from an integer.
    Dyadic(const Integer& p);
    //! \brief Explicit construction from a built-in double-precision value.
    //! \details Tests to ensure that the number is not 'accidentally' created from a rounded version of a string literal,
    //! by comparing the input with it's single-precision approximation.
    explicit Dyadic(double x);
    //! \brief Explicit construction from a floating-point value.
    explicit Dyadic(const RawFloat64& x);
    //! \brief Explicit construction from a floating-point value.
    explicit Dyadic(const RawFloatMP& x);
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
    Dyadic hlf(Dyadic const&);
    //! \brief Convert a floating-point literal to Dyadic i.e. long binary format.
    friend OutputStream& operator<<(OutputStream& os, Dyadic const& x);
};

template<class M, EnableIf<And<IsIntegral<M>,IsUnsigned<M>>>> inline Dyadic::Dyadic(M m) : Dyadic(Integer(m)) { }
template<class N, EnableIf<And<IsIntegral<N>,IsSigned<N>>>> inline Dyadic::Dyadic(N n) : Dyadic(Dyadic(n)) { }


inline Dyadic operator"" _bin(long double x) { return Dyadic(static_cast<double>(x)); }



} // namespace Ariadne

#endif
