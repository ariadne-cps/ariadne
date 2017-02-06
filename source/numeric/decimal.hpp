/***************************************************************************
 *            decimal.hpp
 *
 *  Copyright 2014-17  Pieter Collins
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

/*! \file decimal.hpp
 *  \brief ExactTag decimal numbers, useful for user input.
 */

#ifndef ARIADNE_DECIMAL_HPP
#define ARIADNE_DECIMAL_HPP

#include <string>
#include <iostream>
#include "utility/typedefs.hpp"
#include "numeric/number.decl.hpp"
#include "numeric/integer.hpp"

namespace Ariadne {

//! \ingroup NumericModule
//! \related Rational, Real
//! \brief A decimal number.
class Decimal
    : public DefineComparisonOperators<Decimal,Boolean,Boolean>
{
  public:
    static const Integer _ten;
    Integer _p; Nat _q;
  public:
    typedef ExactTag Paradigm;
    //! \brief Default constructor creates the number 0 (zero).
    Decimal();
    //! \brief Construct the number p/10^q.
    explicit Decimal(Integer p, Nat q);
    explicit Decimal(Integer p, Int q) = delete;
    //! \brief Construct from a double-precision floating-point number representation.
    explicit Decimal(double d);
    //! \brief Construct from a string representation.
    explicit Decimal(String const&);
    //! \brief Convert from a dyadic.
    Decimal(Dyadic const&);
    //! \brief Convert to a rational number.
    explicit operator Rational () const;
    //! \brief Convert to an generic number.
    operator ExactNumber () const;
    //! \brief Unary plus of a decimal value.
    friend Decimal operator+(Decimal const& d);
    //! \brief Negation of a decimal value.
    friend Decimal operator-(Decimal const& d);
    //! \brief Addition of two decimal values.
    friend Decimal operator+(Decimal const& d1, Decimal const& d2);
    //! \brief Subtraction of two decimal values.
    friend Decimal operator-(Decimal const& d1, Decimal const& d2);
    //! \brief Multiplication of two decimal values.
    friend Decimal operator*(Decimal const& d1, Decimal const& d2);
    //! \brief Division of two decimal values yields a rational.
    friend Rational operator/(Decimal const& d1, Decimal const& d2);
    //! \brief Absolute value of a decimal.
    friend Decimal abs(Decimal const& d);
    //! \brief Maximum of two decimal values.
    friend Decimal max(Decimal const& d1, Decimal const& d2);
    //! \brief Minimum of two decimal values.
    friend Decimal min(Decimal const& d1, Decimal const& d2);
    //! \brief Comparison of two decimal values.
    friend Comparison cmp(Decimal const& d1, Decimal const& d2);
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, Decimal const& d);
    //! \brief Construct from a floating-point literal.
    friend Decimal operator"" _decimal (long double dbl);
    friend Decimal operator"" _dec (long double dbl);

    void canonicalize();
};
Decimal operator"" _dec (unsigned long long int n);
Decimal operator"" _dec (long double dbl);
Decimal operator"" _decimal (unsigned long long int n);
Decimal operator"" _decimal (long double dbl);


} // namespace Ariadne

#endif
