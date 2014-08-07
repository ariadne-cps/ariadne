/***************************************************************************
 *            decimal.h
 *
 *  Copyright 2014  Pieter Collins
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

/*! \file decimal.h
 *  \brief Exact decimal numbers, useful for user input.
 */

#ifndef ARIADNE_DECIMAL_H
#define ARIADNE_DECIMAL_H

#include <string>

namespace Ariadne {

class Float;
class Interval;
class Rational;

//! \ingroup NumericModule
//! \related Rational, Real
//! \brief A decimal number.
class Decimal {
    std::string _str;
  public:
    //! \brief Default constructor creates the number 0 (zero).
    Decimal() : _str("0.0") { }
    //! \brief Construct from a double-precision floating-point number representation.
    explicit Decimal(double d);
    //! \brief Construct from a string representation.
    explicit Decimal(std::string);
#ifdef HAVE_GMPXX_H
    //! \brief Convert to a rational number.
    explicit operator Rational () const;
#endif // HAVE_GMPXX_H
    //! \brief Convert to a floating-point interval.
    explicit operator Interval () const;
    //! \brief Convert to a floating-point interval.
    explicit operator Float () const;
    friend std::ostream& operator<<(std::ostream& os, Decimal const& d);
};
} // namespace Ariadne

#endif
