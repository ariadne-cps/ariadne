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
 *  \brief ExactTag decimal numbers, useful for user input.
 */

#ifndef ARIADNE_DECIMAL_H
#define ARIADNE_DECIMAL_H

#include <string>
#include <iostream>
#include "utility/typedefs.h"
#include "numeric/number.decl.h"

namespace Ariadne {

//! \ingroup NumericModule
//! \related Rational, Real
//! \brief A decimal number.
class Decimal {
    StringType _str;
  public:
    typedef ExactTag Paradigm;
    //! \brief Default constructor creates the number 0 (zero).
    Decimal() : _str("0.0") { }
    //! \brief Construct from a double-precision floating-point number representation.
    explicit Decimal(double d);
    //! \brief Construct from a string representation.
    explicit Decimal(StringType);
    //! \brief Convert to a rational number.
    explicit operator Rational () const;
    //! \brief Convert to an generic number.
    operator ExactNumber () const;
    friend OutputStream& operator<<(OutputStream& os, Decimal const& d);
    friend Decimal operator"" _dec (long double dbl);
};
Decimal operator"" _dec (long double dbl);
Decimal operator"" _decimal (long double dbl);


} // namespace Ariadne

#endif
