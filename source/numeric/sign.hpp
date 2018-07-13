/***************************************************************************
 *            numeric/sign.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file numeric/sign.hpp
 *  \brief
 */


#ifndef ARIADNE_SIGN_HPP
#define ARIADNE_SIGN_HPP

#include <iosfwd>

namespace Ariadne {

using OutputStream = std::ostream;

//! \brief The sign of a numerical value.
enum class Sign : char { NEGATIVE=-1, ZERO=0, POSITIVE=+1 };
OutputStream& operator<<(OutputStream& os, Sign s);
static const Sign NEGATIVE = Sign::NEGATIVE;
static const Sign ZERO = Sign::ZERO;
static const Sign POSITIVE = Sign::POSITIVE;
inline Sign operator-(Sign s) { return Sign(-static_cast<char>(s)); }
inline Sign operator*(Sign s1, Sign s2) { return Sign(static_cast<char>(s1)*static_cast<char>(s2)); }

//! \brief The result of a comparison operation.
enum class Comparison : char { LESS=-1, EQUAL=0, GREATER=+1, INCOMPARIBLE=char(-128) };
OutputStream& operator<<(OutputStream& os, Comparison c);
static const Comparison LESS = Comparison::LESS;
static const Comparison EQUAL = Comparison::EQUAL;
static const Comparison GREATER = Comparison::GREATER;
static const Comparison INCOMPARIBLE = Comparison::INCOMPARIBLE;

} // namespace Ariadne

#endif /* ARIADNE_SIGN_HPP */
