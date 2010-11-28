/***************************************************************************
 *            rational.h
 *
 *  Copyright 2008-10  Pieter Collins
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

/*! \file rational.h
 *  \brief Rational number class.
 */
#ifndef ARIADNE_RATIONAL_H
#define ARIADNE_RATIONAL_H

#ifdef HAVE_GMPXX_H
#include <gmpxx.h>
#endif // HAVE_GMPXX_H

typedef unsigned int uint;

namespace Ariadne {

template<class X> X inf();

#ifdef DOXYGEN
//! \ingroup NumericModule
//! \brief %Rational numbers with exact arithmetic.
//! (Only available if the Gnu Multiple Precision library (GMP) is installed.)
class Rational { };
#endif // DOXYGEN


#ifdef HAVE_GMPXX_H
typedef mpq_class Rational;
Rational sqr(const Rational& q);
Rational pow(const Rational& q, uint n);
Rational pow(const Rational& q, int n);
template<> inline Rational inf<Rational>() { return Rational(1,-0); }
#else
#endif // HAVE_GMPXX_H

} // namespace Ariadne

#endif
