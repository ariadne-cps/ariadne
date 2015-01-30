/***************************************************************************
 *            float-exact.h
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

/*! \file float-exact.h
 *  \brief Exact floating-point number class, a subset of dyadic numbers.
 */

#ifndef ARIADNE_DYADIC_H
#define ARIADNE_DYADIC_H

#include <iostream> // For std::floor std::ceil etc
#include <iomanip> // For std::setprecision
#include <cmath> // For std::floor std::ceil etc
#include <algorithm> // For std::max, std::min
#include <limits> // For std::numeric_limits<double>

#include "numeric/logical.h"
#include "numeric/rational.h"
#include "numeric/float.h"

namespace Ariadne {

class Dyadic;

//! \ingroup NumericModule
//! \related Float, ExactInterval
//! \brief A floating-point number, which is taken to represent the \em exact value of a real quantity.
class Dyadic {
    RawFloat _x;
  public:
    //! \brief Default constructor creates the number 0 (zero).
    Dyadic() : _x(0) { }
    //! \brief Convert from a built-in positive integer.
    Dyadic(uint n) : _x(n) { }
    //! \brief Convert from a built-in integer.
    Dyadic(int n) : _x(n) { }
    //! \brief Explicit construction from a built-in double-precision value.
    //! \details Tests to ensure that the number is not 'accidentally' created from a rounded version of a string literal,
    //! by comparing the input with it's single-precision approximation.
    explicit Dyadic(double x) : _x(x) { }
    //! \brief Explicit construction from an approximate floating-point value.
    explicit Dyadic(const RawFloat& x) : _x(x) { }
    //! \brief Convert to a rational number.
    explicit operator Rational () const;
    //! \brief The approximate floating-point number with the same value.
    RawFloat value() const { return _x; }
    //! \brief A double-precision approximateion.
    double get_d() const { return _x.get_d(); }
    friend Dyadic operator"" _bin(long double x);
};

inline Dyadic operator"" _bin(long double x) { return Dyadic(static_cast<double>(x)); }

inline Dyadic operator+(const Dyadic& x) { return Dyadic(+x.value()); }
inline Dyadic operator-(const Dyadic& x) { return Dyadic(-x.value()); }
inline Dyadic operator+(const Dyadic& x1,  const Dyadic& x2);
inline Dyadic operator-(const Dyadic& x1,  const Dyadic& x2);
inline Dyadic operator*(const Dyadic& x1,  const Dyadic& x2);
inline Rational operator/(const Dyadic& x1,  const Dyadic& x2);
inline OutputStream& operator<<(OutputStream& os, const Dyadic& x) { return os << std::showpoint << std::setprecision(18) << x.value(); }

inline Bool operator==(const Dyadic& x1, const Dyadic& x2) { return x1.value()==x2.value(); }
inline Bool operator!=(const Dyadic& x1, const Dyadic& x2) { return x1.value()!=x2.value(); }
inline Bool operator<=(const Dyadic& x1, const Dyadic& x2) { return x1.value()<=x2.value(); }
inline Bool operator>=(const Dyadic& x1, const Dyadic& x2) { return x1.value()>=x2.value(); }
inline Bool operator< (const Dyadic& x1, const Dyadic& x2) { return x1.value()< x2.value(); }
inline Bool operator> (const Dyadic& x1, const Dyadic& x2) { return x1.value()> x2.value(); }

inline Bool operator==(const Dyadic& x1, double x2) { return x1.value()==x2; }
inline Bool operator!=(const Dyadic& x1, double x2) { return x1.value()!=x2; }
inline Bool operator<=(const Dyadic& x1, double x2) { return x1.value()<=x2; }
inline Bool operator>=(const Dyadic& x1, double x2) { return x1.value()>=x2; }
inline Bool operator< (const Dyadic& x1, double x2) { return x1.value()< x2; }
inline Bool operator> (const Dyadic& x1, double x2) { return x1.value()> x2; }


#ifdef HAVE_GMPXX_H
inline Dyadic::operator Rational () const { return Rational(this->get_d(),nullptr); }

inline Bool operator==(const Dyadic& x, const Rational& q) { return Rational(x)==q; }
inline Bool operator!=(const Dyadic& x, const Rational& q) { return Rational(x)!=q; }
inline Bool operator<=(const Dyadic& x, const Rational& q) { return Rational(x)<=q; }
inline Bool operator>=(const Dyadic& x, const Rational& q) { return Rational(x)>=q; }
inline Bool operator< (const Dyadic& x, const Rational& q) { return Rational(x)< q; }
inline Bool operator> (const Dyadic& x, const Rational& q) { return Rational(x)> q; }

inline Bool operator==(const Rational& q, const Dyadic& x) { return q==Rational(x); }
inline Bool operator!=(const Rational& q, const Dyadic& x) { return q!=Rational(x); }
inline Bool operator<=(const Rational& q, const Dyadic& x) { return q<=Rational(x); }
inline Bool operator>=(const Rational& q, const Dyadic& x) { return q>=Rational(x); }
inline Bool operator< (const Rational& q, const Dyadic& x) { return q< Rational(x); }
inline Bool operator> (const Rational& q, const Dyadic& x) { return q> Rational(x); }
#endif // HAVE_GMPXX_H




} // namespace Ariadne

#endif
