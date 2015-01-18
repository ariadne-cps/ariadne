/***************************************************************************
 *            numeric/twoexp.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file numeric/twoexp.h
 *  \brief
 */

#ifndef ARIADNE_TWOEXP_H
#define ARIADNE_TWOEXP_H

#include <cmath>

namespace Ariadne {

/************ TwoExp ***************************************************/

//! \ingroup NumericModule
//! \brief A class representing a number of the form  \c 2<sup>\it n</sup> for some \it n.
//! Useful since floating-point numbers can be exactly multiplied and divided by powers of \c 2.
class TwoExp {
    Int _n;
  public:
    explicit TwoExp(Int n) : _n(n) { }
    Int exponent() const { return this->_n; }
    // NOTE: Use std::pow(2.0,_n) not (1<<_n) since latter does not handle very large exponents
    Float raw() const { return Float(this->get_d()); }
    double get_d() const { return std::pow(2.0,this->_n); }
    operator ExactFloat () const;
    operator ExactFloat64 () const;
    operator ErrorFloat64 () const;
    operator MetricFloat64 () const;
    operator BoundedFloat64 () const;
};
inline TwoExp two_exp(Int n) { return TwoExp(n); }


} // namespace Ariadne

#endif
