/***************************************************************************
 *            numeric/interval.h
 *
 *  Copyright 2005-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file numeric/interval.h
 *  \brief Intervals of integer, rational and floating-point number types.
 */
 
#ifndef ARIADNE_NUMERIC_INTERVAL_H
#define ARIADNE_NUMERIC_INTERVAL_H

#include "numeric/interval-integer.h" // For explicit specialization of Integer interval
#include "numeric/interval-rational.h" // For explicit specialization of Rational interval
#include "numeric/interval-float.h" // For explicit specialization of Float interval

namespace Ariadne {

template<class R> class Interval {
 public:
  Interval(const R& x) : _lower(x), _upper(x) { }
  Interval(const R& l, const R& u) : _lower(l), _upper(u) { }
  const R& lower() const { return this->_lower; }
  const R& upper() const { return this->_upper; }
  R _lower,_upper;
};

} // namespace Ariadne

namespace TBLAS {
  template<class R> int iamax_ (const int N, const R *X, const int incX);
  template<class R> int iamax_ (const int N, const Ariadne::Interval<R> *X, const int incX);
}

#include "numeric/interval.inline.h" // General inline functions
#include "numeric/interval.template.h" // General inline functions

#endif /* ARIADNE_NUMERIC_INTERVAL_H */
