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
 *  \brief Intervals of real number types (currently implemented using Boost).
 */
 
#ifndef ARIADNE_NUMERIC_INTERVAL_H
#define ARIADNE_NUMERIC_INTERVAL_H

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cassert>

#include "base/tribool.h"

#include "numeric/traits.h"
#include "numeric/expression.h"
#include "numeric/operators.h"

#include "numeric/interval.class.h"

#include "numeric/interval_integer.h" // For explicit specialization of Integer interval
//#include "numeric/interval_rational.h" // For explicit specialization of Rational interval

namespace TBLAS {
  template<class real> int iamax_ (const int N, const real *X, const int incX);
  template<class real> int iamax_ (const int N, const Ariadne::Numeric::Interval<real> *X, const int incX);
}

#include "numeric/interval.template.h"
#include "numeric/interval.inline.h"

#endif /* ARIADNE_NUMERIC_INTERVAL_H */
