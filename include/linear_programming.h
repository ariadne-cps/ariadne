/***************************************************************************
 *            linear_programming.h
 *
 *  Copyright 2005-8  Alberto Casagrande, Pieter Collins
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
 
/*! \file linear_programming.h
 *  \brief Linear programming.
 */

#ifndef ARIADNE_LINEAR_PROGRAMMING_H
#define ARIADNE_LINEAR_PROGRAMMING_H 

#include "vector.h"
#include "matrix.h"
#include "numeric.h"
#include "tuple.h"

using namespace boost::numeric;

namespace Ariadne {

//! Test if there exists a point \f$e\f$ with \f$l \leq e \leq u\f$ and \f$l\leq Ae=b\f$.
template<class X> tribool feasible(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u);


} // namespace Ariadne

#endif
