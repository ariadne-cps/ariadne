/***************************************************************************
 *            differential.cc
 *
 *  Copyright 2008-15  Pieter Collins
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

#include "numeric/numeric.h"
#include "geometry/interval.h"

#include "algebra/differential.h"

#include "algebra_operations.h"

#include "differential.tcc"

namespace Ariadne {

template class UnivariateDifferential<Float64>;
template class UnivariateDifferential<ApproximateFloat64>;
template class UnivariateDifferential<ValidatedFloat64>;
template class UnivariateDifferential<UpperInterval>;

template class Differential<Float64>;
template class Differential<ValidatedFloat64>;
template class Differential<ApproximateFloat64>;
template class Differential<UpperInterval>;

}
