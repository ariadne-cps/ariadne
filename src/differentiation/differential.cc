/***************************************************************************
 *            differential.cc
 *
 *  Copyright  2007-8  Pieter Collins
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

#include "numeric/rational.h"
#include "numeric/float.h"
#include "numeric/interval.h"

#include "differentiation/differential.h"
#include "differentiation/differential.code.h"

#include "differentiation/differential_vector.h"
#include "differentiation/differential_vector.code.h"

namespace Ariadne {
    
    
    template class Differential<Rational>;
    template class DifferentialVector<Rational>;

#ifdef ENABLE_FLOAT64
    template class Differential<ApproximateFloat64>;
    template class Differential<Interval64>;
    template class DifferentialVector<ApproximateFloat64>;
    template class DifferentialVector<Interval64>;
#endif
    
#ifdef ENABLE_FLOATMP
    template class Differential<ApproximateFloatMP>;
    template class Differential<IntervalMP>;
    template class DifferentialVector<ApproximateFloatMP>;
    template class DifferentialVector<IntervalMP>;
#endif

  
} // namespace Ariadne
