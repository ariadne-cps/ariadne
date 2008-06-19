/***************************************************************************
 *            discretiser.cc
 *
 *  Copyright  2008  Pieter Collins
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

#include "numeric/float.h"

#include "geometry/zonotope.h"
#include "geometry/grid_approximation_scheme.h"

#include "system/numerical_system.h"
#include "system/map.h"

#include "evaluation/discrete_evolver.h"
#include "evaluation/discrete_evolver.code.h"

namespace Ariadne {
  
    
#ifdef ENABLE_FLOAT64
    template class DiscreteEvolver< Map<Float64> , GridApproximationScheme<Float64>, Zonotope<Float64> >;
//    template class DiscreteEvolver< NumericalSystemInterface< Integer,Zonotope<Float64> >, GridApproximationScheme<Float64>, Zonotope<Float64> >;
//    template class DiscreteEvolver< NumericalSystemInterface< Rational,Zonotope<Float64> >, GridApproximationScheme<Float64>, Zonotope<Float64> >;
#endif
  
#ifdef ENABLE_FLOATMP
//    template class DiscreteEvolver< NumericalSystemInterface< Integer,Zonotope<FloatMP> >, GridApproximationScheme<FloatMP>, Zonotope<FloatMP> >;
//    template class DiscreteEvolver< NumericalSystemInterface< Rational,Zonotope<FloatMP> >, GridApproximationScheme<FloatMP>, Zonotope<FloatMP> >;
#endif

  
} // namespace Ariadne
