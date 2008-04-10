/***************************************************************************
 *            evolver.cc
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

#include "system/map.h"
#include "system/vector_field.h"

#include "evaluation/map_evolver.h"
#include "evaluation/map_evolver.code.h"

#include "evaluation/vector_field_evolver.h"
#include "evaluation/vector_field_evolver.code.h"

#include "evaluation/set_based_hybrid_evolver.h"
#include "evaluation/set_based_hybrid_evolver.code.h"

namespace Ariadne {
  namespace Evaluation {
    using namespace Numeric;

#ifdef ENABLE_FLOAT64
    template class MapEvolver< Geometry::Zonotope<Float64> >;
    template class VectorFieldEvolver< Geometry::Zonotope<Float64> >;
    template class SetBasedHybridEvolver< Geometry::Zonotope<Float64> >;
#endif
  
#ifdef ENABLE_FLOATMP
    template class MapEvolver< Geometry::Zonotope<FloatMP> >;
    template class VectorFieldEvolver< Geometry::Zonotope<FloatMP> >;
    template class SetBasedHybridEvolver< Geometry::Zonotope<FloatMP> >;
#endif

  }
}
