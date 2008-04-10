/***************************************************************************
 *            reducer.cc
 *
 *  Copyright  2007  Pieter Collins
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

#include "geometry/rectangle.h"
#include "geometry/zonotope.h"

#include "evaluation/identity_reducer.h"
#include "evaluation/orthogonal_reducer.h"
#include "evaluation/cascade_reducer.h"

namespace Ariadne {
  namespace Evaluation {
    using namespace Numeric;
    using namespace Geometry;

#ifdef ENABLE_FLOAT64
    template class IdentityReducer< Rectangle<Float64> >;
    template class IdentityReducer< Zonotope<Float64> >;
    template class OrthogonalReducer< Zonotope<Float64> >;
    template class CascadeReducer< Zonotope<Float64> >;
#endif
  
#ifdef ENABLE_FLOATMP
    template class IdentityReducer< Rectangle<FloatMP> >;
    template class IdentityReducer< Zonotope<FloatMP> >;
    template class OrthogonalReducer< Zonotope<FloatMP> >;
    template class CascadeReducer< Zonotope<FloatMP> >;
#endif

  }
}
