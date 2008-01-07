/***************************************************************************
 *            set_based_hybrid_orbiter.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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

#include "evaluation/set_based_hybrid_orbiter.h"
#include "evaluation/set_based_hybrid_orbiter.code.h"

namespace Ariadne {
  namespace Evaluation {
    using namespace Numeric;
    using namespace Geometry;

#ifdef ENABLE_FLOAT64
    template class SetBasedHybridOrbiterInterface<Float64>;
    template class SetBasedHybridOrbiter< Zonotope<Float64,Geometry::UniformErrorTag> >;
#endif
  
#ifdef ENABLE_FLOATMP
    template class SetBasedHybridOrbiterInterface<FloatMP>;
    template class SetBasedHybridOrbiter< Zonotope<FloatMP,Geometry::UniformErrorTag> >;
#endif

  }
}
