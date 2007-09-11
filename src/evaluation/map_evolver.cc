/***************************************************************************
 *            map_evolver.cc
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
#include "geometry/zonotope.h"

#include "evaluation/applicator_plugin.h"
#include "evaluation/applicator_plugin.code.h"

#include "evaluation/map_evolver.h"
#include "evaluation/map_evolver.code.h"

namespace Ariadne {
  namespace Evaluation {
    using namespace Numeric;

#ifdef ENABLE_FLOAT64
    typedef Geometry::Zonotope<Interval<Float64>,Float64> IFZonotope64;
    typedef Geometry::Zonotope<Interval<Float64>,Interval<Float64> > IIZonotope64;
    template class ApplicatorPlugin<IFZonotope64>;
    template class ApplicatorPlugin<IIZonotope64>;
    template class MapEvolver<IFZonotope64>;
    template class MapEvolver<IIZonotope64>;
#endif
  
#ifdef ENABLE_FLOATMP
    typedef Geometry::Zonotope<Interval<FloatMP>,FloatMP> IFZonotopeMP;
    typedef Geometry::Zonotope<Interval<FloatMP>,Interval<FloatMP> > IIZonotopeMP;
    template class ApplicatorPlugin<IFZonotopeMP>;
    template class ApplicatorPlugin<IIZonotopeMP>;
    template class MapEvolver<IFZonotopeMP>;
    template class MapEvolver<IIZonotopeMP>;
#endif

  }
}
