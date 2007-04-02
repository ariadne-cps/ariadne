/***************************************************************************
 *            list_set.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"

#include "geometry/list_set.h"
#include "geometry/list_set.code.h"

namespace Ariadne {
  namespace Geometry {
 
    using namespace Numeric;
    
    template class ListSet< Rectangle<Rational> >;
    template class ListSet< Zonotope<Rational> >;
    
#ifdef ENABLE_FLOAT64
    template class ListSet< Rectangle<Float64> >;
    template class ListSet< Parallelotope<Float64> >;
    template class ListSet< Zonotope<Float64> >;
    template class ListSet< Polytope<Float64> >;
    template class ListSet< Zonotope< Interval<Float64> > >;

    // The following are not defined for all types,
    // so we can't instantiate them in ListSet<R,BS>::_instantiate_geometry_operators()
    template  ListSet< Rectangle<Float64> >::operator  ListSet< Zonotope< Interval<Float64> > >() const;
    template  ListSet< Rectangle<Float64> >::operator  ListSet< Zonotope<Float64> >() const;
    
    template tribool disjoint(const ListSet< Rectangle<Float64> >&, const ListSet< Rectangle<Float64> >&);
    template tribool disjoint(const ListSet< Zonotope<Float64> >&, const ListSet< Zonotope<Float64> >&);

    template tribool subset(const ListSet< Rectangle<Float64> >&, const ListSet< Rectangle<Float64> >&);
    template ListSet< Rectangle<Float64> > open_intersection(const ListSet< Rectangle<Float64> >&, const ListSet< Rectangle<Float64> >&);
#endif
   
#ifdef ENABLE_FLOATMP
    template class ListSet< Rectangle<FloatMP> >;
    template class ListSet< Parallelotope<FloatMP> >;
    template class ListSet< Zonotope<FloatMP> >;
    template class ListSet< Polytope<FloatMP> >;
    template class ListSet< Zonotope< Interval<FloatMP> > >;

    // The following are not defined for all types,
    // so we can't instantiate them in ListSet<R,BS>::_instantiate_geometry_operators()
    template  ListSet< Rectangle<FloatMP> >::operator  ListSet< Zonotope< Interval<FloatMP> > >() const;
    template  ListSet< Rectangle<FloatMP> >::operator  ListSet< Zonotope<FloatMP> >() const;
    
    template tribool disjoint(const ListSet< Rectangle<FloatMP> >&, const ListSet< Rectangle<FloatMP> >&);
    template tribool disjoint(const ListSet< Zonotope<FloatMP> >&, const ListSet< Zonotope<FloatMP> >&);

    template tribool subset(const ListSet< Rectangle<FloatMP> >&, const ListSet< Rectangle<FloatMP> >&);
    template ListSet< Rectangle<FloatMP> > open_intersection(const ListSet< Rectangle<FloatMP> >&, const ListSet< Rectangle<FloatMP> >&);
#endif
   
  }
}
