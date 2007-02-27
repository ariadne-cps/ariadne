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

#include "real_typedef.h"

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"

#include "geometry/list_set.h"
#include "geometry/list_set.code.h"

namespace Ariadne {
  namespace Geometry {
 
    template class ListSet< Rectangle<Real> >;
    template class ListSet< Parallelotope<Real> >;
    template class ListSet< Zonotope<Real> >;
    template class ListSet< Polytope<Real> >;
    template class ListSet< Zonotope< Interval<Real> > >;

    template class ListSet< Rectangle<Rational> >;
    template class ListSet< Zonotope<Rational> >;
    
    // The following are not defined for all types,
    // so we can't instantiate them in ListSet<R,BS>::_instantiate_geometry_operators()
    template  ListSet< Rectangle<Real> >::operator  ListSet< Zonotope< Interval<Real> > >() const;
    template  ListSet< Rectangle<Real> >::operator  ListSet< Zonotope<Real> >() const;
    
    template tribool disjoint(const ListSet< Rectangle<Real> >&, const ListSet< Rectangle<Real> >&);
    template tribool disjoint(const ListSet< Zonotope<Real> >&, const ListSet< Zonotope<Real> >&);

    template tribool subset(const ListSet< Rectangle<Real> >&, const ListSet< Rectangle<Real> >&);
    template ListSet< Rectangle<Real> > open_intersection(const ListSet< Rectangle<Real> >&, const ListSet< Rectangle<Real> >&);
   
  }
}
