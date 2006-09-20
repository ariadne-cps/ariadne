/***************************************************************************
 *            polyhedron.cc
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

#include "geometry/polyhedron.h"
#include "geometry/polyhedron.tpl"

#include "real_typedef.h"

namespace Ariadne {
  namespace Geometry {

    template class Polyhedron<Real>;
    template std::ostream& operator<<(std::ostream&, const Polyhedron<Real>&);
    template std::istream& operator>>(std::istream&, Polyhedron<Real>&);

/*
    
    template bool disjoint(const Polyhedron<Real>&, const Polyhedron<Real>&);
    template bool interiors_intersect(const Polyhedron<Real>&, const Polyhedron<Real>&);
    template bool inner_subset(const Polyhedron<Real>&, const Polyhedron<Real>&);
    template bool subset(const Polyhedron<Real>&, const Polyhedron<Real>&);
    
    template Polyhedron<Real> regular_intersection(const Polyhedron<Real>&, const Polyhedron<Real>&);
    template Polyhedron<Real> intersection(const Polyhedron<Real>&, const Polyhedron<Real>&);
    template Polyhedron<Real> convex_hull(const Polyhedron<Real>&, const Polyhedron<Real>&);
    template Polyhedron<Real> minkowski_sum(const Polyhedron<Real>&, const Polyhedron<Real>&);
    
    template std::ostream& operator<<(std::ostream&, const Polyhedron<Real>&);
    template std::istream& operator>>(std::istream&, Polyhedron<Real>&);
*/
    
  }
}
