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

#include "geometry/rectangle.tpl"
#include "geometry/parallelotope.tpl"
#include "geometry/zonotope.tpl"

#include "geometry/list_set.h"
#include "geometry/list_set.tpl"

#include "real_typedef.h"

namespace Ariadne {
  namespace Geometry {

    template class ListSet<Real,Rectangle>;
    template class ListSet<Real,Parallelotope>;
    template class ListSet<Real,Zonotope>;

    template  ListSet<Real,Rectangle>::operator  ListSet<Real,Parallelotope>() const;
    template  ListSet<Real,Rectangle>::operator  ListSet<Real,Zonotope>() const;
    template  ListSet<Real,Parallelotope>::operator  ListSet<Real,Zonotope>() const;
    

    template bool disjoint(const ListSet<Real,Rectangle>&, const ListSet<Real,Rectangle>&);
    template bool interiors_intersect(const ListSet<Real,Rectangle>&, const ListSet<Real,Rectangle>&);
    template bool interiors_intersect(const ListSet<Real,Parallelotope>&, const ListSet<Real,Parallelotope>&);
    template bool interiors_intersect(const ListSet<Real,Zonotope>&, const ListSet<Real,Zonotope>&);
    template bool interiors_intersect(const Parallelotope<Real>&, const ListSet<Real,Parallelotope>&);
    template bool interiors_intersect(const ListSet<Real,Parallelotope>&, const Parallelotope<Real> &);
    template bool interiors_intersect(const Zonotope<Real>&, const ListSet<Real,Zonotope>&);
    template bool interiors_intersect(const ListSet<Real,Zonotope>&, const Zonotope<Real> &);
    template bool inner_subset(const ListSet<Real,Rectangle>&, const ListSet<Real,Rectangle>&);
    template bool subset(const ListSet<Real,Rectangle>&, const ListSet<Real,Rectangle>&);
    template ListSet<Real,Rectangle> regular_intersection(const ListSet<Real,Rectangle>&, const ListSet<Real,Rectangle>&);
   
    
  }
}
