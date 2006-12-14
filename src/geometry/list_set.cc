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

#include "geometry/list_set.h"
#include "geometry/list_set.code.h"

namespace Ariadne {
  namespace Geometry {

    template class ListSet<Real,Rectangle>;
    template class ListSet<Real,Parallelotope>;
    template class ListSet<Real,Zonotope>;
    template class ListSet<Interval<Real>,Parallelotope>;
    template class ListSet<Interval<Real>,Zonotope>;

    template  ListSet<Real,Rectangle>::operator  ListSet<Real,Parallelotope>() const;
    template  ListSet<Real,Rectangle>::operator  ListSet<Real,Zonotope>() const;
    template  ListSet<Real,Parallelotope>::operator  ListSet<Real,Zonotope>() const;
    

    template tribool disjoint(const ListSet<Real,Rectangle>&, const ListSet<Real,Rectangle>&);
    template tribool disjoint(const ListSet<Real,Zonotope>&, const ListSet<Real,Zonotope>&);
    template tribool disjoint(const ListSet<Real,Parallelotope>&, const ListSet<Real,Parallelotope>&);
/*    
    template tribool disjoint(const ListSet<Real,Rectangle>&, const Rectangle<Real> &);
    template tribool disjoint(const ListSet<Real,Zonotope>&, const Zonotope<Real> &);
    template tribool disjoint(const ListSet<Real,Parallelotope>&, const Parallelotope<Real> &);
    template tribool disjoint(const Rectangle<Real>&, const ListSet<Real,Rectangle> &);
    template tribool disjoint(const Zonotope<Real>&, const ListSet<Real,Zonotope> &);
    template tribool disjoint(const Parallelotope<Real>&, const ListSet<Real,Parallelotope> &);
*/
    template tribool subset(const ListSet<Real,Rectangle>&, const ListSet<Real,Rectangle>&);
    template ListSet<Real,Rectangle> open_intersection(const ListSet<Real,Rectangle>&, const ListSet<Real,Rectangle>&);
   
    
  }
}
