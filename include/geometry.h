/***************************************************************************
 *            geometry.h
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
 
/*! \file geometry.h
 *  \brief Top-level header file for the Geometry module.
 */

#ifndef _ARIADNE_GEOMETRY_H
#define _ARIADNE_GEOMETRY_H

#include <gmpxx.h>
#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/io.hpp>

#include <iostream>
#include <iomanip>

namespace Ariadne {
  namespace Geometry {

    template<typename R> class Point;

    template<typename R> class Rectangle;
    template<typename R> class Polyhedron;
    template<typename R> class Parallelopiped;
    template<typename R> class Simplex;
    
    template<typename R, template<typename> class BS> class ListSet;
    template<typename R> class GridMaskSet;
    template<typename R> class GridCellListSet;
    template<typename R> class GridRectangleListSet;
    template<typename R> class PartitionTreeSet;

/* Need partial template function specialization to do this!
    template<typename R, template<typename> class BS> bool disjoint(const BS<R>&, const BS<R>&);
    template<typename R, template<typename> class BS> bool interiors_intersect(const BS<R>&, const BS<R>&);
    template<typename R, template<typename> class BS> bool subset_of_open_cover(const BS<R>&, const ListSet<R,BS>);
      
    template<class Set1, class Set2> bool disjoint(const Set1&, const Set2&);
    template<class Set1, class Set2> bool interiors_intersect(const Set1&, const Set2&);
    template<class Set1, class Set2> bool inner_subset(const Set1&, const Set2&);
    template<class Set1, class Set2> bool intersect(const Set1&, const Set2&);
    template<class Set1, class Set2> bool subset(const Set1&, const Set2&);

    template<class Set> Set regular_intersection(const Set&, const Set&);
    template<class Set> Set intersection(const Set&, const Set&);
    template<class Set> Set join(const Set&, const Set&);
*/
  }
}



#endif /* _ARIADNE_GEOMETRY_H */
