/***************************************************************************
 *            grid_set.cc
 *
 *  15 February 2006
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

#include "geometry/grid_set.h"
#include "geometry/grid_set.tpl"

#include "real_typedef.h"

namespace Ariadne {
  namespace Geometry {

    template class GridCell<Real>;
    template class GridRectangle<Real>;
    template class GridMaskSet<Real>;
    template class GridCellListSet<Real>;
    template class GridRectangleListSet<Real>;
      
    template class GridMaskSetIterator<Real>;
    template class GridCellListSetIterator<Real>;
    template class GridRectangleListSetIterator<Real>;

    template bool interiors_intersect(const Rectangle<Real>&, const GridMaskSet<Real>&);
    template bool interiors_intersect(const GridRectangle<Real>&, const GridMaskSet<Real>&);
    template bool subset(const Rectangle<Real>&, const GridMaskSet<Real>&);
    template bool subset(const GridMaskSet<Real>&, const GridMaskSet<Real>&);
    template GridMaskSet<Real> regular_intersection(const GridMaskSet<Real>&, 
                                                      const GridMaskSet<Real>&);
    template GridMaskSet<Real> join(const GridMaskSet<Real>&, const GridMaskSet<Real>&);
    template GridMaskSet<Real> difference(const GridMaskSet<Real>&, 
                                            const GridMaskSet<Real>&);

    template GridRectangle<Real>
    over_approximation(const Rectangle<Real>& p, const Grid<Real>& g);

    template
    GridCellListSet<Real>
    over_approximation(const Parallelotope<Real>& p, const Grid<Real>& g);

    template
    GridMaskSet<Real>
    over_approximation(const ListSet<Real,Rectangle>& ls, const FiniteGrid<Real>& g); 

    template
    GridMaskSet<Real>
    over_approximation(const ListSet<Real,Parallelotope>& ls, const FiniteGrid<Real>& g); 

    template
    GridRectangle<Real>
    over_approximation_of_intersection(const Rectangle<Real>& r1, 
                                       const Rectangle<Real>& r2,
                                       const Grid<Real>& g);
    
    template
    GridCellListSet<Real>
    over_approximation_of_intersection(const Parallelotope<Real>& p, 
                                       const Rectangle<Real>& r,
                                       const Grid<Real>& g);

    template
    GridMaskSet<Real>
    over_approximation_of_intersection(const ListSet<Real,Rectangle>& ls, 
                                       const Rectangle<Real>& r, 
                                       const FiniteGrid<Real>& g);
    
    template
    GridMaskSet<Real>
    over_approximation_of_intersection(const ListSet<Real,Parallelotope>& ls, 
                                       const Rectangle<Real>& r, 
                                       const FiniteGrid<Real>& g);
    
    template std::ostream& operator<<(std::ostream&, const GridCell<Real>&);
    template std::ostream& operator<<(std::ostream&, const GridRectangle<Real>&);
    template std::ostream& operator<<(std::ostream&, const GridRectangleListSet<Real>&);
    template std::ostream& operator<<(std::ostream&, const GridCellListSet<Real>&);
    template std::ostream& operator<<(std::ostream&, const GridMaskSet<Real>&);
 
  }
}
