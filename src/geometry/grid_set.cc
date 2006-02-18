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

namespace Ariadne {
  namespace Geometry {

    template class Grid<Dyadic>;
    template class FiniteGrid<Dyadic>;
    template class GridCell<Dyadic>;
    template class GridRectangle<Dyadic>;
    template class GridMaskSet<Dyadic>;
    template class GridCellListSet<Dyadic>;
    template class GridRectangleListSet<Dyadic>;
      
    template class GridMaskSetConstIterator<Dyadic>;
    template class GridCellListSetConstIterator<Dyadic>;
    template class GridRectangleListSetConstIterator<Dyadic>;

    template bool interiors_intersect(const Rectangle<Dyadic>&, const GridMaskSet<Dyadic>&);
    template bool interiors_intersect(const GridRectangle<Dyadic>&, const GridMaskSet<Dyadic>&);
    template bool subset(const Rectangle<Dyadic>&, const GridMaskSet<Dyadic>&);
    template bool subset(const GridMaskSet<Dyadic>&, const GridMaskSet<Dyadic>&);
    template GridMaskSet<Dyadic> regular_intersection(const GridMaskSet<Dyadic>&, 
                                                      const GridMaskSet<Dyadic>&);
    template GridMaskSet<Dyadic> join(const GridMaskSet<Dyadic>&, const GridMaskSet<Dyadic>&);
    template GridMaskSet<Dyadic> difference(const GridMaskSet<Dyadic>&, 
                                            const GridMaskSet<Dyadic>&);

    template GridRectangle<Dyadic>
    over_approximation(const Rectangle<Dyadic>& p, const Grid<Dyadic>& g);

    template
    GridCellListSet<Dyadic>
    over_approximation(const Parallelopiped<Dyadic>& p, const Grid<Dyadic>& g);

    template
    GridMaskSet<Dyadic>
    over_approximation(const ListSet<Dyadic,Rectangle>& ls, const FiniteGrid<Dyadic>& g); 

    template
    GridMaskSet<Dyadic>
    over_approximation(const ListSet<Dyadic,Parallelopiped>& ls, const FiniteGrid<Dyadic>& g); 

    template
    GridRectangle<Dyadic>
    over_approximation_of_intersection(const Rectangle<Dyadic>& r1, 
                                       const Rectangle<Dyadic>& r2,
                                       const Grid<Dyadic>& g);
    
    template
    GridCellListSet<Dyadic>
    over_approximation_of_intersection(const Parallelopiped<Dyadic>& p, 
                                       const Rectangle<Dyadic>& r,
                                       const Grid<Dyadic>& g);

    template
    GridMaskSet<Dyadic>
    over_approximation_of_intersection(const ListSet<Dyadic,Rectangle>& ls, 
                                       const Rectangle<Dyadic>& r, 
                                       const FiniteGrid<Dyadic>& g);
    
    template
    GridMaskSet<Dyadic>
    over_approximation_of_intersection(const ListSet<Dyadic,Parallelopiped>& ls, 
                                       const Rectangle<Dyadic>& r, 
                                       const FiniteGrid<Dyadic>& g);
    
    template std::ostream& operator<<(std::ostream&, const Grid<Dyadic>&);
    template std::ostream& operator<<(std::ostream&, const GridCell<Dyadic>&);
    template std::ostream& operator<<(std::ostream&, const GridRectangle<Dyadic>&);
    template std::ostream& operator<<(std::ostream&, const GridRectangleListSet<Dyadic>&);
    template std::ostream& operator<<(std::ostream&, const GridCellListSet<Dyadic>&);
    template std::ostream& operator<<(std::ostream&, const GridMaskSet<Dyadic>&);
 
  }
}
