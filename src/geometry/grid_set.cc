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
    template class GridBlock<Real>;
    template class GridMaskSet<Real>;
    template class GridCellListSet<Real>;
    template class GridBlockListSet<Real>;
      
    template class GridMaskSetIterator<Real>;
    template class GridCellListSetIterator<Real>;
    template class GridBlockListSetIterator<Real>;

    template tribool disjoint(const GridMaskSet<Real>&, const Rectangle<Real>&);
    template tribool disjoint(const Rectangle<Real>&, const GridMaskSet<Real>&);

    template tribool subset(const Rectangle<Real>&, const GridBlock<Real>&);
    template tribool subset(const Rectangle<Real>&, const GridMaskSet<Real>&);

    template tribool overlap(const GridBlock<Real>&, const GridBlock<Real>&);
    template tribool overlap(const GridBlock<Real>&, const GridMaskSet<Real>&);
    template tribool overlap(const GridMaskSet<Real>&, const GridBlock<Real>&);
    template tribool overlap(const GridMaskSet<Real>&, const GridMaskSet<Real>&);
    
    template tribool subset(const GridCell<Real>&, const GridBlock<Real>&);
    template tribool subset(const GridCell<Real>&, const GridMaskSet<Real>&);
    template tribool subset(const GridBlock<Real>&, const GridBlock<Real>&);
    template tribool subset(const GridBlock<Real>&, const GridMaskSet<Real>&);
    template tribool subset(const GridCellListSet<Real>&, const GridBlock<Real>&);
    template tribool subset(const GridCellListSet<Real>&, const GridMaskSet<Real>&);
    template tribool subset(const GridMaskSet<Real>&, const GridMaskSet<Real>&);

    template GridCellListSet<Real> regular_intersection(const GridCellListSet<Real>&, 
                                                     const GridMaskSet<Real>&);
    template GridCellListSet<Real> regular_intersection(const GridMaskSet<Real>&,
                                                     const GridCellListSet<Real>&);
    template GridMaskSet<Real> regular_intersection(const GridBlock<Real>&, 
                                                 const GridMaskSet<Real>&);
    template GridMaskSet<Real> regular_intersection(const GridMaskSet<Real>&, 
                                                 const GridBlock<Real>&);
    template GridMaskSet<Real> regular_intersection(const GridMaskSet<Real>&, 
                                                 const GridMaskSet<Real>&);
    
    template GridCellListSet<Real> difference(const GridCellListSet<Real>&, 
                                              const GridMaskSet<Real>&);
    template GridMaskSet<Real> difference(const GridMaskSet<Real>&, 
                                          const GridMaskSet<Real>&);
    template GridMaskSet<Real> join(const GridMaskSet<Real>&, const GridMaskSet<Real>&);

    template GridBlock<Real>
    over_approximation(const Rectangle<Real>& r, const Grid<Real>& g);

    template
    GridCellListSet<Real>
    over_approximation(const Zonotope<Real>& z, const Grid<Real>& g);

    template
    GridCellListSet<Real>
    over_approximation(const Polytope<Real>& p, const Grid<Real>& g);


    template
    GridMaskSet<Real>
    over_approximation(const ListSet<Real,Rectangle>& ls, 
                       const FiniteGrid<Real>& g);

    template
    GridMaskSet<Real>
    over_approximation(const ListSet<Real,Parallelotope>& ls, 
                       const FiniteGrid<Real>& g); 

    template
    GridMaskSet<Real>
    over_approximation(const ListSet<Real,Zonotope>& ls, 
                       const FiniteGrid<Real>& fg); 

    template
    GridMaskSet<Real>
    over_approximation(const GridMaskSet<Real>& gms, const FiniteGrid<Real>& fg); 
    
    template
    GridMaskSet<Real>
    over_approximation(const PartitionTreeSet<Real>& pts, const FiniteGrid<Real>& fg); 
    
    template
    GridMaskSet<Real>
    over_approximation(const Set<Real>& set, const FiniteGrid<Real>& fg); 
    
    template
    GridMaskSet<Real>
    over_approximation(const Set<Real>& set, const Grid<Real>& g); 
    


    template GridBlock<Real>
    under_approximation(const Rectangle<Real>& r, const Grid<Real>& g);

    template
    GridCellListSet<Real>
    under_approximation(const Parallelotope<Real>& p, const Grid<Real>& g);

    template
    GridCellListSet<Real>
    under_approximation(const Zonotope<Real>& z, const Grid<Real>& g);

    template
    GridCellListSet<Real>
    under_approximation(const Polytope<Real>& p, const Grid<Real>& g);

    
    
    template
    GridMaskSet<Real>
    under_approximation(const ListSet<Real,Rectangle>& ls, 
                        const FiniteGrid<Real>& g); 

    template
    GridMaskSet<Real>
    under_approximation(const GridMaskSet<Real>& gms, const FiniteGrid<Real>& g);
    
    template
    GridMaskSet<Real>
    under_approximation(const PartitionTreeSet<Real>& pts, const FiniteGrid<Real>& g); 
    
    
    
    template std::ostream& operator<<(std::ostream&, const GridCell<Real>&);
    template std::ostream& operator<<(std::ostream&, 
                                const GridBlock<Real>&);
    template std::ostream& operator<<(std::ostream&, 
                                const GridBlockListSet<Real>&);
    template std::ostream& operator<<(std::ostream&,
                                const GridCellListSet<Real>&);
    template std::ostream& operator<<(std::ostream&, 
                                const GridMaskSet<Real>&);
 
  }
}
