/***************************************************************************
 *            grid_approximation.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
/*! \file grid_approximation.h
 *  \brief Approximation of sets on grids. 
 */

#ifndef ARIADNE_GRID_APPROXIMATION_H
#define ARIADNE_GRID_APPROXIMATION_H

#include "base/types.h"
#include "numeric/declarations.h"
#include "geometry/declarations.h"

namespace Ariadne {
  
    
    template<class R> GridBlock<R> outer_approximation(const Point< Interval<R> >& ipt, const Grid<R>& g);

    template<class R> GridBlock<R> over_approximation(const Box<R>& bx, const Grid<R>& g);
    template<class R> GridBlock<R> under_approximation(const Box<R>& bx, const Grid<R>& g);
    template<class R> GridBlock<R> outer_approximation(const Box<R>& bx, const Grid<R>& g);
    template<class R> GridBlock<R> inner_approximation(const Box<R>& bx, const Grid<R>& g);

    template<class R> GridBlock<R> outer_approximation(const Box<Rational>& bx, const Grid<R>& g);
    template<class R> GridBlock<R> inner_approximation(const Box<Rational>& bx, const Grid<R>& g);

    template<class R, class X> GridCellListSet<R> outer_approximation(const Polytope<X>& pltp, const Grid<R>& g);
    template<class R, class X> GridCellListSet<R> inner_approximation(const Polytope<X>& pltp, const Grid<R>& g);
 
    template<class R, class X> GridCellListSet<R> outer_approximation(const Polyhedron<X>& pltp, const Grid<R>& g);
    template<class R, class X> GridCellListSet<R> inner_approximation(const Polyhedron<X>& pltp, const Grid<R>& g);
 
    template<class R> GridCellListSet<R> outer_approximation(const Zonotope<R>& z, const Grid<R>& g);
    template<class R> GridCellListSet<R> inner_approximation(const Zonotope<R>& z, const Grid<R>& g);
 
    template<class R> GridCellListSet<R> outer_approximation(const TaylorSet<R>& ts, const Grid<R>& g);
    template<class R> GridCellListSet<R> inner_approximation(const TaylorSet<R>& ts, const Grid<R>& g);

    template<class R> GridCellListSet<R> outer_approximation(const SetInterface< Box<R> >& set, const Grid<R>& g);
    template<class R> GridCellListSet<R> inner_approximation(const SetInterface< Box<R> >& set, const Grid<R>& g);

    template<class R> GridCellListSet<R> fuzzy_outer_approximation(const Polyhedron<R>& p, const Grid<R>& g);
    template<class R> GridCellListSet<R> fuzzy_outer_approximation(const Polytope<R>& p, const Grid<R>& g);
    template<class R> GridCellListSet<R> fuzzy_outer_approximation(const Zonotope<R>& r, const Grid<R>& g);
  
    template<class R> GridCellListSet<R> outer_approximation(const BoxListSet<R>& ls, const Grid<R>& g);
    template<class R> GridMaskSet<R> outer_approximation(const BoxListSet<R>& ls, const FiniteGrid<R>& fg);

    template<class R, class BS> GridCellListSet<R> outer_approximation(const ListSet<BS>& ls, const Grid<R>& g);
    template<class R, class BS> GridMaskSet<R> outer_approximation(const ListSet<BS>& ls, const FiniteGrid<R>& fg);
 
    template<class R> GridMaskSet<R> outer_approximation(const SetInterface< Box<R> >& set, const FiniteGrid<R>& fg);
    template<class R> GridMaskSet<R> inner_approximation(const SetInterface< Box<R> >& set, const FiniteGrid<R>& fg);

 	/**
	 *!\ brief This method computes an outer approximation for the set \a theSet on the grid \a theGrid.
	 * Note that, the depth is the total number of subdivisions (in all dimensions) of the unit cello of the grid.
	 * This method does the followig:
	 * 1. Computes the smallest Primary cell enclosing \a theSet
	 * 2. Allocates the pacing for this cell
	 * 3. Minces the paving to the level: depth + <the primary cell hight>
	 * 4. Iterates through the enabled leaf nodes of the paving (all the nodes are initially enabled)
	 * 5. Disables the cells that are disjoint with the \a theSet
	 */
	template<class R, class Set> GridPaving<R> outer_approximation(const Grid<R>& theGrid, const Set& theSet, const uint depth);
	template<class R, class Set> GridPaving<R> inner_approximation(const Grid<R>& theGrid, const Set& theSet, const uint depth);
	

    template<class R> BoxListSet<R> lower_approximation(const SetInterface< Box<R> >& set, const Grid<R>& fg);
    template<class R> BoxListSet<R> lower_approximation(const SetInterface< Box<R> >& set, const FiniteGrid<R>& fg);
 
    template<class R> BoxListSet<R> point_approximation(const SetInterface< Box<R> >& set, const Grid<R>& fg);
    template<class R> BoxListSet<R> point_approximation(const SetInterface< Box<R> >& set, const FiniteGrid<R>& fg);
 
    template<class R> void instantiate_grid_approximation();
  
} // namespace Ariadne


#endif /* ARIADNE_GRID_APRROXIMATION_H */
