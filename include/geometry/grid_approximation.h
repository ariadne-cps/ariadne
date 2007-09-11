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

#include "numeric/declarations.h"
#include "geometry/declarations.h"

namespace Ariadne {
  namespace Geometry {
    
    template<class R> GridBlock<R> outer_approximation(const Point< Numeric::Interval<R> >& ipt, const Grid<R>& g);

    template<class R> GridBlock<R> over_approximation(const Rectangle<R>& r, const Grid<R>& g);
    template<class R> GridBlock<R> under_approximation(const Rectangle<R>& r, const Grid<R>& g);

    template<class R> GridBlock<R> outer_approximation(const Rectangle<R>& r, const Grid<R>& g);
    template<class R> GridBlock<R> inner_approximation(const Rectangle<R>& r, const Grid<R>& g);

    template<class R> GridCellListSet<R> outer_approximation(const Polytope<R>& pltp, const Grid<R>& g);
    template<class R> GridCellListSet<R> inner_approximation(const Polytope<R>& pltp, const Grid<R>& g);
 
    template<class R> GridCellListSet<R> outer_approximation(const Polyhedron<R>& pltp, const Grid<R>& g);
    template<class R> GridCellListSet<R> inner_approximation(const Polyhedron<R>& pltp, const Grid<R>& g);
 
    template<class R, class R0, class R1> GridCellListSet<R> outer_approximation(const Zonotope<R0,R1>& z, const Grid<R>& g);
    template<class R, class R0, class R1> GridCellListSet<R> inner_approximation(const Zonotope<R0,R1>& z, const Grid<R>& g);
 
    template<class R> GridCellListSet<R> outer_approximation(const SetInterface<R>& set, const Grid<R>& g);
    template<class R> GridCellListSet<R> inner_approximation(const SetInterface<R>& set, const Grid<R>& g);

    template<class R> GridCellListSet<R> fuzzy_outer_approximation(const Polyhedron<R>& p, const Grid<R>& g);
    template<class R> GridCellListSet<R> fuzzy_outer_approximation(const Polytope<R>& p, const Grid<R>& g);
    template<class R, class R0, class R1> GridCellListSet<R> fuzzy_outer_approximation(const Zonotope<R0,R1>& r, const Grid<R>& g);

    template<class R, class BS> GridMaskSet<R> outer_approximation(const ListSet<BS>& ls, const FiniteGrid<R>& fg);
 
    template<class R> GridMaskSet<R> outer_approximation(const SetInterface<R>& set, const FiniteGrid<R>& fg);
    template<class R> GridMaskSet<R> inner_approximation(const SetInterface<R>& set, const FiniteGrid<R>& fg);

    template<class R> ListSet< Rectangle<R> > lower_approximation(const SetInterface<R>& set, const Grid<R>& fg);
    template<class R> ListSet< Rectangle<R> > lower_approximation(const SetInterface<R>& set, const FiniteGrid<R>& fg);
 
    template<class R> void instantiate_grid_approximation();
  }
}


#endif /* ARIADNE_GRID_APRROXIMATION_H */
