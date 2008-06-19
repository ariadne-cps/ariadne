/***************************************************************************
 *            grid_cell.h
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
 
/*! \file grid_cell.h
 *  \brief Cells in a grid.
 */

#ifndef ARIADNE_GRID_CELL_H
#define ARIADNE_GRID_CELL_H

#include <iosfwd>

#include <boost/iterator/iterator_adaptor.hpp>

#include "base/array.h"
#include "base/tribool.h"
#include "base/pointer.h"

#include "combinatoric/lattice_set.h"

#include "geometry/declarations.h"
#include "geometry/grid.h"

namespace Ariadne {
  
      
    class basic_set_tag;

    /*! \brief A unit cell in a grid.
     *
     *  A %GridCell is defined by mapping a LatticeCell \a lc into \f$\mathbb{R}^d\f$
     *  via a grid \a g. The bounds of the ith coordinate of the cell are given by
     *  g.subdivision_coordinate(i,lc.lower_bound(i)) and
     *  g.subdivision_coordinate(i,lc.upper_bound(i)).
     *
     *  A %GridCell satisfies the requirements of a RectangleExpression.
     *  \ingroup BasicSet
     *  \ingroup Grid
     */
    template<class R>
    class GridCell 
      : public RectangleExpression< GridCell<R> >
    {
      friend class GridBlock<R>;
      friend class GridMaskSet<R>;
      friend class GridCellListSet<R>;
     public:
      /*! \brief A tag describing the type of set. */
      typedef basic_set_tag set_category;
      /*! \brief The type of denotable real number defining the vertices and cells of the grid. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the set. */
      typedef Point<R> state_type;
      
      /*!\brief Construct the unit cell of the grid. */
      GridCell(const Grid<R>& g);

      /*!\brief Construct from a grid and an unit grid cell. */
      GridCell(const Grid<R>& g, const LatticeCell& pos);

      /*!\brief Construct from a grid and an unit grid cell. */
      GridCell(const Grid<R>& g, const IndexArray& pos);

      /*!\brief The grid containing the cell. */
      const Grid<R>& grid() const;

      /*!\brief The dimension of the cell. */
      dimension_type dimension() const;

      /*!\brief The lower bound of the \a i th coordinate. */
      R lower_bound(dimension_type i) const;

      /*!\brief The upper bound of the \a i th coordinate. */
      R upper_bound(dimension_type i) const;

      /*!\brief The position of the cell in the grid. */
      const LatticeCell& lattice_set() const;

      /*! \brief True if the set is bounded. */
      tribool bounded() const;

      /*!\brief A rectangle containing the grid cell. */
      Box<R> bounding_box() const;

      /*!\brief The one-box neighbourhood of the cell. */
      GridBlock<R> neighbourhood() const;

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     private:
      Grid<R> _grid;
      LatticeCell _lattice_set;
    };

    template<class R> std::ostream& operator<<(std::ostream& os, const GridCell<R>& gc);
      
  
} // namespace Ariadne


#include "grid_cell.inline.h"

#endif /* ARIADNE_GRID_SET_H */
