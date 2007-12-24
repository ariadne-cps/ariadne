/***************************************************************************
 *            grid_tree_cell.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  
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

/*! \file grid_tree_cell.h
 *  \brief A cell in a grid tree. 
 */

#ifndef ARIADNE_GRID_TREE_CELL_H
#define ARIADNE_GRID_TREE_CELL_H

#include <iosfwd>

#include "base/iterator.h"
#include "base/tribool.h"

#include "combinatoric/binary_word.h"
#include "combinatoric/binary_tree.h"
#include "combinatoric/subdivision_tree_set.h"

#include "geometry/declarations.h"
#include "geometry/set_interface.h"
#include "geometry/rectangle_expression.h"


namespace Ariadne {
  namespace Geometry {

    template<class R> class Grid;
    template<class R> class GridTreeCell;

    /* External class declarations. */
    template<class R> class Box;



    /*!\ingroup BasicSet
     * \ingroup GridTree
     * \brief A rectangular cell in a grid tree.
     *
     * Defined as a SubdivisionTreeCell within a base cell given as a Box<R>.
     *
     * Satisfies the requirements of a RectangleExpression.
     */
    template<class R>
    class GridTreeCell
      : public RectangleExpression< GridTreeCell<R> >
    {
     public:
      /*! \brief A tag describing the type of set. */
      typedef basic_set_tag set_category;
      /*! \brief The type of denotable real number used for the cell blounds. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by cells. */
      typedef Point<R> state_type;

      /*!\brief Construct from a rectangle, the subdivision_coordinates and a binary word. */
      GridTreeCell(const Grid<R>& r, 
                   const Combinatoric::BinaryWord& w);

      /*!\brief The underlying grid. */
      const Grid<R>& grid() const;

      /*!\brief The dimension of the cell. */
      dimension_type dimension() const;

      /*!\brief The lower bound of the \a i th coordinate. */
      R lower_bound(dimension_type i) const;
      
      /*!\brief The upper bound of the \a i th coordinate. */
      R upper_bound(dimension_type i) const;

      /*!\brief Subdivide into one of two pieces. */
      GridTreeCell subdivide(bool lr) const;

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     private:
      Grid<R>* _grid_ptr;
      Combinatoric::BinaryWord _word;
    };

    
    
    template<class R>
    GridTreeCell<R> over_approximation(const Box<R>& r, const Grid<R>& g, uint depth);
    
    
    template<class R> 
    std::ostream& operator<<(std::ostream& os, const GridTreeCell<R>& ptc);
    
    
  }
}

#include "grid_tree_cell.inline.h"

#endif /* ARIADNE_GRID_TREE_CELL_H */
