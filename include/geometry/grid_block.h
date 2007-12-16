/***************************************************************************
 *            grid_block.h
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
 
/*! \file grid_block.h
 *  \brief A block of cells in a grid.
 */

#ifndef ARIADNE_GRID_BLOCK_H
#define ARIADNE_GRID_BLOCK_H

#include <iosfwd>

#include <boost/iterator/iterator_adaptor.hpp>

#include "base/array.h"
#include "base/tribool.h"
#include "base/pointer.h"

#include "combinatoric/lattice_set.h"

#include "geometry/declarations.h"
#include "geometry/set_interface.h"
#include "geometry/grid.h"
#include "geometry/grid_cell.h"
#include "geometry/grid_set_iterator.h"

/*TODO: Unify bounds in FiniteGrid, and make GridMaskSet use only finite grid bounds*/

namespace Ariadne {
  namespace Geometry {
      
    template<class R> class GridBlockIterator;

    /*! \brief A rectangular block in a grid.
     *
     *  A %GridBlock is defined by mapping a LatticeBlock \a lb into \f$\mathbb{R}^d\f$
     *  via a Grid \a g. The bounds of the ith coordinate of the cell are given by
     *  g.subdivision_coordinate(i,lb.lower_bound(i)) and
     *  g.subdivision_coordinate(i,lb.upper_bound(i)).
     *
     *  A %GridBlock satisfies all the requirements of a RectangleExpression.
     *  \ingroup BasicSet
     *  \ingroup Grid
     */
    template<class R>
    class GridBlock 
      : public RectangleExpression< GridBlock<R> >
    {
      friend class GridCell<R>;
      friend class GridMaskSet<R>;
     public:
      typedef GridSetIterator< Combinatoric::LatticeBlock::const_iterator, GridCell<R> > iterator;
      typedef GridSetIterator< Combinatoric::LatticeBlock::const_iterator, GridCell<R> > const_iterator;

      /*! \brief The type of denotable real number defining the vertices and cells of the grid. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the set. */
      typedef Point<R> state_type;

      /*!\brief Construct an empty rectangle on a grid. */
      GridBlock(const Grid<R>& g);
      /*!\brief Construct from a grid and a bounding block. */
      GridBlock(const Grid<R>& g, const Combinatoric::LatticeBlock& b);
      /*!\brief Construct from a grid and two integer arrays giving the corners. */
      GridBlock(const Grid<R>& g, const IndexArray& l, const IndexArray& u);
      /*!\brief Construct from a grid and an ordinary rectangle. */
      GridBlock(const Grid<R>& g, const Rectangle<R>& r);
      /*!\brief Construct from a GridCell. */
      GridBlock(const GridCell<R>& gc);
      
      /*!\brief Copy constructor. */
      GridBlock(const GridBlock<R>& gb);
      /*!\brief Copy assignment. */
      GridBlock<R>& operator=(const GridBlock<R>& gb);

      /*!\brief The grid containing the rectangle. */
      const Grid<R>& grid() const;

      /*!\brief The dimension of the rectangle. */
      dimension_type dimension() const;

      /*!\brief The lower bound of the \a i th coordinate. */
      R lower_bound(dimension_type i) const;
      
      /*!\brief The upper bound of the \a i th coordinate. */
      R upper_bound(dimension_type i) const;
      
      /*!\brief The position of the rectangle in the grid. */
      const Combinatoric::LatticeBlock& lattice_set() const;

      /*!\brief Tests if the rectangle is empty. */
      tribool empty() const;

      /*! \brief True if the set is bounded. */
      tribool bounded() const;

      /*!\brief A rectangle containing the grid rectangle. */
      Rectangle<R> bounding_box() const;

      /*!\brief The one-box neighbourhood of the block. */
      GridBlock<R> neighbourhood() const;

      /*!\brief A constant iterator to the lower cell in the grid block. */
      const_iterator begin() const;
      /*!\brief A constant iterator to the past-the-end cell of the grid block. */
      const_iterator end() const;
      
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     private: 
      static void _instantiate_geometry_operators();
     private:
      GridBlock(const Grid<R>* gptr, const Combinatoric::LatticeBlock& lc);
     private:
      const Grid<R>* _grid_ptr;
      Combinatoric::LatticeBlock _lattice_set;
    };

    
    template<class R> tribool subset(const Rectangle<R>&, const GridBlock<R>&);
   
    template<class R> tribool overlap(const GridBlock<R>&, const GridBlock<R>&);
    
    template<class R> tribool subset(const GridCell<R>&, const GridBlock<R>&);
    template<class R> tribool subset(const GridBlock<R>&, const GridBlock<R>&);
    


    template<class R> std::ostream& operator<<(std::ostream& os, const GridBlock<R>& gb);
    
    
  }
}


#include "grid_block.inline.h"

#endif /* ARIADNE_GRID_BLOCK_H */
