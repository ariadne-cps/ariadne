/***************************************************************************
 *            irregular_grid.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 
/*! \file irregular_grid_set.h
 *  \brief Denotable sets on irregular grids.
 */

#ifndef ARIADNE_IRREGULAR_GRID_H
#define ARIADNE_IRREGULAR_GRID_H

#include "../combinatoric/lattice_set.h"

namespace Ariadne {
  namespace Geometry {

    template<class R> class Rectangle;
    template<class BS> class ListSet;
  
    /*! \brief A finite, nonuniform grid of rectangles in Euclidean space. 
     *  \ingroup Grid
     */
    template<class R>
    class IrregularGrid {
     public:
      /*! \brief The type of real number defining the vertices and cells of the grid. */
      typedef R real_type;

      /*! \brief Construct from a list of subdivision coordinates in each 
       * dimension. */
      explicit IrregularGrid(const array< std::vector<R> >& sp);

      /*! \brief Construct from a list of rectangles giving the grid points. */
      explicit IrregularGrid(const ListSet< Rectangle<R> >& ls);

      /*! \brief Join two irregular grids. */
      IrregularGrid(const IrregularGrid& g1,IrregularGrid& g2);

      /*! \brief The underlying dimension of the grid. */
      dimension_type dimension() const;
    
      /*! \brief The coordinate of the \a n th subdivision point in 
       * dimension \a d. */
      real_type subdivision_coordinate(dimension_type d, index_type n) const;
  
      /*! \brief The index of the subdivision at \a x. */
      index_type subdivision_index(dimension_type d, const real_type& x) const;

      /*! \brief The index of the subdivision point below \a x. */
      index_type subdivision_lower_index(dimension_type d, const real_type& x) const;

      /*! \brief The index of the subdivision point above \a x. */
      index_type subdivision_upper_index(dimension_type d, const real_type& x) const;

      /*! \brief Tests whether the grid contains the given lattice rectangle 
       * within its extent. */
      bool encloses(const Rectangle<R>& r) const;

      /*! \brief The rectangle bounding the grid. */
      Rectangle<R> extent() const;
      
      /*! \brief The block of cells enclosing a rectangle. */
      Combinatoric::LatticeBlock index_block(const Rectangle<R>& r) const;

      /*! \brief A rectangle corresponding to a block of cells. */
      Rectangle<R> rectangle(const Combinatoric::LatticeBlock& lb) const;


      /*! \brief The block of valid lattice cells. */
      Combinatoric::LatticeBlock block() const;
      /*! \brief The number of subdivision intervals in each dimension. */
      SizeArray sizes() const;
      /*! \brief The total number of cells. */
      size_type capacity() const;
      /*! \brief The number of subdivision intervals in dimension \a d. */
      size_type size(dimension_type i) const;


      /*! \brief Find the rule to translate elements from a grid to a 
       * refinement. */
      static array< std::vector<index_type> > index_translation(const IrregularGrid<R>& from, const IrregularGrid<R>& to);

      /*!\brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
     private:
      void create();
     private:
      array< array<R> > _subdivision_coordinates;
      array< size_type > _centre_positions;
    };

    template<class R>
    std::ostream& operator<<(std::ostream& os, const IrregularGrid<R>& ig) 
    {
      return ig.write(os);
    }


  }
}

#endif
