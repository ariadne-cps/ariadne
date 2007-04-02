/***************************************************************************
 *            irregular_grid_set.h
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

#ifndef ARIADNE_IRREGULAR_GRID_SET_H
#define ARIADNE_IRREGULAR_GRID_SET_H

#include "../combinatoric/lattice_set.h"

#include "irregular_grid.h"

namespace Ariadne {
  namespace Geometry {

    template<class R> class Point;
    template<class R> class Rectangle;
    template<class R> class GridMaskSet;

    /*! \brief A denotable set on an irregular finite grid, defined using a mask over a block of cells.
     *  
     *  An %IrregularGridMaskSet is useful for computing geometric predicates on ListSet<Rectangle>.
     *
     *  \ingroup DenotableSet
     *  \ingroup Grid
     */
    template<class R>
    class IrregularGridMaskSet 
    {
     public:

      /*! \brief The type of denotable real number defining the vertices and cells of the grid. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the set. */
      typedef Point<R> state_type;

      /*!\brief Copy constructor. */
      IrregularGridMaskSet(const IrregularGridMaskSet<R>& gms);

      /*!\brief Copy assignment. */
      IrregularGridMaskSet<R>& operator=(const GridMaskSet<R>& gms);

      /*!\brief Construct from a list set of rectangles. */
      explicit IrregularGridMaskSet(const ListSet< Rectangle<R> >& rls);

      /*!\brief Convert to a list set of rectangles. */
      operator ListSet< Rectangle<R> > () const;

      //@{
      //! \name Set methods
      /*! \brief The space dimension of the set. */
      dimension_type dimension() const;

      /*!\brief Checks if the set includes a point. */
      tribool contains(const Point<R>& p) const;

      /*! \brief The rectangle bounding the region of the mask. */
      Rectangle<R> bounding_box() const;
      //@}

      /*! \brief The underlying grid. */
      const IrregularGrid<R>& grid() const;

      /*! \brief The underlying grid. */
      const Combinatoric::LatticeMaskSet& lattice_set() const;

      /*! \brief The mask giving the cells in the set. */
      const Combinatoric::LatticeBlock& block() const;

      /*! \brief The mask giving the cells in the set. */
      const BooleanArray& mask() const;

      /*! \brief Returns true if the set is empty. */
      tribool empty() const;
      
      /*! \brief The number of cells in the grid. */
      size_type size() const;
      
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;

#ifdef DOXYGEN
      friend tribool subset<> (const Rectangle<R>&, const IrregularGridMaskSet<R>&);
      friend tribool subset<> (const ListSet< Rectangle<R> >&, const IrregularGridMaskSet<R>&);
      friend IrregularGridMaskSet<R> join<> (const IrregularGridMaskSet<R>&, const IrregularGridMaskSet<R>&);
      friend IrregularGridMaskSet<R> regular_intersection<> (const IrregularGridMaskSet<R>&, const IrregularGridMaskSet<R>&);
      friend IrregularGridMaskSet<R> difference<> (const IrregularGridMaskSet<R>&, const IrregularGridMaskSet<R>&);
#endif
     private: 
      static void _instantiate();
     private:
      IrregularGrid<R> _grid;
      Combinatoric::LatticeMaskSet _lattice_set;
    };
 
    template<class R> 
    tribool subset(const Rectangle<R>&, const IrregularGridMaskSet<R>&);

    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const IrregularGridMaskSet<R>& igms)
    {
      return igms.write(os);
    }


  }
}

#endif
