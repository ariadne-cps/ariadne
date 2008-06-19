/***************************************************************************
 *            grid_box.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
/*! \file grid_box.h
 *  \brief Boxs in a grid.
 */

#ifndef ARIADNE_GRID_BOX_H
#define ARIADNE_GRID_BOX_H

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

    /*! \brief A box in a grid with dyadic relative coordinates.
     *
     *  A %GridBox satisfies the requirements of a RectangleExpression.
     *  \ingroup BasicSet
     *  \ingroup Grid
     */
    template<class R>
    class GridBox 
      : public RectangleExpression< GridBox<R> >
    {
     public:
      /*! \brief A tag describing the type of set. */
      typedef double dyadic_type;
      /*! \brief A tag describing the type of set. */
      typedef basic_set_tag set_category;
      /*! \brief The type of denotable real number defining the vertices and Boxs of the grid. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the set. */
      typedef Point<R> state_type;
      
      /*!\brief Construct from a grid and a list of upper and lower bounds. */
      GridBox(const Grid<R>& g, const dyadic_type* ary);

      /*!\brief The grid containing the box. */
      const Grid<R>& grid() const;

      /*!\brief The dimension of the Box. */
      dimension_type dimension() const;

      /*!\brief The lower bound of the \a i<sup>th</sup>dh coordinate. */
      R lower_bound(dimension_type i) const;

      /*!\brief The upper bound of the \a i<sup>th</sup> coordinate. */
      R upper_bound(dimension_type i) const;

      /*!\brief Subdivide along the \a i<sup>th</sup> coordinate. */
      GridBox<R> subdivide(dimension_type i, bool lr);

      /*!\brief Subdivide into one of two smaller pieces. */
      GridBox<R> subdivide(bool lr);

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream&) const;
     private:
      const Grid<R> _grid;
      array<dyadic_type> _coordinates;
    };

    template<class R> std::ostream& operator<<(std::ostream& os, const GridBox<R>& gb);
      
  }
}


#include "grid_box.inline.h"

#endif /* ARIADNE_GRID_BOX_H */
