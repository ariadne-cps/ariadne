/***************************************************************************
 *            grid_approximation_scheme.h
 *
 *  Copyright  2008  Pieter Collins
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
 
/*! \file grid_approximation_scheme.h
 *  \brief Scheme for approximating on grids
 */

#ifndef ARIADNE_GRID_APPROXIMATION_SCHEME_H
#define ARIADNE_GRID_APPROXIMATION_SCHEME_H

namespace Ariadne {
  namespace Geometry {
  
    class EuclideanSpace;
    template<class R> class Box;
    template<class R> class Grid;
    template<class R> class BoxListSet;
    template<class R> class GridCellListSet;
    template<class R> class GridMaskSet;

    template<class R> 
    class GridApproximationScheme {
     public:
      typedef R Real;
      typedef Box<R> BasicSet;
      typedef Grid<R> Paving;
      typedef BoxListSet<R> CoverListSet;
      typedef GridCellListSet<R> PartitionListSet;
      typedef GridMaskSet<R> PartitionTreeSet;

      typedef EuclideanSpace space_type;
      typedef R real_type;
      typedef Box<R> basic_set_type;
    };

  }
}

#endif /* ARIADNE_GRID_APPROXIMATION_SCHEME_H */


