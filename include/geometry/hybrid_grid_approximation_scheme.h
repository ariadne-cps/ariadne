/***************************************************************************
 *            hybrid_grid_approximation_scheme.h
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
 
/*! \file hybrid_grid_approximation_scheme.h
 *  \brief Scheme for approximating on hybrid_grids
 */

#ifndef ARIADNE_HYBRID_GRID_APPROXIMATION_SCHEME_H
#define ARIADNE_HYBRID_GRID_APPROXIMATION_SCHEME_H

namespace Ariadne {
  
  
    class HybridSpace;
    template<class R> class HybridGrid;
    template<class R> class HybridBox;
    template<class R> class HybridBoxListSet;
    template<class R> class HybridGridCellListSet;
    template<class R> class HybridGridMaskSet;

    template<class R> 
    class HybridGridApproximationScheme {
     public:
      typedef R Real;
      typedef HybridBox<R> BasicSet;
      typedef HybridGrid<R> Paving;
      typedef HybridBoxListSet<R> CoverListSet;
      typedef HybridGridCellListSet<R> PartitionListSet;
      typedef HybridGridMaskSet<R> PartitionTreeSet;

      typedef HybridSpace space_type;
      typedef R real_type;
      typedef HybridBox<R> basic_set_type;
    };

  
} // namespace Ariadne

#endif /* ARIADNE_HYBRID_GRID_APPROXIMATION_SCHEME_H */
