/***************************************************************************
 *            hybrid_grid.h
 *
 *  Copyright  2006-8  Alberto Casagrande,  Pieter Collins
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
 
#ifndef ARIADNE_HYBRID_GRID_H
#define ARIADNE_HYBRID_GRID_H

#include <set>
#include <map>

#include <string>
#include <iostream>
#include <sstream>

#include <boost/shared_ptr.hpp>

#include "base/types.h"

#include "geometry/geometrical_traits.h"
#include "geometry/set_interface.h"
#include "geometry/list_set.h"
#include "geometry/box_list_set.h"
#include "geometry/grid_cell_list_set.h"
#include "geometry/grid_mask_set.h"
#include "geometry/hybrid_space.h"
#include "geometry/hybrid_set_iterator.h"

#include "geometry/hybrid_grid_approximation_scheme.h"


namespace Ariadne {
  namespace Geometry {

    /*! \ingroup HybridSet
     *  \brief A hybrid set comprising of a single basic set for in a discrete mode.
     */
    template<class R> 
    class HybridGrid 
    {
     public:
      typedef DiscreteState discrete_state_type;
      typedef typename std::map< discrete_state_type, Grid<R> >::const_iterator locations_const_iterator;

      /*! \brief */
      HybridGrid() : _grids() { }
      /*! \brief */
      HybridGrid(const HybridSpace& loc, const R& gl) : _grids() { 
        for(typename HybridSpace::const_iterator iter=loc.begin(); iter!=loc.end(); ++iter) {
          this->_grids.insert(std::make_pair(iter->discrete_state(),Grid<R>(iter->dimension(),gl))); }
      }
      /*! \brief */
      template<class GS> HybridGrid(const GS& gs) : _grids() { 
        for(typename GS::locations_const_iterator iter=gs.locations_begin(); iter!=gs.locations_end(); ++iter) {
          this->_grids.insert(std::make_pair(iter->first,iter->second.grid())); }
      }
      /*! \brief */
      void new_location(DiscreteState q, const Grid<R>& g) {
        this->_grids.insert(std::make_pair(q,g)); }
      /*! \brief */
      const Grid<R>& operator[](DiscreteState q) const {
        return this->_grids.find(q)->second; }

      HybridSpace locations() const;
      locations_const_iterator locations_begin() const { return this->_grids.begin(); }
      locations_const_iterator locations_end() const { return this->_grids.end(); }
     private:
      std::map< discrete_state_type, Grid<R> > _grids;
    };

    template<class R> 
    std::ostream& operator<<(std::ostream& os, const HybridGrid<R>& hgr);

   
  }
}

#endif /* ARIADNE_HYBRID_GRID_H */
