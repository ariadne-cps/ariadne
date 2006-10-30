/***************************************************************************
 *            hybrid_set.h
 *
 *  Copyright  2006  Alberto Casagrande,  Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
 
#ifndef _ARIADNE_HYBRID_SET_H
#define _ARIADNE_HYBRID_SET_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "../declarations.h"
#include "../geometry/grid_set.h"
#include "../system/hybrid_automaton.h"

namespace Ariadne {  
namespace Geometry {
  
/*! 
 *  \brief A class for representing the reachable set of a hybrid system.
 */
template< class R >
class HybridGridMaskSet
{
 public:
  /*! \brief Construct a set for \a n discrete modes, each based on the same cells in the same grid. */
  HybridGridMaskSet(size_type n, const FiniteGrid<R>& fg) 
    : _component_sets(n,GridMaskSet<R>(fg)) { }
    
 
  /*! \brief Construct a set for \a n discrete modes, based on a list of finite grids. */
  template<class FGC> HybridGridMaskSet(const FGC& fgs) {
    for(typename FGC::const_iterator grid_iter=fgs.begin();
        grid_iter!=fgs.end(); ++grid_iter)
    {
      this->_component_sets.push_pack(GridMaskSet<R>(*grid_iter));
    }
  }
 
  size_type number_of_discrete_components() const { return _component_sets.size(); }
  GridMaskSet<R>& operator[] (const size_type& i) { return _component_sets[i]; }
  const GridMaskSet<R>& operator[] (const size_type& i) const { return _component_sets[i]; }
 private:
  std::vector< GridMaskSet<R> > _component_sets;
};


}
}

#endif /* _ARIADNE_HYBRID_EVOLVER_H */
