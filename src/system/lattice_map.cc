/***************************************************************************
 *            lattice_map.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it,  Pieter.Collins@cwi.nl
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
 
#include <iostream>

#include "debug.h"

#include "geometry/lattice_set.h"
#include "system/lattice_map.h"

#include "base/stlio.h"

Ariadne::dbgstream cdbg(std::cerr,0);

namespace Ariadne {
  namespace System {
    
    void
    LatticeMultiMap::adjoin_to_image(const Geometry::LatticeCell& lc, const Geometry::LatticeCell& img)
    {
      typedef std::map<Geometry::LatticeCell,Geometry::LatticeCellListSet>::iterator map_iterator;
      map_iterator iter=this->_map.find(lc);
      if(iter==this->_map.end()) {
        this->_map.insert(std::make_pair(lc,Geometry::LatticeCellListSet(this->_result_dimension)));
        iter=this->_map.find(lc);
      }
      iter->second.adjoin(img);
    }
    
    void
    LatticeMultiMap::adjoin_to_image(const Geometry::LatticeCell& lc, const Geometry::LatticeCellListSet& img)
    {
      typedef std::map<Geometry::LatticeCell,Geometry::LatticeCellListSet>::iterator map_iterator;
      map_iterator iter=this->_map.find(lc);
      if(iter==this->_map.end()) {
        this->_map.insert(std::make_pair(lc,Geometry::LatticeCellListSet(this->_result_dimension)));
        iter=this->_map.find(lc);
      }
      iter->second.adjoin(img);
    }
    

    Geometry::LatticeCellListSet
    LatticeMultiMap::apply(const Geometry::LatticeCell& lc) const
    {
      cdbg << "LatticeMap::apply(const Geometry::LatticeCell& lc) const\n";
      typedef std::map<Geometry::LatticeCell,Geometry::LatticeCellListSet>::iterator map_iterator;
      map_iterator iter=this->_map.find(lc);
      if(iter==this->_map.end()) {
        Geometry::LatticeCellListSet img(this->_result_dimension);
        this->_map.insert(std::make_pair(lc,img));
        iter=this->_map.find(lc);
        assert(iter!=this->_map.end());
      }
      return iter->second;
    }
    
    Geometry::LatticeCellListSet
    LatticeMultiMap::operator() (const Geometry::LatticeCell& lc) const
    {
      cdbg << "LatticeMap::operator() (const Geometry::LatticeCell& lc) const\n";
      return this->apply(lc);
    }
    
    Geometry::LatticeCellListSet 
    LatticeMultiMap::operator() (const Geometry::LatticeRectangle& lr) const 
    {
      cdbg << "LatticeMap::operator() (const Geometry::LatticeRectangle& lr) const\n";
      Geometry::LatticeCellListSet result(this->_result_dimension);
      for(Geometry::LatticeRectangle::const_iterator cell_iter=lr.begin(); cell_iter!=lr.end(); ++cell_iter) {
        result.adjoin(this->apply(*cell_iter));
      }
      return result;
    }
    
    Geometry::LatticeCellListSet 
    LatticeMultiMap::operator() (const Geometry::LatticeCellListSet& lcls) const 
    {
      cdbg << "LatticeMap::operator() (const Geometry::LatticeCellListSet& lcls) const\n";
      Geometry::LatticeCellListSet result(this->_result_dimension);
      for(Geometry::LatticeCellListSet::const_iterator cell_iter=lcls.begin(); cell_iter!=lcls.end(); ++cell_iter) {
        result.adjoin(this->apply(*cell_iter));
      }
      cdbg << "Returning" << result << "\n";
      return result;
    }
    
    Geometry::LatticeCellListSet 
    LatticeMultiMap::operator() (const Geometry::LatticeRectangleListSet& lrls) const 
    {
      cdbg << "LatticeMap::operator() (const Geometry::LatticeRectangleListSet& lrls) const\n";
      Geometry::LatticeCellListSet result(this->_result_dimension);
      for(Geometry::LatticeRectangleListSet::const_iterator rect_iter=lrls.begin(); rect_iter!=lrls.end(); ++rect_iter) {
        result.adjoin(this->operator()(*rect_iter));
      }
      return result;
    }
    
    Geometry::LatticeCellListSet 
    LatticeMultiMap::operator() (const Geometry::LatticeMaskSet& lms) const 
    {
      Geometry::LatticeCellListSet result(this->_result_dimension);
      for(Geometry::LatticeMaskSet::const_iterator cell_iter=lms.begin(); cell_iter!=lms.end(); ++cell_iter) {
        result.adjoin(this->apply(*cell_iter));
      }
      return result;
    }
    
    std::ostream&
    operator<<(std::ostream& os, const LatticeMultiMap& m) 
    {
      return os << "LatticeMap(" << m._map << ")";
    }
    
  }
}
