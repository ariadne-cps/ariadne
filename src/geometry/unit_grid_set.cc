/***************************************************************************
 *            unit_grid_set.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <iostream>

#include "geometry/unit_grid_set.h"
#include "geometry/grid_operations.h"

namespace Ariadne {
  namespace Geometry {
    
    
    UnitGridMaskSet::UnitGridMaskSet(const UnitGridMaskSet& ms)
      : _bounds(ms._bounds), _mask(ms._mask)
    {
      this->_compute_cached_attributes();
    }
    

    void
    UnitGridMaskSet::_compute_cached_attributes() {
      _lower=_bounds.lower();
      _upper=_bounds.upper();
      _sizes=_bounds.sizes();
      _strides=_bounds.strides();
    }
    
    bool 
    UnitGridRectangle::empty() const
    {
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        if(this->lower_bound(i)>=this->upper_bound(i)) {
          return true;
        }
      }
      return false;
    }
    
    bool 
    subset(const UnitGridRectangle& r1, const UnitGridRectangle& r2) 
    {
      for(dimension_type i=0; i!=r1.dimension(); ++i) {
        if(r1.lower_bound(i)<r2.lower_bound(i)
            || r1.upper_bound(i)>r2.upper_bound(i))
        {
          return false;
        }
      }
      return true;
    }

    bool 
    interiors_intersect(const UnitGridRectangle& r1, const UnitGridRectangle& r2) 
    {
      for(dimension_type i=0; i!=r1.dimension(); ++i) {
        if(r1.upper_bound(i)<=r2.lower_bound(i)
            || r1.upper_bound(i)<=r2.lower_bound(i))
        {
          return false;
        }
      }
      return true;
    }
    
    UnitGridRectangle 
    regular_intersection(const UnitGridRectangle& r1, const UnitGridRectangle& r2) 
    {
      UnitGridRectangle result(r1.dimension());
      for(dimension_type i=0; i!=r1.dimension(); ++i) {
        result.set_lower_bound(i,std::max(r1.lower_bound(i),r2.lower_bound(i)));
        result.set_upper_bound(i,std::min(r1.upper_bound(i),r2.upper_bound(i)));
      }
      return result;
    }


    bool 
    subset(const UnitGridCell& c, const UnitGridMaskSet& ms) 
    {
      //std::cerr << "subset(const UnitGridCell& c, const UnitGridMaskSet& ms)" << std::endl; 
      //std::cerr << "UnitGridCell: " << c << std::endl;
      return ms.mask()[compute_index(c.lower(),ms.lower(),ms.strides())];
    }
    
    bool 
    interiors_intersect(const UnitGridRectangle& r, const UnitGridMaskSet& ms) 
    {
      //std::cerr << "interiors_intersect(const UnitGridRectangle& r, const UnitGridMaskSet& ms)" << std::endl;
      UnitGridRectangle rstr=regular_intersection(r,ms.bounds());
      if(rstr.empty()) {
        return false;
      }
      //std::cerr << "UnitGridRectangle: " << r << std::endl;
      for(UnitGridRectangle::const_iterator i=rstr.begin(); i!=rstr.end(); ++i) {
        if(subset(*i,ms)) {
          return true;
        }
      }
      return false;
    }
    
    bool 
    subset(const UnitGridRectangle& r, const UnitGridMaskSet& ms) 
    {
      //std::cerr << "subset(const UnitGridRectangle& r, const UnitGridMaskSet& ms)" << std::endl; 
      if(r.empty()) {
        return true;
      }
      if(!subset(r,ms.bounds())) {
        return false;
      }
      //std::cerr << "UnitGridRectangle: " << r << std::endl;
      for(UnitGridRectangle::const_iterator i=r.begin(); i!=r.end(); ++i) {
        if(!subset(*i,ms)) {
          return false;
        }
      }
      return true;
    }
    
    std::ostream& 
    operator<<(std::ostream& os, const UnitGridRectangle& r) 
    {
      if(r.empty() || r.dimension()==0) {
        return os<<"Empty";
      }
      os << r[0];
      for(dimension_type i=1; i!=r.dimension(); ++i) {
        os << "x" << r[i];
      }
      return os;
    }
        
  } 
}
