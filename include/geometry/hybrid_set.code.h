/***************************************************************************
 *            hybrid_set.code.h
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

#include "hybrid_set.h"
#include "box.h"
#include "zonotope.h"
#include "list_set.h"

#include "output/textstream.h"

namespace Ariadne { 


template<class S> 
std::ostream& 
Geometry::HybridSet<S>::write(std::ostream& os) const
{ 
  os << "HybridSet( { \n";
  for(locations_const_iterator iter=this->locations_begin(); iter!=this->locations_end(); ++iter)
  {
    discrete_state_type loc=iter->first;
    const S& set=*iter->second;
    os << "  "<<loc<<": " << set << ",\n";
  }
  os << "} )";
  return os;
}



template<class R> 
std::ostream& 
Geometry::HybridGridMaskSet<R>::write(std::ostream& os) const
{ 
  os << "HybridGridMaskSet( { \n";
  for(typename HybridGridMaskSet<R>::locations_const_iterator iter=this->locations_begin();
      iter!=this->locations_end(); ++iter)
  {
    DiscreteState loc=iter->first;
    const GridMaskSet<R>& set=iter->second;
    os << "  "<<loc<<": " << set.summary() << ",\n";
  }
  os << "} )";
  return os;
}


template<class BS> 
std::ostream& 
Geometry::HybridListSet<BS>::write(std::ostream& os) const
{ 
  os << "HybridListSet<"<<Geometry::name<BS>()<<">( { \n";
  for(typename HybridListSet<BS>::locations_const_iterator iter=this->locations_begin();
      iter!=this->locations_end(); ++iter)
  {
    DiscreteState loc=iter->first;
    const ListSet<BS>& set=iter->second;
    os << "  " << loc << ": " << set.summary() << ",\n";
    /*
    os << "  " << loc << ": { size=" << set.size();
    if(!set.empty()) {
      os << ", front=" << set[0];
    }
    os << "},\n";
    */
  }
  os << "} )";
  return os;
}

} // namespace Ariadne
