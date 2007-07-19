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
#include "rectangle.h"
#include "zonotope.h"
#include "list_set.h"

namespace Ariadne { 


template<class S> 
std::ostream& 
Geometry::HybridSetBase<S>::write(std::ostream& os) const
{ 
  os << "HybridSet( { \n";
  for(const_iterator iter=this->begin(); iter!=this->end(); ++iter)
  {
    id_type loc=iter->first;
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
  for(typename HybridGridMaskSet<R>::const_iterator iter=this->begin();
      iter!=this->end(); ++iter)
  {
    id_type loc=iter->first;
    const GridMaskSet<R>& set=*iter->second;
    os << "  "<<loc<<": GridMaskSet( grid=" << set.grid() << ", extent=" << set.extent() << ", block=" << set.block() << ", size=" << set.size() << " capacity=" << set.capacity() << " ),\n";
  }
  os << "} )";
  return os;
}


template<class BS> 
std::ostream& 
Geometry::HybridListSet<BS>::write(std::ostream& os) const
{ 
  os << "HybridListSet<"<<BS::name()<<">( { \n";
  for(typename HybridListSet<BS>::const_iterator iter=this->begin();
      iter!=this->end(); ++iter)
  {
    id_type loc=iter->first;
    const ListSet<BS>& set=*iter->second;
    os << "  " << loc << ": { size=" << set.size();
    if(!set.empty()) {
      os << ", front=" << set[0];
    }
    os << "},\n";
  }
  os << "} )";
  return os;
}

} // namespace Ariadne
