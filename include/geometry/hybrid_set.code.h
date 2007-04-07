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

namespace Ariadne { namespace Geometry {

template<class S> 
std::ostream& 
HybridSetBase<S>::write(std::ostream& os) const
{ 
  return os << "HybridSet"<<_component_sets;
}


template<class R> 
std::ostream& 
HybridGridMaskSet<R>::write(std::ostream& os) const
{ 
  os << "HybridGridMaskSet( { \n";
  for(typename HybridGridMaskSet<R>::const_iterator iter=this->begin();
      iter!=this->end(); ++iter)
  {
    id_type loc=iter->first;
    const GridMaskSet<R>& set=iter->second;
    os << "  "<<loc<<": GridMaskSet( extent=" << set.extent() << ", block=" << set.block() << ", size=" << set.size() << " ),\n";
  }
  os << "} )";
  return os;
}

}}
