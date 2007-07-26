/***************************************************************************
 *            hybrid_set.inline.h
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

#ifndef ARIADNE_HYBRID_SET_ITERATOR_H
#define ARIADNE_HYBRID_SET_ITERATOR_H

#include <boost/iterator.hpp>
#include <boost/iterator_adaptors.hpp>

#include "hybrid_basic_set.h"

namespace Ariadne { 
  namespace Geometry {

  
  template< class DS, class HBS=HybridBasicSet<typename DS::basic_set_type> >
    class HybridDenotableSetIterator
      : public boost::iterator_facade<HybridDenotableSetIterator<DS>,
                                      HBS,
                                      boost::forward_traversal_tag,
                                      HBS const&,
                                      HBS const*
                                     >

    {
     public:
      HybridDenotableSetIterator(const std::map<location_type,DS>&, bool);
      bool equal(const HybridDenotableSetIterator<DS>&) const;
      const HBS& dereference() const;
      void increment();
     private:
      typename std::map< location_type,DS>::const_iterator loc_iter;
      typename std::map< location_type,DS>::const_iterator loc_end;
      typename DS::const_iterator bs_iter;
      HBS set;
    };


  }
}


namespace Ariadne {

template<class DS, class HBS> inline
Geometry::HybridDenotableSetIterator<DS,HBS>::HybridDenotableSetIterator(const std::map<location_type,DS>& map, bool end)
  : loc_iter(map.begin()),
    loc_end(map.end()),
    bs_iter(loc_iter->second.begin()),
    set(loc_iter->first,*bs_iter)
{
  if(end) { loc_iter=loc_end; }
}


template<class DS, class HBS> inline
bool
Geometry::HybridDenotableSetIterator<DS,HBS>::equal(const HybridDenotableSetIterator<DS>& other) const
{
  return this->loc_iter==other.loc_iter && (this->loc_iter==this->loc_end || this->bs_iter==other.bs_iter);
}


template<class DS, class HBS> inline
const HBS&
Geometry::HybridDenotableSetIterator<DS,HBS>::dereference() const
{
  return this->set;
}


template<class DS, class HBS> inline
void
Geometry::HybridDenotableSetIterator<DS,HBS>::increment() 
{
  ++this->bs_iter;
  if(this->bs_iter==loc_iter->second.end()) {
    ++loc_iter;
    if(this->loc_iter!=this->loc_end) {
      this->bs_iter=this->loc_iter->second.begin();
    }
  }
  this->set=HybridBasicSet<typename DS::basic_set_type>(loc_iter->first,*bs_iter);
}

}

#endif // ARIADNE_HYBRID_SET_ITERATOR_H
