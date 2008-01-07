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

  // FIXME: This class doesn't work

  template< class DS, class HBS=HybridBasicSet<typename DS::basic_set_type> >
    class HybridDenotableSetIterator
      : public boost::iterator_facade<HybridDenotableSetIterator<DS>,
                                      HBS,
                                      boost::forward_traversal_tag,
                                      HBS
                                     >

    {
     public:
      HybridDenotableSetIterator(const std::map<DiscreteState,DS>&, bool);
      bool equal(const HybridDenotableSetIterator<DS>&) const;
      HBS dereference() const;
      void increment();
      const id_type& discrete_state() const { return loc_iter->first; }
      const typename DS::basic_set_type& continuous_state_set() const { return *bs_iter; }
     private:
      void increment_loc();
     private:
      typename std::map< DiscreteState,DS>::const_iterator loc_begin;
      typename std::map< DiscreteState,DS>::const_iterator loc_end;
      typename std::map< DiscreteState,DS>::const_iterator loc_iter;
      typename DS::const_iterator bs_iter;
    };


  }
}


namespace Ariadne {

template<class DS, class HBS> inline
Geometry::HybridDenotableSetIterator<DS,HBS>::HybridDenotableSetIterator(const std::map<DiscreteState,DS>& map, bool end)
  : loc_begin(map.begin()),
    loc_end(map.end()),
    loc_iter(end?loc_end:loc_begin),
    bs_iter()
{
  if(loc_iter!=loc_end) {
    bs_iter=loc_iter->second.begin();
    this->increment_loc();
  }
}


template<class DS, class HBS> inline
bool
Geometry::HybridDenotableSetIterator<DS,HBS>::equal(const HybridDenotableSetIterator<DS>& other) const
{
  return this->loc_iter==other.loc_iter && (this->loc_iter==this->loc_end || this->bs_iter==other.bs_iter);
}


template<class DS, class HBS> inline
HBS
Geometry::HybridDenotableSetIterator<DS,HBS>::dereference() const
{
  return HBS(loc_iter->first,*this->bs_iter);
}


template<class DS, class HBS> inline
void
Geometry::HybridDenotableSetIterator<DS,HBS>::increment() 
{
  ++this->bs_iter;
  this->increment_loc();
}

template<class DS, class HBS> inline
void
Geometry::HybridDenotableSetIterator<DS,HBS>::increment_loc() 
{
  while(bs_iter==loc_iter->second.end()) {
    ++loc_iter;
    if(loc_iter==loc_end) { return; } 
    bs_iter=loc_iter->second.begin();
  }
}

}

#endif // ARIADNE_HYBRID_SET_ITERATOR_H
