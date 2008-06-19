/***************************************************************************
 *            list_set.inline.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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

#include "geometry/geometrical_traits.h"
#include "base/exceptions.h"

namespace Ariadne {

template<class BS> inline
ListSet<BS>::ListSet()
  : _dimension(0), _vector()
{
}

template<class BS> inline
ListSet<BS>::ListSet(dimension_type n) 
  : _dimension(n), _vector() 
{
}

template<class BS> inline
ListSet<BS>::ListSet(const BS& bs)
  : _dimension(bs.dimension()), _vector()
{
  _vector.push_back(bs);
}

template<class BS> template<class BST> inline
ListSet<BS>::ListSet(const ListSet<BST>& ls)
  : _dimension(ls.dimension()), _vector() 
{
  _vector.reserve(ls.size());
  for(typename ListSet<BST>::const_iterator bst_iter=ls.begin();
      bst_iter!=ls.end(); ++bst_iter) {
    this->_vector.push_back(BS(*bst_iter));
  }
}

template<class BS> template<class Iter> inline
ListSet<BS>::ListSet(Iter curr, Iter end)
  : _vector(curr,end) 
{
  if(!_vector.empty()) {
    this->_dimension=_vector.front().dimension();
  }
}


template<class BS> inline
ListSet<BS>::~ListSet() 
{
  this->_vector.clear();
}


template<class BS> inline
size_type 
ListSet<BS>::size() const 
{
  return this->_vector.size();
}

template<class BS> inline
BS 
ListSet<BS>::pop() 
{
  if (this->_vector.empty()) { 
    ARIADNE_THROW(LengthError,"void ListSet<BS>::pop_back()"," empty list");
  }
  BS result=this->_vector.back();
  this->_vector.pop_back();
  return result;
}


template<class BS> inline
dimension_type 
ListSet<BS>::dimension() const 
{
  return this->_dimension;
}

template<class BS> inline
const BS& 
ListSet<BS>::operator[](size_type index) const 
{
  ARIADNE_CHECK_ARRAY_INDEX(*this,index,"void ListSet<BS>::operator[](size_type index)");
  return this->_vector[index];
}



template<class BS> inline
void 
ListSet<BS>::clear() 
{ 
  this->_vector.clear();
}


template<class BS> inline
typename ListSet<BS>::const_iterator 
ListSet<BS>::begin() const 
{
  return _vector.begin();
}

template<class BS> inline
typename ListSet<BS>::const_iterator 
ListSet<BS>::end() const 
{
  return _vector.end();
}

template<class BS> inline
void 
ListSet<BS>::adjoin(const ListSet<BS>& ls) 
{
  if(ls.size()==0) { return; }
  if(this->dimension()==0) { this->_dimension=ls.dimension(); }
  ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,ls,"void ListSet<BS>::adjoin(ListSet<BS> ls)");
  this->_vector.reserve(ls.size());
  for(typename ListSet<BS>::const_iterator iter=ls.begin(); iter!=ls.end(); ++iter) {
    this->_vector.push_back(*iter);
  }
}

template<class BS> inline
void 
ListSet<BS>::adjoin(const BS& bs) 
{
  if(this->dimension()==0) { 
    this->_dimension=bs.dimension(); 
  }
  ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,bs,"void ListSet<BS>::adjoin(BS bs)");
  this->_vector.push_back(bs);
}


template<class BS> template<class S> inline
void 
ListSet<BS>::adjoin(const S& s) 
{
  this->adjoin(s, typename geometrical_traits<S>::set_category());
}


template<class BS> template<class S> inline
void 
ListSet<BS>::adjoin(const S& bs, basic_set_tag) 
{
  this->adjoin(static_cast<BS>(bs));
}


template<class BS> template<class DS> inline
void 
ListSet<BS>::adjoin(const DS& ds, denotable_set_tag) 
{
  for(typename DS::const_iterator iter=ds.begin();
      iter!=ds.end(); ++iter)
    {
      this->adjoin(static_cast<BS>(*iter));
    }
}




template<class BS> inline
std::ostream& 
operator<<(std::ostream& os, const ListSet<BS>& ls)
{
  return ls.write(os);
}

template<class BS> inline
std::istream& 
operator>>(std::istream& is, ListSet<BS>& ls)
{
  return ls.read(is);
}


} // namespace Ariadne
