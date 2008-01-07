/***************************************************************************
 *            box_list_set.inline.h
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

template<class R> inline
Geometry::BoxListSet<R>::BoxListSet()
  : _vector()
{
}

template<class R> inline
Geometry::BoxListSet<R>::BoxListSet(dimension_type d)
  : _vector()
{
}



template<class R> inline
size_type 
Geometry::BoxListSet<R>::size() const 
{
  return this->_vector.size();
}

template<class R> inline
void 
Geometry::BoxListSet<R>::push_back(const Box<R>& bx) 
{
  if (this->size()!=0) { 
    ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,bx,"void BoxListSet<R>::push_back(Box bx)");
  }
  this->_vector.push_back(bx);
}

template<class R> inline
void 
Geometry::BoxListSet<R>::pop_back() 
{
  if (this->_vector.empty()) { 
    ARIADNE_THROW(LengthError,"void BoxListSet<R>::pop_back()"," empty list");
  }
  this->_vector.pop_back();
}

template<class R> inline
dimension_type 
Geometry::BoxListSet<R>::dimension() const 
{
  if(this->_vector.empty()) {
    return 0;
  } else {
    return this->_vector.front().dimension();
  }
}

template<class R> inline
const Geometry::Box<R>& 
Geometry::BoxListSet<R>::get(size_type index) const 
{
  ARIADNE_CHECK_ARRAY_INDEX(*this,index,"Box BoxListSet::get(size_type index)");
  return this->_vector[index];
}

template<class R> inline
void 
Geometry::BoxListSet<R>::set(size_type index, const Box<R>& set) 
{
  ARIADNE_CHECK_ARRAY_INDEX(*this,index,"void BoxListSet::set(size_type index, Box set)");
  this->_vector[index]=set;
}

template<class R> inline
const Geometry::Box<R>& 
Geometry::BoxListSet<R>::operator[](size_type index) const 
{
  ARIADNE_CHECK_ARRAY_INDEX(*this,index,"void BoxListSet<R>::operator[](size_type index)");
  return this->_vector[index];
}



template<class R> inline
tribool 
Geometry::BoxListSet<R>::contains(const Point<R>& p) const 
{
  tribool result=false;
  for (typename BoxListSet<R>::const_iterator i=this->begin(); i!=this->end(); ++i) {
    result=result || i->contains(p);
    if(result) { return result; }
  }
  return result;
}

template<class R> inline
tribool 
Geometry::BoxListSet<R>::empty() const 
{
  tribool result=true;
  for (typename BoxListSet<R>::const_iterator i=this->begin(); i!=this->end(); ++i) {
    result = result && i->empty();
    if(!result) { return result; }
  }
  return result;
}

template<class R> inline
tribool 
Geometry::BoxListSet<R>::bounded() const 
{
  tribool result=true;
  for (typename BoxListSet<R>::const_iterator i=this->begin(); i!=this->end(); ++i) {
    result = result && i->bounded();
    if(!result) { return result; }
  }
  return result;
}

template<class R> inline
Geometry::Box<typename Geometry::BoxListSet<R>::real_type> 
Geometry::BoxListSet<R>::bounding_box() const 
{
  if(this->empty()) { return Box<R>(this->dimension()); }
  Box<R> result=(*this)[0];
  for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
    result=rectangular_hull(result,*iter);
  }
  return result;
}

template<class R> inline
void 
Geometry::BoxListSet<R>::clear() 
{ 
  this->_vector.clear();
}

template<class R> inline
typename Geometry::BoxListSet<R>::iterator 
Geometry::BoxListSet<R>::begin() 
{
  return _vector.begin();
}

template<class R> inline
typename Geometry::BoxListSet<R>::iterator 
Geometry::BoxListSet<R>::end() 
{
  return _vector.end();
}

template<class R> inline
typename Geometry::BoxListSet<R>::const_iterator 
Geometry::BoxListSet<R>::begin() const 
{
  return _vector.begin();
}

template<class R> inline
typename Geometry::BoxListSet<R>::const_iterator 
Geometry::BoxListSet<R>::end() const 
{
  return _vector.end();
}

template<class R> inline
void 
Geometry::BoxListSet<R>::adjoin(const BoxListSet<R>& ls) 
{
  if(ls.size()==0) { return; }
  if(this->_vector.empty()) { *this=ls; }
  ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,ls,"void BoxListSet<R>::adjoin(Geometry::BoxListSet<R> ls)");
  this->_vector.reserve(ls.size());
  for(typename BoxListSet<R>::const_iterator iter=ls.begin(); iter!=ls.end(); ++iter) {
    this->_vector.push_back(*iter);
  }
}

template<class R> inline
void 
Geometry::BoxListSet<R>::adjoin(const Box<R>& bx) 
{
  if(!this->_vector.empty()) {
    ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,bx,"void BoxListSet<R>::adjoin(BS bs)");
  }
  this->_vector.push_back(bx);
}






template<class R> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const BoxListSet<R>& ls)
{
  return ls.write(os);
}

template<class R> inline
std::istream& 
Geometry::operator>>(std::istream& is, BoxListSet<R>& ls)
{
  return ls.read(is);
}


} // namespace Geometry
