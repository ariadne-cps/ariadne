/***************************************************************************
 *            point_list.inline.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it   Pieter.Collins@cwi.nl
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

#include "base/iterator.h"

namespace Ariadne {
  namespace Geometry {

    template<class X>
    class PointListIterator
      : public boost::iterator_facade<PointListIterator<X>,
                                      Point<X>,
                                      boost::forward_traversal_tag,
                                      Point<X> const,
                                      Point<X> const*
                                     >
    {
     public:
      PointListIterator(const PointList<X>* ptl, size_type j) : _ptl(ptl), _j(j) { }
      const Point<X> dereference() const { return (*_ptl)[_j]; }
      void increment() { ++_j; }
      bool equal(const PointListIterator& other) const { 
        return (this->_j==other._j) && (this->_ptl==other._ptl); }
     private:
      const PointList<X>* _ptl;
      size_type _j;
    };

    template<class X>
    class PointListReference
    {
     public:
      PointListReference(PointList<X>& pl, size_type j) : _pl(pl), _j(j) { }
      void operator=(const Point<X>& pt) {
        for(dimension_type i=0; i!=_pl.dimension(); ++i) { _pl._pts(i,_j)=pt[i]; } }
      X& operator[] (size_type i) { return _pl._pts(i,_j); }
      operator Point<X>() const { return _pl.get(_j); }
     private:
      PointList<X>& _pl;
      size_type _j;
    };
    
   } // namespace Geometry
} // namespace Ariadne



namespace Ariadne {

template<class X> inline
Geometry::PointList<X>::PointList(dimension_type d)
  : _size(0), _pts(d+1,1) 
{ 
}

template<class X> inline
Geometry::PointList<X>::PointList(dimension_type d, size_type n)
  : _size(n), _pts(d+1,n) 
{ 
}

template<class X> inline
Geometry::PointList<X>::PointList(const LinearAlgebra::Matrix<X>& G)
  : _size(G.number_of_columns()), _pts(G) 
{ 
}

template<class X> template<class XX> inline
Geometry::PointList<X>::PointList(const PointList<XX>& ptl)
  : _size(ptl.size()), _pts(ptl.generators()) 
{ 
}

template<class X> inline
Geometry::PointList<X>::PointList(const PointList<X>& ptl)
  : _size(ptl._size), _pts(ptl.generators()) 
{ 
}

template<class X> inline
Geometry::PointList<X>& 
Geometry::PointList<X>::operator=(const PointList<X>& ptl)
{
  if(this!=&ptl) { 
    this->_size=ptl.size(); 
    this->_pts=ptl.generators(); 
  } 
  return *this; 
}

template<class X> inline
dimension_type 
Geometry::PointList<X>::dimension() const
{
  return _pts.number_of_rows()-1; 
}

template<class X> inline
size_type 
Geometry::PointList<X>::size() const
{
  return _size; 
}

template<class X> inline
size_type 
Geometry::PointList<X>::capacity() const
{
  return _pts.number_of_columns(); 
}


template<class X> inline
Geometry::Point<X> 
Geometry::PointList<X>::operator[] (const size_type j) const
{
  Point<X> result(this->dimension());
  for(size_type i=0; i!=this->dimension(); ++i) {
    result[i]=this->_pts(i,j);
  }
  return result;
}

template<class X> inline
Geometry::PointListReference<X>  
Geometry::PointList<X>::operator[](const size_type j)
{ 
  return PointListReference<X>(*this,j);
}


template<class X> inline
void 
Geometry::PointList<X>::pop_back()
{
  --this->_size; 
}

template<class X> inline
void 
Geometry::PointList<X>::adjoin(const Point<X>& pt)
{
  this->push_back(pt); 
}
template<class X> inline
void 
Geometry::PointList<X>::clear()
{
  this->_size=0; 
}


template<class X> inline
typename Geometry::PointList<X>::const_iterator 
Geometry::PointList<X>::begin() const
{
  return const_iterator(this,0); 
}

template<class X> inline
typename Geometry::PointList<X>::const_iterator 
Geometry::PointList<X>::end() const
{
  return const_iterator(this,this->_size); 
}



template<class X> inline
void 
Geometry::PointList<X>::_set(size_type j, dimension_type i,const X& x)
{
  _pts(i,j)=x; 
}

template<class X> inline
Geometry::Point<X> 
Geometry::PointList<X>::_get(size_type j) const
{
  return Point<X>(_pts.column(j)); 
}

template<class X> inline
X 
Geometry::PointList<X>::_get(size_type j, dimension_type i) const 
{
  return _pts(i,j); 
}



template<class X> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const PointList<X>& pts)
{
  return pts.write(os);
}



} // namespace Ariadne
