/***************************************************************************
 *            point_list.code.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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

#include <iostream>
#include <stdexcept>

#include "point_list.h"

namespace Ariadne {


    
template<class X> 
void
PointList<X>::reserve(size_type n)
{
  if(this->capacity()>=n) { return; }
  Matrix<X> old(this->_pts);
  _pts=Matrix<X>(this->dimension()+1,n);
  for(size_type j=0; j!=this->size(); ++j) {
    for(size_type i=0; i!=this->dimension()+1u; ++i) {
      _pts(i,j)=old(i,j);
    }
  }
}    


template<class X>
const Matrix<X>&
PointList<X>::generators() const 
{ 
  return this->_pts;
}


template<class X>
void
PointList<X>::push_back(const Point<X>& pt) 
{
  if(this->_size==0) {
    _pts.resize(pt.dimension()+1,1);
  }
  ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,pt,"PointList<X>::push_back(Point<X>)");
  if(this->size()==this->capacity()) {
    reserve(this->size()*2);
  }
  for(size_type i=0; i!=pt.dimension(); ++i) {
    _pts(i,this->size())=pt[i];
  }
  _pts(pt.dimension(),this->size())=static_cast<X>(1);
  ++this->_size;
}

template<class X> 
std::ostream& 
PointList<X>::write(std::ostream& os) const
{
  const PointList<X>& ptl=*this;
  if(ptl.size()==0) { os << "[ "; }
  for(size_type j=0; j!=ptl.size(); ++j) {
    os << ((j==0) ? '[' : ',') << ptl[j]; 
  }
  return os << ']';
}


} // namespace Ariadne
