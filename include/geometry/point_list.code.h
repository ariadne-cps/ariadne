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
  namespace Geometry {

    
    template<class R> 
    void
    PointList<R>::reserve(size_type n)
    {
      if(this->capacity()>=n) { return; }
      LinearAlgebra::Matrix<R> old(this->_pts);
      _pts=LinearAlgebra::Matrix<R>(this->dimension()+1,n);
      for(size_type j=0; j!=this->size(); ++j) {
        for(size_type i=0; i!=this->dimension()+1u; ++i) {
          _pts(i,j)=old(i,j);
        }
      }
    }    
    
    
    template<class R>
    const LinearAlgebra::Matrix<R>&
    PointList<R>::generators() const 
    { 
      return this->_pts;
    }
    
    
    template<class R>
    void
    PointList<R>::push_back(const Point<R>& pt) 
    {
      if(this->_size==0) {
        _pts.resize(pt.dimension()+1,1);
      }
      check_equal_dimensions(*this,pt,"PointList<R>::push_back(Point<R>)");
      if(this->size()==this->capacity()) {
        reserve(this->size()*2);
      }
      for(size_type i=0; i!=pt.dimension(); ++i) {
          _pts(i,this->size())=pt[i];
      }
      _pts(pt.dimension(),this->size())=R(1);
      ++this->_size;
    }
    
    template<class R> 
    std::ostream& 
    PointList<R>::write(std::ostream& os) const
    {
      const PointList<R>& ptl=*this;
      if(ptl.size()==0) { os << "[ "; }
      for(size_type j=0; j!=ptl.size(); ++j) {
        os << ((j==0) ? '[' : ',') << ptl[j]; 
      }
      return os << ']';
    }
      
  }
}
