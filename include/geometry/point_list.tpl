/***************************************************************************
 *            point_list.tpl
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

    template<typename R> 
    void
    PointList<R>::reserve(size_type n)
    {
      if(this->capacity()>=n) { return; }
      LinearAlgebra::Matrix<R> old(this->_pts);
      _pts=LinearAlgebra::Matrix<R>(this->dimension(),n);
      for(size_type j=0; j!=this->size(); ++j) {
        for(size_type i=0; i!=this->dimension(); ++i) {
          _pts(i,j)=old(i,j);
        }
      }
    }    
    
    
    template<typename R>
    LinearAlgebra::Matrix<R> 
    PointList<R>::generators() const 
    { 
      LinearAlgebra::Matrix<R> result(this->dimension(),this->size());
      for(size_type j=0; j!=this->size(); ++j) {
        for(size_type i=0; i!=this->dimension(); ++i) {
          result(i,j)=this->_pts(i,j);
        }
      }
      return result;
    }
    
    
    template<typename R>
    void
    PointList<R>::push_back(const Point<R>& pt) 
    {
      if(this->_size==0) {
        _pts.resize(pt.dimension(),1);
      }
      if(pt.dimension()!=this->dimension()) {
        throw std::runtime_error("Cannot insert point into list of different dimension");
      }
      if(this->size()==this->capacity()) {
        reserve(this->size()*2);
      }
      for(size_type i=0; i!=pt.dimension(); ++i) {
          _pts(i,this->size())=pt[i];
      }  
      ++this->_size;
    }
    
     template<typename R> 
    std::ostream& 
    PointList<R>::write(std::ostream& os) const
    {
      const PointList<R>& pts=*this;
      if(pts.size()==0) { os << "[ "; }
      for(size_type j=0; j!=pts.size(); ++j) {
        os << ((j==0) ? '[' : ',') << pts[j]; 
      }
      return os << ']';
    }
      
  }
}
