/***************************************************************************
 *            point_list.h
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

/*! \file point_list.h
 *  \brief A list of points in Euclidean space.
 */

#ifndef _ARIADNE_POINT_LIST_H
#define _ARIADNE_POINT_LIST_H

#include <iostream>
#include <stdexcept>

#include "../declarations.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../geometry/point.h"

namespace Ariadne {
  namespace Geometry {

    template<typename R>
    class PointListIterator
    {
     public:
      PointListIterator(const LinearAlgebra::Matrix<R>& A, size_type j) : _A(A), _j(j) { }
      Point<R> operator*() const { return Point<R>(column(_A,_j)); }
      void operator++() { ++_j; }
      bool operator==(const PointListIterator& other) { 
        return (this->_j==other._j) && (&(this->_A)==&(other._A)); }
      bool operator!=(const PointListIterator& other) { 
        return !(*this==other); }
     private:
      const LinearAlgebra::Matrix<R>& _A;
      size_type _j;
    };
    
    template<typename R>
    class PointListReference
    {
     public:
      PointListReference(PointList<R>& pl, size_type j) : _pl(pl), _j(j) { }
      void operator=(const Point<R>& pt) {
        for(dimension_type i=0; i!=_pl.dimension(); ++i) { _pl._pts(i,_j)=pt[i]; } }
      R& operator[] (size_type i) { return _pl._pts(i,_j); }
      operator Point<R>() const { return _pl.get(_j); }
     private:
      PointList<R>& _pl;
      size_type _j;
    };
    
    template<typename R>
    class PointList
    {
      friend class PointListReference<R>;
      friend class PointListIterator<R>;
     public:
      typedef PointListIterator<R> const_iterator;
      
      PointList() : _size(0), _pts(0,0) { }
      PointList(dimension_type d) : _size(0), _pts(d,1) { }
      PointList(dimension_type d, size_type n) : _size(n), _pts(d,n) { }
      PointList(const LinearAlgebra::Matrix<R>& g) : _size(g.size2()), _pts(g) { }
      PointList(const PointList& ptl) : _size(ptl._size), _pts(ptl.generators()) { }
      PointList& operator=(const PointList& ptl) {
        if(this!=&ptl) { this->_size=ptl.size(); this->_pts=ptl.generators(); } 
        return *this; }
      dimension_type dimension() const { return _pts.size1(); }
      size_type size() const { return _size; }
      size_type capacity() const { return _pts.size2(); }
      LinearAlgebra::Matrix<R> generators() const { 
        LinearAlgebra::Matrix<R> result(this->dimension(),this->size());
        for(size_type j=0; j!=this->size(); ++j) {
          for(size_type i=0; i!=this->dimension(); ++i) {
            result(i,j)=this->_pts(i,j);
          }
        }
        return result;
      }
      LinearAlgebra::Matrix<R> matrix() const { return this->generators(); }
      operator LinearAlgebra::Matrix<R> () const { return this->generators(); }
      Point<R> operator[] (const size_type j) const { 
        return Point<R>(column(_pts,j)); }
      PointListReference<R> operator[] (const size_type j) { 
        return PointListReference<R>(*this,j);}
      void set(size_type j, dimension_type i,const R& x) { _pts(i,j)=x; }
      Point<R> get(size_type j) const { return Point<R>(column(_pts,j)); }
      R get(size_type j, dimension_type i) const { return _pts(i,j); }
      void push_back(const Point<R>& pt) {
        if(_pts.size2()==0) {
          _pts=LinearAlgebra::Matrix<R>(pt.dimension(),1);
        } else {
          if(pt.dimension()!=_pts.size1()) {
            throw std::runtime_error("Cannot insert point into list of different dimension");
          }
        }
        if(this->size()==this->capacity()) {
          reserve(this->size()*2);
        }
        for(size_type i=0; i!=pt.dimension(); ++i) {
            _pts(i,this->size())=pt[i];
        }
        ++this->_size;
      }
      void pop_back() { --this->_size; }
      void reserve(size_type n) {
        if(this->capacity()>=n) { return; }
        LinearAlgebra::Matrix<R> old(this->_pts);
        _pts=LinearAlgebra::Matrix<R>(this->dimension(),n);
        for(size_type j=0; j!=this->size(); ++j) {
          for(size_type i=0; i!=this->dimension(); ++i) {
            _pts(i,j)=old(i,j);
          }
        }
      }
      const_iterator begin() const { return const_iterator(this->_pts,0); }
      const_iterator end() const { return const_iterator(this->_pts,this->_size); }
     private:
      size_type _size;
      LinearAlgebra::Matrix<R> _pts;
    };
    
    template<typename R> std::ostream& operator<<(std::ostream& os, const PointList<R>& pts);
    
    template<typename R> 
    inline
    std::ostream& 
    operator<<(std::ostream& os, const PointList<R>& pts) 
    {
      if(pts.size()==0) { os << "[ "; }
      for(size_type j=0; j!=pts.size(); ++j) {
        os << ((j==0) ? '[' : ',') << pts[j]; 
      }
      return os << ']';
    }
      
  }
}

#endif /* _ARIADNE_POINT_LIST_H */
