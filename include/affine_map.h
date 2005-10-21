/***************************************************************************
 *            affine_map.h
 *
 *  Wed Feb  2 18:52:36 2005
 *  Copyright  2005  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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
 
#ifndef _AFFINE_MAP_H
#define _AFFINE_MAP_H

#include "map.h"
#include "state.h"
#include "linear_algebra.h"

namespace Ariadne {	
namespace Map{

enum AffineKind {
  HOMOGENEOUS,
  NON_HOMOGENEOUS,
  TRASLATION
};
   
template <typename R>
class AffineMap // : public Map<R,Geometry::State> 
{
 public:
  typedef Geometry::State<R> State;
  typedef R Real;
  
  typedef typename boost::numeric::ublas::matrix<Real> Matrix;
  typedef typename boost::numeric::ublas::vector<Real> Vector;
  
  inline explicit AffineMap() {}
  inline explicit AffineMap(const Matrix& A, const Vector& b) : _A(A), _b(b) { }
  inline explicit AffineMap(const Matrix& A) : _A(A), _b(A.columns()) { }
  inline explicit AffineMap(const Vector& b) : _A(b.size(),b.size()), _b(b) { }
  
  inline AffineMap(const AffineMap<Real>& T) : _A(T._A), _b(T._b) { }
  inline AffineMap<Real>& operator=(const AffineMap<Real>& T) {
    this->_A=T._A; this->_b=T._b; return *this; }
  
  State operator() (const State& x) const;
    
  template<template <typename> class BS>
  BS<R> operator() (const BS<R>& A) const;
  
  inline const Matrix& A() const { return _A; }
  inline const Vector& b() const { return _b; }
  
  inline size_t dimension() const {
    return _b.dimension();
  }
  
  /*! Deprecated. */ 
  inline size_t dim() const {
    return this->dimension();
  }
  
  inline bool invertible() const {
    return _A.invertible(); }
 private:
  Matrix _A;
  Vector _b;
};
  
template <typename R>
template <template<typename> class BS>
BS<R>
AffineMap<R>::operator() (const BS<R>& bs) const
{
  return apply(*this,bs);
}
 
}
}


#endif /* _AFFINE_MAP_H */
