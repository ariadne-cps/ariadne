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

#include <map.h>
#include <denotable_set.h>
#include <linear_algebra.h>

namespace Ariadne {	
namespace Map{

enum AffineKind {
  HOMOGENEOUS,
  NON_HOMOGENEOUS,
  TRASLATION
};
   
template <typename BS_MAP>
class AffineMap {
 public:
  typedef BS_MAP BasicSetMap;
  typedef typename BasicSetMap::DenotableSet DenotableSet;
  typedef typename DenotableSet::BasicSet BasicSet;
  typedef typename BasicSet::State State;
  typedef typename State::Real Real;
  
  typedef typename boost::numeric::ublas::matrix<Real> Matrix;
  typedef typename boost::numeric::ublas::vector<Real> Vector;
  
  AffineMap() {}
  
  AffineMap(const AffineMap<BS_MAP> &T):
    _map(T._map), _B(T._B), _inclusion_map(T._inclusion_map){}
  
  AffineMap(const Matrix &A):_map(A), _inclusion_map(false) {}
  
  AffineMap(const Vector &b): _map(b), _inclusion_map(false) {}
  
  AffineMap(const Matrix &A, const Vector &b):
    _map(A,b), _inclusion_map(false) {}
  
  AffineMap(const Matrix &A, const BasicSet &B):
    _map(A), _B(B), _inclusion_map(true) {}
  
  inline AffineMap<BS_MAP> &operator=(const AffineMap<BS_MAP>  &A) {
    
    this->_map=A._map;
    this->_B=A._B;
    this->_inclusion_map=A._inclusion_map;
    
    return *this;
  }
  
  inline BasicSet operator()(const BasicSet &A) const { 
    
    if (this->_inclusion_map) {
      return (_map(A))+this->_B;
    }
    
    return _map(A); 
  }
  
  inline DenotableSet operator()(const DenotableSet &A) const{ 
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    DenotableSet trans_ds(A.dimension());
    
    for (size_t i=0; i< A.size(); i++) {	
      
      trans_ds.inplace_union(_map(A[i]));
    }
    
    if (this->_inclusion_map) {
      trans_ds=trans_ds+this->_B;
    }
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    return trans_ds;
  }
  
  inline const Matrix& A() const { return _map.A(); }
  
  inline const Vector& b() const { return _map.b(); }
  
  inline const BasicSet& B() const { 
    
    if (!this->_inclusion_map) {
      throw "This is not an inclusion map";
    }
    
    return this->_B; 
  }
  
  inline AffineKind map_kind() { return AFFINE; }
  
  inline const BasicSetMap &basicset_map() const{
    return (this->_map);
  }
  
  inline size_t dimension() const {
    return (this->_map).dimension();
  }
  
  /*! Deprecated. */ 
  inline size_t dim() const {
    return (this->_map).dimension();
  }
  
  inline bool invertible() const {
    return ((!this->_inclusion_map)&&((this->_map).invertible()));
  }
  
 private:
  BasicSetMap _map;
  BasicSet _B;
  bool _inclusion_map;
};
  
  
}
}


#endif /* _AFFINE_MAP_H */
