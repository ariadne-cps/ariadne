/***************************************************************************
 *            map.h
 *
 *  Wed Feb  2 18:33:10 2005
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
 
#ifndef _ARIADNE_MAP_H
#define _ARIADNE_MAP_H

#include "point.h"
#include "list_set.h"

namespace Ariadne {
namespace Map{

enum MapKind {
  LINEAR,
  AFFINE,
  GENERAL
};

enum MapResultKind {
  SINGLE_VALUE,
  MULTI_VALUE
};

}
}

namespace Ariadne {
namespace Map {

template <typename R, template<typename> class S>
class Map {
 public:
  typedef R Real;
  typedef S<R> State;
  
  virtual State apply(const State& x) const = 0;

  // I'm not sure if virtual functions are the way to go here; too restrictive? //
  virtual Geometry::Rectangle<R> apply(const Geometry::Rectangle<R>& A) const {
    throw std::invalid_argument("Not implemented."); }
  virtual Geometry::Polyhedron<R> apply(const Geometry::Polyhedron<R>& A) const {
    throw std::invalid_argument("Not implemented."); }

  template<template<typename> class BS>
  inline BS<R> operator() (const BS<R>& A) const { 
    return apply(*this,A);
  }

  template<template<typename> class BS>
  inline Geometry::ListSet<R,BS> operator() (const Geometry::ListSet<R,BS>& A) const { 
    Geometry::ListSet<R,BS> trans_ds(A.dimension());
    for (size_t i=0; i< A.size(); i++) {
      trans_ds.inplace_union(this->apply(A[i]));
    }
    return trans_ds;
  }
  
  virtual size_t dimension() const = 0;
};
  
  
/* WARNING!!!! Is it the same of an inclusion map? Here I can set
 * threshold in inclusion not, but the map seem similar 
 */ 
template <typename MAP>
class ThresholdMap {
  
 public:
  typedef MAP Map;
  typedef typename Map::DenotableSet DenotableSet;
  typedef typename DenotableSet::BasicSet BasicSet;
  typedef typename BasicSet::State State;
  typedef typename State::Real Real;
  
  ThresholdMap(const Map& T, const BasicSet& threshold)
    : _map(T), _threshold(threshold) {}
  
  ThresholdMap(const Map& T, Real threshold = 0)
    : _map(T), _threshold(T.dimension(), threshold) {}
  
  inline BasicSet operator() (const BasicSet& A) const { 
    return _map(A)+this->_threshold; 
  }
  
  inline void set_threshold(const Real& threshold) { 
    BasicSet new_th((this._map).dimension(), threshold);
    this->_threshold=new_th;
  }
  
  inline DenotableSet operator() (const DenotableSet& A) const{ 
    DenotableSet trans_ds(A.dimension());
    for (size_t i=0; i< A.size(); i++) {
      trans_ds.inplace_union(_map(A[i])+this->threshold);
    }
    return trans_ds;
  }
  
  inline ThresholdMap<MAP>& operator=(const ThresholdMap<MAP>& A) {
    this->_map=A._map;
    this->_threshold=A._threshold;
  }
  
  inline size_t dimension() const {
    return (this->_map).dimension();
  }
  
  /*! Deprecated */
  inline size_t dim() const {
    return (this->_map).dimension();
  }
 private:
  Map _map;
  BasicSet _threshold;
};
  
}
}

#endif /* _ARIADNE_MAP_H */
