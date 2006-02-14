/***************************************************************************
 *            map.h
 *
 *  Wed Feb  2 18:33:10 2005
 *  Copyright  2005, 2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file map.h
 *  \brief Base class for maps.
 */

#ifndef _ARIADNE_MAP_H
#define _ARIADNE_MAP_H

#include "geometry/geometry_declarations.h"

#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelopiped.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

namespace Ariadne {
namespace Evaluation {

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
namespace Evaluation {

template <typename R>
class Map {
 public:
  typedef R Real;
  typedef Geometry::Point<R> State;
  typedef LinearAlgebra::vector<R> Vector;
  typedef LinearAlgebra::matrix<R> Matrix;
  typedef LinearAlgebra::matrix< Interval<R> > IntervalMatrix;
  
  virtual ~Map() { }
  
  // I'm not sure if virtual functions are the way to go here; too restrictive? //
  virtual State apply(const State& x) const {
    throw std::invalid_argument("Not implemented."); }

  virtual Geometry::Rectangle<R> apply(const Geometry::Rectangle<R>& A) const {
    throw std::invalid_argument("Not implemented."); }
  virtual Geometry::Polyhedron<R> apply(const Geometry::Polyhedron<R>& A) const {
    throw std::invalid_argument("Not implemented."); }

  virtual Matrix derivative(const State& r) const {
    throw std::invalid_argument("Derivative at point not implemented."); }
  virtual LinearAlgebra::matrix< Interval<R> > derivative(const Geometry::Rectangle<R>& r) const {
    throw std::invalid_argument("Derivative on Rectangle not implemented."); }
    
  template<template<typename> class BS>
  inline Geometry::ListSet<R,BS> apply(const Geometry::ListSet<R,BS>& A) const { 
    Geometry::ListSet<R,BS> trans_ds(A.dimension());
    for (size_t i=0; i< A.size(); i++) {
      trans_ds.inplace_union(this->apply(A[i]));
    }
    return trans_ds;
  }
  
  virtual dimension_type argument_dimension() const = 0;
  virtual dimension_type result_dimension() const = 0;
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
