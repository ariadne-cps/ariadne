/***************************************************************************
 *            map.tpl
 *
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
 
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"
#include "../geometry/parallelotope.h"
#include "../geometry/polyhedron.h"

#include "../evaluation/map.h"

namespace Ariadne {
  namespace Evaluation {

    template<typename R>
    Map<R>::~Map() 
    {
    }
  
    template<typename R>
    Geometry::Point<R> 
    Map<R>::apply(const Geometry::Point<R>& x) const 
    {
      throw std::invalid_argument("Not implemented."); 
    }
    
    template<typename R>
    Geometry::Rectangle<R>
    Map<R>::apply(const Geometry::Rectangle<R>& r) const 
    {
      throw std::invalid_argument("Not implemented."); 
    }
    
    template<typename R>
    Geometry::Parallelotope<R>
    Map<R>::apply(const Geometry::Parallelotope<R>& p) const 
    {
      throw std::invalid_argument("Not implemented."); 
    }
    
    template<typename R>
    Geometry::Polyhedron<R>
    Map<R>::apply(const Geometry::Polyhedron<R>& p) const 
    {
      throw std::invalid_argument("Not implemented."); 
    }
    
    template<typename R>
    LinearAlgebra::Matrix<R> 
    Map<R>::derivative(const Geometry::Point<R>& x) const 
    {
      throw std::invalid_argument("Derivative at point not implemented."); 
    }

    template<typename R>
    LinearAlgebra::IntervalMatrix<R> 
    Map<R>::derivative(const Geometry::Rectangle<R>& r) const 
    {
      throw std::invalid_argument("Derivative on Rectangle not implemented."); 
    }
    
    
    template<typename R>
    template<template<typename> class BS>
    Geometry::ListSet<R,BS> 
    Map<R>::apply(const Geometry::ListSet<R,BS>& A) const 
    { 
      Geometry::ListSet<R,BS> trans_ds(A.dimension());
      for (size_t i=0; i< A.size(); i++) {
        trans_ds.inplace_union(this->apply(A[i]));
      }
      return trans_ds;
    }  
    
  }
}
