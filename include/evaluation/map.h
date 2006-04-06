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
 *  \brief Map interface.
 */

#ifndef _ARIADNE_MAP_H
#define _ARIADNE_MAP_H

#include <string>

#include "../declarations.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_matrix.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief Base class for (differentiable) functions. */
    template<typename R>
    class Map {
     public:
      typedef R real_type;
      typedef Geometry::Point<R> state_type;
      typedef LinearAlgebra::vector<R> vector_type;
      typedef LinearAlgebra::matrix<R> matrix_type;
      typedef LinearAlgebra::interval_matrix<R> interval_matrix_type;
      
      virtual ~Map();
      
      // I'm not sure if virtual functions are the way to go here; too restrictive? //
      virtual state_type apply(const state_type& x) const;
      virtual Geometry::Rectangle<R> apply(const Geometry::Rectangle<R>& A) const;
      virtual Geometry::Parallelotope<R> apply(const Geometry::Parallelotope<R>& A) const;
      virtual Geometry::Polyhedron<R> apply(const Geometry::Polyhedron<R>& A) const;
    
      virtual LinearAlgebra::matrix<R> derivative(const state_type& r) const;
      virtual LinearAlgebra::interval_matrix<R> derivative(const Geometry::Rectangle<R>& r) const;
        
      template<template<typename> class BS>
      inline Geometry::ListSet<R,BS> apply(const Geometry::ListSet<R,BS>& A) const;
      
      virtual dimension_type argument_dimension() const = 0;
      virtual dimension_type result_dimension() const = 0;
    
      virtual std::string name() const = 0;

      state_type operator() (const state_type& x) const {
        return this->apply(x); }
      Geometry::Rectangle<R> operator() (const Geometry::Rectangle<R>& A) const {
        return this->apply(A); }
      Geometry::Polyhedron<R> operator() (const Geometry::Polyhedron<R>& A) const {
        return this->apply(A); }
    
    };
  
    
  }
}

#endif /* _ARIADNE_MAP_H */
