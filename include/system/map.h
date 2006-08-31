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

namespace Ariadne {
  namespace System {

    /*! \brief Abstract base class for (differentiable) functions.
     *  \ingroup System
     *  \ingroup DiscreteTime
     */
    template<typename R>
    class Map {
     public:
      typedef R real_type;
      typedef Geometry::Point<R> state_type;
      typedef LinearAlgebra::Vector<R> vector_type;
      typedef LinearAlgebra::Matrix<R> matrix_type;
      typedef LinearAlgebra::IntervalMatrix<R> interval_matrix_type;
      
      /*! \brief Virtual destructor. */
      virtual ~Map();
      
      /*! \brief The image of a point, if this can be computed exactly. */
      virtual Geometry::Point<R> operator() (const Geometry::Point<R>& x) const;
      /*! \brief An over-approximation to the image of a rectangle. */
      virtual Geometry::Rectangle<R> operator() (const Geometry::Rectangle<R>& A) const;
      /*! \brief An over-approximation to the image of a parallelotope. */
      virtual Geometry::Parallelotope<R> operator() (const Geometry::Parallelotope<R>& A) const;
      /*! \brief An over-approximation to the image of a zonotope. */
      virtual Geometry::Zonotope<R> operator() (const Geometry::Zonotope<R>& A) const;
      /*! \brief An over-approximation to the image of a simplex. */
      virtual Geometry::Simplex<R> operator() (const Geometry::Simplex<R>& A) const;
      /*! \brief An over-approximation to the image of a polyhedron. */
      virtual Geometry::Polyhedron<R> operator() (const Geometry::Polyhedron<R>& A) const;
    
      /*! \brief The derivative at a point, if this can be computed exactly. */
      virtual LinearAlgebra::Matrix<R> derivative(const Geometry::Point<R>& r) const;
      /*! \brief The derivatives over a rectangle. */
      virtual LinearAlgebra::IntervalMatrix<R> derivative(const Geometry::Rectangle<R>& r) const;
        
      /*! \brief The dimension of the domain space. */
      virtual dimension_type argument_dimension() const = 0;
      /*! \brief The dimension of the range space. */
      virtual dimension_type result_dimension() const = 0;
    
      /*! \brief The name of the map. */
      virtual std::string name() const = 0;

    };
  
    
  }
}

#endif /* _ARIADNE_MAP_H */
