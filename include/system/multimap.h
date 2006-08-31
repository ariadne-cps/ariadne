/***************************************************************************
 *            multimap.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file multimap.h
 *  \brief Mulitvalued map interface.
 */

#ifndef _ARIADNE_MULTIMAP_H
#define _ARIADNE_MULTIMAP_H

#include <string>

#include "../declarations.h"

namespace Ariadne {
  namespace System {

    /*! \brief Abstract base class for multivalued functions.
     *  \ingroup System
     *  \ingroup DiscreteTime
     */
     template<typename R>
    class MultiMap {
     public:
      typedef R real_type;
      typedef Geometry::Point<R> state_type;
      typedef LinearAlgebra::Vector<R> vector_type;
      typedef LinearAlgebra::Matrix<R> matrix_type;
      typedef LinearAlgebra::IntervalMatrix<R> interval_matrix_type;
      
      virtual ~MultiMap();
      
      /*! \brief A minimal cover of the image of a rectangle. */
      virtual Geometry::ListSet<R,Geometry::Rectangle> operator() (const Geometry::Rectangle<R>& A) const;
      /*! \brief A minimal cover of the image of a parallelotope. */
      virtual Geometry::ListSet<R,Geometry::Parallelotope> operator() (const Geometry::Parallelotope<R>& A) const;
      /*! \brief A minimal cover of the image of a zonotope. */
      virtual Geometry::ListSet<R,Geometry::Zonotope> operator() (const Geometry::Zonotope<R>& A) const;
      /*! \brief A minimal cover of the image of a simplex. */
      virtual Geometry::ListSet<R,Geometry::Simplex> operator() (const Geometry::Simplex<R>& A) const;
      /*! \brief A minimal cover of the image of a polyhedron. */
      virtual Geometry::ListSet<R,Geometry::Polyhedron> operator() (const Geometry::Polyhedron<R>& A) const;
    
      /*! \brief The dimension of the domain space. */
      virtual dimension_type argument_dimension() const = 0;
      /*! \brief The dimension of the range space. */
      virtual dimension_type result_dimension() const = 0;
    
      /*! \brief The name of the map. */
      virtual std::string name() const = 0;
    };
  
    
  }
}

#endif /* _ARIADNE_MULTIMAP_H */
