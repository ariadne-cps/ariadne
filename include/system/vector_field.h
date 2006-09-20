/***************************************************************************
 *            vector_field.h
 *
 *  Thu Feb  3 21:06:54 2005
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
 
/*! \file vector_field.h
 *  \brief Vector_type field interface.
 */
 
#ifndef _ARIADNE_VECTOR_FIELD_H
#define _ARIADNE_VECTOR_FIELD_H

#include "../declarations.h"

namespace Ariadne {
  namespace System {

    /*! \brief Abstract base class for (differentiable) vector fields.
     *  \ingroup System
     *  \ingroup ContinuousTime
     */
    template <typename R>
    class VectorField {
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type of denotable state the system acts on. */
      typedef Geometry::Point<R> state_type;
      
      /*! \brief Virtual destructor. */
      virtual ~VectorField();
     
      /*! \brief An approximation to the vector field at a point. */
      virtual LinearAlgebra::Vector<R> operator() (const Geometry::Point<R>& x) const;
      /*! \brief A bound for the vector field over a rectangle. */
      virtual LinearAlgebra::Vector< Interval<R> > operator() (const Geometry::Rectangle<R>& A) const;

      /*! \brief An approximation to the Jacobian derivative at a point. */
      virtual LinearAlgebra::Matrix<R> derivative(const Geometry::Point<R>& x) const;
      /*! \brief A bound for the Jacobian derivative over a rectangle. */
      virtual LinearAlgebra::Matrix< Interval<R> > derivative(const Geometry::Rectangle<R>& A) const;
    
      /*! \brief The dimension of the space the vector field lives in. */
      virtual dimension_type dimension() const = 0;

      /*! \brief The name of the system. */
      virtual std::string name() const = 0;
    };
   
  }
}

#endif /* _ARIADNE_VECTOR_FIELD_H */
