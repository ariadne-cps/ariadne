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
 *  \brief vector_type field interface.
 */
 
#ifndef _ARIADNE_VECTOR_FIELD_H
#define _ARIADNE_VECTOR_FIELD_H

#include "../declarations.h"

#include "../numeric/interval.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_vector.h"
#include "../linear_algebra/interval_matrix.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief Abstract base class for (differentiable) vector fields. */
    template <typename R>
    class VectorField {
     public:
      typedef R real_type;
      typedef Geometry::Point<R> state_type;
      
      virtual ~VectorField();
     
      virtual LinearAlgebra::vector<R> apply(const Geometry::Point<R>& x) const;
      virtual LinearAlgebra::interval_vector<R> apply(const Geometry::Rectangle<R>& A) const;

      virtual LinearAlgebra::matrix<R> derivative(const Geometry::Point<R>& x) const;
      virtual LinearAlgebra::interval_matrix<R> derivative(const Geometry::Rectangle<R>& A) const;
    
      virtual dimension_type dimension() const = 0;

      virtual std::string name() const = 0;

      LinearAlgebra::vector<R> operator() (const Geometry::Point<R>& x) const  {
        return this->apply(x); }

      LinearAlgebra::interval_vector<R> operator() (const Geometry::Rectangle<R>& r) const {
        return this->apply(r); }
     };
    
  }
}

#endif /* _ARIADNE_VECTOR_FIELD_H */
