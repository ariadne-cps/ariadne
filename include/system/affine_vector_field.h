/***************************************************************************
 *            affine_vector_field.h
 *
 *  Fri Feb  4 08:57:39 2005
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
 
 /*! \file affine_vector_field.h
 *  \brief Vector_type fields of affine form of the form \f$\dot{x}=Ax+b\f$.
 */

#ifndef ARIADNE_AFFINE_VECTOR_FIELD_H
#define ARIADNE_AFFINE_VECTOR_FIELD_H

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "function/affine_function.h"
#include "system/vector_field.h"

namespace Ariadne {
  

    /*!\ingroup ContinuousTime
     * \brief An affine vector field in Euclidean space, given by \f$f(x)=Ax+b\f$.
     */
    template<class R>
    class AffineVectorField
      : public VectorField<R> 
    {
      typedef typename traits<R>::interval_type F;
     public:
      /*! \brief Construct from the matrix \a A and the vector \a b.. */
      AffineVectorField(const Matrix<R> &A, const Vector<R> &b)
        : VectorField<R>(AffineFunction<R>(A,b)) { }
      Matrix<F> A() const { return static_cast<AffineFunction<R>&>(this->function()).A(); }
      Vector<F> b() const { return static_cast<AffineFunction<R>&>(this->function()).b(); }
    };
 
    
  
} // namespace Ariadne

#endif /* ARIADNE_AFFINE_VECTOR_FIELD_H */
