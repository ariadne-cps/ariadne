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

#ifndef _AFFINE_VECTOR_FIELD_H
#define _AFFINE_VECTOR_FIELD_H

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_vector.h"
#include "../linear_algebra/interval_matrix.h"

#include "../system/vector_field.h"
#include "../system/affine_map.h"

namespace Ariadne {
  namespace System {

    /*! \brief An affine vector field in Euclidean space.
     *  \ingroup ContinuousTime
     */
    template <typename R>
    class AffineVectorField : public VectorField<R> 
    {
      typedef typename Geometry::Polyhedron<R> Polyhedron;
      typedef typename Geometry::Rectangle<R> Rectangle;
    
     public:
      typedef typename Geometry::Point<R> state_type;
      
      typedef LinearAlgebra::Matrix<R> matrix_type;
      typedef LinearAlgebra::Vector<R> mector_type;
    
      typedef LinearAlgebra::IntervalVector<R> IntervalVector_type;
      typedef LinearAlgebra::IntervalMatrix<R> IntervalMatrix_type;
    
      virtual ~AffineVectorField();
      
      AffineVectorField(const AffineVectorField<R>& F) : _A(F.A()), _b(F.b()) { }
      AffineVectorField(const LinearAlgebra::Matrix<R> &A, const LinearAlgebra::Vector<R> &b) : _A(A), _b(b) { }
    
      LinearAlgebra::Vector<R> operator() (const Geometry::Point<R>& s) const;
      LinearAlgebra::IntervalVector<R> operator() (const Geometry::Rectangle<R>& r) const;
    
      LinearAlgebra::Matrix<R> derivative(const Geometry::Point<R>& x) const;
      LinearAlgebra::IntervalMatrix<R> derivative(const Geometry::Rectangle<R>& r) const;
      
      const LinearAlgebra::Matrix<R>& A() const { return this->_A; }
      const LinearAlgebra::Vector<R>& b() const { return this->_b; }
      
      dimension_type dimension() const {
        return this->_b.size();
      }
      
      std::string name() const { return "AffineVectorField"; }
      
     private:
      LinearAlgebra::Matrix<R> _A;
      LinearAlgebra::Vector<R> _b;
    };
 
    template<typename R> std::ostream& operator<<(std::ostream& os, const AffineVectorField<R>& vf);
    
  }
}

namespace Ariadne {
  namespace LinearAlgebra {
    template <typename R>
    Matrix<R> 
    exp_Ah_approx(const Matrix<R>& A, 
                  const R& h, 
                  const R& e); 

    /*! \brief Compute \f$A^{-1}(e^{Ah}-I) = h\sum_{n=0}^{\infty} \frac{{(Ah)}^{n}}{(n+1)!}\f$. */
    template <typename R> 
    Matrix<R> 
    exp_Ah_sub_id_div_A_approx(const Matrix<R>& A, 
                               const R& h, 
                               const R& e); 

  }
}

#endif /* _AFFINE_VECTOR_FIELD_H */
