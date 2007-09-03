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

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../system/vector_field.h"

namespace Ariadne {
  namespace System {

    /*!\ingroup ContinuousTime
     * \brief An affine vector field in Euclidean space, given by \f$f(x)=Ax+b\f$.
     */
    template<class R>
    class AffineVectorField : public VectorFieldInterface<R> 
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef typename Numeric::traits<R>::interval_type I;
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type of denotable state the system acts on. */
      typedef Geometry::Point<R> state_type;
      
      /*! \brief The type of vector used to represent the affine transformation. */
      typedef LinearAlgebra::Matrix<R> matrix_type;
      /*! \brief The type of matrix used to represent the affine transformation. */
      typedef LinearAlgebra::Vector<R> mector_type;
    
      /*! \brief Virtual destructor. */
      virtual ~AffineVectorField();
      
      /*! \brief Copy constructor. */
      AffineVectorField(const AffineVectorField<R>& avf) : _a(avf.A()), _b(avf.b()) { }
      /*! \brief Construct from the matrix \a A and the vector \a b.. */
      AffineVectorField(const LinearAlgebra::Matrix<R> &A, const LinearAlgebra::Vector<R> &b) : _a(A), _b(b) { }
      /*! \brief Make a copy (clone) of the vector field. */
      AffineVectorField<R>* clone() const { return new AffineVectorField<R>(this->A(),this->b()); }
      
      /*! \brief An approximation to the vector field at a point. */
      virtual LinearAlgebra::Vector<F> image(const Geometry::Point<F>& x) const;

      /*! \brief An approximation to the Jacobian derivative at a point. */
      virtual LinearAlgebra::Matrix<F> jacobian(const Geometry::Point<F>& x) const;
      
      /*! \brief The matrix \f$A\f$. */
      const LinearAlgebra::Matrix<R>& A() const { return this->_a; }
      /*! \brief The vector \f$b\f$. */
      const LinearAlgebra::Vector<R>& b() const { return this->_b; }
      
      /*! \brief The dimension of the result. */
      smoothness_type smoothness() const { return std::numeric_limits<smoothness_type>::max(); }
 
      /*! \brief The dimension of the vector field is given by the size of \f$b\f$. */
      dimension_type dimension() const { return this->_b.size(); }
      
      /*! \brief  The name of the system. */
      std::string name() const { return "AffineVectorField"; }
      
      /*! \brief  The name of the system. */
      virtual std::ostream& write(std::ostream& os) const;
     private:
      LinearAlgebra::Matrix<R> _a;
      LinearAlgebra::Vector<R> _b;
    };
 
    
  }
}

namespace Ariadne {
  namespace LinearAlgebra {
    /*! \brief Compute \f$e^{Ah}\f$. */
    template<class R>
    Matrix<typename Numeric::traits<R>::arithmetic_type> 
    exp_Ah_approx(const Matrix<R>& A, 
                  const R& h); 

    /*! \brief Compute \f$A^{-1}(e^{Ah}-I) = h\sum_{n=0}^{\infty} \frac{{(Ah)}^{n}}{(n+1)!}\f$. */
    template<class R>
    Matrix<typename Numeric::traits<R>::arithmetic_type> 
    exp_Ah_sub_id_div_A_approx(const Matrix<R>& A, 
                               const R& h); 

  }
}

#endif /* ARIADNE_AFFINE_VECTOR_FIELD_H */
