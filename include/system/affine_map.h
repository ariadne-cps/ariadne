/***************************************************************************
 *            affine_map.h
 *
 *  Wed Feb  2 18:52:36 2005
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
 
/*! \file affine_map.h
 *  \brief Maps of affine form \f$x\rightarrow Ax+b\f$.
 */

#ifndef ARIADNE_AFFINE_MAP_H
#define ARIADNE_AFFINE_MAP_H

#include "../declarations.h"

#include "../exceptions.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../system/map.h"


namespace Ariadne {
  namespace System {

    /*! \brief An affine map \f$f(x)=Ax+b\f$ on Euclidean space. 
     *  \ingroup DiscreteTime
     */
    template<class R>
    class AffineMap : public Map<R> 
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef typename Map<R>::I I;
     public:
      /*! \brief The type of denotable real number used to describe the system. */
      typedef R real_type;
      /*! \brief The type of denotable state the system acts on. */
      typedef Geometry::Point<R> state_type;
      
      /*! \brief Default constructor constructs a map on a zero-dimensional space. */
      explicit AffineMap() {}
      /*! \brief Construct from the matrix \f$A\f$ and the vector \f$b\f$. */
      explicit AffineMap(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& b)
        : _a(A), _b(b) { }
      /*! \brief Construct a linear map from the matrix \f$A\f$. */
      explicit AffineMap(const LinearAlgebra::Matrix<R>& A)
        : _a(A), _b(A.number_of_rows()) { }
      /*! \brief Construct a translation from the vector \f$b\f$. */
      explicit AffineMap(const LinearAlgebra::Vector<R>& b)
        : _a(LinearAlgebra::Matrix<R>::identity(b.size())), _b(b) { }
      
      /*! \brief Copy constructor. */
      AffineMap(const AffineMap<R>& T) : _a(T._a), _b(T._b) { }
      /*! \brief Assignment operator. */
      AffineMap<R>& operator=(const AffineMap<R>& T) {
        this->_a=T._a; this->_b=T._b; return *this; }
      /*! \brief Returns a pointer to a dynamically-allocated copy of the map. */
      virtual AffineMap<R>* clone() const { return new AffineMap<R>(*this); }

      
      /*! \brief  An approximation to the image of an approximate point. */
      Geometry::Point<F> image(const Geometry::Point<F>& A) const;
      
      /*! \brief  The map applied to a zonotope. */
      Geometry::Zonotope<F> image(const Geometry::Zonotope<F>& A) const;

      /*! \brief  The map applied to a polytope. */
      Geometry::Polytope<F> image(const Geometry::Polytope<F>& A) const;

      /*! \brief  The map applied to a zonotope basic set. */
      Geometry::Zonotope<F> operator() (const Geometry::Zonotope<F>& A) const {
        return this->image(A); }
              
      /*! \brief  The map applied to a polytopic basic set. */
      Geometry::Polytope<F> operator() (const Geometry::Polytope<R>& A) const{
        return this->image(A); };
              
      /*! \brief  The linear transformation of the map. */
      const LinearAlgebra::Matrix<R>& A() const { return _a; }
      /*! \brief  The offset vector of the map. */
      const LinearAlgebra::Vector<R>& b() const { return _b; }
      
      /*! \brief  The dimension of the argument. */
      virtual dimension_type argument_dimension() const {
        return _a.number_of_columns();
      }
      
      /*! \brief The dimension of the result. */
      virtual size_type smoothness() const { return (size_type) -1; }
      
      /*! \brief The dimension of the result. */
      virtual dimension_type result_dimension() const {
        return _b.size();
      }
      
      /*! \brief The Jacobian derivative matrix at a point. */
      virtual LinearAlgebra::Matrix<F> jacobian(const Geometry::Point<F>& pt) const;

      /*! \brief True if the map is invertible, which is equivalent to invertiblity of
       *  the matrix A. */
      bool invertible() const { throw NotImplemented("bool AffineMap<R>::invertible() const"); }
  
      /*! \brief  The name of the system. */
      std::string name() const { return "AffineMap"; }
      
      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;
     protected:
      LinearAlgebra::Matrix<R> _a;
      LinearAlgebra::Vector<R> _b;
    };


  }
}


#endif /* ARIADNE_AFFINE_MAP_H */
