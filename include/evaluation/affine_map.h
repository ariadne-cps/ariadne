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

#ifndef _ARIADNE_AFFINE_MAP_H
#define _ARIADNE_AFFINE_MAP_H

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../evaluation/map.h"

#include "../geometry/point.h"
#include "../geometry/simplex.h"
#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/polyhedron.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief An affine map on Euclidean space. */
    template <typename R>
    class AffineMap // : public Map<R,Geometry::Point> 
    {
     public:
      typedef R Real;
      typedef Geometry::Point<Real> State;
      
      typedef Ariadne::LinearAlgebra::matrix<Real> Matrix;
      typedef Ariadne::LinearAlgebra::vector<Real> Vector;
      
      inline explicit AffineMap() {}
      inline explicit AffineMap(const Matrix& A, const Vector& b) : _A(A), _b(b) { }
      inline explicit AffineMap(const Matrix& A) : _A(A), _b(A.size2()) { }
      inline explicit AffineMap(const Vector& b) : _A(b.size(),b.size()), _b(b) { }
      
      inline AffineMap(const AffineMap<Real>& T) : _A(T._A), _b(T._b) { }
      inline AffineMap<Real>& operator=(const AffineMap<Real>& T) {
        this->_A=T._A; this->_b=T._b; return *this; }
      
      /*! \brief  The map applied to a state. */
      State operator() (const State& x) const;
        
      /*! \brief  The map applied to a simplex basic set. */
      Geometry::Simplex<R> operator() (const Geometry::Simplex<R>& A) const;
      
      /*! \brief  The map applied to a rectangle basic set. */
      Geometry::Rectangle<R> operator() (const Geometry::Rectangle<R>& A) const;
      
      /*! \brief  The map applied to a parallelotope basic set. */
      Geometry::Parallelotope<R> operator() (const Geometry::Parallelotope<R>& A) const;
      
      /*! \brief  The map applied to a zonotope basic set. */
      Geometry::Zonotope<R> operator() (const Geometry::Zonotope<R>& A) const;
      
      /*! \brief  The map applied to a polyhedron basic set. */
      Geometry::Polyhedron<R> operator() (const Geometry::Polyhedron<R>& A) const;
      
      /*! \brief  The linear transformation of the map. */
      inline const Matrix& A() const { return _A; }
      /*! \brief  The offset vector of the map. */
      inline const Vector& b() const { return _b; }
      
      /*! \brief  The dimension of the argument. */
      inline size_type argument_dimension() const {
        return _A.size2();
      }
      
      /*! \brief The dimension of the result. */
      inline size_type result_dimension() const {
        return _b.size();
      }
      
      inline bool invertible() const { assert(false); return false; }
     private:
      Matrix _A;
      Vector _b;
    };
      
    
  }
}


#endif /* _ARIADNE_AFFINE_MAP_H */
