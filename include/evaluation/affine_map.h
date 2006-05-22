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
#include "../linear_algebra/interval_matrix.h"

#include "../geometry/point.h"
#include "../geometry/simplex.h"
#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/polyhedron.h"
#include "../geometry/list_set.h"
#include "../geometry/grid_set.h"

#include "../evaluation/map.h"


namespace Ariadne {
  namespace Evaluation {

    /*! \brief An affine map on Euclidean space. */
    template <typename R>
    class AffineMap : public Map<R> 
    {
     public:
      typedef R real_type;
      typedef Geometry::Point<R> state_type;
      
      typedef Ariadne::LinearAlgebra::Matrix<R> Matrix_type;
      typedef Ariadne::LinearAlgebra::Vector<R> Vector_type;
      
      explicit AffineMap() {}
      explicit AffineMap(const Matrix_type& A, const Vector_type& b) : _A(A), _b(b) { }
      explicit AffineMap(const Matrix_type& A) : _A(A), _b(A.size2()) { }
      explicit AffineMap(const Vector_type& b) : _A(b.size(),b.size()), _b(b) { }
      
      AffineMap(const AffineMap<real_type>& T) : _A(T._A), _b(T._b) { }
      AffineMap<real_type>& operator=(const AffineMap<real_type>& T) {
        this->_A=T._A; this->_b=T._b; return *this; }
      
      /*! \brief  The map applied to a state. */
      state_type operator() (const state_type& x) const;
        
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
      
      /*! \brief  The map applied to a list of basic sets. */
      template <template<class> class BS>
      Geometry::ListSet<R,BS> operator() (const Geometry::ListSet<R,BS>& A) const {
        Geometry::ListSet<R,BS> output(A.dimension());

        for (size_type i=0; i<A.size(); i++) { 
          output.push_back((*this)(A[i]));
        }

        return output;
      }

      /*! \brief  The map applied to a gridmaskset. */
      Geometry::ListSet<R,Geometry::Parallelotope> operator() (const Geometry::GridMaskSet<R>& ) const;
      
      /*! \brief  The linear transformation of the map. */
      const LinearAlgebra::Matrix<R>& A() const { return _A; }
      /*! \brief  The offset vector of the map. */
      const LinearAlgebra::Vector<R>& b() const { return _b; }
      
      /*! \brief  The dimension of the argument. */
      dimension_type argument_dimension() const {
        return _A.size2();
      }
      
      /*! \brief The dimension of the result. */
      dimension_type result_dimension() const {
        return _b.size();
      }
      
      bool invertible() const { assert(false); return false; }
  
      std::string name() const { return "AffineMap"; }
     private:
      Matrix_type _A;
      Vector_type _b;
    };
      
    
  }
}


#endif /* _ARIADNE_AFFINE_MAP_H */
