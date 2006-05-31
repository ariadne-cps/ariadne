/***************************************************************************
 *            affine_multimap.h
 *
 *  Thr May 30 18:52:36 2005
 *  Copyright  2006  Alberto Casagrande
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
 
/*! \file affine_multimap.h
 *  \brief Maps of affine form \f$x\rightarrow Ax \oplus S \f$.
 */

#ifndef _ARIADNE_AFFINE_MULTIMAP_H
#define _ARIADNE_AFFINE_MULTIMAP_H

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

#include "../evaluation/affine_map.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief An affine multimap on Euclidean space. */
    template <typename R, template<typename> class BS>
    class AffineMultiMap : public AffineMap<R> 
    {
     public:
      typedef R real_type;
      typedef Geometry::Point<R> state_type;
      
      typedef Ariadne::LinearAlgebra::Matrix<R> Matrix_type;
      typedef Ariadne::LinearAlgebra::Vector<R> Vector_type;
      typedef BS<R> Set_type;

      explicit AffineMultiMap() {}
      explicit AffineMultiMap(const Matrix_type& A, const Vector_type& b):
                               AffineMap<R>(A,b), _S(A.size1()) { }

      explicit AffineMultiMap(const Matrix_type& A, const Set_type& S):
                               AffineMap<R>(A), _S(S) { }

      explicit AffineMultiMap(const Matrix_type& A) : 
                               AffineMap<R>(A), _S(A.size1()) { }

      explicit AffineMultiMap(const Set_type& S) :  _S(S) { 

        this->_A=Matrix_type(S.dimension(),S.dimension());
        this->_b=Vector_type(S.dimension());
      }
      
      AffineMultiMap(const AffineMultiMap<R,BS>& T) : 
                               AffineMap<R>((AffineMap<R> &)T), _S(T._S) { }

      AffineMultiMap<R,BS>& 
      operator=(const AffineMultiMap<R,BS>& T) 
      {
        this->_A=T._A; 
	this->_b=T._b; 
	this->_S=T._S; 
	return *this; 
      }
      
      /*! \brief  The map applied to a state. */
      Set_type operator() (const state_type& x) const;
        
      /*! \brief  The map applied to a basic set. */
      template<template<typename> class BS2>
      inline
      BS2<R> operator() (const BS2<R>& A) const {
         AffineMap<real_type> *amap = (AffineMap<real_type> *)this;

         if (Geometry::is_a<BS,BS2>())
	   return Geometry::minkowski_sum((*amap)(A), (BS2<R>)(this->_S));
         
         throw std::invalid_argument("AffineMultiMap<R,BS>::operator() (const BS2<R>& A)");

	 return BS<R>(A.dimension());
      }
           
      /*! \brief  The map applied to a list of basic sets. */
      template <template<class> class BS2>
      Geometry::ListSet<R,BS2> operator() (const Geometry::ListSet<R,BS2>& A) const {
        Geometry::ListSet<R,BS2> output(A.dimension());

        for (size_type i=0; i<A.size(); i++) { 
          output.push_back((*this)(A[i]));
        }

        return output;
      }

      /*! \brief  The map applied to a gridmaskset. */
      Geometry::ListSet<R,Geometry::Zonotope> operator() (const Geometry::GridMaskSet<R>& ) const;
      
      /*! \brief  The offset set of the map. */
      const BS<R>& S() const { return _S; }
      
      bool invertible() const { return false; }

      std::string name() const { return "AffineMultiMap"; }
     private:
      Set_type _S;
    };
      
    
  }
}


#endif /* _ARIADNE_AFFINE_MULTIMAP_H */

