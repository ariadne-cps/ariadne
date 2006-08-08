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

#include "../system/affine_map.h"

namespace Ariadne {
  namespace System {

    /*! \brief An affine multimap on Euclidean space, with the form \f$x\mapsto Ax+S \f$ for some basic set \f$S\f$. 
     */
    template <typename R,template <typename> class BS>
    class AffineMultiMap
    {
     public:
      typedef R real_type;
      typedef Geometry::Point<R> state_type;
      typedef LinearAlgebra::Matrix<R> matrix_type;
      typedef LinearAlgebra::Vector<R> vector_type;
      typedef BS<R> set_type;

      explicit AffineMultiMap(const matrix_type& A, const set_type& S)
        : _A(A), _S(S) { }
        
      AffineMultiMap<R,BS>& operator=(const AffineMultiMap<R,BS>& T) {
        this->_A=T._A; this->_S=T._S; return *this; 
      }
      
      /*! \brief  The map applied to a state. */
      BS<R> operator() (const Geometry::Point<R>& x) const;
        
      /*! \brief  The map applied to a BS<R>. */
      BS<R> operator() (const BS<R>& z) const;
           
      /*! \brief  The map applied to a grid mask set. */
      Geometry::ListSet<R,BS> operator() (const Geometry::GridMaskSet<R>& ) const;
      
      /*! \brief  The matrix of the map. */
      const LinearAlgebra::Matrix<R>& A() const { return _A; }
      
      /*! \brief  The offset set of the map. */
      const BS<R>& S() const { return _S; }
      
      /*! \brief  The dimension of the argument. */
      dimension_type argument_dimension() const {
        return _A.size2();
      }
      
      /*! \brief The dimension of the result. */
      dimension_type result_dimension() const {
        return _A.size1();
      }
      

      std::string name() const { return "AffineMultiMap"; }
     private:
      matrix_type _A;
      set_type _S;
    };
      
    
  }
}


#endif /* _ARIADNE_AFFINE_MULTIMAP_H */
