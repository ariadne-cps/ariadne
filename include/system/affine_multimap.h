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

#ifndef ARIADNE_AFFINE_MULTIMAP_H
#define ARIADNE_AFFINE_MULTIMAP_H

#include "multimap.h"

#include "linear_algebra/declarations.h"
#include "geometry/declarations.h"

namespace Ariadne {
  namespace System {

    /*! \brief An affine multimap on Euclidean space, with the form \f$x\mapsto Ax+S \f$ for some basic set \f$S\f$. 
     *  \ingroup DiscreteTime
    */
    template<class R,template<class> class BS>
    class AffineMultiMap
    {
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type of denotable state the system acts on. */
      typedef Geometry::Point<R> state_type;
      /*! \brief The type of set describing the image of a point. */
      typedef BS<R> set_type;
      
      /*! \brief The type of vector used to represent the affine transformation. */
      typedef LinearAlgebra::Vector<R> vector_type;
      /*! \brief The type of matrix used to represent the affine transformation. */
      typedef LinearAlgebra::Matrix<R> matrix_type;

      /*! \brief Construct from a matrix describing a linear transformation
       *  and a set decribing the offset. */
      explicit AffineMultiMap(const LinearAlgebra::Matrix<R>& A, 
                              const BS<R>& S)
        : _mx(A), _set(S) { }
        
      /*! \brief Copy constructor. */
      AffineMultiMap(const AffineMultiMap<R,BS>& T) : _mx(T._mx), _set(T._set) { }
      
      /*! \brief Assignment operator. */
      AffineMultiMap<R,BS>& operator=(const AffineMultiMap<R,BS>& T) {
        this->_mx=T._mx; this->_set=T._set; return *this; 
      }
      
      /*! \brief  The map applied to a state. */
      BS<R> operator() (const Geometry::Point<R>& x) const;
        
      /*! \brief  The map applied to a BS<R>. */
      BS<R> operator() (const BS<R>& z) const;
           
      /*! \brief  The matrix of the map. */
      const LinearAlgebra::Matrix<R>& A() const { return _mx; }
      
      /*! \brief  The offset set of the map. */
      const BS<R>& S() const { return _set; }
      
      /*! \brief  The dimension of the argument. */
      dimension_type argument_dimension() const {
        return _mx.number_of_columns();
      }
      
      /*! \brief The dimension of the result. */
      dimension_type result_dimension() const {
        return _mx.number_of_rows();
      }
      

      /*! \brief  The name of the system. */
      std::string name() const { return "AffineMultiMap"; }
     private:
      matrix_type _mx;
      set_type _set;
    };
      
    
  }
}


#endif /* ARIADNE_AFFINE_MULTIMAP_H */
