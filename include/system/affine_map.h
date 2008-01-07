/***************************************************************************
 *            affine_map.h
 *
 *  Copyright  2005-7  Alberto Casagrande Pieter Collins
 *
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
 *  \brief MapInterfaces of affine form \f$x\rightarrow Ax+b\f$.
 */

#ifndef ARIADNE_AFFINE_MAP_H
#define ARIADNE_AFFINE_MAP_H

#include "base/types.h"

#include "numeric/declarations.h"
#include "numeric/traits.h"

#include "linear_algebra/declarations.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "geometry/declarations.h"

#include "system/map.h"


namespace Ariadne {
  namespace System {

    /*! \brief An affine map \f$f(x)=Ax+b\f$ on Euclidean space. 
     *  \ingroup DiscreteTime
     */
    template<class R>
    class AffineMap
      : public Map<R> 
    {
     public:
      /*! \brief Construct from the matrix \f$A\f$ and the vector \f$b\f$. */
      explicit AffineMap(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& b)
        : Map<R>(Function::AffineFunction<R>(A,b)) { }
    };


  }
}


#endif /* ARIADNE_AFFINE_MAP_H */
