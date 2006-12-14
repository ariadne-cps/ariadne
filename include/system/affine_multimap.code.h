/***************************************************************************
 *            affine_multimap.code.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
 
#include "affine_multimap.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/simplex.h"
#include "../geometry/polyhedron.h"

#include "../system/affine_map.h"

namespace Ariadne {
  namespace System {

    template<class R, template<class> class BS>
    BS<R>
    AffineMultiMap<R,BS>::operator() (const Geometry::Point<R>& pt) const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);

      using namespace Ariadne::LinearAlgebra;
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef typename Numeric::traits<R>::interval_type I;
     
      check_argument_dimension(*this,pt,__PRETTY_FUNCTION__);
      Vector<I> iv=this->A()*Vector<F>(pt.position_vector());
//      return over_approximation(minkowski_sum(this->S(),BS<R>(Geometry::Rectangle<R>(iv))));
    }
    
    template<class R, template<class> class BS>
    BS<R>
    AffineMultiMap<R,BS>::operator() (const BS<R>& bs) const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);

      using namespace Ariadne::LinearAlgebra;
      using namespace Ariadne::Geometry;

      check_argument_dimension(*this,bs,__PRETTY_FUNCTION__);
//      return over_approximation(minkowski_sum(AffineMap<R>(this->A()).image(bs),this->S()));
    }
    
  }
}
