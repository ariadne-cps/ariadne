/***************************************************************************
 *            vector_field.code.h
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
#include "vector_field.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"



namespace Ariadne {
  namespace System {

    template<class R>
    VectorField<R>::~VectorField() 
    {
    }
     
    template<class R>
    typename VectorField<R>::F 
    VectorField<R>::derivative(const Geometry::Point<F>& x, const size_type& i, const multi_index_type& j) const
    {
      throw DeferredImplementation(__PRETTY_FUNCTION__);
    }
   

    template<class R>
    LinearAlgebra::Matrix<typename VectorField<R>::F> 
    VectorField<R>::jacobian(const Geometry::Point<F>& x) const 
    {
      throw DeferredImplementation(__PRETTY_FUNCTION__);
    }
    
    template<class R>
    std::ostream&
    VectorField<R>::write(std::ostream& os) const 
    {
      return os << "VectorField()";
    }
    
    
  }
}
