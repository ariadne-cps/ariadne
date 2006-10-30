/***************************************************************************
 *            map.tpl
 *
 *  Copyright  2005, 2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it,  Pieter.Collins@cwi.nl
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
 
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"

#include "../system/map.h"

namespace Ariadne {
  namespace System {

    template<class R>
    Map<R>::~Map() 
    {
    }
    
    template<class R>
    typename Map<R>::result_type
    Map<R>::image(const Geometry::Point<R>& r) const 
    {
      throw DeferredImplementation(__PRETTY_FUNCTION__);
    }
    
    template<class R>
    Geometry::Rectangle<R>
    Map<R>::image(const Geometry::Rectangle<R>& r) const 
    {
      throw DeferredImplementation(__PRETTY_FUNCTION__);
    }
    
    template<class R>
    typename Map<R>::F
    Map<R>::derivative(const Geometry::Point<R>& r, const size_type& i, const multi_index_type& j) const 
    {
      throw DeferredImplementation(__PRETTY_FUNCTION__);
    }
    
    template<class R>
    typename Map<R>::I
    Map<R>::derivative(const Geometry::Rectangle<R>& r, const size_type& i, const multi_index_type& j) const 
    {
      throw DeferredImplementation(__PRETTY_FUNCTION__);
    }
    
    template<class R>
    LinearAlgebra::Matrix<typename Map<R>::F>
    Map<R>::jacobian(const Geometry::Point<R>& r) const 
    {
      throw DeferredImplementation(__PRETTY_FUNCTION__);
    }
    
    template<class R>
    LinearAlgebra::Matrix<typename Map<R>::I> 
    Map<R>::jacobian(const Geometry::Rectangle<R>& r) const 
    {
      throw DeferredImplementation(__PRETTY_FUNCTION__);
    }
    
    
   
  }
}
