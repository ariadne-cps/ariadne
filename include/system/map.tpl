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

    template<typename R>
    Map<R>::~Map() 
    {
    }
  
    template<typename R>
    size_type
    Map<R>::smoothness() const
    {
      return (size_type)(-1);
    }
  
    template<typename R>
    typename Map<R>::result_type
    Map<R>::operator() (const Geometry::Point<R>& r) const 
    {
      throw std::invalid_argument(this->name()+"::operator() (Rectangle) not implemented."); 
    }
    
    template<typename R>
    Geometry::Rectangle<R>
    Map<R>::operator() (const Geometry::Rectangle<R>& r) const 
    {
      throw std::invalid_argument(this->name()+"::operator() (Rectangle) not implemented."); 
    }
    
    template<typename R>
    LinearAlgebra::Matrix<typename Map<R>::F>
    Map<R>::jacobian(const Geometry::Point<R>& r) const 
    {
      throw std::invalid_argument(this->name()+"::jacobian(Point) not implemented."); 
    }
    
    template<typename R>
    LinearAlgebra::Matrix< Interval<R> > 
    Map<R>::jacobian(const Geometry::Rectangle<R>& r) const 
    {
      throw std::invalid_argument(this->name()+"::jacobian(Rectangle) not implemented."); 
    }
    
    
   
  }
}
