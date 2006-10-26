/***************************************************************************
 *            multimap.tpl
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/simplex.h"
#include "../geometry/polytope.h"
#include "../geometry/list_set.h"

#include "../system/multimap.h"

namespace Ariadne {
  namespace System {

    template<class R>
    MultiMap<R>::~MultiMap() 
    {
    }
  
    template<class R>
    Geometry::ListSet<R,Geometry::Rectangle>
    MultiMap<R>::operator() (const Geometry::Rectangle<R>& r) const 
    {
      throw std::invalid_argument(this->name()+"::operator() (Rectangle) not implemented."); 
    }
    
    template<class R>
    Geometry::ListSet<R,Geometry::Parallelotope>
    MultiMap<R>::operator() (const Geometry::Parallelotope<R>& p) const 
    {
      throw std::invalid_argument(this->name()+"::operator() (Parallelotope) not implemented."); 
    }
    
    template<class R>
    Geometry::ListSet<R,Geometry::Zonotope>
    MultiMap<R>::operator() (const Geometry::Zonotope<R>& p) const 
    {
      throw std::invalid_argument(this->name()+"::operator() (Zonotope) not implemented."); 
    }
    
    template<class R>
    Geometry::ListSet<R,Geometry::Simplex>
    MultiMap<R>::operator() (const Geometry::Simplex<R>& p) const 
    {
      throw std::invalid_argument(this->name()+"::operator() (Simplex) not implemented."); 
    }
    
    template<class R>
    Geometry::ListSet<R,Geometry::Polytope>
    MultiMap<R>::operator() (const Geometry::Polytope<R>& p) const 
    {
      throw std::invalid_argument(this->name()+"::operator() (Polytope) not implemented."); 
    }
    
    
  }
}
