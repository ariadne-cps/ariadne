/***************************************************************************
 *            discrete_time_system.code.h
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
 
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/list_set.h"

#include "system/discrete_time_system.h"

namespace Ariadne {
  

    template<class R>
    DiscreteTimeSystem<R>::~DiscreteTimeSystem() 
    {
    }
    
    
    
    
    template<class R>
    Matrix<F> 
    DiscreteTimeSystem<R>::jacobian(const Point<F>& x,
                                    const Point<F>& u,
                                    const Point<F>& v) const
    {
      throw DeferredImplementation(__PRETTY_FUNCTION__);
    }
