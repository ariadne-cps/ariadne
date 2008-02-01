/***************************************************************************
 *            taylor_set.code.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
#include <iostream>
#include <vector>
#include <algorithm>

#include "taylor_set.h"

#include "exceptions.h"
#include "base/array.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "geometry/point.h"
#include "geometry/box.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"

#include "output/logging.h"

namespace Ariadne {
    
extern int Geometry::verbosity; 
   
template<class R>
void
Geometry::TaylorSet<R>::_instantiate_geometry_operators() 
{
  Zonotope<R>* z=0;
  TaylorSet<R>* ts=0;
  *z=over_approximation(*ts);
}
    
template<class R> 
Geometry::Zonotope<R> 
Geometry::over_approximation(const TaylorSet<R>&) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



} // namespace Ariadne
