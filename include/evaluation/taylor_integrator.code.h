/***************************************************************************
 *            taylor_integrator.code.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 
//#define DEBUG

#include <iosfwd>
#include <string>
#include <sstream>
#include <algorithm>
#include <typeinfo>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "taylor_integrator.h"

#include "base/array.h"
#include "base/exceptions.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"

#include "linear_algebra/matrix.h"
#include "linear_algebra/matrix_function.h"

#include "geometry/box.h"
#include "geometry/zonotope.h"

#include "system/vector_field.h"
#include "system/affine_vector_field.h"

#include "output/logging.h"


namespace Ariadne { 

 static int& verbosity = integrator_verbosity; }



template<class R>
TaylorIntegrator<R>*
TaylorIntegrator<R>::clone() const
{
  return new TaylorIntegrator<R>();
}


template<class R>
TaylorIntegrator<R>::TaylorIntegrator()
{
}



template<class R>
TaylorFlow<typename TaylorIntegrator<R>::I>
TaylorIntegrator<R>::flow(const VectorField<R>& vector_field, 
                                      const Box<R>& bounding_box) const
{
}



template<class R>
Point<typename TaylorIntegrator<R>::I>
TaylorIntegrator<R>::flow_step(const VectorField<R>& vector_field, 
                                           const Point<I>& initial_point, 
                                           const Interval<R>& step_size, 
                                           const Box<R>& bounding_box) const
{
}




template<class R>
Zonotope<typename TaylorIntegrator<R>::I>
TaylorIntegrator<R>::integration_step(const VectorField<R>& vector_field, 
                                                  const Zonotope<I,I>& initial_set, 
                                                  const Interval<R>& step_size, 
                                                  const Box<R>& bounding_box) const
{
}


template<class R>
Zonotope<typename TaylorIntegrator<R>::I>
TaylorIntegrator<R>::reachability_step(const VectorField<R>& vector_field, 
                                                   const Zonotope<I,I>& initial_set, 
                                                   const Interval<R>& step_size, 
                                                   const Box<R>& bounding_box) const
{
}



template<class R>
std::ostream&
TaylorIntegrator<R>::write(std::ostream& os) const
{
  return os << "TaylorIntegrator( )";
}



}
