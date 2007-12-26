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
#include "numeric/arithmetic.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"

#include "linear_algebra/matrix.h"
#include "linear_algebra/matrix_function.h"

#include "geometry/box.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"

#include "system/vector_field.h"
#include "system/affine_vector_field.h"

#include "output/logging.h"


namespace Ariadne { 

namespace Evaluation { static int& verbosity = integrator_verbosity; }



template<class R>
Evaluation::TaylorIntegrator<R>*
Evaluation::TaylorIntegrator<R>::clone() const
{
  return new TaylorIntegrator<R>();
}


template<class R>
Evaluation::TaylorIntegrator<R>::TaylorIntegrator()
{
}



template<class R>
System::TaylorFlow<typename Evaluation::TaylorIntegrator<R>::I>
Evaluation::TaylorIntegrator<R>::flow(const System::VectorFieldInterface<R>& vector_field, 
                                      const Geometry::Box<R>& bounding_box) const
{
}



template<class R>
Geometry::Point<typename Evaluation::TaylorIntegrator<R>::I>
Evaluation::TaylorIntegrator<R>::flow_step(const System::VectorFieldInterface<R>& vector_field, 
                                           const Geometry::Point<I>& initial_point, 
                                           const Numeric::Interval<R>& step_size, 
                                           const Geometry::Box<R>& bounding_box) const
{
}




template<class R>
Geometry::Zonotope<typename Evaluation::TaylorIntegrator<R>::I>
Evaluation::TaylorIntegrator<R>::integration_step(const System::VectorFieldInterface<R>& vector_field, 
                                                  const Geometry::Zonotope<I,I>& initial_set, 
                                                  const Numeric::Interval<R>& step_size, 
                                                  const Geometry::Box<R>& bounding_box) const
{
}


template<class R>
Geometry::Zonotope<typename Evaluation::TaylorIntegrator<R>::I>
Evaluation::TaylorIntegrator<R>::reachability_step(const System::VectorFieldInterface<R>& vector_field, 
                                                   const Geometry::Zonotope<I,I>& initial_set, 
                                                   const Numeric::Interval<R>& step_size, 
                                                   const Geometry::Box<R>& bounding_box) const
{
}



template<class R>
std::ostream&
Evaluation::TaylorIntegrator<R>::write(std::ostream& os) const
{
  return os << "TaylorIntegrator( )";
}



}
