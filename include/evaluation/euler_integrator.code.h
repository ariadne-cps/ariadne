/***************************************************************************
 *            euler_integrator.code.h
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

#include "euler_integrator.h"

#include "base/array.h"
#include "numeric/arithmetic.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "geometry/rectangle.h"

#include "system/vector_field.h"

#include "output/logging.h"

namespace Ariadne {
    

namespace Evaluation { static int& verbosity = integrator_verbosity; }



template<class R>
Evaluation::EulerIntegrator<R>::EulerIntegrator()
{
}



template<class R>
Evaluation::EulerIntegrator<R>*
Evaluation::EulerIntegrator<R>::clone() const
{
  return new EulerIntegrator<R>();
}



template<class R>
Geometry::Rectangle<R> 
Evaluation::EulerIntegrator<R>::integration_step(const System::VectorFieldInterface<R>& vector_field, 
                                                         const Geometry::Rectangle<R>& initial_set, 
                                                         const Numeric::Interval<R>& step_size, 
                                                         const Geometry::Box<R>& bounding_set) const
{
  ARIADNE_LOG(6,"EulerIntegrator::integration_step(VectorFieldInterface,Rectangle,Interval,Box) const\n");
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"EulerIntegrator::integration_step(VectorFieldInterface,Rectangle,Interval,Box) const");

  return initial_set + I(step_size) * vector_field(bounding_set);
}



template<class R>
Geometry::Rectangle<R> 
Evaluation::EulerIntegrator<R>::reachability_step(const System::VectorFieldInterface<R>& vector_field, 
                                                          const Geometry::Rectangle<R>& initial_set, 
                                                          const Numeric::Interval<R>& step_size, 
                                                          const Geometry::Box<R>& bounding_set) const
{
  ARIADNE_LOG(6,"EulerIntegrator::reachability_step(VectorFieldInterface,Rectangle,Numeric::Interval<R>) const\n");
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set(),"EulerIntegrator::reachability_step(VectorFieldInterface,Rectangle,Numeric::Interval<R>) const");
  
  return initial_set + I(0,step_size.upper()) * vector_field(bounding_set);
}


template<class R>
std::ostream&
Evaluation::EulerIntegrator<R>::write(std::ostream& os) const
{
  return os << "EulerIntegrator( )";
}


}

