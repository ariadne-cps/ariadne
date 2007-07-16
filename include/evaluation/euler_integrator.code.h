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

#include "../base/array.h"
#include "../numeric/arithmetic.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/rectangle.h"

#include "../system/vector_field.h"

#include "../evaluation/integrator.h"

#include "../output/logging.h"

namespace Ariadne {
    

namespace Evaluation { static int& verbosity = integrator_verbosity; }



template<class R>
Evaluation::EulerIntegrator<R>::EulerIntegrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_basic_set_radius)
  : Base_(maximum_step_size,lock_to_grid_time,maximum_basic_set_radius)
{
}



template<class R>
Geometry::Point< Numeric::Interval<R> >
Evaluation::EulerIntegrator<R>::bounded_flow(const System::VectorFieldInterface<R>& vector_field, 
                                             const Geometry::Point<I>& initial_point, 
                                             const Geometry::Rectangle<R>& bounding_set, 
                                             const time_type& step_size) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  if(verbosity>6) { std::clog << "EulerIntegrator::integration_step(VectorFieldInterface,Rectangle,time_type) const" << std::endl; }
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_point,"EulerIntegrator::integration_step(VectorFieldInterface,Point,Rectangle,time_type) const");
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,bounding_set,"EulerIntegrator::integration_step(VectorFieldInterface,Point,Rectangle,time_type) const");
  
  return initial_point + ( I(step_size) * vector_field(Point<I>(bounding_set)) );
}


template<class R>
LinearAlgebra::Matrix< Numeric::Interval<R> > 
Evaluation::EulerIntegrator<R>::bounded_flow_jacobian(const System::VectorFieldInterface<R>& vector_field, 
                                                      const Geometry::Point<I>& initial_point, 
                                                      const Geometry::Rectangle<R>& bounding_set, 
                                                      const time_type& step_size) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}




template<class R>
Geometry::Rectangle<R> 
Evaluation::EulerIntegrator<R>::bounded_integration_step(const System::VectorFieldInterface<R>& vector_field, 
                                                         const Geometry::Rectangle<R>& initial_set, 
                                                         const Geometry::Rectangle<R>& bounding_set, 
                                                         const time_type& step_size) const
{
  if(verbosity>6) { std::clog << "EulerIntegrator::integration_step(VectorFieldInterface,Rectangle,Rectangle,time_type) const" << std::endl; }
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"EulerIntegrator::integration_step(VectorFieldInterface,Rectangle,Rectangle,time_type) const");

  return initial_set + I(step_size) * vector_field(bounding_set);
}



template<class R>
Geometry::Rectangle<R> 
Evaluation::EulerIntegrator<R>::bounded_reachability_step(const System::VectorFieldInterface<R>& vector_field, 
                                                          const Geometry::Rectangle<R>& initial_set, 
                                                          const Geometry::Rectangle<R>& bounding_set, 
                                                          const time_type& step_size) const
{
  if(verbosity>6) { std::clog << "EulerIntegrator::reachability_step(VectorFieldInterface,Rectangle,time_type) const" << std::endl; }
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set(),"EulerIntegrator::reachability_step(VectorFieldInterface,Rectangle,time_type) const");
  
  return initial_set + I(0,step_size) * vector_field(bounding_set);
}


template<class R>
Geometry::Rectangle<R> 
Evaluation::EulerIntegrator<R>::integration_step(const System::VectorFieldInterface<R>& vector_field, 
                                                 const Geometry::Rectangle<R>& initial_set, 
                                                 time_type& step_size) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  if(verbosity>6) { std::clog << "EulerIntegrator::integration_step(VectorFieldInterface,Rectangle,time_type) const" << std::endl; }
  
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"EulerIntegrator::integration_step(VectorFieldInterface,Rectangle,time_type) const");
  
  const VectorFieldInterface<R>& vf(vector_field);
  Rectangle<R> r=initial_set;
  Rectangle<R> q=estimate_flow_bounds(vf,r,step_size);
  
  Interval<R> h=step_size;      
  Vector< Interval<R> > fq=vf(q);
  r=r+(h*fq);
  
  if(verbosity>0) {
    std::clog << "suggested stepsize=" << step_size << std::endl;
    
    std::clog << "stepsize=" << h << std::endl;
    std::clog << "bound=" << q << std::endl;
    
    std::clog << "derivative=" << fq << std::endl;
    
    std::clog << "position=" << r << std::endl;
  }
  return r;
}



template<class R>
Geometry::Rectangle<R> 
Evaluation::EulerIntegrator<R>::reachability_step(const System::VectorFieldInterface<R>& vector_field, 
                                                  const Geometry::Rectangle<R>& initial_set, 
                                                  time_type& step_size) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  if(verbosity>6) { std::clog << "EulerIntegrator::reachability_step(VectorFieldInterface,Rectangle,time_type) const" << std::endl; }
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set(),"EulerIntegrator::reachability_step(VectorFieldInterface,Rectangle,time_type) const");
  
  const VectorFieldInterface<R>& vf(vector_field);
  Rectangle<R> r=initial_set;
  time_type& h=step_size;
  
  Rectangle<R> q=estimate_flow_bounds(vf,r,h);
  
  Vector< Interval<R> > fq=vf(q);
  
  r=r+Vector< Interval<R> >(Interval<R>(R(0),h)*fq);
  
  if(verbosity>1) { 
    std::clog << "suggested stepsize=" << step_size << std::endl;
    
    std::clog << "stepsize=" << h << std::endl;
    std::clog << "bound=" << q << std::endl;
    
    std::clog << "derivative=" << fq << std::endl;
    std::clog << "position=" << r << std::endl;
  }
  
  return r;
}



}

