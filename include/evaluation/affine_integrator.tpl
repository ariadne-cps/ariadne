/***************************************************************************
 *            affine_integrator.tpl
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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

#include "affine_integrator.h"

#include "../base/array.h"
#include "../debug.h"
#include "../exceptions.h"
#include "../numeric/arithmetic.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"

#include "../linear_algebra/matrix.h"
#include "../linear_algebra/matrix_function.h"

#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/list_set.h"
#include "../geometry/grid.h"
#include "../geometry/grid_set.h"

#include "../system/vector_field.h"
#include "../system/affine_vector_field.h"

#include "../evaluation/integrator.h"

namespace Ariadne { namespace Evaluation { extern int verbosity; } }

template<class R>
Ariadne::Evaluation::AffineIntegrator<R>::AffineIntegrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_basic_set_radius)
  : Integrator<R>(maximum_step_size,lock_to_grid_time,maximum_basic_set_radius)
{
}



template<class R>
Ariadne::Geometry::Parallelotope<R> 
Ariadne::Evaluation::AffineIntegrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                      const Geometry::Parallelotope<R>& initial_set, 
                                      time_type& step_size) const
{
  if(verbosity>6) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }

  Geometry::Zonotope<R> phiz=
    this->integration_step(vector_field,static_cast<const Geometry::Zonotope<R>&>(initial_set),step_size);
  return Geometry::Parallelotope<R>(phiz.centre(),phiz.generators());
}



template<class R>
Ariadne::Geometry::Zonotope<R> 
Ariadne::Evaluation::AffineIntegrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                      const Geometry::Zonotope<R>& initial_set, 
                                      time_type& step_size) const
{
  if(verbosity>6) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }

  const System::AffineVectorField<R>& vf=*dynamic_cast<const System::AffineVectorField<R>*>(&vector_field);
  if(!&vf) {
    throw std::runtime_error(std::string(__FUNCTION__)+": dynamic_cast to AffineVectorField failed");
  }
  
  Geometry::Zonotope<R> z=initial_set;
  Interval<R> h=step_size;
 
  R max_error=(norm(z.generators())/Interval<R>(65536)).upper();
  assert(max_error>0);
  
  if(verbosity>7) { 
    std::cerr << "zonotope generators=" << z.generators() << std::endl;
    std::cerr << "maximum allowed error=" << max_error << std::endl;
  
    std::cerr << "jacobian=" << vf.A() << std::endl;
    std::cerr << "step size=" << h << std::endl;
  }
  

  LinearAlgebra::Matrix< Interval<R> > D=LinearAlgebra::exp_Ah_approx(vf.A(),h.upper(),max_error);
  if(verbosity>7) { std::cerr << "approximate derivative=" << D << std::endl; }
  LinearAlgebra::Matrix< Interval<R> > P=LinearAlgebra::exp_Ah_sub_id_div_A_approx(vf.A(),h.upper(),max_error);
  if(verbosity>7) { std::cerr << "twist=" << P << std::endl; }
  
  LinearAlgebra::Vector< Interval<R> > ib=vf.b();
  LinearAlgebra::Matrix< Interval<R> > iD=D;
  if(verbosity>7) { std::cerr << "approximating derivative=" << iD << std::endl; }
  LinearAlgebra::Matrix< Interval<R> > iP=P;
  if(verbosity>7) { std::cerr << "approximating twist=" << iP << std::endl; }
  //LinearAlgebra::Vector< Interval<R> > iC=iD*z.centre().position_vector()+iP*vf.b();
  LinearAlgebra::Vector< Interval<R> > iv1=iD*z.centre().position_vector();
  if(verbosity>7) { std::cerr << "iv1=" << iv1 << std::endl; }
   LinearAlgebra::Vector< Interval<R> > iv2=iP*ib;
  if(verbosity>7) { std::cerr << "iv2=" << iv2 << std::endl; }
  LinearAlgebra::Vector< Interval<R> > iCv=iv1+iv2;
  Geometry::Point< Interval<R> > iC(iCv);
  
  if(verbosity>7) { std::cerr << "interval centre=" << iC << std::endl; }
  
  z=Geometry::over_approximation(Geometry::Zonotope< Interval<R> >(iC,iD*z.generators()));
  if(verbosity>7) { std::cerr << "zonotope=" << z << std::endl; }
  //IntervalZonotope<R> img(iD*z.centre().position_vector()+iP*vf.b(),iD*z.generators());
  
  return z;      
}

template<class R>
Ariadne::Geometry::Zonotope<R> 
Ariadne::Evaluation::AffineIntegrator<R>::reachability_step(const System::VectorField<R>& vector_field, 
                                       const Geometry::Zonotope<R>& initial_set, 
                                       time_type& step_size) const
{
  if(verbosity>6) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }
  throw NotImplemented(__PRETTY_FUNCTION__);

}
