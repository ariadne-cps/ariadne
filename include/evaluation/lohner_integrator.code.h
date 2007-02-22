/***************************************************************************
 *            lohner_integrator.code.h
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

#include "lohner_integrator.h"

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
Ariadne::Evaluation::LohnerIntegrator<R>::LohnerIntegrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_basic_set_radius)
  : Integrator<R>(maximum_step_size,lock_to_grid_time,maximum_basic_set_radius)
{
}


template<class R>
Ariadne::Geometry::Zonotope< Ariadne::Interval<R> >
Ariadne::Evaluation::LohnerIntegrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                      const Geometry::Zonotope< Interval<R> >& initial_set, 
                                      time_type& step_size) const
{
  if(verbosity>6) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }
  if(verbosity>6) { std::cerr << "  step_size=" << conv_approx<double>(step_size) << "  initial_set=" << initial_set << std::endl; }
  const Geometry::Zonotope<I>& z=initial_set;
  const Geometry::Point<I>& c=z.centre();
  const LinearAlgebra::Matrix<I>& G=z.generators();
  time_type suggested_step_size=step_size;
  
  Geometry::Rectangle<R> bbox=z.bounding_box();
  bbox=this->estimate_flow_bounds(vector_field,bbox,step_size);
  if(verbosity>4) { if(suggested_step_size!=step_size) { std::cerr << "  using step_size=" << conv_approx<double>(step_size) << std::endl; } }
  const System::VectorField<R>& vf=vector_field;
  const size_type n=vf.dimension();
  Interval<R> h=step_size;
  const LinearAlgebra::Matrix<I> id=LinearAlgebra::Matrix<I>::identity(n);
  
  LinearAlgebra::Vector<I> f=vf(bbox);
  LinearAlgebra::Matrix<I> df=vf.jacobian(bbox);
  LinearAlgebra::Matrix<I> hdf=h*df;
  LinearAlgebra::Matrix<I> dphi=exp(hdf);
  
  Geometry::Rectangle<R> cbbox(c);
  cbbox=refine_flow_bounds(vf,cbbox,bbox,step_size);
  LinearAlgebra::Vector<I> fc=vf.image(cbbox);
  LinearAlgebra::Matrix<I> dfc=vf.jacobian(cbbox);
  Geometry::Point<I> phic=c+h*fc;
  LinearAlgebra::Matrix<I> phiG=dphi*G;

  if(verbosity>7) {
    std::cerr << "  flow_bounds=" << bbox << std::endl; 
    std::cerr << "  centre_flow_bounds=" << cbbox << std::endl; 
    std::cerr << "  f_for_centre=" << fc << std::endl; 
    std::cerr << "  f_for_set=" << f << std::endl; 
    std::cerr << "  df_for_set=" << df << std::endl; 
    std::cerr << "  hdf_for_set=" << hdf << std::endl; 
    std::cerr << "  exp_hdf_for_set=" << dphi << std::endl; 

    std::cerr << "  new_centre=" << phic << std::endl;
    std::cerr << "  new_generators=" << phiG << std::endl;
  }
  return Geometry::Zonotope<I>(phic,phiG);
}





template<class R>
Ariadne::Geometry::Zonotope<R> 
Ariadne::Evaluation::LohnerIntegrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                      const Geometry::Zonotope<R>& initial_set, 
                                      time_type& step_size) const
{
  if(verbosity>6) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }
  check_equal_dimensions(vector_field,initial_set);

  Geometry::Zonotope< Interval<R> > fuzzy_zonotopic_initial_set(initial_set);
  Geometry::Zonotope< Interval<R> > fuzzy_zonotopic_final_set=this->integration_step(vector_field,fuzzy_zonotopic_initial_set,step_size);
  
  Geometry::Zonotope<R> final_set=Geometry::over_approximation(fuzzy_zonotopic_final_set);
  return final_set;
}



template<class R>
Ariadne::Geometry::Zonotope<R> 
Ariadne::Evaluation::LohnerIntegrator<R>::reachability_step(const System::VectorField<R>& vector_field, 
                                       const Geometry::Zonotope<R>& initial_set, 
                                       time_type& step_size) const
{
  if(verbosity>6) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }
  if(verbosity>6) { std::cerr << "step_size=" << conv_approx<double>(step_size) << "  initial_set=" << initial_set << std::endl; }
  
  typedef typename traits<R>::arithmetic_type F;

  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;

  check_equal_dimensions(vector_field,initial_set);

  const VectorField<R>& vf(vector_field);
  Zonotope<R> z=initial_set;
  const size_type n=z.dimension();
  const Matrix<R> id=Matrix<R>::identity(n);
  
  /* Throws exception if we can't find flow bounds for given stepsize. */
  time_type input_step_size=step_size;
  Rectangle<R> b=estimate_flow_bounds(vf,z.bounding_box(),step_size,4);
  if(verbosity>4) { if(step_size<input_step_size) { std::cerr << "  using step_size=" << conv_approx<double>(step_size); } }
  Interval<R> h=step_size;
  Interval<R> hh=h/R(2);
  
  Vector< Interval<R> > f=vf(b);
  Matrix< Interval<R> > df=vf.jacobian(b);
  
  Matrix< Interval<R> > dphi=id+Interval<R>(0,step_size)*df;
  
  Point<R> c=z.centre();
  Rectangle<R> rc=Geometry::Rectangle<R>(c,c);
  Rectangle<R> cbbox=refine_flow_bounds(vf,rc,b,step_size/2);
  Point< Interval<R> > phic(c.position_vector()+hh*vf(cbbox));
  
  Vector< Interval<R> > fh=hh*f;
  
  Matrix<R> zfh=symmetrize(fh);
  Matrix< Interval<R> > imdf=dphi*z.generators();
  Matrix<R> mdf(imdf.number_of_rows(),imdf.number_of_columns());
  for(size_type i=0; i!=imdf.number_of_rows(); ++i) {
    for(size_type j=0; j!=imdf.number_of_rows(); ++j) {
      mdf(i,j)=(imdf(i,j)>=R(0)) ? imdf(i,j).upper() : imdf(i,j).lower();
    }
  }
  
  z=Geometry::over_approximation(Geometry::Zonotope< Interval<R> >(phic,zfh,mdf));
  if(verbosity>7) {
    std::cerr << "suggested stepsize=" << step_size << std::endl;
    
    std::cerr << "stepsize=" << h << std::endl;
    std::cerr << "half_stepsize=" << h << std::endl;
    std::cerr << "bound=" << b << std::endl;
    
    std::cerr << "flow=" << f << "=" << std::endl;
    std::cerr << "jacobian=" << df << std::endl;
    std::cerr << "flow derivative=" << dphi << std::endl;
      
    std::cerr << "centre=" << c << std::endl;
    std::cerr << "flowed on centre=" << phic << std::endl;
    
    std::cerr << "flow times stepsize=" << fh << std::endl;
    std::cerr << "symmetrised flow=" << zfh << std::endl;
    std::cerr << "over approximating Matrix=" << mdf << std::endl;
    
    std::cerr << "approximating zonotope " << z;
  }

  return z;
}



template<class R>
Ariadne::Geometry::Zonotope< Ariadne::Interval<R> > 
Ariadne::Evaluation::LohnerIntegrator<R>::reachability_step(const System::VectorField<R>& vector_field, 
                                       const Geometry::Zonotope< Interval<R> >& initial_set, 
                                       time_type& step_size) const
{
  if(verbosity>6) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }

  typedef Interval<R> I;

  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;

  check_equal_dimensions(vector_field,initial_set);

  const VectorField<R>& vf(vector_field);
  Zonotope<I> z=initial_set;
  const size_type n=z.dimension();
  const Matrix<R> id=Matrix<R>::identity(n);
  
  /* Throws exception if we can't find flow bounds for given stepsize. */
  Rectangle<R> bbox=estimate_flow_bounds(vf,z.bounding_box(),step_size,256);
  Interval<R> hi(0,step_size);
  Interval<R> hh(step_size/2);
  
  Vector<I> f=vf(bbox);
  Matrix<I> df=vf.jacobian(bbox);
  Matrix<I> dphi=id+hi*df;
  
  Point<I> c=z.centre();
  Point<I> phic=c+Vector<I>(hh*f);
  
  Vector<I> fh=(hh*f);
  
  Matrix<I> zfh=symmetrize(fh);
  Matrix<I> mdf=dphi*z.generators();
  
  z=Zonotope<I>(phic,mdf,zfh);

  if(verbosity>7) {
    std::cerr << "suggested stepsize=" << step_size << std::endl;
      
    std::cerr << "stepsize=" << conv_approx<double>(step_size) << std::endl;
    std::cerr << "bound=" << bbox << std::endl;
    
    std::cerr << "flow=" << f << "=" << std::endl;
    std::cerr << "jacobian=" << df << std::endl;
    std::cerr << "flow derivative=" << dphi << std::endl;
      
    std::cerr << "centre=" << c << std::endl;
    std::cerr << "bounds on centre=" << phic << std::endl;
    
    std::cerr << "flow times stepsize=" << fh << std::endl;
    std::cerr << "symmetrised flow=" << zfh << std::endl;
    std::cerr << "over approximating Matrix=" << mdf << std::endl;
    
    std::cerr << "approximating zonotope " << z;
  }

  return z;
}
