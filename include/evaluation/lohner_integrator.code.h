/***************************************************************************
 *            lohner_integrator.code.h
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

#include "lohner_integrator.h"

#include "../base/array.h"
#include "../base/exceptions.h"
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

#include "../output/logging.h"

namespace Ariadne { 
namespace Evaluation { 

extern int verbosity;

template<class R>
LinearAlgebra::Matrix<R>
symmetrize(const LinearAlgebra::Vector< Numeric::Interval<R> >& iv)
{
  LinearAlgebra::Matrix<R> A(iv.size(),iv.size()+1);
  for(size_type i=0; i!=A.number_of_rows(); ++i) {
    A(i,i)=iv(i).radius();
    A(i,iv.size())=iv(i).centre();
  }
  return A;
}

}
}


    
template<class R>
Ariadne::Evaluation::LohnerIntegrator<R>::LohnerIntegrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_basic_set_radius)
  : Base_(maximum_step_size,lock_to_grid_time,maximum_basic_set_radius)
{
}


template<class R>
Ariadne::Geometry::Zonotope<typename Ariadne::Evaluation::LohnerIntegrator<R>::I>
Ariadne::Evaluation::LohnerIntegrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                      const Geometry::Zonotope<I>& initial_set, 
                                      time_type& step_size) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  if(verbosity>6) { std::clog << "LohnerIntegrator::integration_step(VectorField,Zonotope<Interval>,time_type) const" << std::endl; }
  if(verbosity>6) { std::clog << "  step_size=" << conv_approx<double>(step_size) << "  initial_set=" << initial_set << std::endl; }
  const Zonotope<I>& z=initial_set;
  const Point<I>& c=z.centre();
  const Matrix<I>& G=z.generators();
  time_type suggested_step_size=step_size;
  
  Rectangle<R> bbox=z.bounding_box();
  bbox=this->estimate_flow_bounds(vector_field,bbox,step_size);
  if(verbosity>4) { if(suggested_step_size!=step_size) { std::clog << "  using step_size=" << conv_approx<double>(step_size) << std::endl; } }
  const VectorField<R>& vf=vector_field;
  const size_type n=vf.dimension();
  Interval<R> h=step_size;
  const Matrix<I> id=LinearAlgebra::Matrix<I>::identity(n);
  
  Vector<I> f=vf(bbox);
  Matrix<I> df=vf.jacobian(bbox);
  Matrix<I> hdf=h*df;
  Matrix<I> dphi=exp(hdf);
  
  Rectangle<R> cbbox(c);
  cbbox=refine_flow_bounds(vf,cbbox,bbox,step_size);
  Vector<I> fc=vf.image(cbbox);
  Matrix<I> dfc=vf.jacobian(cbbox);
  Point<I> phic=c+h*fc;
  Matrix<I> phiG=dphi*G;

  if(verbosity>7) {
    std::clog << "  flow_bounds=" << bbox << std::endl; 
    std::clog << "  centre_flow_bounds=" << cbbox << std::endl; 
    std::clog << "  f_for_centre=" << fc << std::endl; 
    std::clog << "  f_for_set=" << f << std::endl; 
    std::clog << "  df_for_set=" << df << std::endl; 
    std::clog << "  hdf_for_set=" << hdf << std::endl; 
    std::clog << "  exp_hdf_for_set=" << dphi << std::endl; 

    std::clog << "  new_centre=" << phic << std::endl;
    std::clog << "  new_generators=" << phiG << std::endl;
  }
  return Zonotope<I>(phic,phiG);
}








template<class R>
Ariadne::Geometry::Zonotope<typename Ariadne::Evaluation::LohnerIntegrator<R>::I> 
Ariadne::Evaluation::LohnerIntegrator<R>::reachability_step(const System::VectorField<R>& vector_field, 
                                                            const Geometry::Zonotope<I>& initial_set, 
                                                            time_type& step_size) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  if(verbosity>6) { std::clog << "LohnerIntegrator::reachability_step(VectorField,Zonotope<Interval>,time_type) const" << std::endl; }

  check_equal_dimensions(vector_field,initial_set);

  const VectorField<R>& vf(vector_field);
  Zonotope<I> z=initial_set;
  const size_type n=z.dimension();
  const Matrix<R> id=Matrix<R>::identity(n);
  
  /* Throws exception if we can't find flow bounds for given stepsize. */
  Rectangle<R> bbox=estimate_flow_bounds(vf,z.bounding_box(),step_size,256);
  Interval<R> hi(0,step_size);
  // FIXME: There should be no need to convert to time_type / Rational
  Interval<R> hh(time_type(step_size/2));
  
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
    std::clog << "suggested stepsize=" << step_size << std::endl;
      
    std::clog << "stepsize=" << conv_approx<double>(step_size) << std::endl;
    std::clog << "bound=" << bbox << std::endl;
    
    std::clog << "flow=" << f << "=" << std::endl;
    std::clog << "jacobian=" << df << std::endl;
    std::clog << "flow derivative=" << dphi << std::endl;
      
    std::clog << "centre=" << c << std::endl;
    std::clog << "bounds on centre=" << phic << std::endl;
    
    std::clog << "flow times stepsize=" << fh << std::endl;
    std::clog << "symmetrised flow=" << zfh << std::endl;
    std::clog << "over approximating Matrix=" << mdf << std::endl;
    
    std::clog << "approximating zonotope " << z;
  }

  return z;
}
