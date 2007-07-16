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

namespace {

using namespace Ariadne;

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



namespace Ariadne { 

namespace Evaluation { static int& verbosity = integrator_verbosity; }



template<class R>
Evaluation::LohnerIntegrator<R>::LohnerIntegrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_basic_set_radius)
  : Base_(maximum_step_size,lock_to_grid_time,maximum_basic_set_radius)
{
}


template<class R>
Geometry::Point<typename Evaluation::LohnerIntegrator<R>::I>
Evaluation::LohnerIntegrator<R>::bounded_flow(const System::VectorFieldInterface<R>& vector_field, 
                                              const Geometry::Point<I>& initial_point, 
                                              const Geometry::Rectangle<R>& bounding_box, 
                                              const time_type& step_size) const
{
  // Use second order formula \f$ \Phi(t,p) = p + tf(p) + t^2/2 Df(B)f(B) \f$
  const System::VectorFieldInterface<R>& vf=vector_field;
  const Geometry::Point<I>& p=initial_point;
  Geometry::Point<I> b=bounding_box;
  I h=step_size;

  return p + h * ( vf(p) + (h/2) * ( vf.jacobian(b) * vf(b) ) );
}


template<class R>
LinearAlgebra::Matrix<typename Evaluation::LohnerIntegrator<R>::I>
Evaluation::LohnerIntegrator<R>::bounded_flow_jacobian(const System::VectorFieldInterface<R>& vector_field, 
                                                       const Geometry::Point<I>& initial_point, 
                                                       const Geometry::Rectangle<R>& bounding_box, 
                                                       const time_type& step_size) const
{
  // Use first order formula \f$ D\Phi(t,p) = I + t Df(B) W \f$ where W is a bound for D\Phi([0,h],p)
  // Use ||W-I|| < e^{Lh}-1, where L is the  norm of Df
  const dimension_type d=vector_field.dimension();
  const System::VectorFieldInterface<R>& vf=vector_field;
  const Geometry::Point<I>& p=initial_point;
  Geometry::Point<I> b=bounding_box;
  I h=step_size;

  LinearAlgebra::Matrix<I> Id = LinearAlgebra::Matrix<I>::identity(d);

  LinearAlgebra::Matrix<I> Df = vf.jacobian(b);
  R l = LinearAlgebra::norm(Df).upper();
  I e = Numeric::sub_up(Numeric::exp_up(mul_up(h.upper(),l)),R(1))*I(-1,1);

  LinearAlgebra::Matrix<I> W(d,d,&e,0,0);
  for(uint i=0; i!=d; ++i) { W(i,i) = W(i,i)+1; }
  verbosity=7;
  ARIADNE_LOG(7,"Df="<<Df<<"\nW="<<W<<"\n");
  return Id + h * (Df * W);
}


template<class R>
Geometry::Zonotope<typename Evaluation::LohnerIntegrator<R>::I>
Evaluation::LohnerIntegrator<R>::bounded_integration_step(const System::VectorFieldInterface<R>& vector_field, 
                                                          const Geometry::Zonotope<I>& initial_set, 
                                                          const Geometry::Rectangle<R>& bounding_box, 
                                                          const time_type& step_size) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  const Zonotope<I>& z=initial_set;
  const Point<I>& c=z.centre();
  const Matrix<I>& G=z.generators();

  const VectorFieldInterface<R>& vf=vector_field;
  const Rectangle<R>& bbox=bounding_box;
  const size_type n=vf.dimension();
  Interval<R> h=step_size;
  const Matrix<I> id=LinearAlgebra::Matrix<I>::identity(n);
  
  Vector<I> f=vf(bounding_box);
  Matrix<I> df=vf.jacobian(bounding_box);
  Matrix<I> hdf=h*df;
  Matrix<I> dphi=exp(hdf);
  
  Rectangle<R> cbbox(c);
  cbbox=refine_flow_bounds(vf,cbbox,bounding_box,step_size);
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
Geometry::Zonotope<typename Evaluation::LohnerIntegrator<R>::I> 
Evaluation::LohnerIntegrator<R>::bounded_reachability_step(const System::VectorFieldInterface<R>& vector_field, 
                                                           const Geometry::Zonotope<I>& initial_set,
                                                           const Geometry::Rectangle<R>& bounding_box,
                                                           const time_type& step_size) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  if(verbosity>6) { std::clog << "LohnerIntegrator::reachability_step(VectorFieldInterface,Zonotope<Interval>,time_type) const" << std::endl; }
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"LohnerIntegrator::reachability_step(VectorFieldInterface,Zonotope<Interval>,time_type)");
  
  const VectorFieldInterface<R>& vf(vector_field);
  Zonotope<I> z=initial_set;
  const Rectangle<R>& bbox=bounding_box;
  const size_type n=z.dimension();
  const Matrix<R> id=Matrix<R>::identity(n);
  
  Interval<R> hi(0,step_size);
  // FIXME: There should be no need to convert to time_type / Rational
  Interval<R> hh(time_type(step_size/2));
  
  Vector<I> f=vf(bbox);
  Matrix<I> df=vf.jacobian(bbox);
  Matrix<I> dphi=id+hi*df;
  
  Point<I> c=z.centre();
  Point<I> phic=c+Vector<I>(hh*f);
  
  Vector<I> fh=(hh*f);
  
  Matrix<I> mdf=dphi*z.generators();
  
  Matrix<I> zfh=symmetrize(fh);
  //z=Zonotope<I>(phic,mdf,zfh);

  //FIXME: Is the below formula correct?
  z=Zonotope<I>(phic,mdf,fh);
  
  if(verbosity>7) {
    std::clog << "suggested stepsize=" << step_size << std::endl;
    
    std::clog << "stepsize=" << conv_approx<double>(step_size) << std::endl;
    std::clog << "bound=" << bbox << std::endl;
    
    std::clog << "flow=" << f << std::endl;
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



}
