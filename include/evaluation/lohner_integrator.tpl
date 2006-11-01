/***************************************************************************
 *            lohner_integrator.tpl
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

namespace Ariadne { namespace Evaluation { int verbosity=8; } }

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
  const LinearAlgebra::Matrix<I> id=LinearAlgebra::identity_matrix<I>(n);
  
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
    std::cerr << "  exp_hdf_for_set=" << dphi << std::endl; 

    std::cerr << "  new_centre=" << phic << std::endl;
    std::cerr << "  new _generators=" << phiG << std::endl;
  }
  return Geometry::Zonotope<I>(phic,phiG);
}



template<class R>
Ariadne::Geometry::Parallelotope< Ariadne::Interval<R> > 
Ariadne::Evaluation::LohnerIntegrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                      const Geometry::Parallelotope< Interval<R> >& initial_set, 
                                      time_type& step_size) const
{
  if(verbosity>0) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }
  
  const Geometry::Zonotope< Interval<R> >& zonotopic_initial_set=initial_set;
  Geometry::Zonotope< Interval<R> > zonotopic_final_set=this->integration_step(vector_field,zonotopic_initial_set,step_size);
  return Geometry::Parallelotope< Interval<R> >(zonotopic_final_set);
}



template<class R>
Ariadne::Geometry::Parallelotope<R> 
Ariadne::Evaluation::LohnerIntegrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                      const Geometry::Parallelotope<R>& initial_set, 
                                      time_type& step_size) const
{
  if(verbosity>0) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }
  check_equal_dimensions(vector_field,initial_set);

  Geometry::Zonotope< Interval<R> > fuzzy_zonotopic_initial_set(initial_set);
  Geometry::Zonotope< Interval<R> > fuzzy_zonotopic_final_set=this->integration_step(vector_field,fuzzy_zonotopic_initial_set,step_size);
  Geometry::Parallelotope< Interval<R> > fuzzy_parallelotopic_final_set=fuzzy_zonotopic_final_set;
  Geometry::Parallelotope<R> final_set(Geometry::over_approximation(fuzzy_parallelotopic_final_set));
  return final_set;
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

  typedef typename traits<R>::arithmetic_type F;

  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;

  check_equal_dimensions(vector_field,initial_set);

  const VectorField<R>& vf(vector_field);
  Zonotope<R> z=initial_set;
  const size_type n=z.dimension();
  const Matrix<R> id=identity_matrix<R>(n);
  
  /* Throws exception if we can't find flow bounds for given stepsize. */
  Rectangle<R> b=estimate_flow_bounds(vf,z.bounding_box(),step_size,256);
  Interval<R> h=step_size;
  
  Vector< Interval<R> > f=vf(b);
  Matrix< Interval<R> > df=vf.jacobian(b);
  
  Matrix< Interval<R> > dphi=id+Interval<R>(0,step_size)*df;
  
  Point<R> c=z.centre();
  Rectangle<R> rc=Geometry::Rectangle<R>(c,c);
  //Rectangle<R> phic=refine_flow_bounds(vf,rc,b,step_size/2);
  // FIXME: Compute centre
  Point< Interval<R> > phic=c;
  
  Vector< Interval<R> > fh=(h/R(2)*f);
  
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
    std::cerr << "bound=" << b << std::endl;
    
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
  const Matrix<R> id=identity_matrix<R>(n);
  
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

template<class R>
Ariadne::Geometry::Parallelotope<R> 
Ariadne::Evaluation::LohnerIntegrator<R>::integration_step(const System::AffineVectorField<R>& vector_field, 
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
Ariadne::Evaluation::LohnerIntegrator<R>::integration_step(const System::AffineVectorField<R>& vector_field, 
                                      const Geometry::Zonotope<R>& initial_set, 
                                      time_type& step_size) const
{
  if(verbosity>6) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }

  const System::AffineVectorField<R>& vf=vector_field;
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
Ariadne::Evaluation::LohnerIntegrator<R>::reachability_step(const System::AffineVectorField<R>& vector_field, 
                                       const Geometry::Zonotope<R>& initial_set, 
                                       time_type& step_size) const
{
  if(verbosity>6) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }
  throw NotImplemented(__PRETTY_FUNCTION__);

}
