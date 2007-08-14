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
    A(i,i)=radius(iv(i));
    A(i,iv.size())=midpoint(iv(i));
  }
  return A;
}

}



namespace Ariadne { 

namespace Evaluation { static int& verbosity = integrator_verbosity; }



template<class R>
Evaluation::LohnerIntegrator<R>*
Evaluation::LohnerIntegrator<R>::clone() const
{
  return new LohnerIntegrator<R>(this->maximum_step_size(),this->lock_to_grid_time(),this->maximum_basic_set_radius());
}


template<class R>
Evaluation::LohnerIntegrator<R>::LohnerIntegrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_basic_set_radius)
  : Base_(maximum_step_size,lock_to_grid_time,maximum_basic_set_radius)
{
}


template<class R>
Geometry::ListSet< Geometry::Zonotope<typename Evaluation::LohnerIntegrator<R>::I> >
Evaluation::LohnerIntegrator<R>::subdivide(const Geometry::Zonotope<I>& basic_set) const
{
  Geometry::Zonotope<I> orthogonal_zonotope(Geometry::orthogonal_over_approximation(basic_set));
  return orthogonal_zonotope.subdivide();
}


template<class R>
Geometry::Point<typename Evaluation::LohnerIntegrator<R>::I>
Evaluation::LohnerIntegrator<R>::bounded_flow(const System::VectorFieldInterface<R>& vector_field, 
                                              const Geometry::Point<I>& initial_point, 
                                              const Geometry::Rectangle<R>& bounding_box, 
                                              const Numeric::Interval<R>& step_size) const
{
  // Use second order formula \f$ \Phi(t,p) = p + tf(p) + t^2/2 Df(B)f(B) \f$
  const System::VectorFieldInterface<R>& vf=vector_field;
  const Geometry::Point<I>& p=initial_point;
  Geometry::Point<I> b=this->refine_flow_bounds(vector_field,Geometry::Rectangle<R>(initial_point),bounding_box,step_size.upper());
  I h=step_size;

  
  return p + h * ( vf(p) + (h/2) * ( vf.jacobian(b) * vf(b) ) );
}


template<class R>
LinearAlgebra::Matrix<typename Evaluation::LohnerIntegrator<R>::I>
Evaluation::LohnerIntegrator<R>::bounded_flow_jacobian(const System::VectorFieldInterface<R>& vector_field, 
                                                       const Geometry::Point<I>& initial_point, 
                                                       const Geometry::Rectangle<R>& bounding_box, 
                                                       const Numeric::Interval<R>& step_size) const
{
  // Don't implement since this is a C^0 integrator.
  throw NotImplemented(__PRETTY_FUNCTION__);
}


/*! Use the formula \f$ y(c+Ge,t) \in  c + tf(c) + \frac{t^2}{2} Df(B) f(B) + ( I + t Df(X) ) G e\f$  */
template<class R>
Geometry::Zonotope<typename Evaluation::LohnerIntegrator<R>::I>
Evaluation::LohnerIntegrator<R>::bounded_integration_step(const System::VectorFieldInterface<R>& vector_field, 
                                                          const Geometry::Zonotope<I>& initial_set, 
                                                          const Geometry::Rectangle<R>& bounding_box, 
                                                          const Numeric::Interval<R>& step_size) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  const Zonotope<I>& z=initial_set;
  const Point<I>& c=z.centre();
  const Matrix<I>& G=z.generators();

  const VectorFieldInterface<R>& vf=vector_field;
  const Rectangle<R>& bb=bounding_box;
  const Interval<R> h=step_size;
  const Rectangle<R> r=z.bounding_box();
  
  Vector<I> fc=vf.image(c);
  Vector<I> f=vf(bb);

  Matrix<I> df=vf.jacobian(bb);
  Matrix<I> dfc=vf.jacobian(c);
  Matrix<I> df0=vf.jacobian(r);


  Point<I> phic=c+h*(fc+(h/2)*(df*f));
  Matrix<I> phiG=G+h*df0*G;

  ARIADNE_LOG(7,"  c="<<c<<"\n")
  ARIADNE_LOG(7,"  G="<<G<<"\n")
  ARIADNE_LOG(7,"  B="<<bb<<"\n")

  ARIADNE_LOG(7,"  f(c)="<<fc<<"\n")
  ARIADNE_LOG(7,"  f(B)="<<f<<"\n")
  ARIADNE_LOG(7,"  df(c) for centre="<<dfc<<"\n")
  ARIADNE_LOG(7,"  df(X)="<<df0<<"\n")
  ARIADNE_LOG(7,"  df(B)="<<df<<"\n")

  ARIADNE_LOG(7,"  new_centre="<<phic<<"\n")
  ARIADNE_LOG(7,"  new_generators="<<phiG<<"\n")

  return Zonotope<I>(phic,phiG);

  
}




template<class R>
Geometry::Zonotope<typename Evaluation::LohnerIntegrator<R>::I> 
Evaluation::LohnerIntegrator<R>::bounded_reachability_step(const System::VectorFieldInterface<R>& vector_field, 
                                                           const Geometry::Zonotope<I>& initial_set,
                                                           const Geometry::Rectangle<R>& bounding_box,
                                                           const Numeric::Interval<R>& step_size) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  ARIADNE_LOG(7,"LohnerIntegrator::reachability_step(VectorFieldInterface,Zonotope<Interval>,Numeric::Interval<R>) const\n")
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"LohnerIntegrator::reachability_step(VectorFieldInterface,Zonotope<Interval>,Numeric::Interval<R>)");
  
  const VectorFieldInterface<R>& vf(vector_field);
  Zonotope<I> z=initial_set;
  const Rectangle<R>& bbox=bounding_box;
  const size_type n=z.dimension();
  const Matrix<R> id=Matrix<R>::identity(n);
  
  Interval<R> hi(0,step_size.upper());
  Interval<R> hh(step_size/2);
  
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
  
  ARIADNE_LOG(7,"stepsize="<<step_size<<"\n");
  ARIADNE_LOG(7,"bound="<<bbox<<"\n");
    
  ARIADNE_LOG(7,"flow="<<f<<"\n");
  ARIADNE_LOG(7,"jacobian="<<df<<"\n");
  ARIADNE_LOG(7,"flow derivative="<<dphi<<"\n");
    
  ARIADNE_LOG(7,"centre="<<c<<"\n");
  ARIADNE_LOG(7,"bounds on centre="<<phic<<"\n");
    
  ARIADNE_LOG(7,"flow times stepsize="<<fh<<"\n");
  ARIADNE_LOG(7,"symmetrised flow="<<zfh<<"\n");
  ARIADNE_LOG(7,"over approximating Matrix="<<mdf<<"\n");
    
  ARIADNE_LOG(7,"approximating zonotope "<<z<<"\n");

  return z;
}







template<class R>
Evaluation::C1LohnerIntegrator<R>*
Evaluation::C1LohnerIntegrator<R>::clone() const
{
  return new C1LohnerIntegrator<R>(this->maximum_step_size(),this->lock_to_grid_time(),this->maximum_basic_set_radius());
}


template<class R>
Evaluation::C1LohnerIntegrator<R>::C1LohnerIntegrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_basic_set_radius)
  : Base_(maximum_step_size,lock_to_grid_time,maximum_basic_set_radius)
{
}


template<class R>
Geometry::ListSet< Geometry::Zonotope<typename Evaluation::C1LohnerIntegrator<R>::I> >
Evaluation::C1LohnerIntegrator<R>::subdivide(const Geometry::Zonotope<I>& basic_set) const
{
  Geometry::Zonotope<I> orthogonal_zonotope(Geometry::orthogonal_over_approximation(basic_set));
  return orthogonal_zonotope.subdivide();
}


template<class R>
Geometry::Point<typename Evaluation::C1LohnerIntegrator<R>::I>
Evaluation::C1LohnerIntegrator<R>::bounded_flow(const System::VectorFieldInterface<R>& vector_field, 
                                              const Geometry::Point<I>& initial_point, 
                                              const Geometry::Rectangle<R>& bounding_box, 
                                              const Numeric::Interval<R>& step_size) const
{
  // Use second order formula \f$ \Phi(t,p) = p + tf(p) + t^2/2 Df(B)f(B) \f$
  const System::VectorFieldInterface<R>& vf=vector_field;
  const Geometry::Point<I>& pt=initial_point;
  Geometry::Rectangle<R> r(pt);
  Geometry::Rectangle<R> bb=bounding_box;
  I h=step_size;

  bb=this->refine_flow_bounds(vf,r,bb,step_size.upper());
  Geometry::Point<I> b=bb;
  
  return pt + h * ( vf(pt) + (h/2) * ( vf.jacobian(b) * vf(b) ) );
}


template<class R>
LinearAlgebra::Matrix<typename Evaluation::C1LohnerIntegrator<R>::I>
Evaluation::C1LohnerIntegrator<R>::bounded_flow_jacobian(const System::VectorFieldInterface<R>& vector_field, 
                                                       const Geometry::Point<I>& initial_point, 
                                                       const Geometry::Rectangle<R>& bounding_box, 
                                                       const Numeric::Interval<R>& step_size) const
{
  // Use first order formula \f$ D\Phi(t,p) = I + t Df(B) W \f$ where W is a bound for D\Phi([0,h],p)
  // Use ||W-I|| < e^{Lh}-1, where L is the  norm of Df
  const dimension_type d=vector_field.dimension();
  const System::VectorFieldInterface<R>& vf=vector_field;
  //const Geometry::Point<I>& pr=initial_point;
  Geometry::Point<I> b=bounding_box;
  I h=step_size;

  LinearAlgebra::Matrix<I> Id = LinearAlgebra::Matrix<I>::identity(d);

  LinearAlgebra::Matrix<I> Df = vf.jacobian(b);
  R l = LinearAlgebra::norm(Df).upper();
  I e = Numeric::sub_up(Numeric::exp_up(mul_up(h.upper(),l)),R(1))*I(-1,1);

  LinearAlgebra::Matrix<I> W = LinearAlgebra::Matrix<I>::identity(d)+e*LinearAlgebra::Matrix<I>::one(d,d);

  ARIADNE_LOG(7,"Df="<<Df<<"\nW="<<W<<"\n");

  // Perform a couple of steps
  ARIADNE_LOG(7,"  W0="<<W<<"\n");
  W=Id + I(0,h.upper()) * (Df * W);
  ARIADNE_LOG(7,"  W1="<<W<<"\n");
  W=Id + h * (Df * W);
  ARIADNE_LOG(7,"  DPhi="<<W<<"\n");
  return W;
}


template<class R>
Geometry::Zonotope<typename Evaluation::C1LohnerIntegrator<R>::I>
Evaluation::C1LohnerIntegrator<R>::bounded_integration_step(const System::VectorFieldInterface<R>& vector_field, 
                                                          const Geometry::Zonotope<I>& initial_set, 
                                                          const Geometry::Rectangle<R>& bounding_box, 
                                                          const Numeric::Interval<R>& step_size) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  const Zonotope<I>& z=initial_set;
  const Point<I>& c=z.centre();
  const Matrix<I>& G=z.generators();

  const VectorFieldInterface<R>& vf=vector_field;
  const Rectangle<R>& bb=bounding_box;
  const Interval<R>& h=step_size;
  
  Matrix<I> Dphi=bounded_flow_jacobian(vf,c,bb,h);
  Point<I> phic=bounded_flow(vf,c,bb,h);
  Matrix<I> phiG=Dphi*G;
  ARIADNE_LOG(7,"  flow_jacobian="<<Dphi<<"\n");
  ARIADNE_LOG(7,"  new_centre="<<phic<<"\n  new_generators="<<phiG<<"\n");
  return Zonotope<I>(phic,phiG);
}




template<class R>
Geometry::Zonotope<typename Evaluation::C1LohnerIntegrator<R>::I> 
Evaluation::C1LohnerIntegrator<R>::bounded_reachability_step(const System::VectorFieldInterface<R>& vector_field, 
                                                           const Geometry::Zonotope<I>& initial_set,
                                                           const Geometry::Rectangle<R>& bounding_box,
                                                           const Numeric::Interval<R>& step_size) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  ARIADNE_LOG(6,"C1LohnerIntegrator::reachability_step(VectorFieldInterface,Zonotope<Interval>,Interval) const\n");
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"C1LohnerIntegrator::reachability_step(VectorFieldInterface,Zonotope<Interval>,Interval)");
  
  const VectorFieldInterface<R>& vf(vector_field);
  const Zonotope<I>& z=initial_set;
  const Point<I>& c=z.centre();
  const Matrix<I>& G=z.generators();
  const Rectangle<R>& bb=bounding_box;
  const Interval<R>& h=step_size;
  
  const dimension_type d=vf.dimension();
  const Matrix<I> id=Matrix<I>::identity(d);

  Point<I> phic=this->bounded_flow(vf,c,bb,h/2);
  Matrix<I> Dphi=this->bounded_flow_jacobian(vf,c,bb,h/2);
  Vector<I> hhf=(h/2)*vf(bb);

  Zonotope<I> result(phic,Dphi*G,hhf);
  
  
  ARIADNE_LOG(7,"h="<<h.midpoint()<<"\n");
  ARIADNE_LOG(7,"B="<<bb<<"\n");
    
  ARIADNE_LOG(7,"c="<< c<<"\n");
  ARIADNE_LOG(7,"Phi(c,h/2)="<<phic<<"\n");
  ARIADNE_LOG(7,"DPhi(c,h/2)="<<Dphi<<"\n");
  ARIADNE_LOG(7,"(h/2)f(B)="<<hhf<<"\n");
    
  ARIADNE_LOG(7,"z="<<result);

  return result;
}



}
