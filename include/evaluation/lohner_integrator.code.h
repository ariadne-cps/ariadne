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
  return new LohnerIntegrator<R>();
}


template<class R>
Evaluation::LohnerIntegrator<R>::LohnerIntegrator()
{
}



template<class R>
Geometry::Point<typename Evaluation::LohnerIntegrator<R>::I>
Evaluation::LohnerIntegrator<R>::flow_step(const System::VectorFieldInterface<R>& vector_field, 
                                           const Geometry::Point<I>& initial_point, 
                                           const Numeric::Interval<R>& step_size, 
                                           const Geometry::Rectangle<R>& bounding_box) const
{
  // Use second order formula \f$ \Phi(t,p) = p + tf(p) + t^2/2 Df(B)f(B) \f$
  const System::VectorFieldInterface<R>& vf=vector_field;
  const Geometry::Point<I>& p=initial_point;
  Geometry::Point<I> b=bounding_box;
  I h=step_size;

  
  return p + h * ( vf(p) + (h/2) * ( vf.jacobian(b) * vf(b) ) );
}



/*! Use the formula \f$ y(c+Ge,t) \in  c + tf(c) + \frac{t^2}{2} Df(B) f(B) + ( I + t Df(X) ) G e\f$  */
template<class R>
Geometry::Zonotope<typename Evaluation::LohnerIntegrator<R>::I>
Evaluation::LohnerIntegrator<R>::integration_step(const System::VectorFieldInterface<R>& vector_field, 
                                                  const Geometry::Zonotope<I,I>& initial_set, 
                                                  const Numeric::Interval<R>& step_size, 
                                                  const Geometry::Rectangle<R>& bounding_box) const
{
  ARIADNE_LOG(2,"LohnerIntegrator::integration_step(VectorField,Zonotope<I,I>,Time,Rectangle)\n");

  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  const Zonotope<I,I>& z=initial_set;
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
  Zonotope<I> iz(phic,phiG);

  ARIADNE_LOG(7,"  c="<<c<<"\n");
  ARIADNE_LOG(7,"  G="<<G<<"\n");
  ARIADNE_LOG(7,"  B="<<bb<<"\n");

  ARIADNE_LOG(7,"  f(c)="<<fc<<"\n");
  ARIADNE_LOG(7,"  f(B)="<<f<<"\n");
  ARIADNE_LOG(7,"  df(c) for centre="<<dfc<<"\n");
  ARIADNE_LOG(7,"  df(X)="<<df0<<"\n");
  ARIADNE_LOG(7,"  df(B)="<<df<<"\n");

  ARIADNE_LOG(7,"  new_centre="<<phic<<"\n");
  ARIADNE_LOG(7,"  new_generators="<<phiG<<"\n");

  // Zonotope<I,I> oaiz(Geometry::orthogonal_over_approximation(iz));
  // ARIADNE_LOG(7,"  new_regularised_centre="<<oaiz.centre()<<"\n");
  // ARIADNE_LOG(7,"  new_regularised_generators="<<oaiz.generators()<<"\n");

  return iz;
  
}




template<class R>
Geometry::Zonotope<typename Evaluation::LohnerIntegrator<R>::I>
Evaluation::LohnerIntegrator<R>::reachability_step(const System::VectorFieldInterface<R>& vector_field, 
                                                   const Geometry::Zonotope<I,I>& initial_set,
                                                   const Numeric::Interval<R>& step_size,
                                                   const Geometry::Rectangle<R>& bounding_box) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  ARIADNE_LOG(7,"LohnerIntegrator::reachability_step(VectorFieldInterface,Zonotope<Interval>,Numeric::Interval<R>) const\n")
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"LohnerIntegrator::reachability_step(VectorFieldInterface,Zonotope<Interval>,Numeric::Interval<R>)");
  
  const VectorFieldInterface<R>& vf(vector_field);
  Zonotope<I,I> z=initial_set;
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
  Zonotope<I,I> iz=Zonotope<I>(phic,mdf,fh);
  
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

  return iz;
}



template<class R>
Geometry::Zonotope<typename Evaluation::LohnerIntegrator<R>::I,R>
Evaluation::LohnerIntegrator<R>::integration_step(const System::VectorFieldInterface<R>& vector_field, 
                                                  const Geometry::Zonotope<I,R>& initial_set, 
                                                  const Numeric::Interval<R>& step_size, 
                                                  const Geometry::Rectangle<R>& bounding_box) const
{
  ARIADNE_LOG(2,"LohnerIntegrator::integration_step(VectorField,Zonotope<I,R>,Time,Rectangle)\n");
  return Geometry::over_approximation(this->integration_step(vector_field,Geometry::Zonotope<I,I>(initial_set),step_size,bounding_box));
}

template<class R>
Geometry::Zonotope<typename Evaluation::LohnerIntegrator<R>::I,R>
Evaluation::LohnerIntegrator<R>::reachability_step(const System::VectorFieldInterface<R>& vector_field, 
                                                  const Geometry::Zonotope<I,R>& initial_set, 
                                                  const Numeric::Interval<R>& step_size, 
                                                  const Geometry::Rectangle<R>& bounding_box) const
{
  return Geometry::over_approximation(this->reachability_step(vector_field,Geometry::Zonotope<I,I>(initial_set),step_size,bounding_box));
}



template<class R>
std::ostream&
Evaluation::LohnerIntegrator<R>::write(std::ostream& os) const
{
  return os << "LohnerIntegrator( )";
}






template<class R>
Evaluation::C1LohnerIntegrator<R>*
Evaluation::C1LohnerIntegrator<R>::clone() const
{
  return new C1LohnerIntegrator<R>();
}


template<class R>
Evaluation::C1LohnerIntegrator<R>::C1LohnerIntegrator()
{
}


template<class R>
Geometry::Point<typename Evaluation::C1LohnerIntegrator<R>::I>
Evaluation::C1LohnerIntegrator<R>::flow_step(const System::VectorFieldInterface<R>& vector_field, 
                                             const Geometry::Point<I>& initial_point, 
                                             const Numeric::Interval<R>& step_size, 
                                             const Geometry::Rectangle<R>& bounding_box) const
{
  // Use second order formula \f$ \Phi(t,p) = p + tf(p) + t^2/2 Df(B)f(B) \f$
  const System::VectorFieldInterface<R>& vf=vector_field;
  const Geometry::Point<I>& pt=initial_point;
  Geometry::Rectangle<R> r(pt);
  Geometry::Rectangle<R> bb=bounding_box;
  I h=step_size;

  Geometry::Point<I> b=bb;
  
  return pt + h * ( vf(pt) + (h/2) * ( vf.jacobian(b) * vf(b) ) );
}


template<class R>
LinearAlgebra::Matrix<typename Evaluation::C1LohnerIntegrator<R>::I>
Evaluation::C1LohnerIntegrator<R>::flow_step_jacobian(const System::VectorFieldInterface<R>& vector_field, 
                                                      const Geometry::Point<I>& initial_point, 
                                                      const Numeric::Interval<R>& step_size, 
                                                      const Geometry::Rectangle<R>& bounding_box) const
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
Evaluation::C1LohnerIntegrator<R>::integration_step(const System::VectorFieldInterface<R>& vector_field, 
                                                    const Geometry::Zonotope<I,I>& initial_set, 
                                                    const Numeric::Interval<R>& step_size, 
                                                    const Geometry::Rectangle<R>& bounding_box) const
{
  ARIADNE_LOG(2,"C1LohnerIntegrator::integration_step(VectorField,Zonotope<I,I>,Time,Rectangle)\n");
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
  
  Matrix<I> Dphi=flow_step_jacobian(vf,c,h,bb);
  Point<I> phic=flow_step(vf,c,h,bb);
  Matrix<I> phiG=Dphi*G;
  ARIADNE_LOG(7,"  flow_jacobian="<<Dphi<<"\n");
  ARIADNE_LOG(7,"  new_centre="<<phic<<"\n  new_generators="<<phiG<<"\n");
  Zonotope<I> iz(phic,phiG);
  return iz;
  //Zonotope<I> oaiz(Geometry::orthogonal_over_approximation(iz));
  //return oaiz;
}




template<class R>
Geometry::Zonotope<typename Evaluation::C1LohnerIntegrator<R>::I> 
Evaluation::C1LohnerIntegrator<R>::reachability_step(const System::VectorFieldInterface<R>& vector_field, 
                                                     const Geometry::Zonotope<I,I>& initial_set,
                                                     const Numeric::Interval<R>& step_size,
                                                     const Geometry::Rectangle<R>& bounding_box) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  ARIADNE_LOG(6,"C1LohnerIntegrator::reachability_step(VectorFieldInterface,Zonotope<Interval>,Interval,Rectangle) const\n");
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"C1LohnerIntegrator::reachability_step(VectorFieldInterface,Zonotope<Interval>,Interval,Rectangle)");
  
  const VectorFieldInterface<R>& vf(vector_field);
  const Zonotope<I>& z=initial_set;
  const Point<I>& c=z.centre();
  const Matrix<I>& G=z.generators();
  const Rectangle<R>& bb=bounding_box;
  const Interval<R>& h=step_size;
  
  const dimension_type d=vf.dimension();
  const Matrix<I> id=Matrix<I>::identity(d);

  Point<I> phic=this->flow_step(vf,c,h/2,bb);
  Matrix<I> Dphi=this->flow_step_jacobian(vf,c,h/2,bb);
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



template<class R>
std::ostream&
Evaluation::C1LohnerIntegrator<R>::write(std::ostream& os) const
{
  return os << "C1LohnerIntegrator( )";
}



}
