/***************************************************************************
 *            standard_integrator.code.h
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

#include "standard_flower.h"

#include "base/array.h"
#include "base/exceptions.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/covector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/matrix_function.h"

#include "function/taylor_series.h"
#include "function/taylor_variable.h"
#include "function/taylor_derivative.h"
#include "function/taylor_model.h"

#include "function/affine_variable.h"
#include "function/affine_model.h"

#include "geometry/box.h"

#include "system/vector_field.h"

#include "evaluation/bounder_interface.h"

#include "output/logging.h"

#include "function/taylor_series.code.h"
#include "function/taylor_variable.code.h"

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





template<class R>
Geometry::Point< Numeric::Interval<R> >
Evaluation::StandardFlower<R>::flow_step(const System::VectorField<R>& vector_field, 
                                         const Geometry::Point<I>& initial_point, 
                                         const Numeric::Rational& step_size, 
                                         const Geometry::Box<R>& bounding_box) const
{
  typedef Numeric::Interval<R> I;
  
  // Use second order formula \f$ \Phi(t,p) = p + tf(p) + t^2/2 Df(B)f(B) \f$
  const System::VectorField<R>& vf=vector_field;
  const Geometry::Point<I>& p=initial_point;
  Geometry::Point<I> b=bounding_box;
  I h=step_size;
  
  return p + h * ( vf(p) + (h/2) * ( vf.jacobian(b) * vf(b) ) );
}

template<class R>
LinearAlgebra::Matrix< Numeric::Interval<R> >
Evaluation::StandardFlower<R>::flow_step_jacobian(const System::VectorField<R>& vector_field, 
                                                  const Geometry::Point<I>& initial_point, 
                                                  const Numeric::Rational& step_size, 
                                                  const Geometry::Box<R>& bounding_box) const
{
  // Use first order formula \f$ D\Phi(t,p) = I + t Df(B) W \f$ where W is a bound for D\Phi([0,h],p)
  // Use ||W-I|| < e^{Lh}-1, where L is the  norm of Df
  typedef Numeric::Interval<R> I;
  const dimension_type d=vector_field.dimension();
  const System::VectorField<R>& vf=vector_field;
  //const Geometry::Point<I>& pr=initial_point;
  Geometry::Point<I> b=bounding_box;
  I h=step_size;

  LinearAlgebra::Matrix<I> Id = LinearAlgebra::Matrix<I>::identity(d);

  LinearAlgebra::Matrix<I> Df = vf.jacobian(b);
  R l = LinearAlgebra::norm(Df).upper();
  I e = Numeric::sub_up(Numeric::exp_up(mul_up(h.upper(),l)),R(1))*I(-1,1);

  LinearAlgebra::Matrix<I> W = LinearAlgebra::Matrix<I>::identity(d)+e*LinearAlgebra::Matrix<I>::one(d,d);

  // Perform a couple of steps
  W=Id + I(0,h.upper()) * (Df * W);
  W=Id + h * (Df * W);
  return W;
}






namespace Function {

template<class X>
class TaylorSeriesAffineVariable
  : public TaylorSeries< AffineVariable<X> >
{
 public:
  TaylorSeriesAffineVariable() : TaylorSeries< AffineVariable<X> >() { }
  TaylorSeriesAffineVariable(const TaylorSeries< AffineVariable<X> >& x)
    : TaylorSeries< AffineVariable<X> >(x) { }
  static TaylorSeriesAffineVariable<X> constant(uint n, uint d, X v) {
    TaylorSeries< AffineVariable<X> > x(d);
    x[0]=AffineVariable<X>::constant(n,v);
    for(uint j=1; j<=d; ++j) {
      x[j]=AffineVariable<X>::constant(n,0.0); 
    }
    return x;
  }
  static TaylorSeriesAffineVariable<X> constant_variable(uint n, uint d, X v, uint i) {
    TaylorSeries< AffineVariable<X> > x(d);
    x[0]=AffineVariable<X>::variable(n,v,i);
    for(uint j=1; j<=d; ++j) {
      x[j]=AffineVariable<X>::constant(n,0.0); 
    }
    return x;
  }
  static TaylorSeriesAffineVariable<X> variable(uint n, uint d, X v, uint i) {
    ARIADNE_ASSERT(d>=1);
    TaylorSeries< AffineVariable<X> > x(d);
    x[0]=AffineVariable<X>::variable(n,v,i);
    x[1]=AffineVariable<X>::variable(n,1.0,i);
    for(uint j=2; j<=d; ++j) {
      x[j]=AffineVariable<X>::constant(n,0.0); 
    }
    return x;
  }
};


template<class X>
class TaylorSeriesTaylorVariable
  : public TaylorSeries< TaylorVariable<X> >
{
 public:
  TaylorSeriesTaylorVariable() : TaylorSeries< TaylorVariable<X> >() { }
  TaylorSeriesTaylorVariable(const TaylorSeries< TaylorVariable<X> >& x)
    : TaylorSeries< TaylorVariable<X> >(x) { }
  static TaylorSeriesTaylorVariable<X> constant_variable(uint n, uint ot, uint ox, X v, uint i) {
    TaylorSeries< TaylorVariable<X> > x(ot);
    x[0]=TaylorVariable<X>::variable(n,ox,v,i);
    for(uint j=1; j<=ot; ++j) {
      x[j]=TaylorVariable<X>::constant(n,ox,0.0); 
    }
    return x;
  }
  static TaylorSeriesTaylorVariable<X> variable(uint n, uint ot, uint ox, X v, uint i) {
    TaylorSeries< TaylorVariable<X> > x(ot);
    x[0]=TaylorVariable<X>::variable(n,ox,v,i);
    x[1]=TaylorVariable<X>::variable(n,ox,1.0,i);
    for(uint j=2; j<=ot; ++j) {
      x[j]=TaylorVariable<X>::constant(n,ox,0.0); 
    }
    return x;
  }
};

template<class R>
AffineVariable<R> midpoint(const Function::AffineVariable< Numeric::Interval<R> >& iav)
{
  R v=midpoint(iav.value());
  LinearAlgebra::Covector<R> cv=midpoint(iav.derivative());
  return AffineVariable<R>(v,cv);
}



template<class X> 
array< TaylorSeries< AffineVariable<X> > >
integrate(const TaylorDerivative<X>& vf, const Geometry::Point<X> x)
{
  ARIADNE_ASSERT(vf.result_size()==vf.argument_size());
  ARIADNE_ASSERT(vf.argument_size()==x.dimension());
  dimension_type n=x.dimension();
  smoothness_type d=vf.degree();
  array< TaylorSeries< AffineVariable<X> > > y(n);
  array< TaylorSeries< AffineVariable<X> > > yp(n);
  for(size_type i=0; i!=n; ++i) {
    y[i]=TaylorSeries< TaylorVariable<X> >(0);
    y[i][0]=AffineVariable<X>::variable(n,x[i],i);
  }
  for(uint j=0; j<d; ++j) {
    yp=evaluate(vf,y);
    //cout << "j="<<j<<"\n y="<<y<<"\n yp="<<yp<<endl;
    for(uint i=0; i!=n; ++i) {  
      y[i]=antiderivative(yp[i],y[i][0]);
    }
  } 
  return y;
}


template<class X> 
array< TaylorSeries< TaylorVariable<X> > >
integrate(const TaylorDerivative<X>& vf, const Geometry::Point<X> x, smoothness_type ox)
{
  ARIADNE_ASSERT(vf.result_size()==vf.argument_size());
  ARIADNE_ASSERT(vf.argument_size()==x.dimension());
  dimension_type n=x.dimension();
  smoothness_type ot=vf.degree();
  array< TaylorSeries< TaylorVariable<X> > > y(n);
  array< TaylorSeries< TaylorVariable<X> > > yp(n);
  for(size_type i=0; i!=n; ++i) {
    y[i]=TaylorSeries< TaylorVariable<X> >(0);
    y[i][0]=TaylorVariable<X>::variable(n,ox,x[i],i);
  }
  for(uint j=0; j<ot; ++j) {
    yp=evaluate(vf,y);
    //cout << "j="<<j<<"\n y="<<y<<"\n yp="<<yp<<endl;
    for(uint i=0; i!=n; ++i) {  
      y[i]=antiderivative(yp[i],y[i][0]);
    }
  } 
  return y;
}








} // namespace Function







template<class R>
Function::AffineModel<R>
Evaluation::StandardFlower<R>::affine_flow_model(const System::VectorField<R>& vector_field, 
                                                 const Geometry::Point<R>& initial_point, 
                                                 const Geometry::Box<R>& initial_domain, 
                                                 const Numeric::Rational& step_size, 
                                                 const Geometry::Box<R>& bounding_box) const
{
  // Convert from Taylor flow model
  using namespace Function;
  typedef Numeric::Interval<R> I;

  TaylorModel<R> taylor_flow_model=this->taylor_flow_model(vector_field,initial_point,initial_domain,step_size,bounding_box);
  return AffineModel<R>(Geometry::Box<R>(taylor_flow_model.domain()),
                        Geometry::Point<R>(taylor_flow_model.centre()),
                        Geometry::Point<I>(taylor_flow_model.evaluate(taylor_flow_model.centre())),
                        LinearAlgebra::Matrix<I>(taylor_flow_model.jacobian(taylor_flow_model.domain().position_vectors())));

  uint to=this->temporal_order();
  dimension_type n=initial_point.dimension();
  I h=step_size;
  
  // Make dvf contain the vector field derivatives at the centre of the initial set,
  // except for the highest-order term, which contains the derivatives over the entire set.
  TaylorDerivative<I> dvf=vector_field.derivative(Geometry::Point<I>(bounding_box),to);
  TaylorDerivative<I> cvf=vector_field.derivative(initial_point,to-1);
  for(uint i=0; i!=n; ++i) {
    dvf[i].assign(cvf[i]);
  }
 
  // Set up array of flow derivative values
  // Each component is a constant in time and a variable in space.
  array< TaylorSeriesAffineVariable<I> > y(n);
  for(size_type i=0; i!=n; ++i) {
    y[i]=TaylorSeriesAffineVariable<I>::constant_variable(n,0,initial_point[i],i);
  }

  // Compute the Taylor series of the state and first variation
  array< TaylorSeriesAffineVariable<I> > yp(n);
  for(uint j=0; j<to; ++j) {
    yp=evaluate(dvf,y);
    for(uint i=0; i!=n; ++i) {  
      y[i]=antiderivative(yp[i],y[i][0]);
    }
  } 

  //for(uint j=0; j<=to; ++j) { cout << "y["<<j<<"]=\n"; for(uint i=0; i!=n; ++i) { cout << " " << midpoint(y[i][j]) << endl; } }

  // Compute the state and first variation at the final time
  array< AffineVariable<I> > r(n);
  for(uint i=0; i!=n; ++i) {
    r[i]=y[i][0];
  }
  I c=1;
  for(uint j=1; j<=to; ++j) {
    c*=h; c/=j;
    for(uint i=0; i!=n; ++i) {
      r[i] += y[i][j]*c;
    }
  }

  return AffineModel<R>(bounding_box,initial_point,r);
}



template<class R>
Function::TaylorModel<R>
Evaluation::StandardFlower<R>::taylor_flow_model(const System::VectorField<R>& vector_field, 
                                                 const Geometry::Point<R>& initial_point, 
                                                 const Geometry::Box<R>& initial_domain, 
                                                 const Numeric::Rational& step_size, 
                                                 const Geometry::Box<R>& bounding_box) const
{
  uint verbosity=0;
  ARIADNE_LOG(6,"taylor_flow_model(...)\n");
  using namespace Function;
  dimension_type n=initial_point.dimension();
  smoothness_type ot=this->temporal_order();
  smoothness_type ox=this->spacial_order();
  Numeric::Interval<R> h=step_size;

  TaylorDerivative<I> vfc=vector_field.derivative(initial_point,ot);
  TaylorDerivative<I> vfb=vector_field.derivative(bounding_box,ot);

  const array<R>& x=initial_point.data();
  array< TaylorSeries< TaylorVariable<I> > > y(n);
  for(uint i=0; i!=n; ++i) {
    // y[i][0]=TaylorVariable<I>::variable(n,ox,x[i],i);
    y[i][0]=TaylorVariable<I>::variable(n,ox,0.0,i);
  }
  ARIADNE_LOG(7,"y="<<y<<"\n");
  array< TaylorSeries< TaylorVariable<I> > > yp(n);
  for(uint j=0; j<=ot; ++j) {
    yp=Function::evaluate(vfc,y);
    for(uint i=0; i!=n; ++i) {  
      y[i]=antiderivative(yp[i],y[i][0]);
    }
    ARIADNE_LOG(7,"yp="<<yp<<"\ny="<<y<<"\n");
  }

  for(uint i=0; i!=n; ++i) {  
    y[i][0].value()=x[i];
  }
  ARIADNE_LOG(7,"\ny="<<y<<"\n\n");

  TaylorDerivative<I> phi(n,n,ox);
  for(uint j=0; j<=ot; ++j) {
    for(uint i=0; i!=n; ++i) {
      phi[i]+=y[i][j]*Numeric::pow(h,j);
    }
    ARIADNE_LOG(7,"phi="<<phi<<"\n");
  }

  //FIXME: Put rigorous error bounds in flow model
  return Function::TaylorModel<R>(initial_domain,initial_point,phi,phi);
}




}
