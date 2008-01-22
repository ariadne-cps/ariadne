/***************************************************************************
 *            integrator.code.h
 *
 *  Copyright  2007  Pieter Collins
 *
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

#include "integrator.h"

#include "base/array.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "function/taylor_series.h"
#include "function/affine_variable.h"
#include "function/affine_model.h"
#include "function/taylor_variable.h"
#include "function/taylor_model.h"

#include "geometry/point.h"
#include "geometry/box.h"

#include "system/vector_field.h"

#include "evaluation/standard_integrator.h"

#include "output/logging.h"

#include "function/taylor_series.code.h"
#include "function/affine_variable.code.h"


using namespace std;

namespace Ariadne { 

namespace Evaluation { static int& verbosity = integrator_verbosity; }

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



template<class R> inline
std::pair< Numeric::Rational, Geometry::Box<R> >
Evaluation::IntegratorBase<R>::flow_bounds(const System::VectorField<R>& vf, 
                                             const Geometry::Box<R>& bx,
                                             const Numeric::Rational& t) const
{
  return Evaluation::standard_flow_bounds(vf,bx,t);
}

template<class R> inline
std::pair< Numeric::Rational, Function::TaylorDerivative<Numeric::Interval<R> > >
Evaluation::IntegratorBase<R>::variation_flow_bounds(const System::VectorField<R>& vf, 
                                                     const Geometry::Box<R>& bx,
                                                     const Numeric::Rational& t,
                                                     smoothness_type o) const
{
  Numeric::Rational h=t;
  I hi=I(t)*I(0,1);
  Function::TaylorDerivative<I> d=Function::TaylorDerivative<I>::variable(bx.position_vectors(),o);
  while(false) {   // expand flow bounds
    Function::TaylorDerivative<I> vfd=vf.derivative(Geometry::Point<I>(d.value()),1);
    Function::TaylorDerivative<I> nd=bx+compose(vfd,d)*I(2*hi);
    d=nd;
  }
  while(true) {
    d=bx+compose(vf.derivative(Geometry::Point<I>(d.value()),1),d)*I(2*hi);
  }
  return make_pair(h,d);
}



template<class R>
Function::AffineModel<R>
Evaluation::IntegratorBase<R>::affine_flow_model(const System::VectorField<R>& vector_field, 
                                                 const Geometry::Point<R>& initial_point, 
                                                 const Numeric::Rational& step_size, 
                                                 const Geometry::Box<R>& bounding_box) const
{
  using namespace Function;
  typedef Numeric::Interval<R> I;

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
Evaluation::IntegratorBase<R>::taylor_flow_model(const System::VectorField<R>& vector_field, 
                                                 const Geometry::Point<R>& initial_point, 
                                                 const Numeric::Rational& step_size, 
                                                 const Geometry::Box<R>& bounding_box) const
{
  using namespace Function;
  dimension_type n=initial_point.dimension();
  smoothness_type ot=this->temporal_order();
  smoothness_type ox=this->spacial_order();
  Numeric::Interval<R> h=step_size;

  const array<R>& x=initial_point.data();
  array< TaylorSeriesTaylorVariable<I> > y(n);
  for(uint i=0; i!=n; ++i) {
    y[i]=TaylorSeriesTaylorVariable<I>::variable(n,ot,ox,x[i],i);
  }
  std::cout << "y="<<y<<std::endl;
  array< TaylorSeriesTaylorVariable<I> > yp(n);
  for(uint j=0; j<ot; ++j) {
    vector_field.compute(yp.begin(),y.begin());
    for(uint i=0; i!=n; ++i) {  
      y[i]=antiderivative(yp[i],y[i][0]);
    }
  }

  Geometry::Point<I> r(n);
  I t=step_size;
  I c=1;
  for(uint j=0; j<=ot; ++j) {
    I hpj=pow(h,j);
    for(uint i=0; i!=n; ++i) {
      r[i]+=y[i][j].value()*c;
    }
    c*=t;
  }
  std::cout << r << std::endl;

  assert(false); // Not implemented
  // Return function model
}


}
