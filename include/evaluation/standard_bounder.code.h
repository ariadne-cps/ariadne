/***************************************************************************
 *            standard_bounder.code.h
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
 
#include "base/tuple.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "function/taylor_derivative.h"
#include "geometry/box.h"
#include "system/vector_field_interface.h"
#include "output/logging.h"

#include "bounder_interface.h"
#include "standard_bounder.h"

namespace Ariadne {

static const double DEFAULT_MAXIMUM_STEP_SIZE=0.125;






template<class R>
StandardBounder<R>::StandardBounder()
  : _maximum_step_size(DEFAULT_MAXIMUM_STEP_SIZE)
{
}

template<class R>
StandardBounder<R>::StandardBounder(const Rational& mss)
  : _maximum_step_size(mss)
{
}


template<class R>
StandardBounder<R>*
StandardBounder<R>::clone() const
{
  return new StandardBounder<R>(*this);
}


template<class R>
Rational
StandardBounder<R>::maximum_step_size() const
{
  return this->_maximum_step_size;
}


template<class R>
std::pair<Rational, Box<R> >
StandardBounder<R>::flow_bounds(const VectorField<R>& vf,
                                            const Box<R>& r) const
{
  return flow_bounds(vf,r,this->maximum_step_size());
}

template<class R>
std::pair<Rational, Box<R> >
StandardBounder<R>::flow_bounds(const VectorField<R>& vf,
                                            const Box<R>& r,
                                            const Rational& hmax) const
{
  // Try to find a time h and a set b such that subset(r+Interval<R>(0,h)*vf(b),b) holds
  ARIADNE_LOG(6,"flow_bounds(VectorField,Box,Time hmax)\n");
  ARIADNE_LOG(7,"  r="<<r<<" hmax="<<hmax<<"\n");
  
  ARIADNE_ASSERT(vf.dimension()==r.dimension());

  // Set up constants of the method.
  // TODO: Better estimates of constants
  const R INITIAL_MULTIPLIER=2;
  const R MULTIPLIER=1.125;
  const R BOX_RADIUS_MULTIPLIER=1.03125;
  const uint EXPANSION_STEPS=8;
  const uint REDUCTION_STEPS=8;
  const uint REFINEMENT_STEPS=4;
  Box<R> b,nb;
  Vector<I> eps(r.dimension(),I(-Ariadne::eps<R>(),Ariadne::eps<R>()));
  Vector<I> delta=r.position_vectors()-r.centre().position_vector();
  
  Rational h=hmax;
  Rational hmin=hmax/(1<<REDUCTION_STEPS);
  bool success=false;
  while(!success) {
    ARIADNE_ASSERT(h>hmin);
    Interval<R> ih(0,h);
    b=r+INITIAL_MULTIPLIER*ih*vf(r)+delta;
    for(uint i=0; i!=EXPANSION_STEPS; ++i) {
      Vector<I> df=vf(b);
      nb=r+ih*df;
      ARIADNE_LOG(9,"  h="<<h<<" b="<<b<<" vf="<<vf(b)<<" nb="<<nb<<"\n");
      if(possibly(subset(nb,b))) {
        if(subset(nb,b)) { } else { std::cerr<<"WARNING: bounding box is not strict subset"<<std::endl; }
        success=true;
        break;
      } else {
        b=r+MULTIPLIER*ih*df+delta;
      }
    }
    if(!success) {
      h/=2;
    }
  }

  ARIADNE_ASSERT(possibly(subset(nb,b)));
  b=nb;
  
  Interval<R> ih(0,h);
  for(uint i=0; i!=REFINEMENT_STEPS; ++i) {
     b=r+ih*vf(b);
  }
  
  // Check result of operation
  // We use "possibly" here since the bound may touch 
  ARIADNE_ASSERT(possibly(subset(r+ih*vf(b),b)));
  
  ARIADNE_LOG(7,"  h="<<h<<" b="<<b<<" r+[0,h]*f(b)="<<r+ih*vf(b)<<"\n");

  return std::make_pair(h,b);
}



  
template<class R>
bool
StandardBounder<R>::check_flow_bounds(const VectorField<R>& vf,
                                                  const Box<R>& r,
                                                  const Rational& h,
                                                  const Box<R>& b) const
{
  ARIADNE_LOG(6,"StandardBounder::check_flow_bounds");
  
  
  return subset(r+Interval<R>(0,h)*vf(b),b);
}






template<class R>
Box<R>
StandardBounder<R>::refine_flow_bounds(const VectorField<R>& vector_field,
                                              const Box<R>& initial_set,
                                              const Box<R>& estimated_bounds,
                                              const Rational& step_size) const
{
  ARIADNE_LOG(6,"StandardBounder::refine_flow_bounds(VectorField vf, Box r, Box bb, Time t)\n");
  ARIADNE_LOG(7,"  h="<<step_size.get_d()<<", r="<<initial_set<<", bb="<<estimated_bounds<<", ");
  
  
  const VectorField<R>& vf=vector_field;
  Box<R> rx=initial_set;
  Box<R> b=estimated_bounds;
  Interval<R> h=Interval<R>(0,step_size);
  
  Box<R> xb=rx+h*vf(b);
  Box<R> xxb=rx+h*vf(xb);
  
  ARIADNE_LOG(7,"  nbb="<<xxb<<"\n");

  
  return xb;
}

 


template<class R> inline
std::pair< Rational, TaylorDerivative<Interval<R> > >
StandardBounder<R>::variation_flow_bounds(const VectorField<R>& vf, 
                                                      const Box<R>& bx,
                                                      const Rational& t,
                                                      smoothness_type o) const
{
  Rational h=t;
  I hi=I(t)*I(0,1);
  TaylorDerivative<I> d=TaylorDerivative<I>::variable(bx.position_vectors(),o);
  while(false) {   // expand flow bounds
    TaylorDerivative<I> vfd=vf.derivative(Point<I>(d.value()),1);
    TaylorDerivative<I> nd=bx+compose(vfd,d)*I(2*hi);
    d=nd;
  }
  while(true) {
    d=bx+compose(vf.derivative(Point<I>(d.value()),1),d)*I(2*hi);
  }
  return make_pair(h,d);
}





template<class R>
Matrix< Interval<R> >
StandardBounder<R>::estimate_flow_jacobian_bounds(const VectorField<R>& vf,
                                                              const Box<R>& b,
                                                              const Rational& h) const
{
  dimension_type d=vf.dimension();
  Matrix<I> Df = vf.jacobian(b);
  R l = norm(Df).upper();
  I e = sub_up(exp_up(mul_up(Interval<R>(h).upper(),l)),R(1))*I(-1,1);

  Matrix<I> W = Matrix<I>::identity(d)+e*Matrix<I>::one(d,d);
  return W;
}

  
} // namespace Ariadne
