/***************************************************************************
 *            bounder.code.h
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
 
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../geometry/rectangle.h"
#include "../system/vector_field_interface.h"
#include "../output/logging.h"

#include "bounder_interface.h"
#include "bounder.h"

namespace Ariadne {

namespace Evaluation { static int& verbosity = integrator_verbosity; }

template<class R>
Evaluation::Bounder<R>::Bounder() 
{
}


template<class R>
Evaluation::Bounder<R>*
Evaluation::Bounder<R>::clone() const
{
  return new Bounder<R>();
}


template<class R>
Geometry::Rectangle<R>
Evaluation::Bounder<R>::flow_bounds(const System::VectorFieldInterface<R>& vf,
                                    const Geometry::Rectangle<R>& r,
                                    Numeric::Rational& h) const
{
  Numeric::Rational oldh=h;
  ARIADNE_LOG(5,"Bounder::flow_bounds(VectorField vf, Recangle R, Time h)\n");
  ARIADNE_LOG(6,"  h="<<h.get_d()<<", r="<<r<<"\n");
  Geometry::Rectangle<R> bb=this->estimate_flow_bounds(vf,r,h);
  bb=this->refine_flow_bounds(vf,r,bb,h);
  bb=this->refine_flow_bounds(vf,r,bb,h);
  ARIADNE_LOG(6,"  h="<<h.get_d()<<", bb="<<bb<<"\n");
  assert(h==oldh);
  return bb;
}

template<class R>
Geometry::Rectangle<R>
Evaluation::Bounder<R>::flow_bounds(const System::VectorFieldInterface<R>& vf,
                                    const Geometry::Rectangle<R>& r,
                                    const Numeric::Rational& h) const
{
  Geometry::Rectangle<R> bb=this->estimate_flow_bounds(vf,r,h);
  bb=this->refine_flow_bounds(vf,r,bb,h);
  bb=this->refine_flow_bounds(vf,r,bb,h);
  return bb;
}

  
template<class R>
bool
Evaluation::Bounder<R>::check_flow_bounds(const System::VectorFieldInterface<R>& vf,
                                             const Geometry::Rectangle<R>& r,
                                             const Geometry::Rectangle<R>& b,
                                             const Numeric::Rational& h) const
{
  if(verbosity>6) { std::clog << "Bounder::check_flow_bounds" << std::endl; }
  using namespace Geometry;
  using namespace Numeric;
  return subset(r+Interval<R>(0,h)*vf(b),b);
}



template<class R>
Geometry::Rectangle<R>
Evaluation::Bounder<R>::estimate_flow_bounds(const System::VectorFieldInterface<R>& vf,
                                                   const Geometry::Rectangle<R>& r,
                                                   const Numeric::Rational& h) const
{
  return this->estimate_flow_bounds(vf,r,h,12);
}

  
template<class R>
Geometry::Rectangle<R>
Evaluation::Bounder<R>::estimate_flow_bounds(const System::VectorFieldInterface<R>& vf,
                                                   const Geometry::Rectangle<R>& r,
                                                   const Numeric::Rational& h,
                                                   const unsigned int& maximum_iterations) const
{
  using namespace Geometry;
  using namespace Numeric;
  
  ARIADNE_LOG(8,"Bounder::estimate_flow_bounds" << " (VectorField vf, Rectangle r, Time t, int n)\n");
  ARIADNE_LOG(9,"  h="<<h<<", r="<<r<<", n="<<maximum_iterations<<"\n");
  
  typedef typename Numeric::traits<R>::arithmetic_type F;
  uint iteration=0;
  R multiplier=1.125;
  time_type t=h;
  Rectangle<R> reach(vf.dimension());
  Rectangle<R> bounds(vf.dimension());
  reach=r;
  
  while(t>0) {
    bounds=reach+multiplier*Numeric::Interval<R>(0,h)*vf(reach);
    LinearAlgebra::Vector< Interval<R> > df=vf(bounds);
    
    time_type dt=t;
    for(dimension_type i=0; i!=vf.dimension(); ++i) {
      if(df(i).upper()>0) {
        dt=min(dt,time_type(div_up(sub_up(bounds[i].upper(),reach[i].upper()),df(i).upper())));
      }
      if(df(i).lower()<0) {
        dt=min(dt,time_type(div_up(sub_up(bounds[i].lower(),reach[i].lower()),df(i).lower())));
      }
    }
    reach=bounds;
    t-=dt;
    
    ++iteration;
    if(iteration==maximum_iterations) {
      throw std::runtime_error(std::string(__FUNCTION__)+": Cannot find bounding box for flow");
    }
  }
  ARIADNE_LOG(9,"  bounds="<<bounds<<"\n");
  return reach;
}



template<class R>
Geometry::Rectangle<R>
Evaluation::Bounder<R>::estimate_flow_bounds(const System::VectorFieldInterface<R>& vf,
                                                const Geometry::Rectangle<R>& r,
                                                Numeric::Rational& h) const
{
  using namespace Geometry;
  using namespace Numeric;
  
  if(verbosity>6) { std::clog << "Bounder::estimate_flow_bounds" << std::endl; }
  
  static const unsigned int max_tries=12;
  
  unsigned int max_iterations=12;
  unsigned int remaining_tries=max_tries;
  
  Rectangle<R> bounds(vf.dimension());
  while(bounds.empty()) {
    try {
      bounds=estimate_flow_bounds(vf,r,h,max_iterations);
    }
    catch(std::runtime_error) { 
      h/=2;
      max_iterations+=1;
      --remaining_tries;
      if(remaining_tries==0) {
        throw std::runtime_error(std::string(__FUNCTION__)+": cannnot find bounding box for flow");
      }
    }
  }
  
  if(verbosity>7) { std::clog << "  h=" << conv_approx<double>(h) << "  b=" << bounds << std::endl; }
  
  return bounds;
}



template<class R>
Geometry::Rectangle<R>
Evaluation::Bounder<R>::refine_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                              const Geometry::Rectangle<R>& initial_set,
                                              const Geometry::Rectangle<R>& estimated_bounds,
                                              const Numeric::Rational& step_size) const
{
  ARIADNE_LOG(6,"Bounder::refine_flow_bounds(VectorField vf, Rectangle r, Rectangle bb, Time t)\n");
  ARIADNE_LOG(7,"  h="<<step_size.get_d()<<", r="<<initial_set<<", bb="<<estimated_bounds<<", ");
  
  using namespace System;
  using namespace Geometry;
  using namespace LinearAlgebra;
  using namespace Numeric;
  const VectorFieldInterface<R>& vf=vector_field;
  Rectangle<R> rx=initial_set;
  Rectangle<R> b=estimated_bounds;
  Interval<R> h=Interval<R>(0,step_size);
  
  Rectangle<R> xb=rx+h*vf(b);
  Rectangle<R> xxb=rx+h*vf(xb);
  
  ARIADNE_LOG(7,"  nbb="<<xxb<<"\n");

  
  return xb;
}


template<class R>
Geometry::Rectangle<R>
Evaluation::Bounder<R>::refine_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                              const Geometry::Point<I>& initial_point,
                                              const Geometry::Rectangle<R>& estimated_bounds,
                                              const Numeric::Rational& step_size) const
{
  if(verbosity>6) { std::clog << "Bounder::refine_flow_bounds(VectorField,Point,Rectangle,Time)" << std::endl; }
  
  using namespace System;
  using namespace Geometry;
  using namespace LinearAlgebra;
  using namespace Numeric;
  const VectorFieldInterface<R>& vf=vector_field;
  Rectangle<R> rx(initial_point);
  Rectangle<R> b=estimated_bounds;
  Interval<R> h=Interval<R>(0,step_size);
  
  Rectangle<R> xb=rx+h*vf(b);
  Rectangle<R> xxb=rx+h*vf(xb);
  
  if(verbosity>7) { std::clog << "new_bounds " << xxb << "," << xb << " vs old_bounds " << b << "  " << subset(xb,b) << std::endl; }
  
  return xb;
}


template<class R>
Geometry::Rectangle<R>
Evaluation::Bounder<R>::estimate_interval_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                                         const Geometry::Rectangle<R>& initial_set,
                                                         Numeric::Interval<R>& step_size) const
{
  if(verbosity>6) { std::clog << "Bounder::estimate_flow_bounds(VectorField,Point,TimeInterval)" << std::endl; }
  
  Numeric::Rational backwards_step_size=step_size.lower();
  Numeric::Rational forwards_step_size=step_size.upper();
  
  assert(backwards_step_size==0);
  Geometry::Rectangle<R> estimated_bounds=estimate_flow_bounds(vector_field,initial_set,forwards_step_size);

  // FIXME: need to round inwards here...
  step_size=Numeric::Interval<R>(backwards_step_size,forwards_step_size);
  assert(step_size.lower()>=backwards_step_size);
  assert(step_size.upper()>=forwards_step_size);

  return estimated_bounds;
}


template<class R>
Geometry::Rectangle<R>
Evaluation::Bounder<R>::refine_interval_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                                       const Geometry::Rectangle<R>& initial_set,
                                                       const Geometry::Rectangle<R>& estimated_bounds,
                                                       const Numeric::Interval<R>& step_size) const
{
  if(verbosity>6) { std::clog << "Bounder::refine_flow_bounds(VectorField,Point,Rectangle,TimeInterval)" << std::endl; }
  
  using namespace System;
  using namespace Geometry;
  using namespace LinearAlgebra;
  using namespace Numeric;
  const VectorFieldInterface<R>& vf=vector_field;
  const Rectangle<R>& rx(initial_set);
  Rectangle<R> b=estimated_bounds;
  const Interval<R>& h=step_size;
  
  Rectangle<R> xb=rx+h*vf(b);
  Rectangle<R> xxb=rx+h*vf(xb);
  
  if(verbosity>7) { std::clog << "new_bounds " << xxb << "," << xb << " vs old_bounds " << b << "  " << subset(xb,b) << std::endl; }
  
  return xb;
}





template<class R>
LinearAlgebra::Matrix< Numeric::Interval<R> >
Evaluation::Bounder<R>::estimate_flow_jacobian_bounds(const System::VectorFieldInterface<R>& vf,
                                                      const Geometry::Rectangle<R>& b,
                                                      const Numeric::Rational& h) const
{
  dimension_type d=vf.dimension();
  LinearAlgebra::Matrix<I> Df = vf.jacobian(b);
  R l = LinearAlgebra::norm(Df).upper();
  I e = Numeric::sub_up(Numeric::exp_up(mul_up(Numeric::Interval<R>(h).upper(),l)),R(1))*I(-1,1);

  LinearAlgebra::Matrix<I> W = LinearAlgebra::Matrix<I>::identity(d)+e*LinearAlgebra::Matrix<I>::one(d,d);
  return W;
}

  
} // namespace Ariadne
