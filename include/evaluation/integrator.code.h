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
  static TaylorSeriesAffineVariable<X> variable(uint n, uint d, X v, uint i) {
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


} // namespace Function



template<class R> inline
std::pair< Numeric::Rational, Geometry::Box<R> >
Evaluation::Integrator<R>::flow_bounds(const System::VectorField<R>& vf, 
                                             const Geometry::Box<R>& bx,
                                             const Numeric::Rational& t) const
{
  return Evaluation::standard_flow_bounds(vf,bx,t);
}


template<class R>
Function::AffineModel<R>
Evaluation::Integrator<R>::affine_flow_model(const System::VectorField<R>& vector_field, 
                                             const Geometry::Point<R>& initial_point, 
                                             const Numeric::Rational& step_size, 
                                             const Geometry::Box<R>& bounding_box) const
{
  using namespace Function;
  dimension_type n=initial_point.dimension();
  smoothness_type d=this->temporal_order();
  const array<R>& x=initial_point.data();
  array< TaylorSeriesAffineVariable<I> > y(n);
  for(uint i=0; i!=n; ++i) {
    y[i]=TaylorSeriesAffineVariable<I>::variable(n,d,x[i],i);
  }
  std::cout << "y="<<y<<std::endl;
  array< TaylorSeriesAffineVariable<I> > yp(n);
  for(uint j=0; j<d; ++j) {
    vector_field.compute(yp.begin(),y.begin());
    for(uint i=0; i!=n; ++i) {  
      y[i]=antiderivative(yp[i],y[i][0]);
    }
  }

  Geometry::Point<I> r(n);
  I t=step_size;
  I c=1;
  for(uint j=0; j<=d; ++j) {
    for(uint i=0; i!=n; ++i) {
      r[i]+=y[i][j].value()*c;
    }
    c*=t;
  }
  std::cout << r << std::endl;

  // Return function model
}



template<class R>
Function::TaylorModel<R>
Evaluation::Integrator<R>::flow_model(const System::VectorField<R>& vector_field, 
                                      const Geometry::Point<R>& initial_point, 
                                      const Numeric::Rational& step_size, 
                                      const Geometry::Box<R>& bounding_box) const
{
  using namespace Function;
  dimension_type n=initial_point.dimension();
  smoothness_type ot=this->temporal_order();
  smoothness_type ox=this->spacial_order();

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
    for(uint i=0; i!=n; ++i) {
      r[i]+=y[i][j].value()*c;
    }
    c*=t;
  }
  std::cout << r << std::endl;

  // Return function model
}


}
