/***************************************************************************
 *            affine_integrator.code.h
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
#include <cassert>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "affine_integrator.h"

#include "../base/array.h"

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

#include "../system/exceptions.h"
#include "../system/vector_field.h"
#include "../system/affine_vector_field.h"

#include "../evaluation/integrator.h"

#include "../output/logging.h"

    
template<class R> 
R
Ariadne::Evaluation::gexp_up(const R& x, uint k)
{
  using namespace Numeric;
  
  R result=div_up(static_cast<R>(1),static_cast<R>(factorial(k)));
  uint n=k;
  R term=result;
  while(add_down(result,term)==result) {
    n=n+1;
    term=div_up(mul_up(result,term),static_cast<R>(n));
    result=add_up(result,term);
  }
  return result;
}

template<class R> 
Ariadne::LinearAlgebra::Vector< Ariadne::Numeric::Interval<R> > 
Ariadne::Evaluation::gexp(
    const LinearAlgebra::Matrix<R>& A, 
    const LinearAlgebra::Vector<R>& b, 
    const time_type& qt, 
    const uint& k)
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  
  //FIXME: Make number of steps depend on precision
  static const int MAX_STEPS=12;
  typedef Interval<R> I;

  int ns=MAX_STEPS;

  Matrix<R> id=LinearAlgebra::Matrix<R>::identity(A.number_of_rows());
  
  Vector<I> result=b/static_cast<R>(factorial(k));
  Vector<I> term=result;
  I t(qt);
  
  for(uint n=k+1; n!=k+ns; ++n) {
    term=(A*term)*(t/static_cast<R>(n));
    result=result+term;
  }
  
  R err=mul_up(pow_up((A.norm()*t).upper(),ns),gexp_up((A.norm()*t).upper(),k+ns));
  I ierr=err*I(-1,1);
  result+=Vector<I>(result.size(),ierr);

  if(Evaluation::verbosity>7) { 
    std::clog << "A=" << A << ",  t=" << qt << ",  k=" << k << "\n" 
              << "gexp(A,t,k)=" << result << ", err=" << err << std::endl; 
  }

  return result;
}

template<class R> 
Ariadne::LinearAlgebra::Matrix< Ariadne::Numeric::Interval<R> > 
Ariadne::Evaluation::gexp(
    const LinearAlgebra::Matrix<R>& A, 
    const time_type& qt, 
    const uint& k)
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  
  //FIXME: Make number of steps depend on precision
  static const int MAX_STEPS=12;
  typedef Interval<R> I;

  int ns=MAX_STEPS;

  Matrix<R> id=LinearAlgebra::Matrix<R>::identity(A.number_of_rows());
  
  Matrix<I> result=id/static_cast<R>(factorial(k));
  Matrix<I> term=result;
  I t(qt);
  
  for(uint n=k+1; n!=k+ns; ++n) {
    term=(A*term)*(t/static_cast<R>(n));
    result=result+term;
  }
  
  R err=mul_up(pow_up((A.norm()*t).upper(),ns),gexp_up((A.norm()*t).upper(),k+ns));
  I ierr=err*I(-1,1);
  result+=Matrix<I>(result.number_of_rows(),result.number_of_columns(),&ierr,0,0);
  
  if(Evaluation::verbosity>7) { 
    std::clog << "A=" << A << ",  t=" << qt << ",  k=" << k << "\n" 
              << "gexp(A,t,k)=" << result << ", err=" << err << std::endl; 
  }
  
  return result;
}

    
template<class R>
Ariadne::Evaluation::AffineIntegrator<R>::AffineIntegrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_basic_set_radius)
  : Base_(maximum_step_size,lock_to_grid_time,maximum_basic_set_radius)
{
}



template<class R>
Ariadne::Geometry::Zonotope< Ariadne::Numeric::Interval<R> > 
Ariadne::Evaluation::AffineIntegrator<R>::integration_step(
    const System::VectorField<R>& vector_field, 
    const Geometry::Zonotope<I>& initial_set, 
    time_type& step_size) const
{
  if(verbosity>6) { std::clog << "AffineIntegrator::integration_step(VectorField,Zonotope<Interval>,time_type) const" << std::endl; }

  //std::type_info info(&vector_field);
  //std::clog << "Vector field type is:" << info.name() << std::endl;
  const System::AffineVectorField<R>* affine_vector_field_ptr=dynamic_cast<const System::AffineVectorField<R>*>(&vector_field);
  if(!affine_vector_field_ptr) {
    throw std::runtime_error(std::string(__FUNCTION__)+": dynamic_cast to AffineVectorField failed");
  }
  return integration_step(*affine_vector_field_ptr,initial_set,step_size);
}



template<class R>
Ariadne::Geometry::Zonotope< Ariadne::Numeric::Interval<R> > 
Ariadne::Evaluation::AffineIntegrator<R>::reachability_step(
    const System::VectorField<R>& vector_field, 
    const Geometry::Zonotope< Numeric::Interval<R> >& initial_set, 
    time_type& step_size) const
{
  if(verbosity>6) { std::clog << "AffineIntegrator::reachability_step(VectorField,Zonotope<Interval>,time_type) const" << std::endl; }

  const System::AffineVectorField<R>* affine_vector_field_ptr=dynamic_cast<const System::AffineVectorField<R>*>(&vector_field);
  if(!affine_vector_field_ptr) {
    throw std::runtime_error(std::string(__FUNCTION__)+": dynamic_cast to AffineVectorField failed");
  }
  
  return reachability_step(*affine_vector_field_ptr,initial_set,step_size);
}



template<class R>
Ariadne::Geometry::Zonotope< Ariadne::Numeric::Interval<R> > 
Ariadne::Evaluation::AffineIntegrator<R>::integration_step(
    const System::AffineVectorField<R>& affine_vector_field, 
    const Geometry::Zonotope< Numeric::Interval<R> >& initial_set, 
    time_type& step_size) const
{
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  

  if(verbosity>6) { std::clog << "AffineIntegrator::integration_step(AffineVectorField,Zonotope<Interval>,time_type) const" << std::endl; }

  const AffineVectorField<R>& vf=affine_vector_field;
  Zonotope<I> z=initial_set;
  I h=step_size;
 
  if(verbosity>7) { 
    std::clog << "zonotope generators=" << z.generators() << std::endl;
  
    std::clog << "jacobian=" << vf.A() << std::endl;
    std::clog << "step size=" << h << std::endl;
  }
  
  // Use the formula x(t) = x0 + h * P * ( A * x0 + b ) 
  //                      = D * x0 + h * P * b
  // where P = gexp(A*h,1) = sum (Ah)^n / (n+1)!
  const Matrix<R>& A=vf.A();
  const Vector<R>& b=vf.b();
  
  Matrix<I> iP=gexp(A,h.upper(),1);
  Matrix<I> iD=iP*(h*A)+Matrix<R>::identity(vf.dimension());

  if(verbosity>7) { std::clog << "approximating derivative=" << iD << std::endl; }
  if(verbosity>7) { std::clog << "approximating twist=" << iP << std::endl; }

  Vector<I> iv1=(iD*z.centre().position_vector());
  if(verbosity>7) { std::clog << "iv1=" << iv1 << std::endl; }
  Vector<I> iv2=h*(iP*b);
  if(verbosity>7) { std::clog << "iv2=" << iv2 << std::endl; }
  Vector<I> icv=iv1+iv2;
  Point<I> ic(icv);
  
  if(verbosity>7) { std::clog << "interval centre=" << ic << std::endl; }
  
  z=Zonotope<I>(ic,iD*z.generators());
  return z;
}



template<class R>
Ariadne::Geometry::Zonotope< Ariadne::Numeric::Interval<R> > 
Ariadne::Evaluation::AffineIntegrator<R>::reachability_step(
    const System::AffineVectorField<R>& vector_field, 
    const Geometry::Zonotope< Numeric::Interval<R> >& initial_set, 
    time_type& step_size) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  if(verbosity>6) { std::clog << "AffineIntegrator::reachability_step(AffineVectorField,Zonotope<Interval>,time_type) const" << std::endl; }


  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"AffineIntegrator::reachability_step(AffineVectorField,Zonotope<Interval>,Time)");

  const AffineVectorField<R>& avf(vector_field);
  Zonotope<I> iz=initial_set;
  const size_type n=vector_field.dimension();
  const Matrix<R> id=Matrix<R>::identity(n);
  time_type hc=step_size/2;
  // FIXME: There should be no need to convert to time_type / Rational
  Interval<R> hh(time_type(step_size/2));
  
  const Matrix<R>& A=avf.A();
  const Vector<R>& b=avf.b();
  
  /* Throws exception if we can't find flow bounds for given stepsize. */
  iz=this->integration_step(avf,iz,hc);
  
  /* Use centre c, generators G and t*(Ac+b), and error (exp(At)-I)Ge+inv(A)(exp(At)-I-At)(Ac+b)
   * 
   */
  const Point<I>& c=iz.centre();
  const Matrix<I>& G=iz.generators();
  Vector<I> Acpb=A*c.position_vector()+b;
  //Acpb=Acpb+b;
  R nrmA=norm(A).upper();
  R nrmAh=(nrmA*hh).upper();
  R err=(gexp_up(nrmAh,1)*nrmA*hh*norm(G)+gexp_up(nrmAh,2)*hh*hh*nrmA*norm(Acpb).upper()).upper();
  Vector<I> errv=Vector<I>(avf.dimension(),Interval<R>(-err,err));
  iz=Zonotope<I>(c+errv,G,hh*Acpb);
  
  return iz;
}
