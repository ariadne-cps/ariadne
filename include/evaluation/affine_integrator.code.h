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

#include "../output/logging.h"

namespace Ariadne {

namespace Evaluation { static int& verbosity = integrator_verbosity; }



template<class R> 
R
Evaluation::gexp_up(const R& x, uint k)
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
LinearAlgebra::Vector< Numeric::Interval<R> > 
Evaluation::gexp(
                 const LinearAlgebra::Matrix<R>& A, 
                 const LinearAlgebra::Vector<R>& b, 
                 const Numeric::Interval<R>& t, 
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
  
  for(uint n=k+1; n!=k+ns; ++n) {
    term=(A*term)*(t/static_cast<R>(n));
    result=result+term;
  }
  
  R err=mul_up(pow_up(I(norm(A)*t).upper(),ns),gexp_up(I(norm(A)*t).upper(),k+ns));
  I ierr=err*I(-1,1);
  result+=Vector<I>(result.size(),ierr);
  
  if(Evaluation::verbosity>7) { 
    std::clog << "A=" << A << ",  t=" << t << ",  k=" << k << "\n" 
              << "gexp(A,t,k)=" << result << ", err=" << err << std::endl; 
  }
  
  return result;
}

template<class R> 
LinearAlgebra::Matrix< Numeric::Interval<R> > 
Evaluation::gexp(
                 const LinearAlgebra::Matrix<R>& A, 
                 const Numeric::Interval<R>& t, 
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
  
  for(uint n=k+1; n!=k+ns; ++n) {
    term=(A*term)*(t/static_cast<R>(n));
    result=result+term;
  }
  
  R err=mul_up(pow_up(I(norm(A)*t).upper(),ns),gexp_up(I(norm(A)*t).upper(),k+ns));
  I ierr=err*I(-1,1);
  result+=Matrix<I>(result.number_of_rows(),result.number_of_columns(),&ierr,0,0);
  
  if(Evaluation::verbosity>7) { 
    std::clog << "A=" << A << ",  t=" << t << ",  k=" << k << "\n" 
              << "gexp(A,t,k)=" << result << ", err=" << err << std::endl; 
  }
  
  return result;
}


template<class R>
Evaluation::AffineIntegrator<R>::AffineIntegrator()
{
}



template<class R>
Evaluation::AffineIntegrator<R>*
Evaluation::AffineIntegrator<R>::clone() const
{
  return new AffineIntegrator<R>();
}



template<class R> 
Geometry::Point< Numeric::Interval<R> > 
Evaluation::AffineIntegrator<R>::flow_step(const System::VectorFieldInterface<R>& vf,
                                           const Geometry::Point<I>& p,
                                           const Numeric::Interval<R>& h,
                                           const Geometry::Rectangle<R>& bb) const
{
  const System::AffineVectorField<R> avf=dynamic_cast<const System::AffineVectorField<R>&>(vf);
  if(!&avf) {
    ARIADNE_THROW(std::runtime_error,"AffineIntegrator::flow_step(...)","vector_field is not affine.");
  }
  return this->flow_step(avf,p,h);
}
     
 
template<class R> 
LinearAlgebra::Matrix< Numeric::Interval<R> > 
Evaluation::AffineIntegrator<R>::flow_step_jacobian(const System::VectorFieldInterface<R>& vf,
                                                    const Geometry::Point<I>& p,
                                                    const Numeric::Interval<R>& h,
                                                    const Geometry::Rectangle<R>& bb) const
{
  const System::AffineVectorField<R> avf=dynamic_cast<const System::AffineVectorField<R>&>(vf);
  if(!&avf) {
    ARIADNE_THROW(std::runtime_error,"AffineIntegrator::flow_step_jacobian(...)","vector_field is not affine.");
  }
  return this->flow_step_jacobian(avf,p,h);
}

     



template<class R>
Geometry::Zonotope< Numeric::Interval<R>,R > 
Evaluation::AffineIntegrator<R>::integration_step(const System::VectorFieldInterface<R>& vector_field, 
                                                  const Geometry::Zonotope<I,R>& initial_set, 
                                                  const Numeric::Interval<R>& step_size, 
                                                  const Geometry::Rectangle<R>& bounding_set) const
{
  return Geometry::over_approximation(this->integration_step(vector_field,Geometry::Zonotope<I>(initial_set),step_size,bounding_set));
}



template<class R>
Geometry::Zonotope< Numeric::Interval<R>,R > 
Evaluation::AffineIntegrator<R>::reachability_step(const System::VectorFieldInterface<R>& vector_field, 
                                                   const Geometry::Zonotope<I,R>& initial_set, 
                                                   const Numeric::Interval<R>& step_size, 
                                                   const Geometry::Rectangle<R>& bounding_set) const
{
  return Geometry::over_approximation(this->reachability_step(vector_field,Geometry::Zonotope<I>(initial_set),step_size,bounding_set));
}


template<class R>
Geometry::Zonotope< Numeric::Interval<R> > 
Evaluation::AffineIntegrator<R>::integration_step(const System::VectorFieldInterface<R>& vector_field, 
                                                  const Geometry::Zonotope<I,I>& initial_set, 
                                                  const Numeric::Interval<R>& step_size, 
                                                  const Geometry::Rectangle<R>& bounding_set) const
{
  const System::AffineVectorField<R> affine_vector_field=dynamic_cast<const System::AffineVectorField<R>&>(vector_field);
  if(!&affine_vector_field) {
    ARIADNE_THROW(std::runtime_error,"AffineIntegrator::integration_step(...)","vector_field is not affine.");
  }
  return this->integration_step(affine_vector_field,initial_set,step_size);
}



template<class R>
Geometry::Zonotope< Numeric::Interval<R> > 
Evaluation::AffineIntegrator<R>::reachability_step(const System::VectorFieldInterface<R>& vector_field, 
                                                   const Geometry::Zonotope<I,I>& initial_set, 
                                                   const Numeric::Interval<R>& step_size, 
                                                   const Geometry::Rectangle<R>& bounding_set) const
{
  const System::AffineVectorField<R> affine_vector_field=dynamic_cast<const System::AffineVectorField<R>&>(vector_field);
  if(!&affine_vector_field) {
    ARIADNE_THROW(std::runtime_error,"AffineIntegrator::reachability_step(...)","vector_field is not affine.");
  }
  return this->reachability_step(affine_vector_field,initial_set,step_size);
}






template<class R>
Geometry::Zonotope<typename Evaluation::AffineIntegrator<R>::I,R>
Evaluation::AffineIntegrator<R>::integration_step(const System::AffineVectorField<R>& affine_vector_field, 
                                                  const Geometry::Zonotope<I,R>& initial_set, 
                                                  const Numeric::Interval<R>& step_size) const
{
  return Geometry::over_approximation(this->integration_step(affine_vector_field,Geometry::Zonotope<I,I>(initial_set),step_size));
}



template<class R>
Geometry::Zonotope<typename Evaluation::AffineIntegrator<R>::I,R>
Evaluation::AffineIntegrator<R>::reachability_step(const System::AffineVectorField<R>& affine_vector_field, 
                                                   const Geometry::Zonotope<I,R>& initial_set, 
                                                   const Numeric::Interval<R>& step_size) const
{
  return Geometry::over_approximation(this->reachability_step(affine_vector_field,Geometry::Zonotope<I,I>(initial_set),step_size));
}



template<class R> 
Geometry::Point< Numeric::Interval<R> > 
Evaluation::AffineIntegrator<R>::flow_step(const System::AffineVectorField<R>& avf,
                                           const Geometry::Point<I>& p,
                                           const Numeric::Interval<R>& h) const
{
  ARIADNE_LOG(6,"AffineIntegrator::flow_step(AffineVectorField,Point<Interval>,Time) const\n");
  const LinearAlgebra::Matrix<R>& A=avf.A();
  const LinearAlgebra::Vector<R>& b=avf.b();
  return Geometry::Point<I>(gexp(A,h,0)*p.position_vector() + gexp(A,h,1)*b);
}
     
 
template<class R> 
LinearAlgebra::Matrix< Numeric::Interval<R> > 
Evaluation::AffineIntegrator<R>::flow_step_jacobian(const System::AffineVectorField<R>& avf,
                                                     const Geometry::Point<I>& p,
                                                     const Numeric::Interval<R>& h) const
{
  ARIADNE_LOG(6,"AffineIntegrator::flow_step_jacobian(AffineVectorField,Point<Interval>,Time) const\n");
  const LinearAlgebra::Matrix<R>& A=avf.A();
  return gexp(A,h,0);
}


template<class R>
Geometry::Zonotope< Numeric::Interval<R> > 
Evaluation::AffineIntegrator<R>::integration_step(const System::AffineVectorField<R>& affine_vector_field, 
                                                  const Geometry::Zonotope<I,I>& initial_set, 
                                                  const Numeric::Interval<R>& step_size) const
{
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  
  ARIADNE_LOG(6,"AffineIntegrator::integration_step(AffineVectorField,Zonotope<Interval>,Time) const\n");
  
  const AffineVectorField<R>& vf=affine_vector_field;
  Zonotope<I> z=initial_set;
  I h=step_size;
  
  ARIADNE_LOG(9,"vf="<<vf<<"\n");
  ARIADNE_LOG(7,"z="<<z<<"\n");
  ARIADNE_LOG(7,"h="<<h<<"\n");
  ARIADNE_LOG(7,"vf(c)="<<vf(z.centre())<<", h="<<h<<"\n");
  
  // Use the formula x(t) = x0 + h * P * ( A * x0 + b ) 
  //                      = D * x0 + h * P * b
  // where P = gexp(A*h,1) = sum (Ah)^n / (n+1)!
  const Matrix<R>& A=vf.A();
  const Vector<R>& b=vf.b();
  
  Matrix<I> iP=gexp(A,h,1);
  Matrix<I> iD=iP*(h*A)+Matrix<R>::identity(vf.dimension());
  
  ARIADNE_LOG(9,"approximating derivative=" << iD << "\n");
  ARIADNE_LOG(9,"approximating twist=" << iP << "\n");
  
  Vector<I> iv1=(iD*z.centre().position_vector());
  ARIADNE_LOG(9,"iv1="<<iv1<<"\n");
  Vector<I> iv2=h*(iP*b);
  ARIADNE_LOG(9,"iv2="<<iv2<<"\n");
  Vector<I> icv=iv1+iv2;
  Point<I> ic(icv);
  ARIADNE_LOG(9,"ic="<<ic<<"\n");
  
  z=Zonotope<I>(ic,iD*z.generators());
  ARIADNE_LOG(7,"result="<<z<<"\n");
  return z;
}



template<class R>
Geometry::Zonotope< Numeric::Interval<R> > 
Evaluation::AffineIntegrator<R>::reachability_step(const System::AffineVectorField<R>& vector_field, 
                                                   const Geometry::Zonotope<I,I>& initial_set, 
                                                   const Numeric::Interval<R>& step_size) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  ARIADNE_LOG(6,"AffineIntegrator::reachability_step(AffineVectorField,Zonotope<Interval>,Time)\n");
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"AffineIntegrator::reachability_step(AffineVectorField,Zonotope<Interval>,Time)");
  
  const AffineVectorField<R>& avf(vector_field);
  Zonotope<I> iz=initial_set;
  const size_type n=vector_field.dimension();
  const Matrix<R> id=Matrix<R>::identity(n);
  const I hc=step_size/2;
  // FIXME: There should be no need to convert to time_type / Rational
  Interval<R> hh(step_size/2);
  
  const Matrix<R>& A=avf.A();
  const Vector<R>& b=avf.b();
  
  // No change of step size for affine integrator
  iz=this->integration_step(avf,iz,hc);
  
  /* Use centre c, generators G and t*(Ac+b), and error (exp(At)-I)Ge+inv(A)(exp(At)-I-At)(Ac+b)
   * 
   */
  const Point<I>& c=iz.centre();
  const Matrix<I>& G=iz.generators();
  Vector<I> Acpb=A*c.position_vector()+b;
  //Acpb=Acpb+b;
  R nrmA=norm(A).upper();
  R nrmAh=I(nrmA*hh).upper();
  R err=I(gexp_up(nrmAh,1)*nrmA*hh*norm(G)+gexp_up(nrmAh,2)*hh*hh*nrmA*norm(Acpb).upper()).upper();
  Vector<I> errv=Vector<I>(avf.dimension(),Interval<R>(-err,err));
  iz=Zonotope<I>(c+errv,G,hh*Acpb);
  
  return iz;
}


template<class R>
std::ostream&
Evaluation::AffineIntegrator<R>::write(std::ostream& os) const
{
  return os << "AffineIntegrator( )";
}


}
