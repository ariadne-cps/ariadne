/***************************************************************************
 *            affine_integrator.code.h
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
#include <cassert>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "affine_integrator.h"

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

namespace Ariadne { namespace Evaluation { extern int verbosity; static bool warning=true; } }

template<class R> 
R
Ariadne::Evaluation::gexp_up(const R& x, uint k)
{
  R result=div_up(static_cast<R>(1),static_cast<R>(k));
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
    const uint& k, 
    const R& err)
{
  check_size(b,A.number_of_rows());
  check_size(b,A.number_of_columns());
  
  LinearAlgebra::Vector< Interval<R> > result=b/static_cast<R>(factorial(k));
  
  LinearAlgebra::Vector< Interval<R> > term=result;
  Numeric::Interval<R> t(qt);
  R nrm=(A.norm()*t).upper();
  
  uint n=k;
  R e=(n<=nrm) ? static_cast<R>(0) : (term.norm() / (static_cast<R>(1)-nrm/n) ).upper();
  while( ( n<=nrm ) || (e>= err) ) {
    n=n+1;
    term=(A*term)*(t/static_cast<R>(n));
    result=result+term;
    e=(n<=nrm) ? static_cast<R>(0) : (term.norm() / (static_cast<R>(1)-nrm/n) ).upper();
    //std::cerr << n << " " << term << " " << result << " " << e << std::endl;
  }
  
  return result+LinearAlgebra::Vector< Interval<R> >(result.size(),Interval<R>(-e,e));
}

    
template<class R>
Ariadne::Evaluation::AffineIntegrator<R>::AffineIntegrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_basic_set_radius)
  : Integrator<R>(maximum_step_size,lock_to_grid_time,maximum_basic_set_radius)
{
}



template<class R>
Ariadne::Geometry::Parallelotope<R> 
Ariadne::Evaluation::AffineIntegrator<R>::integration_step(
    const System::VectorField<R>& vector_field, 
    const Geometry::Parallelotope<R>& initial_set, 
    time_type& step_size) const
{
  if(verbosity>6) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }

  const Geometry::Zonotope<R>& initial_zonotope=static_cast<const Geometry::Zonotope<R>&>(initial_set);
  const System::AffineVectorField<R>* affine_vector_field_ptr=dynamic_cast<const System::AffineVectorField<R>*>(&vector_field);
  if(!affine_vector_field_ptr) {
    // FIXME: dynamic_cast doesn't work with boost python interface
    if(warning) {
      warning=false;
      std::cerr << "\nWarning: using static_cast to AffineVectorField\n" << std::endl;
    }
    affine_vector_field_ptr=static_cast<const System::AffineVectorField<R>*>(&vector_field);
    //throw std::runtime_error(std::string(__FUNCTION__)+": dynamic_cast to AffineVectorField failed");
  }
  Geometry::Zonotope<R> phiz=
    this->integration_step(*affine_vector_field_ptr,initial_zonotope,step_size);
  return Geometry::Parallelotope<R>(phiz.centre(),phiz.generators());
}



template<class R>
Ariadne::Geometry::Zonotope<R> 
Ariadne::Evaluation::AffineIntegrator<R>::integration_step(
    const System::VectorField<R>& vector_field, 
    const Geometry::Zonotope<R>& initial_set, 
    time_type& step_size) const
{
  if(verbosity>6) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }

  //std::type_info info(&vector_field);
  //std::cerr << "Vector field type is:" << info.name() << std::endl;
  const System::AffineVectorField<R>* affine_vector_field_ptr=dynamic_cast<const System::AffineVectorField<R>*>(&vector_field);
  if(!affine_vector_field_ptr) {
    // FIXME: dynamic_cast doesn't work with boost python interface
    if(warning) {
      warning=false;
      std::cerr << "\nWarning: using static_cast to AffineVectorField\n" << std::endl;
    }
    affine_vector_field_ptr=static_cast<const System::AffineVectorField<R>*>(&vector_field);
    //throw std::runtime_error(std::string(__FUNCTION__)+": dynamic_cast to AffineVectorField failed");
  }
  return integration_step(*affine_vector_field_ptr,initial_set,step_size);
}



template<class R>
Ariadne::Geometry::Zonotope<R> 
Ariadne::Evaluation::AffineIntegrator<R>::reachability_step(
    const System::VectorField<R>& vector_field, 
    const Geometry::Zonotope<R>& initial_set, 
    time_type& step_size) const
{
  if(verbosity>6) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }

  const System::AffineVectorField<R>* affine_vector_field_ptr=dynamic_cast<const System::AffineVectorField<R>*>(&vector_field);
  if(!affine_vector_field_ptr) {
    // FIXME: dynamic_cast doesn't work with boost python interface
    if(warning) {
      warning=false;
      std::cerr << "\nWarning: using static_cast to AffineVectorField\n" << std::endl;
    }
    affine_vector_field_ptr=static_cast<const System::AffineVectorField<R>*>(&vector_field);
    //throw std::runtime_error(std::string(__FUNCTION__)+": dynamic_cast to AffineVectorField failed");
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
  
  if(verbosity>6) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }

  const AffineVectorField<R>& vf=affine_vector_field;
  Zonotope< Interval<R> > z=initial_set;
  Interval<R> h=step_size;
 
  R max_error=(norm(z.generators())/Interval<R>(65536)).upper();
  assert(max_error>0);
  
  if(verbosity>7) { 
    std::cerr << "zonotope generators=" << z.generators() << std::endl;
    std::cerr << "maximum allowed error=" << max_error << std::endl;
  
    std::cerr << "jacobian=" << vf.A() << std::endl;
    std::cerr << "step size=" << h << std::endl;
  }
  

  Matrix< Interval<R> > D=exp_Ah_approx(vf.A(),h.upper(),max_error);
  if(verbosity>7) { std::cerr << "approximate derivative=" << D << std::endl; }
  Matrix< Interval<R> > P=exp_Ah_sub_id_div_A_approx(vf.A(),h.upper(),max_error);
  if(verbosity>7) { std::cerr << "twist=" << P << std::endl; }
  
  Vector< Interval<R> > ib=vf.b();
  Matrix< Interval<R> > iD=D;
  if(verbosity>7) { std::cerr << "approximating derivative=" << iD << std::endl; }
  Matrix< Interval<R> > iP=P;
  if(verbosity>7) { std::cerr << "approximating twist=" << iP << std::endl; }
  //Vector< Interval<R> > iC=iD*z.centre().position_vector()+iP*vf.b();
  Vector< Interval<R> > iv1=iD*z.centre().position_vector();
  if(verbosity>7) { std::cerr << "iv1=" << iv1 << std::endl; }
   Vector< Interval<R> > iv2=iP*ib;
  if(verbosity>7) { std::cerr << "iv2=" << iv2 << std::endl; }
  Vector< Interval<R> > iCv=iv1+iv2;
  Point< Interval<R> > iC(iCv);
  
  if(verbosity>7) { std::cerr << "interval centre=" << iC << std::endl; }
  
  z=Zonotope< Interval<R> >(iC,iD*z.generators());
  return z;
}

template<class R>
Ariadne::Geometry::Zonotope<R> 
Ariadne::Evaluation::AffineIntegrator<R>::integration_step(
    const System::AffineVectorField<R>& affine_vector_field, 
    const Geometry::Zonotope<R>& initial_set, 
    time_type& step_size) const
{
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  if(verbosity>6) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }

  const AffineVectorField<R>& vf=affine_vector_field;
  Zonotope<R> z=initial_set;
  Interval<R> h=step_size;
 
  R max_error=(norm(z.generators())/Interval<R>(65536)).upper();
  assert(max_error>0);
  
  if(verbosity>7) { 
    std::cerr << "zonotope generators=" << z.generators() << std::endl;
    std::cerr << "maximum allowed error=" << max_error << std::endl;
  
    std::cerr << "jacobian=" << vf.A() << std::endl;
    std::cerr << "step size=" << h << std::endl;
  }
  

  Matrix< Interval<R> > D=exp_Ah_approx(vf.A(),h.upper(),max_error);
  if(verbosity>7) { std::cerr << "approximate derivative=" << D << std::endl; }
  Matrix< Interval<R> > P=exp_Ah_sub_id_div_A_approx(vf.A(),h.upper(),max_error);
  if(verbosity>7) { std::cerr << "twist=" << P << std::endl; }
  
  Vector< Interval<R> > ib=vf.b();
  Matrix< Interval<R> > iD=D;
  if(verbosity>7) { std::cerr << "approximating derivative=" << iD << std::endl; }
  Matrix< Interval<R> > iP=P;
  if(verbosity>7) { std::cerr << "approximating twist=" << iP << std::endl; }
  //Vector< Interval<R> > iC=iD*z.centre().position_vector()+iP*vf.b();
  Vector< Interval<R> > iv1=iD*z.centre().position_vector();
  if(verbosity>7) { std::cerr << "iv1=" << iv1 << std::endl; }
   Vector< Interval<R> > iv2=iP*ib;
  if(verbosity>7) { std::cerr << "iv2=" << iv2 << std::endl; }
  Vector< Interval<R> > iCv=iv1+iv2;
  Point< Interval<R> > iC(iCv);
  
  if(verbosity>7) { std::cerr << "interval centre=" << iC << std::endl; }
  
  z=over_approximation(Zonotope< Interval<R> >(iC,iD*z.generators()));
  if(verbosity>7) { std::cerr << "zonotope=" << z << std::endl; }
  //IntervalZonotope<R> img(iD*z.centre().position_vector()+iP*vf.b(),iD*z.generators());
  
  return z;      
}

template<class R>
Ariadne::Geometry::Zonotope<R> 
Ariadne::Evaluation::AffineIntegrator<R>::reachability_step(
    const System::AffineVectorField<R>& vector_field, 
    const Geometry::Zonotope<R>& initial_set, 
    time_type& step_size) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  if(verbosity>6) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }

  typedef Interval<R> I;

  check_equal_dimensions(vector_field,initial_set);

  const AffineVectorField<R>& avf(vector_field);
  const VectorField<R>& vf(vector_field);
  Zonotope<I> iz=initial_set;
  const size_type n=vector_field.dimension();
  const Matrix<R> id=Matrix<R>::identity(n);
  time_type hc=step_size/2;
  Interval<R> hh(step_size/2);
  
  const Matrix<R>& A=avf.A();
  const Vector<R>& b=avf.b();
  
  /* Throws exception if we can't find flow bounds for given stepsize. */
  iz=this->integration_step(avf,iz,hc);
  
  /* Use centre c, generators G and t*(Ac+b), and error (exp(At)-I)Ge+inv(A)(exp(At)-I-At)(Ac+b)
   * 
   */
  const Point< Interval<R> >& c=iz.centre();
  const Matrix< Interval<R> >& G=iz.generators();
  Vector< Interval<R> > Acpb=A*c.position_vector()+b;
  //Acpb=Acpb+b;
  R nrmA=norm(A).upper();
  R nrmAh=(nrmA*hh).upper();
  R err=(gexp_up(nrmAh,1)*nrmA*hh*norm(G)+gexp_up(nrmAh,2)*hh*hh*nrmA*norm(Acpb).upper()).upper();
  Vector<I> errv=Vector<I>(avf.dimension(),Interval<R>(-err,err));
  iz=Zonotope<I>(c+errv,G,hh*Acpb);
  Zonotope<R> z=over_approximation(iz);
  
  return z;
}
