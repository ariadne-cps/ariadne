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

#include "base/array.h"
#include "base/exceptions.h"
#include "numeric/arithmetic.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "function/affine_model.h"

#include "geometry/box.h"
#include "geometry/zonotope.h"

#include "system/vector_field.h"

#include "evaluation/standard_integrator.h"

#include "output/logging.h"

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



template<class R> inline
std::pair< Numeric::Rational, Geometry::Box<R> >
Evaluation::LohnerIntegrator<R>::flow_bounds(const System::VectorField<R>& vf, 
                                             const Geometry::Box<R>& bx,
                                             const Numeric::Rational& t) const
{
  return Evaluation::standard_flow_bounds(vf,bx,t);
}


template<class R>
Geometry::Point<typename Evaluation::LohnerIntegrator<R>::I>
Evaluation::LohnerIntegrator<R>::flow_step(const System::VectorField<R>& vector_field, 
                                           const Geometry::Point<I>& initial_point, 
                                           const Numeric::Interval<R>& step_size, 
                                           const Geometry::Box<R>& bounding_box) const
{
  // Use second order formula \f$ \Phi(t,p) = p + tf(p) + t^2/2 Df(B)f(B) \f$
  const System::VectorField<R>& vf=vector_field;
  const Geometry::Point<I>& p=initial_point;
  Geometry::Point<I> b=bounding_box;
  I h=step_size;

  
  return p + h * ( vf(p) + (h/2) * ( vf.jacobian(b) * vf(b) ) );
}

template<class R>
LinearAlgebra::Matrix<typename Evaluation::LohnerIntegrator<R>::I>
Evaluation::LohnerIntegrator<R>::flow_step_jacobian(const System::VectorField<R>& vector_field, 
                                                      const Geometry::Point<I>& initial_point, 
                                                      const Numeric::Interval<R>& step_size, 
                                                      const Geometry::Box<R>& bounding_box) const
{
  // Use first order formula \f$ D\Phi(t,p) = I + t Df(B) W \f$ where W is a bound for D\Phi([0,h],p)
  // Use ||W-I|| < e^{Lh}-1, where L is the  norm of Df
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
Geometry::Zonotope<R>
Evaluation::LohnerIntegrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                                    const Geometry::Zonotope<R>& initial_set, 
                                                    const Numeric::Interval<R>& step_size, 
                                                    const Geometry::Box<R>& bounding_box) const
{
  ARIADNE_LOG(2,"LohnerIntegrator::integration_step(VectorField,Zonotope,Time,Box)\n");
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  const Zonotope<R>& z=initial_set;
  const Point<R>& c=z.centre();

  const VectorField<R>& vf=vector_field;
  const Box<R>& bb=bounding_box;
  const Interval<R>& h=step_size;
  
  Point<I> phic=flow_step(vf,z.centre(),h,bb);
  Matrix<I> Dphi=flow_step_jacobian(vf,z.bounding_box(),h,bb);

  Function::AffineModel<R> model(z.bounding_box(),z.centre(),phic,Dphi);
  
  return orthogonal_over_approximation(Geometry::apply(model,z));
}




template<class R>
Geometry::Zonotope<R> 
Evaluation::LohnerIntegrator<R>::reachability_step(const System::VectorField<R>& vector_field, 
                                                   const Geometry::Zonotope<R>& initial_set,
                                                   const Numeric::Interval<R>& step_size,
                                                   const Geometry::Box<R>& bounding_box) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;
  
  ARIADNE_LOG(6,"LohnerIntegrator::reachability_step(VectorField,Zonotope<Interval>,Interval,Box) const\n");
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"LohnerIntegrator::reachability_step(VectorField,Zonotope<Interval>,Interval,Box)");
  
  const VectorField<R>& vf(vector_field);
  const Zonotope<R>& z=initial_set;
  const Box<R>& bb=bounding_box;
  const Interval<R>& h=step_size;
  
  const dimension_type d=vf.dimension();
  const Matrix<I> id=Matrix<I>::identity(d);

  Point<I> phic=this->flow_step(vf,z.centre(),h/2,bb);
  Matrix<I> Dphi=this->flow_step_jacobian(vf,z.bounding_box(),h/2,bb);
  Matrix<I> gen=Dphi*z.generators();
  Vector<I> hhf=(h/2)*vf(bb);
  Vector<I> err=Dphi*(I(-1,1)*z.error());

  Zonotope<R> result(phic+err,concatenate_columns(gen,hhf));
  
  
  ARIADNE_LOG(7,"h="<<h.midpoint()<<"\n");
  ARIADNE_LOG(7,"B="<<bb<<"\n");
    
  ARIADNE_LOG(7,"c="<< z.centre() <<"\n");
  ARIADNE_LOG(7,"Phi(c,h/2)="<<phic<<"\n");
  ARIADNE_LOG(7,"DPhi(c,h/2)="<<Dphi<<"\n");
  ARIADNE_LOG(7,"(h/2)f(B)="<<hhf<<"\n");
    
  ARIADNE_LOG(7,"z="<<result);

  return result;
}



template<class R>
std::ostream&
Evaluation::LohnerIntegrator<R>::write(std::ostream& os) const
{
  return os << "LohnerIntegrator( )";
}



}
