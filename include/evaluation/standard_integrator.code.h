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

#include "lohner_integrator.h"

#include "base/array.h"
#include "base/exceptions.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"

#include "linear_algebra/matrix.h"
#include "linear_algebra/matrix_function.h"

#include "geometry/box.h"
#include "geometry/zonotope.h"

#include "system/vector_field.h"

#include "evaluation/standard_bounder.h"

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

// Provide as a function for convenience
template<class R> 
std::pair< Numeric::Rational,Geometry::Box<R> >
Evaluation::standard_flow_bounds(const System::VectorField<R>& vf, 
                                 const Geometry::Box<R>& ibb, 
                                 const Numeric::Rational& max_h)
{
  return StandardBounder<R>().flow_bounds(vf,ibb,max_h);
}





template<class R>
Geometry::Point< Numeric::Interval<R> >
Evaluation::standard_flow_step(const System::VectorField<R>& vector_field, 
                               const Geometry::Point< Numeric::Interval<R> >& initial_point, 
                               const Numeric::Interval<R>& step_size, 
                               const Geometry::Box<R>& bounding_box) 
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
Evaluation::standard_flow_step_jacobian(const System::VectorField<R>& vector_field, 
                                        const Geometry::Point< Numeric::Interval<R> >& initial_point, 
                                        const Numeric::Interval<R>& step_size, 
                                        const Geometry::Box<R>& bounding_box) 
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


}
