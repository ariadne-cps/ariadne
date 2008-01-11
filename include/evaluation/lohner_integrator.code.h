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

#include "base/stlio.h"
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

#include "evaluation/integrator.h"
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
  return new LohnerIntegrator<R>(*this);
}


template<class R>
Evaluation::LohnerIntegrator<R>::LohnerIntegrator(smoothness_type temporal_order)
  : _integrator(new IntegratorBase<R>(temporal_order, 1u))
{
}



template<class R> inline
std::pair< Numeric::Rational, Geometry::Box<R> >
Evaluation::LohnerIntegrator<R>::flow_bounds(const System::VectorField<R>& vf, 
                                             const Geometry::Box<R>& bx,
                                             const Numeric::Rational& t) const
{
  return this->_integrator->flow_bounds(vf,bx,t);
}


template<class R>
Geometry::Zonotope<R>
Evaluation::LohnerIntegrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                                  const Geometry::Zonotope<R>& initial_set, 
                                                  const Numeric::Rational& step_size, 
                                                  const Geometry::Box<R>& bounding_box) const
{
  ARIADNE_LOG(2,"LohnerIntegrator::integration_step(VectorField,Zonotope,Time,Box)\n");
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Function;
  using namespace Geometry;
  using namespace System;
  
  ARIADNE_LOG(5,"LohnerIntegrator::integration_step(VectorField,Zonotope,Time,Box)\n");
  ARIADNE_LOG(6,"bounding_box="<<bounding_box<<"\n");
  ARIADNE_LOG(6,"initial_set="<<initial_set<<"\n");
  AffineModel<R> flow_model=this->_integrator->affine_flow_model(vector_field,initial_set.centre(),step_size,bounding_box);
  ARIADNE_LOG(6,"flow_model="<<flow_model<<"\n");
  Zonotope<R> flow_set=Geometry::apply(flow_model,initial_set);
  ARIADNE_LOG(6,"flow_set="<<flow_set<<"\n");
  return Geometry::orthogonal_over_approximation(flow_set);
}




template<class R>
Geometry::Zonotope<R> 
Evaluation::LohnerIntegrator<R>::reachability_step(const System::VectorField<R>& vector_field, 
                                                   const Geometry::Zonotope<R>& initial_set,
                                                   const Numeric::Rational& step_size,
                                                   const Geometry::Box<R>& bounding_box) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Function;
  using namespace Geometry;
  using namespace System;
  
  ARIADNE_LOG(6,"LohnerIntegrator::reachability_step(VectorField,Zonotope<Interval>,Interval,Box) const\n");
  Rational half_step_size=step_size/2;

  AffineModel<R> flow_model=this->_integrator->affine_flow_model(vector_field,initial_set.centre(),half_step_size,bounding_box);
  Point<I> phic=flow_model.value();
  Matrix<I> Dphi=flow_model.jacobian();
  Matrix<I> gen=Dphi*initial_set.generators();
  Vector<I> hhf=I(half_step_size)*vector_field(bounding_box);
  Vector<I> err=Dphi*(I(-1,1)*initial_set.error());

  Zonotope<R> result(phic+err,concatenate_columns(gen,hhf));

  return result;
}



template<class R>
std::ostream&
Evaluation::LohnerIntegrator<R>::write(std::ostream& os) const
{
  return os << "LohnerIntegrator( temporal_order=" << this->_integrator->temporal_order() <<" )";
}



}
