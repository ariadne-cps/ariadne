/***************************************************************************
 *            kuhn_integrator.code.h
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

#include "kuhn_integrator.h"

#include "base/array.h"
#include "base/exceptions.h"

#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "geometry/box.h"
#include "geometry/zonotope.h"

#include "system/vector_field.h"

#include "evaluation/standard_integrator.h"
#include "evaluation/standard_integrator.code.h"

#include "output/logging.h"


namespace Ariadne { 

namespace Evaluation { static int& verbosity = integrator_verbosity; }



template<class R> inline
std::pair< Numeric::Rational, Geometry::Box<R> >
Evaluation::KuhnIntegrator<R>::flow_bounds(const System::VectorField<R>& vf, 
                                             const Geometry::Box<R>& bx,
                                             const Numeric::Rational& t) const
{
  return Evaluation::standard_flow_bounds(vf,bx,t);
}



template<class R>
Geometry::Zonotope<R>
Evaluation::KuhnIntegrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                                const Geometry::Zonotope<R>& initial_set, 
                                                const Numeric::Rational& step_size, 
                                                const Geometry::Box<R>& flow_bounding_box) const
{
  ARIADNE_LOG(2,"KuhnIntegrator::integration_step(VectorField,Zonotope,Time,Box)\n");

  using namespace LinearAlgebra;
  using namespace Geometry;

  Point<I> nic=Evaluation::standard_flow_step(vector_field,Point<I>(initial_set.centre()),I(step_size),flow_bounding_box);
  Matrix<I> dPhi=Evaluation::standard_flow_step_jacobian(vector_field,Point<I>(initial_set.bounding_box()),I(step_size),flow_bounding_box);
  Matrix<I> niG=dPhi*initial_set.generators();
  
  dimension_type d=initial_set.dimension();
  size_type ng=initial_set.number_of_generators();
  Point<R> nc=midpoint(nic);
  Matrix<R> nG(d,ng+d);
  nG(slice(0,d),slice(0,ng))=midpoint(niG);

  for(size_type i=0; i!=d; ++i) {
    R& err=nG(i,ng+i);
    err=radius(nic[i]);
    for(size_type j=0; j!=ng; ++j) {
      err=add_up(err,midpoint(niG(i,j)));
    }
  }
  
  return Geometry::cascade_over_approximation(Zonotope<R>(nc,nG),this->_cascade_size);
}

template<class R>
Geometry::Zonotope<R>
Evaluation::KuhnIntegrator<R>::reachability_step(const System::VectorField<R>& vector_field, 
                                                 const Geometry::Zonotope<R>& initial_set, 
                                                 const Numeric::Rational& step_size, 
                                                 const Geometry::Box<R>& flow_bounding_box) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;

  const VectorField<R>& vf=vector_field;
  const Zonotope<R>& z=initial_set;
  const Point<R>& c=z.centre();
  const Point<I> ic(c);
  const Matrix<R>& G=z.generators();
  const Point<I> ibpt=z.bounding_box();
  const Box<R>& fbb=flow_bounding_box;
  const Rational& h=step_size;
  const Rational hh=h/2;
  const Interval<R> ihh=hh;

  Point<I> nic=Evaluation::standard_flow_step(vf,ic,ihh,fbb);
  Matrix<I> dPhi=Evaluation::standard_flow_step_jacobian(vf,ibpt,ihh,fbb);
  Vector<I> nif=ihh*vf.evaluate(fbb);
  Matrix<I> niG=dPhi*G;

  dimension_type d=initial_set.dimension();
  size_type ng=initial_set.number_of_generators();
  Point<R> nc=midpoint(nic);
  Matrix<R> nG(d,ng+1u+d);
  nG(slice(0,d),slice(0,ng))=midpoint(niG);
  nG.column(ng)=midpoint(nif);

  for(size_type i=0; i!=d; ++i) {
    R& err=nG(i,ng+1u+i);
    err=radius(nic[i]);
    for(size_type j=0; j!=ng; ++j) {
      err=add_up(err,midpoint(niG(i,j)));
    }
    err=add_up(err,midpoint(nif(i)));
  }
  
  return Zonotope<R>(nc,nG);
  
}



template<class R>
std::ostream&
Evaluation::KuhnIntegrator<R>::write(std::ostream& os) const
{
  return os << "KuhnIntegrator( cascade_size="<<this->_cascade_size<<" )";
}



}
