/***************************************************************************
 *            standard_applicator.code.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 
#include "standard_applicator.h"

#include <iosfwd>
#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/diagonal_matrix.h"

#include "combinatoric/lattice_set.h"

#include "geometry/box.h"
#include "geometry/zonotope.h"
#include "geometry/basic_set_adaptor.h"
#include "geometry/list_set.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "geometry/rectangular_set.h"

#include "system/grid_multimap.h"


#include "system/map.h"
#include "system/discrete_time_system.h"

#include "output/logging.h"

namespace Ariadne {


using namespace Numeric;
using namespace LinearAlgebra;
using namespace Geometry;
using namespace System;

static int& verbosity = Evaluation::applicator_verbosity; 


template<class R>
Rectangle<R>
Evaluation::StandardApplicator<R>::apply(const Map<R>& f, const Rectangle<R>& r) const
{
  static int& verbosity = Evaluation::applicator_verbosity; 
  ARIADNE_LOG(6,"Rectangle<Float> StandardApplicator::apply(Map f, Rectangle<Float> r)\n");
  ARIADNE_LOG(7,"  r="<<r<<"\n");
  ARIADNE_LOG(8,"  f(r)="<<Rectangle<R>(f.image(Point< Interval<R> >(r)))<<"\n");
  return Rectangle<R>(f.image(Point< Interval<R> >(r)));
}



template<class R>
Zonotope<R>
Evaluation::StandardApplicator<R>::apply(const Map<R>& f, const Zonotope<R>& z) const
{
  static int& verbosity = Evaluation::applicator_verbosity; 
  ARIADNE_LOG(6,"Zontope<UniformErrorTag> StandardApplicator::apply(Map f, Zonotope<UniformErrorTag> z)\n");
  ARIADNE_LOG(7,"  z="<<z<<"\n");
  typedef Interval<R> I;
  
  dimension_type d=z.dimension();
  dimension_type m=z.number_of_generators();
  Matrix<I> df = f.jacobian(Point<I>(z.bounding_box()));
  Point<I> ic=f(z.centre());
  Matrix<I> iG = df*z.generators();
  
  Point<R> nc=midpoint(ic);
  Matrix<R> nG=midpoint(iG);
  Vector<R> ne(d);
  for(dimension_type i=0; i!=d; ++i) {
    R& err=ne[i];
    err=add_up(err,ic[i].radius());
    for(size_type j=0; j!=m; ++j) {
      err=add_up(err,iG(i,j).radius());
    }
  }

  if(!(z.error()==0)) {
    for(dimension_type i=0; i!=d; ++i) {
      R& err=ne[i];
      for(size_type j=0; j!=d; ++j) {
        err=add_up(err,mul_up(abs(df(i,j)).upper(),z.error()[j]));
      }
    }
  }

  return Zonotope<R>(nc,nG,ne);
}


}
