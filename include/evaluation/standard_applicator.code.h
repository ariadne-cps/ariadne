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
Evaluation::StandardApplicator<R>::apply(const MapInterface<R>& f, const Rectangle<R>& r) const
{
  static int& verbosity = Evaluation::applicator_verbosity; 
  ARIADNE_LOG(6,"Rectangle<Float> StandardApplicator::apply(MapInterface f, Rectangle<Float> r)\n");
  ARIADNE_LOG(7,"  r="<<r<<"\n");
  ARIADNE_LOG(8,"  f(r)="<<Rectangle<R>(f.image(Point< Interval<R> >(r)))<<"\n");
  return Rectangle<R>(f.image(Point< Interval<R> >(r)));
}


template<class R>
Zonotope<R,ExactTag>
Evaluation::StandardApplicator<R>::apply(const MapInterface<R>& f, const Zonotope<R,ExactTag>& z) const
{
  typedef Numeric::Interval<R> I;
  static int& verbosity = Evaluation::applicator_verbosity; 
  ARIADNE_LOG(6,"Zonotope StandardApplicator::apply(MapInterface f, Zonotope z)\n");
  ARIADNE_LOG(7,"  z="<<z<<"\n");
  
  const dimension_type d=z.dimension();
  const size_type ng=z.number_of_generators();
  
  const Point<R>& c=z.centre();
  const Matrix<R>& G=z.generators();
  
  Point<I> ic=f(c);
  Matrix<I> Df = f.jacobian(z.bounding_box());
  Matrix<I> iG = Df*G;
  
  Point<R> nc=midpoint(ic);
  Matrix<R> nG(d,ng+d);
  nG(slice(0,d),slice(0,ng))=midpoint(iG);
  for(size_type i=0; i!=d; ++i) {
    R& err=nG(i,ng+i);
    err=ic[i].radius();
    for(size_type j=0; j!=ng; ++j) {
      err=add_up(err,iG(i,j).radius());
    }
  }

  Zonotope<R> result(nc,nG);
  ARIADNE_LOG(8,"  f(z)="<<result<<"\n");
  return result;
}





template<class R>
Zonotope<R,UniformErrorTag>
Evaluation::StandardApplicator<R>::apply(const MapInterface<R>& f, const Zonotope<R,UniformErrorTag>& z) const
{
  static int& verbosity = Evaluation::applicator_verbosity; 
  ARIADNE_LOG(6,"Zontope<UniformErrorTag> StandardApplicator::apply(MapInterface f, Zonotope<UniformErrorTag> z)\n");
  ARIADNE_LOG(7,"  z="<<z<<"\n");
  typedef Interval<R> I;
  
  Point<I> img_centre=f(z.centre());
  Matrix<I> df_on_set = f.jacobian(Point<I>(z.bounding_box()));
  Matrix<I> img_generators = df_on_set*z.generators();
  Zonotope<R,IntervalTag> interval_zonotope(img_centre,img_generators);
  Zonotope<R,UniformErrorTag> result;
  Geometry::over_approximate(result,interval_zonotope);
  ARIADNE_LOG(8,"  f(z)="<<result<<"\n");
  return result;
}


}
