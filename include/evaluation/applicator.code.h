/***************************************************************************
 *            applicator.code.h
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
 
#include "applicator.h"

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

#include "geometry/rectangle.h"
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
Evaluation::apply(const MapInterface<R>& f, const Rectangle<R>& r) 
{
  static int& verbosity = Evaluation::applicator_verbosity; 
  ARIADNE_LOG(6,"Rectangle<Float> apply(MapInterface f, Rectangle<Float> r)\n");
  ARIADNE_LOG(7,"  r="<<r<<"\n");
  ARIADNE_LOG(8,"  f(r)="<<Rectangle<R>(f.image(Point< Interval<R> >(r)))<<"\n");
  return Rectangle<R>(f.image(Point< Interval<R> >(r)));
}


template<class R>
Zonotope<R> 
Evaluation::apply(const MapInterface<R>& f, const Zonotope<R>& z)  
{
  static int& verbosity = Evaluation::applicator_verbosity; 
  ARIADNE_LOG(6,"Zonotope<Float> MapEvolver::apply(MapInterface f, Zonotope<Float,Float> z)\n");
  ARIADNE_LOG(7,"  z="<<z<<"\n");
  typedef typename traits<R>::arithmetic_type F;
  
  const size_type m=z.dimension();
  const size_type n=z.dimension();
  
  Vector< Interval<R> > cuboid_vector(m);
  const Interval<R> unit_interval(-1,1);
  for(size_type i=0; i!=cuboid_vector.size(); ++i) {
    cuboid_vector(i)=Interval<R>(-1,1);
  }
  
  const Point<R>& c=z.centre();
  const Matrix<R>& g=z.generators();
  
  Point< Interval<R> > img_centre=f(c);
  Matrix< Interval<R> > df_on_set = f.jacobian(z.bounding_box());
  Matrix< Interval<R> > df_at_centre = f.jacobian(c);
  
  Matrix< Interval<R> > img_generators = df_at_centre*g;
  
  Matrix< Interval<R> > img_generators_inverse = inverse(Matrix< Interval<R> >(img_generators));
  
  Matrix< Interval<R> > img_generators_on_set = df_on_set * g;
  Matrix< Interval<R> > cuboid_transform = img_generators_inverse * img_generators_on_set;
  
  Vector< Interval<R> > new_cuboid = cuboid_transform * cuboid_vector;
  
  R new_cuboid_sup(0);
  for(size_type j=0; j!=n; ++j) {
    new_cuboid_sup=std::max( new_cuboid_sup, R(abs(new_cuboid(j).lower())) );
    new_cuboid_sup=std::max( new_cuboid_sup, R(abs(new_cuboid(j).upper())) );
  }
  
  // FIXME: This is incorrect; need over-approximations
  Point<R> nc=midpoint(img_centre);
  Matrix<R> ng=midpoint(img_generators);
  
  Zonotope<R> result(nc,ng);
  ARIADNE_LOG(8,"  f(z)="<<result<<"\n");
  return result;
}





template<class R>
Zonotope<R,UniformErrorTag> 
Evaluation::apply(const MapInterface<R>& f, const Zonotope<R,UniformErrorTag>& z)  
{
  static int& verbosity = Evaluation::applicator_verbosity; 
  ARIADNE_LOG(6,"Zontope<R,UniformErrorTag> MapEvolver::apply(MapInterface f, Zonotope<R,UniformErrorTag> z)\n");
  ARIADNE_LOG(7,"  z="<<z<<"\n");
  typedef Interval<R> I;
  
  Point<I> img_centre=f(z.centre());
  Matrix<I> df_on_set = f.jacobian(over_approximation(z.bounding_box()));
  Matrix<I> img_generators = df_on_set*z.generators();
  Zonotope<R,IntervalTag> interval_zonotope(img_centre,img_generators);
  Zonotope<R,UniformErrorTag> result;
  Geometry::over_approximate(result,interval_zonotope);
  ARIADNE_LOG(8,"  f(z)="<<result<<"\n");
  return result;
}





template<class BS>
Evaluation::Applicator<BS>::Applicator() 
{
}


template<class R>
Evaluation::Applicator<R>*
Evaluation::Applicator<R>::clone() const 
{
  return new Applicator<R>();
}



template<class R>
Rectangle<R>
Evaluation::Applicator<R>::apply(const MapInterface<R>& f, const Rectangle<R>& r) const
{
  return Evaluation::apply(f,r);
}



template<class R>
Zonotope<R,ExactTag>
Evaluation::Applicator<R>::apply(const MapInterface<R>& f, const Zonotope<R,ExactTag>& z) const
{
  return Evaluation::apply(f,z);
}

template<class R>
Zonotope<R,UniformErrorTag>
Evaluation::Applicator<R>::apply(const MapInterface<R>& f, const Zonotope<R,UniformErrorTag>& z) const
{
  return Evaluation::apply(f,z);
}

template<class R>
void
Evaluation::Applicator<R>::instantiate()
{
  /*
  * bs=0;
  MapInterface<R>* f=0;
  Evaluation::apply(*f,*bs);
  */
}



}
