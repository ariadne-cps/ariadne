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

#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../combinatoric/lattice_set.h"

#include "../geometry/rectangle.h"
#include "../geometry/zonotope.h"
#include "../geometry/basic_set_adaptor.h"
#include "../geometry/list_set.h"
#include "../geometry/grid.h"
#include "../geometry/grid_set.h"
#include "../geometry/rectangular_set.h"

#include "../system/grid_multimap.h"


#include "../system/map.h"
#include "../system/discrete_time_system.h"

#include "../output/logging.h"

namespace Ariadne {

static int& verbosity = Evaluation::applicator_verbosity; 

template<class R>
Geometry::Rectangle<R> 
Evaluation::apply(const System::MapInterface<R>& f, const Geometry::Rectangle<R>& r) 
{
  ARIADNE_LOG(6,"Rectangle<Float> apply(MapInterface f, Rectangle<Float> r)\n");
  ARIADNE_LOG(7,"  r="<<r<<"\n");
  ARIADNE_LOG(8,"  f(r)="<<Geometry::Rectangle<R>(f.image(Geometry::Point< Numeric::Interval<R> >(r)))<<"\n");
  return Geometry::Rectangle<R>(f.image(Geometry::Point< Numeric::Interval<R> >(r)));
}


template<class R>
Geometry::Zonotope<R> 
Evaluation::apply(const System::MapInterface<R>& f, const Geometry::Zonotope<R>& z)  
{
  ARIADNE_LOG(6,"Zonotope<Float> MapEvolver::apply(MapInterface f, Zonotope<Float,Float> z)\n");
  ARIADNE_LOG(7,"  z="<<z<<"\n");
  typedef typename Numeric::traits<R>::arithmetic_type F;
  
  const size_type m=z.dimension();
  const size_type n=z.dimension();
  
  LinearAlgebra::Vector< Numeric::Interval<R> > cuboid_vector(m);
  const Numeric::Interval<R> unit_interval(-1,1);
  for(size_type i=0; i!=cuboid_vector.size(); ++i) {
    cuboid_vector(i)=Numeric::Interval<R>(-1,1);
  }
  
  const Geometry::Point<R>& c=z.centre();
  const LinearAlgebra::Matrix<R>& g=z.generators();
  
  Geometry::Point< Numeric::Interval<R> > img_centre=f(c);
  LinearAlgebra::Matrix< Numeric::Interval<R> > df_on_set = f.jacobian(z.bounding_box());
  LinearAlgebra::Matrix< Numeric::Interval<R> > df_at_centre = f.jacobian(c);
  
  LinearAlgebra::Matrix< Numeric::Interval<R> > img_generators = df_at_centre*g;
  
  LinearAlgebra::Matrix< Numeric::Interval<R> > img_generators_inverse = LinearAlgebra::inverse(LinearAlgebra::Matrix< Numeric::Interval<R> >(img_generators));
  
  LinearAlgebra::Matrix< Numeric::Interval<R> > img_generators_on_set = df_on_set * g;
  LinearAlgebra::Matrix< Numeric::Interval<R> > cuboid_transform = img_generators_inverse * img_generators_on_set;
  
  LinearAlgebra::Vector< Numeric::Interval<R> > new_cuboid = cuboid_transform * cuboid_vector;
  
  R new_cuboid_sup(0);
  for(size_type j=0; j!=n; ++j) {
    new_cuboid_sup=std::max( new_cuboid_sup, R(abs(new_cuboid(j).lower())) );
    new_cuboid_sup=std::max( new_cuboid_sup, R(abs(new_cuboid(j).upper())) );
  }
  
  // FIXME: This is incorrect; need over-approximations
  Geometry::Point<R> nc=Geometry::midpoint(img_centre);
  LinearAlgebra::Matrix<R> ng=midpoint(img_generators);
  
  Geometry::Zonotope<R> result(nc,ng);
  ARIADNE_LOG(8,"  f(z)="<<result<<"\n");
  return result;
}





template<class R>
Geometry::Zonotope<Numeric::Interval<R>,R> 
Evaluation::apply(const System::MapInterface<R>& f, const Geometry::Zonotope<Numeric::Interval<R>,R>& z)  
{
  ARIADNE_LOG(6,"Zontope<Interval,Float> MapEvolver::apply(MapInterface f, Zonotope<Interval,Float> z)\n");
  ARIADNE_LOG(7,"  z="<<z<<"\n");
  typedef Numeric::Interval<R> I;
  
  Geometry::Point<I> img_centre=f(z.centre());
  LinearAlgebra::Matrix<I> df_on_set = f.jacobian(over_approximation(z.bounding_box()));
  LinearAlgebra::Matrix<I> img_generators = df_on_set*z.generators();
  Geometry::Zonotope<I,I> interval_zonotope(img_centre,img_generators);
  Geometry::Zonotope<I,R> result(over_approximation(interval_zonotope));
  ARIADNE_LOG(8,"  f(z)="<<result<<"\n");
  return result;
}



template<class R>
Geometry::Zonotope< Numeric::Interval<R> > 
Evaluation::apply(const System::MapInterface<R>& f, const Geometry::Zonotope< Numeric::Interval<R> >& z)  
{
  ARIADNE_LOG(6,"Zontope<Interval,Interval> MapEvolver::apply(MapInterface f, Zonotope<Interval,Interval> z)\n");
  ARIADNE_LOG(7,"  z="<<z<<"\n");
  typedef Numeric::Interval<R> I;
  
  Geometry::Point<I> img_centre=f(z.centre());
  LinearAlgebra::Matrix<I> df_on_set = f.jacobian(over_approximation(z.bounding_box()));
  LinearAlgebra::Matrix<I> img_generators = df_on_set*z.generators();
  
  Geometry::Zonotope<I> result(img_centre,img_generators);
  ARIADNE_LOG(8,"  f(z)="<<result<<"\n");
  return result;
}



template<class BS>
Evaluation::Applicator<BS>::Applicator() 
{
}


template<class BS>
Evaluation::Applicator<BS>*
Evaluation::Applicator<BS>::clone() const 
{
  return new Applicator<BS>();
}



template<class BS>
BS
Evaluation::Applicator<BS>::apply(const System::MapInterface<R>& f, const BS& bs) const
{
  return Evaluation::apply(f,bs);
}

template<class BS>
void
Evaluation::Applicator<BS>::instantiate()
{
  BS* bs=0;
  System::MapInterface<R>* f=0;
  Evaluation::apply(*f,*bs);
}



}
