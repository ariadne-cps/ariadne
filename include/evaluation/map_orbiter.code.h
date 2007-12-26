/***************************************************************************
 *            map_orbiter.code.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 
#include "map_orbiter.h"

#include <iosfwd>
#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "base/tuple.h"

#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "combinatoric/lattice_set.h"

#include "geometry/box.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"
#include "geometry/grid_approximation.h"
#include "geometry/rectangular_set.h"
#include "geometry/orbit.h"

#include "system/grid_multimap.h"


#include "system/map.h"
#include "system/discrete_time_system.h"

#include "evaluation/standard_applicator.h"

#include "output/logging.h"

namespace Ariadne {
  
namespace Evaluation { static int& verbosity = applicator_verbosity; }






template<class BS>
Evaluation::MapOrbiter<BS>::MapOrbiter(const EvolutionParameters<R>& parameters, 
                                       const ApplicatorInterface<BS>& applicator, 
                                       const ApproximatorInterface<BS>& approximator)
  : _applicator(applicator.clone()),
    _approximator(approximator.clone())
{
}


template<class BS>
Evaluation::MapOrbiter<BS>::MapOrbiter(const MapOrbiter<BS>& orbiter)
  : _applicator(orbiter._applicator->clone()),
    _approximator(orbiter._approximator->clone())
{
}


template<class BS>
Evaluation::MapOrbiter<BS>*
Evaluation::MapOrbiter<BS>::clone() const 
{
  return new MapOrbiter<BS>(*this);
}


template<class BS>
typename BS::real_type
Evaluation::MapOrbiter<BS>::maximum_basic_set_radius() const 
{
  return 1.0;
}

template<class BS>
Geometry::Grid<typename BS::real_type>
Evaluation::MapOrbiter<BS>::grid(dimension_type d) const 
{
  return Geometry::Grid<R>(d,0.125);
}


template<class BS>
Geometry::GridCellListSet<typename BS::real_type> 
Evaluation::MapOrbiter<BS>::upper_evolve(const System::MapInterface<R>& f, 
                                         const Geometry::Box<R>& bx, 
                                         const Numeric::Integer& n) const
{
  const Geometry::Grid<R> grid=this->grid(bx.dimension());
  BS bs=this->_approximator->over_approximation(bx);
  for(Numeric::Integer i=0; i!=n; ++i) {
    bs=this->_applicator->apply(f,bs);
  }
  return this->_approximator->outer_approximation(bs,grid);
}

template<class BS>
Geometry::GridCellListSet<typename BS::real_type> 
Evaluation::MapOrbiter<BS>::upper_reach(const System::MapInterface<R>& f, 
                                        const Geometry::Box<R>& bx, 
                                        const Numeric::Integer& n) const
{
  const Geometry::Grid<R> grid=this->grid(bx.dimension());
  Geometry::GridCellListSet<R> gcls(grid);
  BS bs=this->_approximator->over_approximation(bx);
  gcls.adjoin(this->_approximator->outer_approximation(bs,grid));
  for(Numeric::Integer i=0; i!=n; ++i) {
    bs=this->_applicator->apply(f,bs);
    gcls.adjoin(this->_approximator->outer_approximation(bs,grid));
  }
  return gcls;
}


template<class BS>
Geometry::ListSet< Geometry::Box<typename BS::real_type> >
Evaluation::MapOrbiter<BS>::lower_evolve(const System::MapInterface<R>& f, 
                                             const Geometry::Box<R>& bx, 
                                             const Numeric::Integer& n) const
{
  Geometry::ListSet< Geometry::Box<R> > result;
  std::vector< std::pair<T, BS> > stack;
  T t=0;
  BS bs=this->_approximator->over_approximation(bx);
  stack.push_back(std::make_pair(t,bs));
  while(!stack.empty()) {
    make_lpair(t,bs)=stack.back();
    stack.pop_back();
    do {
      bs=this->_applicator->apply(f,bs);
      t+=1;
    } while(t!=n && bs.radius()<this->maximum_basic_set_radius());
    if(t==n) {
      result.adjoin(this->_approximator->bounding_box(bs));
    } else {
      std::pair<BS,BS> split=this->_approximator->subdivide(bs);
      stack.push_back(std::make_pair(t,split.first));
      stack.push_back(std::make_pair(t,split.second));
    }
  }
  return result;
}

template<class BS>
Geometry::ListSet< Geometry::Box<typename BS::real_type> >
Evaluation::MapOrbiter<BS>::lower_reach(const System::MapInterface<R>& f, 
                                            const Geometry::Box<R>& s, 
                                            const Numeric::Integer& n) const
{
  Geometry::ListSet< Geometry::Box<R> > ls;
  BS bs=this->_approximator->over_approximation(s);
  ls.adjoin(this->_approximator->bounding_box(bs));
  for(Numeric::Integer i=0; i!=n; ++i) {
    bs=this->_applicator->apply(f,bs);
    ls.adjoin(this->_approximator->bounding_box(bs));
  }
  return ls;
}





template<class BS>
Geometry::Orbit< Numeric::Integer, BS >
Evaluation::MapOrbiter<BS>::orbit(const System::MapInterface<R>& f, const BS& bs, const Numeric::Integer& n) const
{
  assert(n>=0);
  Geometry::Orbit<Numeric::Integer, BS >  orbit(bs);
  BS es=bs;
  for(Numeric::Integer i=0; i!=n; ++i) {
    es=this->_applicator->apply(f,es);
    orbit.push_back(i,es.bounding_box());
  }
  return orbit;
}


template<class BS>
Geometry::Orbit< Numeric::Integer, BS>*
Evaluation::MapOrbiter<BS>::orbit(const System::MapInterface<R>& f, const Geometry::Box<R>& bx, const Numeric::Integer& n) const
{
  BS bs=this->_approximator->over_approximation(bx);
  return new Geometry::Orbit<T,BS>(this->orbit(f,bs,n));
}






}
