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
#include "geometry/box_list_set.h"
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
  : _maximum_basic_set_radius(parameters.maximum_basic_set_radius()),
    _applicator(applicator.clone()),
    _approximator(approximator.clone())
{
}


template<class BS>
Evaluation::MapOrbiter<BS>::MapOrbiter(const MapOrbiter<BS>& orbiter)
  : _maximum_basic_set_radius(orbiter._maximum_basic_set_radius),
    _applicator(orbiter._applicator->clone()),
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
  return this->_maximum_basic_set_radius;
}



template<class BS>
Geometry::GridCellListSet<typename BS::real_type> 
Evaluation::MapOrbiter<BS>::upper_evolve(const System::Map<R>& f, 
                                         const Geometry::GridCell<R>& gc, 
                                         const Numeric::Integer& n) const
{
  const Geometry::Grid<R>& grid=gc.grid();
  BS bs=this->over_approximation(Geometry::Box<R>(gc));
  std::vector<BSL> orbit=this->_orbit(f,bs,n,upper_semantics);
  GCLS result=this->outer_approximation(orbit[n],grid);
  result.unique_sort();
  return result;
}


template<class BS>
Geometry::GridCellListSet<typename BS::real_type> 
Evaluation::MapOrbiter<BS>::upper_reach(const System::Map<R>& f, 
                                        const Geometry::GridCell<R>& gc, 
                                        const Numeric::Integer& n) const
{
  using namespace std;
  const Geometry::Grid<R>& grid=gc.grid();
  BS bs=this->over_approximation(Geometry::Box<R>(gc));
  std::vector<BSL> orbit=this->_orbit(f,bs,n,upper_semantics);
  GCLS result(grid);
  for(size_type i=0; i!=orbit.size(); ++i) {
    cout << "i="<<i<<endl;
    cout << "grid="<<grid<<endl;
    cout << "result="<<result<<endl;
    cout << "orbit[i]="<<orbit[i]<<endl;
    GCLS goa=Geometry::outer_approximation(orbit[i],grid);
    cout << "os="<<endl;
    GCLS oa(this->outer_approximation(orbit[i],grid));
    cout << "os="<<oa<<endl;
    result.adjoin(this->outer_approximation(orbit[i],grid));
  }
  result.unique_sort();
  return result;
}


template<class BS>
std::pair< Numeric::Integer, Geometry::Box<typename BS::real_type> >
Evaluation::MapOrbiter<BS>::lower_evolve(const System::Map<R>& f, 
                                         const Geometry::Box<R>& bx, 
                                         const Numeric::Integer& t) const
{
  BS bs=this->over_approximation(bx);
  Numeric::Integer n=0;
  if(bs.radius()>this->maximum_basic_set_radius()) {
    return std::make_pair(n,this->bounding_box(bs));
  }
  for(n=0; n!=t; ++n) {
    BS rs=this->apply(f,bs);
    if(rs.radius()>this->maximum_basic_set_radius()) {
      break;
    }
    bs=rs;
  }
  return std::make_pair(n,this->bounding_box(bs));
}


template<class BS>
Geometry::BoxListSet<typename BS::real_type>
Evaluation::MapOrbiter<BS>::lower_reach(const System::Map<R>& f, 
                                        const Geometry::Box<R>& bx, 
                                        const Numeric::Integer& t) const
{
  Geometry::BoxListSet<R> result;
  BS bs=this->over_approximation(bx);
  for(uint n=0; n<=t; ++n) {
    bs=this->apply(f,bs);
    if(bs.radius()>this->maximum_basic_set_radius()) {
      break;
    }
    result.push_back(this->bounding_box(bs));
  }
  return result;
} 




template<class BS>
Geometry::Orbit< Numeric::Integer, BS >
Evaluation::MapOrbiter<BS>::orbit(const System::Map<R>& f, const BS& bs, const Numeric::Integer& n) const
{
  assert(n>=0);
  Geometry::Orbit<Numeric::Integer, BS >  orbit(bs);
  BS es=bs;
  for(Numeric::Integer i=0; i!=n; ++i) {
    es=this->apply(f,es);
    orbit.push_back(i,es.bounding_box());
  }
  return orbit;
}


template<class BS>
Geometry::Orbit< Numeric::Integer, BS>*
Evaluation::MapOrbiter<BS>::orbit(const System::Map<R>& f, const Geometry::Box<R>& bx, const Numeric::Integer& n) const
{
  BS bs=this->over_approximation(bx);
  return new Geometry::Orbit<T,BS>(this->orbit(f,bs,n));
}


template<class BS>
std::ostream&
Evaluation::MapOrbiter<BS>::write(std::ostream& os) const
{
  return os << "MapOrbiter<"<<Geometry::name<BS>()<<">"
            << "( maximum_basic_set_radius="<<this->maximum_basic_set_radius()<<" )";
}






}
