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

#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../combinatoric/lattice_set.h"

#include "../geometry/rectangle.h"
#include "../geometry/zonotope.h"
#include "../geometry/list_set.h"
#include "../geometry/grid.h"
#include "../geometry/grid_set.h"
#include "../geometry/partition_tree_set.h"
#include "../geometry/grid_approximation.h"
#include "../geometry/rectangular_set.h"
#include "../geometry/orbit.h"

#include "../system/grid_multimap.h"


#include "../system/map.h"
#include "../system/discrete_time_system.h"

#include "../evaluation/applicator.h"

#include "../output/logging.h"

namespace Ariadne {
  
namespace Evaluation { static int& verbosity = applicator_verbosity; }



template<class BS> inline
BS
Evaluation::MapOrbiter<BS>::apply(const System::MapInterface<R>& f, const BS& bs) const
{
  return this->_applicator->apply(f,bs);
}




template<class BS>
Evaluation::MapOrbiter<BS>::MapOrbiter(const EvolutionParameters<R>& parameters, const ApplicatorInterface<BS>& applicator)
  : _applicator(applicator.clone())
{
}


template<class BS>
Evaluation::MapOrbiter<BS>::MapOrbiter(const MapOrbiter<BS>& orbiter)
  : _applicator(orbiter._applicator->clone())
{
}


template<class BS>
Evaluation::MapOrbiter<BS>*
Evaluation::MapOrbiter<BS>::clone() const 
{
  return new MapOrbiter<BS>(*this);
}


template<class BS>
Geometry::Rectangle<typename BS::real_type>
Evaluation::MapOrbiter<BS>::apply(const System::MapInterface<R>& f, const Geometry::Rectangle<R>& r) const 
{
  ARIADNE_LOG(4,"Rectangle MapOrbiter::apply(MapInterface,Rectangle)\n");
  return this->_applicator->apply(f,BS(r)).bounding_box();
}


template<class BS>
Geometry::GridCellListSet<typename BS::real_type>
Evaluation::MapOrbiter<BS>::apply(const System::MapInterface<R>& f, const Geometry::GridCell<R>& gc) const
{
  ARIADNE_LOG(4,"GridCellListSet MapOrbiter::apply(MapInterface,GridCell)\n");
  return fuzzy_outer_approximation(this->_applicator->apply(f,BS(gc)),gc.grid());
}


template<class BS>
Geometry::GridCellListSet<typename BS::real_type>
Evaluation::MapOrbiter<BS>::apply(const System::MapInterface<R>& f, const Geometry::GridCell<R>& gc, const Geometry::Grid<R>& g) const
{
  ARIADNE_LOG(4,"GridCellListSet MapOrbiter::apply(MapInterface,GridCell)\n");
  return fuzzy_outer_approximation(this->_applicator->apply(f,BS(gc)),g);
}



template<class BS>
Geometry::DiscreteTimeOrbit< Numeric::Integer, Geometry::Rectangle<typename BS::real_type> >
Evaluation::MapOrbiter<BS>::orbit(const System::MapInterface<R>& f, const Geometry::Rectangle<R>& r, const Numeric::Integer& n) const
{
  ARIADNE_LOG(4,"DiscreteTimeOrbit<Integer,Rectangle> MapOrbiter::orbit(MapInterface,RectangleRectangle)\n");
  assert(n>=0);
  Geometry::DiscreteTimeOrbit<Numeric::Integer, Geometry::Rectangle<R> > orbit(r);
  BS bs(r);
  for(Numeric::Integer i=0; i!=n; ++i) {
    bs=this->_applicator->apply(f,bs);
    orbit.push_back(i,bs.bounding_box());
  }
  return orbit;
}


template<class BS>
Geometry::DiscreteTimeOrbit< Numeric::Integer, Geometry::Rectangle<typename BS::real_type> >
Evaluation::MapOrbiter<BS>::orbit(const System::MapInterface<R>& f, const Geometry::Rectangle<R>& r, const Numeric::Integer& n, const R& mbsr) const
{
  ARIADNE_LOG(4,"DiscreteTimeOrbit<Integer,Rectangle> MapOrbiter::orbit(MapInterface,Rectangle,Integer,Float)\n");
  assert(n>=0);
  Geometry::DiscreteTimeOrbit<Numeric::Integer, Geometry::Rectangle<R> > orbit(r);
  BS bs(r);
  size_type i=0;
  while(i!=n && orbit.final().set().radius()<mbsr) {
    bs=this->_applicator->apply(f,bs);
    orbit.push_back(i,bs.bounding_box());
  }
  return orbit;
}


template<class BS>
Geometry::DiscreteTimeOrbit< Numeric::Integer, Geometry::GridCellListSet<typename BS::real_type> >
Evaluation::MapOrbiter<BS>::orbit(const System::MapInterface<R>& f, const Geometry::GridCell<R>& gc, const Numeric::Integer& n) const
{
  ARIADNE_LOG(4,"DiscreteTimeOrbit<Integer,GridCellListSet> MapOrbiter::orbit(MapInterface,GridCell,Integer)\n");
  assert(n>=0);
  Geometry::GridCellListSet<R> gcls(gc.grid());
  gcls.adjoin(gc);
  Geometry::DiscreteTimeOrbit<Numeric::Integer, Geometry::GridCellListSet<R> > orbit(gcls);
  const Geometry::Grid<R>& g(gc.grid());
  BS bs(gc);
  for(int i=0; i!=n; ++i) {
    bs=this->_applicator->apply(f,bs);
    gcls=fuzzy_outer_approximation(bs,g);
    orbit.push_back(i,gcls);
  }
  return orbit;
}





}
