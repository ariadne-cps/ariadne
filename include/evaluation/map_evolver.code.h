/***************************************************************************
 *            map_evolver.code.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
#include "map_evolver.h"
#include "map_orbiter.h"

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
#include "geometry/list_set.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"
#include "geometry/grid_approximation.h"
#include "geometry/orbit.h"

#include "system/grid_multimap.h"


#include "system/map_interface.h"
#include "system/transition_system.h"

#include "evaluation/standard_applicator.h"
#include "evaluation/standard_approximator.h"
#include "evaluation/model_checker.h"

#include "output/logging.h"

namespace Ariadne {
  
namespace Evaluation { 
static int& verbosity = applicator_verbosity; 
static const double DEFAULT_MAXIMUM_BASIC_SET_RADIUS=0.25;
static const double DEFAULT_GRID_LENGTH=0.125;
}



template<class R> 
Evaluation::EvolutionParameters<R>*
Evaluation::MapEvolver<R>::default_parameters() 
{
  return new EvolutionParameters<R>();
}

template<class R> 
Evaluation::MapOrbiterInterface<R>*
Evaluation::MapEvolver<R>::default_orbiter() 
{
  typedef Geometry::Zonotope<R,Geometry::UniformErrorTag> BS;
  const EvolutionParameters<R>& parameters=*this->_parameters;
  StandardApplicator<R> applicator;
  StandardApproximator<BS> approximator;
  return new MapOrbiter<BS>(parameters,applicator,approximator);
}

template<class R> 
Evaluation::ModelChecker<R>
Evaluation::MapEvolver<R>::model_checker() const
{
  return ModelChecker<R>();
}

template<class R> 
System::TransitionSystem<R>
Evaluation::MapEvolver<R>::discrete_map(const System::MapInterface<R>& map) const
{
  return System::TransitionSystem<R>(map,*this->_orbiter,1u);
}


template<class R> 
Geometry::GridMaskSet<R>
Evaluation::MapEvolver<R>::outer_approximation(const Geometry::SetInterface<R>& set) const
{
  Geometry::FiniteGrid<R> grid=this->_parameters->finite_grid(set.dimension());
  return Geometry::outer_approximation(set,grid);
}

template<class R> 
Geometry::GridMaskSet<R>
Evaluation::MapEvolver<R>::inner_approximation(const Geometry::SetInterface<R>& set) const
{
  Geometry::FiniteGrid<R> grid(this->_parameters->finite_grid(set.dimension()));
  return Geometry::inner_approximation(set,grid);
}

template<class R> 
Geometry::ListSet< Geometry::Box<R> >
Evaluation::MapEvolver<R>::lower_approximation(const Geometry::SetInterface<R>& set) const
{
  return Geometry::point_approximation(set,this->_parameters->grid(set.dimension()));
}

template<class R>
Evaluation::MapEvolver<R>::MapEvolver() 
  : _parameters(default_parameters()),
    _orbiter(default_orbiter())
{
  _parameters->set_maximum_basic_set_radius(DEFAULT_MAXIMUM_BASIC_SET_RADIUS);
  _parameters->set_grid_length(DEFAULT_GRID_LENGTH);
}


template<class R>
Evaluation::MapEvolver<R>::MapEvolver(const EvolutionParameters<R>& parameters)
  : _parameters(new EvolutionParameters<R>(parameters)),
    _orbiter(default_orbiter())
{
}



template<class R>
Evaluation::MapEvolver<R>::MapEvolver(const MapEvolver<R>& other) 
  : _parameters(new EvolutionParameters<R>(other.parameters())),
    _orbiter(other._orbiter->clone())
{
}



template<class R>
Evaluation::MapEvolver<R>::~MapEvolver() 
{
  delete this->_parameters;
  delete this->_orbiter;
}


template<class R>
Evaluation::MapEvolver<R>*
Evaluation::MapEvolver<R>::clone() const 
{
  return new MapEvolver<R>(*this);
}




template<class R>
const Evaluation::EvolutionParameters<R>&
Evaluation::MapEvolver<R>::parameters() const
{
  return *this->_parameters;
}


template<class R>
Evaluation::EvolutionParameters<R>&
Evaluation::MapEvolver<R>::parameters() 
{
  return *this->_parameters;
}









template<class R>
Geometry::OrbitInterface< Numeric::Integer >*
Evaluation::MapEvolver<R>::orbit(const System::MapInterface<R>& map, const Geometry::Box<R>& set, const Numeric::Integer& time) const
{
  ARIADNE_LOG(4,"Orbit<Integer,Box> MapEvolver::orbit(MapInterface,BoxBox)\n");
  return this->_orbiter->orbit(map,set,time);
}



template<class R>
Geometry::SetInterface<R>*
Evaluation::MapEvolver<R>::lower_evolve(const System::MapInterface<R>& map, 
                                        const Geometry::SetInterface<R>& initial_set,
                                        const Numeric::Integer& steps) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
  //return this->model_checker().evolve(this->discrete_map(map),lower_approximation(initial_set),steps);
}

template<class R>
Geometry::SetInterface<R>*
Evaluation::MapEvolver<R>::lower_reach(const System::MapInterface<R>& map, 
                                       const Geometry::SetInterface<R>& initial_set,
                                       const Numeric::Integer& steps) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
  //return this->model_checker().lower_reach(this->discrete_map(map),lower_approximation(initial_set),steps);
}

template<class R>
Geometry::SetInterface<R>*
Evaluation::MapEvolver<R>::upper_evolve(const System::MapInterface<R>& map, 
                                        const Geometry::SetInterface<R>& initial_set,
                                        const Numeric::Integer& steps) const
{
  return this->model_checker().evolve(this->discrete_map(map),outer_approximation(initial_set),steps).clone();
}

template<class R>
Geometry::SetInterface<R>*
Evaluation::MapEvolver<R>::upper_reach(const System::MapInterface<R>& map, 
                                       const Geometry::SetInterface<R>& initial_set,
                                       const Numeric::Integer& steps) const
{
  return this->model_checker().reach(this->discrete_map(map),outer_approximation(initial_set),steps).clone();
}



template<class R>
Geometry::SetInterface<R>*
Evaluation::MapEvolver<R>::chainreach(const System::MapInterface<R>& map, 
                                      const Geometry::SetInterface<R>& initial_set) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
Geometry::SetInterface<R>*
Evaluation::MapEvolver<R>::chainreach(const System::MapInterface<R>& map, 
                                      const Geometry::SetInterface<R>& initial_set, 
                                      const Geometry::Box<R>& bounding_box) const
{
  return this->model_checker().chainreach(this->discrete_map(map),this->outer_approximation(initial_set),bounding_box).clone();
}




template<class R>
Geometry::SetInterface<R>*
Evaluation::MapEvolver<R>::viable(const System::MapInterface<R>& map, 
                                  const Geometry::SetInterface<R>& bounding_set) const
{
  return this->model_checker().viable(this->discrete_map(map),this->inner_approximation(bounding_set)).clone();
}



template<class R>
tribool
Evaluation::MapEvolver<R>::verify(const System::MapInterface<R>& map, 
                                  const Geometry::SetInterface<R>& initial_set, 
                                  const Geometry::SetInterface<R>& safe_set) const
{
  return this->model_checker().verify(this->discrete_map(map),this->outer_approximation(initial_set),this->inner_approximation(safe_set));
}




template<class R>
System::GridMultiMap<R> 
Evaluation::MapEvolver<R>::discretize(const System::MapInterface<R>& f, 
                                      const Geometry::GridMaskSet<R>& domain,
                                      const Geometry::Grid<R>& range_grid) const
{
  ARIADNE_LOG(2,"GridMultiMap*Evaluation::MapEvolver::discretize(MapInterface map, GridMaskSet domain, Grid range_grid)\n");
  ARIADNE_LOG(3,"domain="<<domain<<"\n"<<"range_grid="<<range_grid);
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  System::GridMultiMap<R> result(domain.grid(),range_grid);
  for(typename GridMaskSet<R>::const_iterator dom_iter=domain.begin();
      dom_iter!=domain.end(); ++dom_iter)
    {
      const GridCell<R>& gc=*dom_iter;
      GridCellListSet<R> gcls=this->_orbiter->upper_evolve(f,gc,1);
      result.adjoin_to_image(gc,gcls);
    }
  return result;
}





}
