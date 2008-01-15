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
#include "geometry/box_list_set.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"
#include "geometry/grid_approximation.h"


#include "system/map.h"

#include "evaluation/standard_applicator.h"
#include "evaluation/standard_approximator.h"
#include "evaluation/standard_subdivider.h"

#include "output/logging.h"

namespace Ariadne {
  
namespace Evaluation { 
static const double DEFAULT_MAXIMUM_BASIC_SET_RADIUS=0.25;
static const double DEFAULT_GRID_LENGTH=0.125;
}


template<class BS>
Evaluation::MapEvolver<BS>::MapEvolver(const EvolutionParameters<R>& parameters)
  : _parameters(parameters.clone()),
    _applicator(new StandardApplicator<R>),
    _approximator(new StandardApproximator<BS>),
    _subdivider(new StandardSubdivider<BS>)
{
}

template<class BS> 
Evaluation::MapEvolver<BS>::MapEvolver(const EvolutionParameters<R>& parameters,
                                       const ApplicatorInterface<BS>& applicator)
  : _parameters(parameters.clone()),
    _applicator(applicator.clone()),
    _approximator(new StandardApproximator<BS>),
    _subdivider(new StandardSubdivider<BS>)
{
}

template<class BS> 
Evaluation::MapEvolver<BS>::MapEvolver(const EvolutionParameters<R>& parameters,
                                       const ApplicatorInterface<BS>& applicator,
                                       const ApproximatorInterface<BS>& approximator,
                                       const SubdividerInterface<BS>& subdivider)
  : _parameters(parameters.clone()),
    _applicator(applicator.clone()),
    _approximator(approximator.clone()),
    _subdivider(subdivider.clone())
{
}



template<class BS>
const Evaluation::EvolutionParameters<typename BS::real_type>&
Evaluation::MapEvolver<BS>::parameters() const
{
  return *this->_parameters;
}


template<class BS>
Evaluation::EvolutionParameters<typename BS::real_type>&
Evaluation::MapEvolver<BS>::parameters() 
{
  return *this->_parameters;
}







template<class BS>
void
Evaluation::MapEvolver<BS>::_step(BSL& evolve,
                                  BSL& reach, 
                                  TBSL& working,
                                  const Mp& map,
                                  const T& time,
                                  Semantics semantics) const
{
  uint verbosity=this->_parameters->verbosity();
  ARIADNE_LOG(5,"  working.size()="<<working.size()<<"\n");
  TBS tbs=working.pop();
  ARIADNE_LOG(5,"  tbs="<<tbs<<", r="<<this->radius(tbs)<<"\n");
  assert(tbs.time()<20);
  if(this->radius(tbs) > this->maximum_basic_set_radius()) {
    if(semantics==upper_semantics) {
      ARIADNE_LOG(5,"    subdivide...\n");
      //ARIADNE_LOG(7,"      "<<this->subdivide(tbs)<<"\n");
      working.adjoin(this->subdivide(tbs));
    } else if(semantics==lower_semantics) {
      ARIADNE_LOG(5,"    blocking...\n");
    }
  } else if(tbs.time()==time) {
      ARIADNE_LOG(5,"    end...\n");
    evolve.adjoin(tbs.set());
  } else {
    ARIADNE_LOG(5,"  apply... ");
    tbs=this->apply(map,tbs);
    ARIADNE_LOG(7,tbs);
    ARIADNE_LOG(5,"\n");
    reach.adjoin(tbs.set());
    working.adjoin(tbs);
  }
}
  


template<class BS>
Geometry::GridCellListSet<typename BS::real_type>
Evaluation::MapEvolver<BS>::_upper_evolve(const System::Map<R>& map, 
                                          const Geometry::GridCellListSet<R>& initial_set,
                                          const Numeric::Integer& time) const
{
  BSL reach, evolve;
  reach=this->basic_set_list(initial_set);
  TBSL working=this->timed_basic_set_list(initial_set);
  while(working.size()!=0) { 
    this->_step(evolve,reach,working,map,time,upper_semantics);
  }
  return this->outer_approximation(evolve,initial_set.grid());
}

template<class BS>
Geometry::GridCellListSet<typename BS::real_type>
Evaluation::MapEvolver<BS>::_upper_reach(const System::Map<R>& map, 
                                         const Geometry::GridCellListSet<R>& initial_set,
                                         const Numeric::Integer& time) const
{
  BSL reach, evolve;
  TBSL working=this->timed_basic_set_list(initial_set);
  while(working.size()!=0) { 
    this->_step(evolve,reach,working,map,time,upper_semantics);
  }
  return this->outer_approximation(reach,initial_set.grid());
}

    

template<class BS>
Geometry::SetInterface<typename BS::real_type>*
Evaluation::MapEvolver<BS>::lower_evolve(const System::Map<R>& map, 
                                         const Geometry::SetInterface<R>& initial_set,
                                         const Numeric::Interval<Numeric::Integer>& steps) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class BS>
Geometry::SetInterface<typename BS::real_type>*
Evaluation::MapEvolver<BS>::upper_evolve(const System::Map<R>& map, 
                                         const Geometry::SetInterface<R>& initial_set,
                                         const Numeric::Interval<Numeric::Integer>& steps) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class BS>
Geometry::SetInterface<typename BS::real_type>*
Evaluation::MapEvolver<BS>::lower_evolve(const System::Map<R>& map, 
                                         const Geometry::SetInterface<R>& initial_set,
                                         const Numeric::Integer& steps) const
{
  BSL reach, evolve;
  TBSL working=this->timed_basic_set_list(this->lower_approximation(initial_set));
  while(working.size()!=0) { 
    this->_step(evolve,reach,working,map,steps,lower_semantics);
  }
  return new GMS(this->outer_approximation(reach,this->grid(map.dimension())));
}

template<class BS>
Geometry::SetInterface<typename BS::real_type>*
Evaluation::MapEvolver<BS>::lower_reach(const System::Map<R>& map, 
                                       const Geometry::SetInterface<R>& initial_set,
                                       const Numeric::Integer& steps) const
{
  uint verbosity=this->_parameters->verbosity();
  BSL reach, evolve;
  ARIADNE_LOG(2,"MapEvolver::lower_reach(map,set,steps)\n");
  ARIADNE_LOG(3,"intial_set="<<initial_set<<"\n");
  BxLS initial_cells=this->lower_approximation(initial_set);
  ARIADNE_LOG(3,"intial_cells="<<initial_cells<<"\n");
  reach=this->basic_set_list(initial_cells);
  TBSL working=this->timed_basic_set_list(initial_cells);
  ARIADNE_LOG(3,"intial_working_sets="<<working<<"\n");
  while(working.size()!=0) { 
    this->_step(evolve,reach,working,map,steps,lower_semantics);
  }
  ARIADNE_LOG(3,"evolve="<<evolve<<"\n");
  ARIADNE_LOG(3,"reach="<<reach<<"\n");
  //GMS* result=new GMS(this->lower_approximation(evolve,this->grid(map.dimension())));
  Gr grid=this->grid(map.dimension());
  ARIADNE_LOG(3,"grid="<<grid<<"\n");
  GCLS approx=this->outer_approximation(reach,grid);
  ARIADNE_LOG(3,"approx="<<approx.summary()<<"\n");
  GMS* result=new GMS(approx);
  ARIADNE_LOG(3,"result="<<result->summary()<<"\n");
  return result;
}

template<class BS>
Geometry::SetInterface<typename BS::real_type>*
Evaluation::MapEvolver<BS>::upper_evolve(const System::Map<R>& map, 
                                         const Geometry::SetInterface<R>& initial_set,
                                         const Numeric::Integer& steps) const
{
  GCLS evolve=this->outer_approximation(initial_set);
  T quotient=quot(steps,this->lock_to_grid_steps());
  T remainder=rem(steps,this->lock_to_grid_steps());
  for(size_type i=0; i!=quotient; ++i) {
    evolve=this->_upper_evolve(map,evolve,this->lock_to_grid_steps());
  }
  evolve=this->_upper_evolve(map,evolve,remainder);
  return new GMS(evolve);
}


template<class BS>
Geometry::SetInterface<typename BS::real_type>*
Evaluation::MapEvolver<BS>::upper_reach(const System::Map<R>& map, 
                                       const Geometry::SetInterface<R>& initial_set,
                                       const Numeric::Integer& steps) const
{
  uint verbosity=this->_parameters->verbosity();
  ARIADNE_LOG(2,"MapEvolver::upper_reach(map,set,steps)\n");
  ARIADNE_LOG(3,"intial_set="<<initial_set<<"\n");
  GCLS evolve=this->outer_approximation(initial_set);
  GCLS reach=evolve;
  T quotient=quot(steps,this->lock_to_grid_steps());
  T remainder=rem(steps,this->lock_to_grid_steps());
  for(size_type i=0; i!=quotient; ++i) {
    evolve=this->_upper_evolve(map,evolve,this->lock_to_grid_steps());
    reach.adjoin(this->_upper_reach(map,evolve,this->lock_to_grid_steps()));
  }
  reach.adjoin(this->_upper_reach(map,evolve,remainder));
  return new GMS(reach);
}




template<class BS>
Geometry::SetInterface<typename BS::real_type>*
Evaluation::MapEvolver<BS>::chainreach(const System::Map<R>& map, 
                                       const Geometry::SetInterface<R>& initial_set) const
{
  Bx bb=this->bounding_domain(map);
  Gr grid=this->grid(map.dimension());
  T time=this->lock_to_grid_steps();
  GMS* result=new GMS(grid,bb);
  GB bounds=result->bounds();
  GCLS found=this->outer_approximation(initial_set);
  found=this->_upper_reach(map,found,time);
  while(!found.empty()) {
    result->adjoin(found);
    found=this->_upper_evolve(map,found,time);
    found.remove(*result);
    found.restrict(bounds);
  }
  return result;
}




template<class BS>
Geometry::SetInterface<typename BS::real_type>*
Evaluation::MapEvolver<BS>::viable(const System::Map<R>& map, 
                                   const Geometry::SetInterface<R>& bounding_set) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



template<class BS>
tribool
Evaluation::MapEvolver<BS>::verify(const System::Map<R>& map, 
                                   const Geometry::SetInterface<R>& initial_set, 
                                   const Geometry::SetInterface<R>& safe_set) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}









}
