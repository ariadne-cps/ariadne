/***************************************************************************
 *            reachability_analyser.code.h
 *
 *  Copyright  2006-8  Alberto Casagrande, Pieter Collins
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
 
#include "reachability_analyser.h"

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

#include "evaluation/evolution_parameters.h"
#include "evaluation/standard_applicator.h"
#include "evaluation/standard_approximator.h"
#include "evaluation/standard_subdivider.h"

#include "output/logging.h"

namespace {
using namespace Ariadne;
template<class R> void compute_lock_to_grid_time(Integer& t, const EvolutionParameters<R>& parameters) {
  t=parameters.lock_to_grid_steps(); }
template<class R> void compute_lock_to_grid_time(Rational& t, const EvolutionParameters<R>& parameters) {
  t=parameters.lock_to_grid_time(); }
}



namespace Ariadne {
  
 
static const double DEFAULT_MAXIMUM_ENCLOSURE_RADIUS=0.25;
static const double DEFAULT_GRID_LENGTH=0.125;



template<class T, class Aprx>
ReachabilityAnalyser<TransitionSystemInterface<T,Aprx>,Aprx>::
ReachabilityAnalyser(const EvolutionParameters<R>& parameters)
  : _parameters(parameters.clone())
{
}




template<class T, class Aprx>
const EvolutionParameters<typename Aprx::real_type>&
ReachabilityAnalyser<TransitionSystemInterface<T,Aprx>,Aprx>::parameters() const
{
  return *this->_parameters;
}


template<class T, class Aprx>
EvolutionParameters<typename Aprx::real_type>&
ReachabilityAnalyser<TransitionSystemInterface<T,Aprx>,Aprx>::parameters() 
{
  return *this->_parameters;
}



template<class T, class Aprx>
T
ReachabilityAnalyser<TransitionSystemInterface<T,Aprx>,Aprx>::lock_to_grid_time() const
{
  T t; compute_lock_to_grid_time(t,*this->_parameters); return t;
}









template<class T, class Aprx>
typename ReachabilityAnalyser<TransitionSystemInterface<T,Aprx>,Aprx>::SetPointer
ReachabilityAnalyser<TransitionSystemInterface<T,Aprx>,Aprx>::
lower_evolve(const System& system, 
             const Set& initial_set,
             const Time& time) const
{
  BxLS initial=this->_lower_approximation(initial_set);
  BxLS final;
  for(typename BxLS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
    final.adjoin(this->_lower_evolve(system,*bs_iter,time));
  }
  return new GMS(this->_outer_approximation(final));
}


template<class T, class Aprx>
typename ReachabilityAnalyser<TransitionSystemInterface<T,Aprx>,Aprx>::SetPointer
ReachabilityAnalyser<TransitionSystemInterface<T,Aprx>,Aprx>::
lower_reach(const System& system, 
            const Set& initial_set,
            const Time& time) const
{
  BxLS initial=this->_lower_approximation(initial_set);
  BxLS reach;
  for(typename BxLS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
    reach.adjoin(this->_lower_reach(system,*bs_iter,time));
  }
  return new GMS(this->_outer_approximation(reach));
}


template<class T, class Aprx>
typename ReachabilityAnalyser<TransitionSystemInterface<T,Aprx>,Aprx>::SetPointer
ReachabilityAnalyser<TransitionSystemInterface<T,Aprx>,Aprx>::
upper_evolve(const System& system, 
             const Set& initial_set,
             const Time& time) const
{
  GCLS evolve=this->_outer_approximation(initial_set);
  Time lock_to_grid_time=this->lock_to_grid_time();
  Integer time_steps=quot(time,lock_to_grid_time);
  Time remainder_time=rem(time,lock_to_grid_time);
  for(size_type i=0; i!=time_steps; ++i) {
    evolve=this->_upper_evolve(system,evolve,lock_to_grid_time);
  }
  evolve=this->_upper_evolve(system,evolve,remainder_time);
  return new GMS(evolve);
}


template<class T, class Aprx>
typename ReachabilityAnalyser<TransitionSystemInterface<T,Aprx>,Aprx>::SetPointer
ReachabilityAnalyser<TransitionSystemInterface<T,Aprx>,Aprx>::
upper_reach(const System& system, 
            const Set& initial_set,
            const Time& time) const
{
  uint verbosity=this->_parameters->verbosity();
  ARIADNE_LOG(2,"ReachabilityAnalyser::upper_reach(system,set,time)\n");
  ARIADNE_LOG(3,"intial_set="<<initial_set<<"\n");
  GCLS evolve=this->_outer_approximation(initial_set);
  GCLS reach=evolve;
  Time lock_to_grid_time = this->lock_to_grid_time();
  Integer time_steps=quot(time,lock_to_grid_time);
  Time remainder_time=rem(time,lock_to_grid_time);
  for(size_type i=0; i!=time_steps; ++i) {
    evolve=this->_upper_evolve(system,evolve,lock_to_grid_time);
    reach.adjoin(this->_upper_reach(system,evolve,lock_to_grid_time));
  }
  reach.adjoin(this->_upper_reach(system,evolve,remainder_time));
  return new GMS(reach);
}




template<class T, class Aprx>
typename ReachabilityAnalyser<TransitionSystemInterface<T,Aprx>,Aprx>::SetPointer
ReachabilityAnalyser<TransitionSystemInterface<T,Aprx>,Aprx>::
chain_reach(const System& system, 
            const Set& initial_set) const
{
  Bx bb=this->bounding_domain(system);
  Gr grid=this->grid(system.state_space());
  T time=this->lock_to_grid_time();
  GMS* result=new GMS(grid,bb);
  GB bounds=result->bounds();
  GCLS found=this->_outer_approximation(initial_set);
  found=this->_upper_reach(system,found,time);
  while(!found.empty()) {
    result->adjoin(found);
    found=this->_upper_evolve(system,found,time);
    found.remove(*result);
    found.restrict(bounds);
  }
  return result;
}




template<class T, class Aprx>
typename ReachabilityAnalyser<TransitionSystemInterface<T,Aprx>,Aprx>::SetPointer
ReachabilityAnalyser<TransitionSystemInterface<T,Aprx>,Aprx>::
viable(const System& system, 
       const Set& bounding_set) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class T, class Aprx>
tribool
ReachabilityAnalyser<TransitionSystemInterface<T,Aprx>,Aprx>::
verify(const System& system, 
       const Set& initial_set, 
       const Set& safe_set) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}












} // namespace Ariadne
