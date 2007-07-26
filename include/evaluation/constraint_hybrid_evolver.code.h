/***************************************************************************
 *            constraint_hybrid_evolver.code.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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
 
#include <vector>

#include "../geometry/rectangle_expression.h"
#include "../geometry/set_interface.h"
#include "../geometry/hybrid_set.h"
#include "../geometry/timed_set.h"
#include "../system/map.h"
#include "../system/vector_field.h"
#include "../system/constraint_hybrid_automaton.h"
#include "../evaluation/applicator.h"
#include "../evaluation/integrator.h"

#include "../evaluation/lohner_integrator.h"
#include "../evaluation/constraint_hybrid_evolver_plugin.h"

#include "../output/epsfstream.h"
#include "../output/logging.h"

#include "constraint_hybrid_evolver.h"


namespace {

using namespace Ariadne;

template<class V, class Iter> inline
void append(V& v, Iter begin, Iter end) 
{
  for(;begin!=end;++begin) {
    v.push_back(*begin);
  }
}

template<class V, class C> inline
void append(V& v, const C& c) 
{
  for(typename C::const_iterator iter=c.begin(); iter!=c.end(); ++iter) {
    v.push_back(*iter);
  }
}


} 



namespace Ariadne {
 
namespace Evaluation { static int& verbosity = hybrid_evolver_verbosity; }

using namespace Numeric;
using namespace LinearAlgebra;
using namespace Geometry;
using namespace System;

template<class R>
Evaluation::ConstraintHybridEvolver<R>::~ConstraintHybridEvolver()
{
  delete this->_plugin;
}


template<class R>
Evaluation::ConstraintHybridEvolver<R>::ConstraintHybridEvolver(Applicator<R>& a, Integrator<R>& i)
  : _plugin(new ConstraintHybridEvolverPlugin<R>(a,i))
{
}

template<class R>
Evaluation::ConstraintHybridEvolver<R>::ConstraintHybridEvolver(const ConstraintHybridEvolver<R>& evolver)
  :_plugin(evolver._plugin)
{
}



template<class R>
const std::vector<typename Evaluation::ConstraintHybridEvolver<R>::timed_set_type>&
Evaluation::ConstraintHybridEvolver<R>::trace() const
{
  return this->_trace;
}


template<class R>
time_type
Evaluation::ConstraintHybridEvolver<R>::maximum_step_size() const
{
  return this->_plugin->_integrator->maximum_step_size();
}


template<class R>
R
Evaluation::ConstraintHybridEvolver<R>::maximum_basic_set_radius() const
{
  return this->_plugin->_integrator->maximum_basic_set_radius();
}


template<class R>
time_type
Evaluation::ConstraintHybridEvolver<R>::lock_to_grid_time() const
{
  return this->_plugin->_integrator->lock_to_grid_time();
}



template<class R> 
typename Evaluation::ConstraintHybridEvolver<R>::working_sets_type
Evaluation::ConstraintHybridEvolver<R>::_compute_working_sets(const hybrid_list_set_type& set) const
{
  working_sets_type working_sets;
  for(typename hybrid_list_set_type::const_iterator loc_iter=set.begin();
      loc_iter!=set.end(); ++loc_iter)
  {
    id_type loc_id=loc_iter->discrete_state();
    const continuous_basic_set_type& bs=loc_iter->continuous_state_set();
    working_sets.push_back(timed_set_type(0,0,loc_id,bs));
  }
  return working_sets;
}


template<class R> inline
typename Evaluation::ConstraintHybridEvolver<R>::hybrid_list_set_type
Evaluation::ConstraintHybridEvolver<R>::_compute_list_set(const working_sets_type& working_sets, const HybridSpace& locations) const
{
  hybrid_list_set_type result_set(locations);
  for(typename working_sets_type::const_iterator iter=working_sets.begin();
      iter!=working_sets.end(); ++iter)
  {
    result_set.adjoin(hybrid_basic_set_type(iter->discrete_state(),iter->continuous_state_set()));
  }

    return result_set;
}



template<class R>
Geometry::HybridListSet<typename Evaluation::ConstraintHybridEvolver<R>::continuous_basic_set_type> 
Evaluation::ConstraintHybridEvolver<R>::upper_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                     const Geometry::HybridListSet<continuous_basic_set_type>& initial_set, 
                                                     time_type evolution_time,
                                                     size_type maximum_number_of_events) const
{
  verbosity=8;
  ARIADNE_LOG(2,"HybridEvolver::upper_evolve(HybridAutomaton automaton, ListSet initial_set, Time time, Integer maximum_number_of_events)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<", evolution_time="<<evolution_time<<"\n");
  
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  std::vector< timed_set_type > working_sets=this->_compute_working_sets(initial_set);
  std::vector< timed_set_type > final_sets;
  std::vector< timed_set_type >& trace_sets=this->_trace;
  trace_sets=working_sets;

  while(!working_sets.empty()) {
    ARIADNE_LOG(9,"working_sets.size()="<<working_sets.size()<<"\n");
    timed_set_type working_set=working_sets.back(); working_sets.pop_back();
    ARIADNE_LOG(9,"working_set="<<working_set<<"\n");
    assert(working_set.time()<=evolution_time);
    
    if(working_set.time()==evolution_time) {
      ARIADNE_LOG(9,"  reached evolution time\n");
      final_sets.push_back(working_set);
    } else if(working_set.continuous_state_set().radius()>this->maximum_basic_set_radius()) {
      ::append(working_sets,this->_plugin->subdivide(working_set));
    } else {
      ::append(working_sets,this->_plugin->upper_evolution_step(automaton,working_set,evolution_time));
      ::append(this->_trace,this->_plugin->upper_evolution_step(automaton,working_set,evolution_time));
    }
  }

  return this->_compute_list_set(final_sets,automaton.locations());
}


template<class R>
Geometry::HybridListSet<typename Evaluation::ConstraintHybridEvolver<R>::continuous_basic_set_type> 
Evaluation::ConstraintHybridEvolver<R>::upper_reach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridListSet<continuous_basic_set_type>& initial_set, 
                                                    time_type evolution_time, 
                                                    size_type maximum_number_of_events) const
{
  ARIADNE_LOG(2,"HybridEvolver::upper_evolve(HybridAutomaton automaton, ListSet initial_set, Time time, Integer maximum_number_of_events)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<", evolution_time="<<evolution_time<<"\n");
  ARIADNE_LOG(7,"maximum_step_size="<<this->maximum_step_size()<<", maximum_basic_set_radius="<<this->maximum_basic_set_radius()<<"\n");

  
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  std::vector< timed_set_type > working_sets=this->_compute_working_sets(initial_set);
  std::vector< timed_set_type > reached_sets;

  while(!working_sets.empty()) {
    timed_set_type working_set=working_sets.back(); working_sets.pop_back();
    if(working_set.time()==evolution_time) {
    } else if(working_set.continuous_state_set().radius()>this->maximum_basic_set_radius()) {
      ::append(working_sets,this->_plugin->subdivide(working_set));
    } else {
      ::append(working_sets,this->_plugin->upper_evolution_step(automaton,working_set,evolution_time));
      ::append(reached_sets,this->_plugin->upper_reachability_step(automaton,working_set,evolution_time));
      if(!working_sets.empty()) { std::cout << "evolution step: " << working_sets.back() << std::endl; }
      if(!reached_sets.empty()) { std::cout << "reach step: " << reached_sets.back() << std::endl; }
    }
  }

  return this->_compute_list_set(reached_sets,automaton.locations());
}




template<class R>
Geometry::HybridGridCellListSet<R>
Evaluation::ConstraintHybridEvolver<R>::upper_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                     const Geometry::HybridGridCell<R>& initial_set, 
                                                     const Geometry::HybridGrid<R>& hybrid_grid, 
                                                     time_type evolution_time, 
                                                     size_type maximum_number_of_events) const
{
  HybridGridCellListSet<R> final_set(hybrid_grid);
  
  HybridListSet< Zonotope<I> > list_set(automaton.locations());
  list_set.adjoin(initial_set);
  
  list_set=this->upper_evolve(automaton,list_set,evolution_time,maximum_number_of_events);
  
  final_set.adjoin_outer_approximation(list_set);

  return final_set;
}


template<class R>
Geometry::HybridGridCellListSet<R>
Evaluation::ConstraintHybridEvolver<R>::upper_reach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridGridCell<R>& initial_set, 
                                                    const Geometry::HybridGrid<R>& hybrid_grid, 
                                                    time_type evolution_time, 
                                                    size_type maximum_number_of_events) const
{
  HybridGridCellListSet<R> final_set(hybrid_grid);
  
  HybridListSet< Zonotope<I> > list_set(automaton.locations());
  list_set.adjoin(initial_set);
  
  list_set=this->upper_reach(automaton,list_set,evolution_time,maximum_number_of_events);
  
  final_set.adjoin_outer_approximation(list_set);

  return final_set;
}



template<class R>
Geometry::HybridGridMaskSet<R>
Evaluation::ConstraintHybridEvolver<R>::upper_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                     const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                     time_type evolution_time, 
                                                     size_type maximum_number_of_events) const
{
  Geometry::HybridGrid<R> hybrid_grid=initial_set.grid();
  Geometry::HybridGridMaskSet<R> current_set=initial_set;
  Geometry::HybridGridMaskSet<R> next_set=initial_set;
  next_set.clear();

  Geometry::HybridGridMultiMap<R> map(initial_set.grid(),initial_set.grid());
  
  time_type lock_to_grid_time = this->lock_to_grid_time();
  size_type evolution_steps = int_up<int>(time_type(evolution_time/lock_to_grid_time));
  lock_to_grid_time = evolution_time/evolution_steps;

  for(uint i=0; i!=evolution_steps; ++i) {
    for(typename HybridGridMaskSet<R>::const_iterator cell_iter=current_set.begin();
        cell_iter!=current_set.end(); ++cell_iter)
    {
      //const HybridBasicSet< GridCell<R> >& cell=*cell_iter;
      HybridGridCell<R> cell=*cell_iter;
      if(!map.has_key(cell)) {
        map.set_image(cell,this->upper_evolve(automaton,cell,hybrid_grid,lock_to_grid_time,maximum_number_of_events));
      }
      next_set.adjoin(map.image(cell));
    }
    current_set=next_set;
    next_set.clear();
  }

  return current_set;
}



template<class R>
Geometry::HybridGridMaskSet<R>
Evaluation::ConstraintHybridEvolver<R>::upper_reach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                    time_type evolution_time, 
                                                    size_type maximum_number_of_events) const
{
  Geometry::HybridGrid<R> hybrid_grid=initial_set.grid();
  Geometry::HybridGridMaskSet<R> current_set=initial_set;
  Geometry::HybridGridMaskSet<R> next_set(current_set);
  next_set.clear();

  Geometry::HybridGridMultiMap<R> map(initial_set.grid(),initial_set.grid());

  size_type evolution_steps = int_up<int>(time_type(evolution_time/this->lock_to_grid_time()));
  time_type lock_to_grid_time = evolution_time/evolution_steps;

  for(typename HybridGridMaskSet<R>::const_iterator cell_iter=current_set.begin();
      cell_iter!=current_set.end(); ++cell_iter)
  {
    const HybridGridCell<R>& cell=*cell_iter;
    next_set.adjoin(this->upper_reach(automaton,cell,hybrid_grid,lock_to_grid_time,maximum_number_of_events));
  }
  current_set=next_set;
  next_set.clear();

  for(uint i=0; i!=evolution_steps-1; ++i) {
    for(typename HybridGridMaskSet<R>::const_iterator cell_iter=current_set.begin();
        cell_iter!=current_set.end(); ++cell_iter)
    {
      const HybridGridCell<R>& cell=*cell_iter;
      if(!map.has_key(cell)) {
        map.set_image(cell,this->upper_evolve(automaton,cell,hybrid_grid,lock_to_grid_time,maximum_number_of_events));
      }
      next_set.adjoin(map.image(cell));
    }
    current_set.adjoin(next_set);
    next_set.clear();
  }

  return current_set;
}


template<class R>
Geometry::HybridGridMaskSet<R>
Evaluation::ConstraintHybridEvolver<R>::viable(const System::ConstraintHybridAutomaton<R>& automaton, 
                                               const Geometry::HybridGridMaskSet<R>& bounding_set) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
Geometry::HybridGridMaskSet<R>
Evaluation::ConstraintHybridEvolver<R>::chainreach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                   const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                   const Geometry::HybridGridMaskSet<R>& bounding_set) const
{
  Geometry::HybridGrid<R> hybrid_grid=initial_set.grid();
  Geometry::HybridGridMaskSet<R> reach_set=initial_set;
  Geometry::HybridGridMaskSet<R> current_set=initial_set;
  Geometry::HybridGridMaskSet<R> found_set=initial_set;

  time_type lock_to_grid_time = this->lock_to_grid_time();
  size_type maximum_number_of_events = (size_type) -1;
  
  for(typename HybridGridMaskSet<R>::const_iterator cell_iter=reach_set.begin();
      cell_iter!=current_set.end(); ++cell_iter)
  {
    const HybridGridCell<R>& cell=*cell_iter;
    found_set.adjoin(this->upper_reach(automaton,cell,hybrid_grid,lock_to_grid_time,maximum_number_of_events));
  }

  reach_set=found_set;
  found_set.remove(reach_set);
  while(!found_set.empty()) {
    current_set=found_set;
    found_set.clear();

    for(typename HybridGridMaskSet<R>::const_iterator cell_iter=current_set.begin();
        cell_iter!=current_set.end(); ++cell_iter)
    {
      const HybridGridCell<R>& cell=*cell_iter;
      found_set.adjoin(this->upper_evolve(automaton,cell,hybrid_grid,lock_to_grid_time,maximum_number_of_events));
    }

    found_set.remove(reach_set);
    reach_set.adjoin(found_set);
  }

  return reach_set;
}



} // namespace Ariadne
