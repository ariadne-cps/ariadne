/***************************************************************************
 *            constraint_based_hybrid_evolver.code.h
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

#include "geometry/rectangle_expression.h"
#include "geometry/set_interface.h"
#include "geometry/hybrid_set.h"
#include "geometry/timed_set.h"
#include "system/map.h"
#include "system/vector_field.h"
#include "system/constraint_based_hybrid_automaton.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/applicator_interface.h"
#include "evaluation/integrator_interface.h"
#include "evaluation/detector_interface.h"
#include "evaluation/applicator.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/detector.h"

#include "evaluation/constraint_based_hybrid_scheduler.h"

#include "output/epsstream.h"
#include "output/logging.h"

#include "constraint_based_hybrid_evolver.h"


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
Evaluation::ConstraintBasedHybridEvolver<R>::~ConstraintBasedHybridEvolver()
{
  delete this->_parameters;
  delete this->_scheduler;
}


template<class R>
Evaluation::ConstraintBasedHybridEvolver<R>::ConstraintBasedHybridEvolver(const EvolutionParameters<R>& p)
  : _parameters(new EvolutionParameters<R>(p)),
    _scheduler(new ConstraintBasedHybridScheduler<R>(Applicator<R>(),LohnerIntegrator<R>(),Detector<R>()))
{
}

template<class R>
Evaluation::ConstraintBasedHybridEvolver<R>::ConstraintBasedHybridEvolver(const EvolutionParameters<R>& p, const ApplicatorInterface<BS>& a, const IntegratorInterface<BS>& i)
  : _parameters(new EvolutionParameters<R>(p)),
    _scheduler(new ConstraintBasedHybridScheduler<R>(a,i,Detector<R>()))
{
}

template<class R>
Evaluation::ConstraintBasedHybridEvolver<R>::ConstraintBasedHybridEvolver(const EvolutionParameters<R>& p, const ApplicatorInterface<BS>& a, const IntegratorInterface<BS>& i, const DetectorInterface<R>& d)
  : _parameters(new EvolutionParameters<R>(p)),
    _scheduler(new ConstraintBasedHybridScheduler<R>(a,i,d))
{
}

template<class R>
Evaluation::ConstraintBasedHybridEvolver<R>::ConstraintBasedHybridEvolver(const ConstraintBasedHybridEvolver<R>& evolver)
  : _parameters(new EvolutionParameters<R>(*evolver._parameters)),
    _scheduler(new ConstraintBasedHybridScheduler<R>(*evolver._scheduler))
{
}



template<class R>
const std::vector<typename Evaluation::ConstraintBasedHybridEvolver<R>::timed_set_type>&
Evaluation::ConstraintBasedHybridEvolver<R>::trace() const
{
  return this->_trace;
}


template<class R>
Evaluation::EvolutionParameters<R>&
Evaluation::ConstraintBasedHybridEvolver<R>::parameters() 
{
  return *this->_parameters;
}


template<class R>
const Evaluation::EvolutionParameters<R>&
Evaluation::ConstraintBasedHybridEvolver<R>::parameters() const
{
  return *this->_parameters;
}


template<class R>
Numeric::Rational
Evaluation::ConstraintBasedHybridEvolver<R>::maximum_step_size() const
{
  return this->_parameters->maximum_step_size();
}


template<class R>
R
Evaluation::ConstraintBasedHybridEvolver<R>::maximum_basic_set_radius() const
{
  return this->_parameters->maximum_basic_set_radius();
}


template<class R>
Numeric::Rational
Evaluation::ConstraintBasedHybridEvolver<R>::lock_to_grid_time() const
{
  return this->_parameters->lock_to_grid_time();
}



template<class R> 
typename Evaluation::ConstraintBasedHybridEvolver<R>::working_sets_type
Evaluation::ConstraintBasedHybridEvolver<R>::_compute_working_sets(const hybrid_list_set_type& set) const
{
  working_sets_type working_sets;
  /*
  for(typename hybrid_list_set_type::const_iterator iter=set.begin();
      iter!=set.end(); ++iter)
  {
    id_type loc_id=iter->discrete_state();
    std::cerr << "loc_id=" << loc_id << std::endl;
    const continuous_basic_set_type& bs=iter->continuous_state_set();
    std::cerr << "bs=" << bs << std::endl;
    timed_set_type timed_set(Rational(0),Integer(0),loc_id,bs);
    std::cerr << "ts=" << timed_set << std::endl;
    working_sets.push_back(timed_set);
    std::cerr << "working_sets.size()=" << working_sets.size() << std::endl;
  }
  */
  for(typename hybrid_list_set_type::locations_const_iterator loc_iter=set.locations_begin();
      loc_iter!=set.locations_end(); ++loc_iter)
  {
    discrete_state_type loc_id=loc_iter->first;
    for(typename ListSet<BS>::const_iterator set_iter=loc_iter->second.begin();
        set_iter!=loc_iter->second.end(); ++set_iter)
    {
      const continuous_basic_set_type& bs=*set_iter;
      timed_set_type timed_set(Rational(0),Integer(0),loc_id,bs);
      working_sets.push_back(timed_set);
    }
  }

  return working_sets;
}


template<class R> inline
typename Evaluation::ConstraintBasedHybridEvolver<R>::hybrid_list_set_type
Evaluation::ConstraintBasedHybridEvolver<R>::_compute_list_set(const working_sets_type& working_sets, const HybridSpace& locations) const
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
Geometry::HybridSet<R>
Evaluation::ConstraintBasedHybridEvolver<R>::discrete_step(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                            const Geometry::HybridSet<R>& initial_set) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
Geometry::HybridSet<R>
Evaluation::ConstraintBasedHybridEvolver<R>::continuous_chainreach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridSet<R>& initial_set, 
                                                    const Geometry::HybridSet<R>& bounding_set) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
Geometry::HybridSet<R> 
Evaluation::ConstraintBasedHybridEvolver<R>::lower_evolve(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                           const Geometry::HybridSet<R>&,
                                           time_type evolution_time,
                                           size_type maximum_number_of_events) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}
      

template<class R>
Geometry::HybridSet<R> 
Evaluation::ConstraintBasedHybridEvolver<R>::upper_evolve(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                           const Geometry::HybridSet<R>&,
                                           time_type evolution_time,
                                           size_type maximum_number_of_events) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}
      

template<class R>
Geometry::HybridSet<R> 
Evaluation::ConstraintBasedHybridEvolver<R>::lower_reach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                          const Geometry::HybridSet<R>&,
                                          time_type initial_time, 
                                          time_type final_time, 
                                          size_type maximum_number_of_events) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}
      

template<class R>
Geometry::HybridSet<R> 
Evaluation::ConstraintBasedHybridEvolver<R>::upper_reach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                          const Geometry::HybridSet<R>&,
                                          time_type initial_time, 
                                          time_type final_time, 
                                          size_type maximum_number_of_events) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}
      

template<class R>
Geometry::HybridSet<R>
Evaluation::ConstraintBasedHybridEvolver<R>::chainreach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                         const Geometry::HybridSet<R>& initial_set, 
                                         const Geometry::HybridSet<R>& bounding_set) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
Geometry::HybridSet<R>
Evaluation::ConstraintBasedHybridEvolver<R>::viable(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                     const Geometry::HybridSet<R>& bounding_set) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
tribool
Evaluation::ConstraintBasedHybridEvolver<R>::verify(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                     const Geometry::HybridSet<R>& initial_set, 
                                     const Geometry::HybridSet<R>& sage_set) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}





template<class R>
Geometry::HybridListSet<typename Evaluation::ConstraintBasedHybridEvolver<R>::continuous_basic_set_type> 
Evaluation::ConstraintBasedHybridEvolver<R>::upper_evolve(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                     const Geometry::HybridListSet<continuous_basic_set_type>& initial_set, 
                                                     time_type evolution_time,
                                                     size_type maximum_number_of_events) const
{
  ARIADNE_LOG(2,"ConstraintBasedHybridEvolver::upper_evolve(ConstraintBasedHybridAutomaton automaton, ListSet initial_set, Time time, Integer maximum_number_of_events)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<", evolution_time="<<evolution_time<<"\n");
  
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  std::vector< timed_set_type > working_sets=this->_compute_working_sets(initial_set);
  ARIADNE_LOG(5,"working_sets="<<working_sets<<"\n");
  std::vector< timed_set_type > final_sets;
  std::vector< timed_set_type >& trace_sets=this->_trace;
  trace_sets=working_sets;

  while(!working_sets.empty()) {
    ARIADNE_LOG(9,"working_sets.size()="<<working_sets.size()<<"\n");
    timed_set_type working_set=working_sets.back(); working_sets.pop_back();
    ARIADNE_LOG(9,"unregularised_working_set="<<working_set<<"\n");
    working_set=this->_scheduler->regularize(working_set);
    ARIADNE_LOG(9,"regularised_working_set="<<working_set<<"\n");
    if(possibly(working_set.time()>evolution_time)) {
      ARIADNE_LOG(3,"\n  working_set.time()="<<working_set.time());
      ARIADNE_LOG(3,", evolution_time="<<evolution_time<<"\n\n");
      ARIADNE_LOG(3,"Warning: evolution time may be exceeded\n\n");
      final_sets.push_back(working_set);
    } else {
      if(working_set.time()==evolution_time) {
        ARIADNE_LOG(9,"  reached evolution time\n");
        final_sets.push_back(working_set);
      } else if(working_set.continuous_state_set().radius()>this->maximum_basic_set_radius()) {
        ::append(working_sets,this->_scheduler->subdivide(working_set));
      } else {
        ::append(working_sets,this->_scheduler->upper_integration_step(automaton,working_set,evolution_time,this->maximum_step_size()));
        ::append(this->_trace,this->_scheduler->trace);
      }
    }
  }

  return this->_compute_list_set(final_sets,automaton.locations());
}


template<class R>
Geometry::HybridListSet<typename Evaluation::ConstraintBasedHybridEvolver<R>::continuous_basic_set_type> 
Evaluation::ConstraintBasedHybridEvolver<R>::upper_reach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridListSet<continuous_basic_set_type>& initial_set, 
                                                    time_type evolution_time, 
                                                    size_type maximum_number_of_events) const
{
  ARIADNE_LOG(2,"ConstraintBasedHybridEvolver::upper_reach(ConstraintBasedHybridAutomaton automaton, ListSet initial_set, Time time, Integer maximum_number_of_events)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<", evolution_time="<<evolution_time<<"\n");
  ARIADNE_LOG(7,"maximum_step_size="<<this->maximum_step_size()<<", maximum_basic_set_radius="<<this->maximum_basic_set_radius()<<"\n");

  
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  std::vector< timed_set_type > working_sets=this->_compute_working_sets(initial_set);
  std::vector< timed_set_type > reached_sets;

  while(!working_sets.empty()) {
    timed_set_type working_set=working_sets.back(); working_sets.pop_back();
    if(working_set.time()>=evolution_time) {
    } else if(working_set.continuous_state_set().radius()>this->maximum_basic_set_radius()) {
      ::append(working_sets,this->_scheduler->subdivide(working_set));
    } else {
      ::append(working_sets,this->_scheduler->upper_reachability_step(automaton,working_set,evolution_time,this->maximum_step_size()));
      // A reachability step returns a list whose final element is the continuous reached set
      reached_sets.push_back(working_sets.back()); working_sets.pop_back();
      if(!working_sets.empty()) { std::cout << "evolution step: " << working_sets.back() << std::endl; }
      if(!reached_sets.empty()) { std::cout << "reach step: " << reached_sets.back() << std::endl; }
    }
  }

  return this->_compute_list_set(reached_sets,automaton.locations());
}




template<class R>
Geometry::HybridGridCellListSet<R>
Evaluation::ConstraintBasedHybridEvolver<R>::upper_evolve(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                     const Geometry::HybridGridCell<R>& initial_set, 
                                                     const Geometry::HybridGrid<R>& hybrid_grid, 
                                                     time_type evolution_time, 
                                                     size_type maximum_number_of_events) const
{
  HybridGridCellListSet<R> final_set(hybrid_grid);
  
  HybridListSet<BS> list_set(automaton.locations());
  list_set.adjoin(initial_set);
  
  list_set=this->upper_evolve(automaton,list_set,evolution_time,maximum_number_of_events);
  
  final_set.adjoin_outer_approximation(list_set);

  return final_set;
}


template<class R>
Geometry::HybridGridCellListSet<R>
Evaluation::ConstraintBasedHybridEvolver<R>::upper_reach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridGridCell<R>& initial_set, 
                                                    const Geometry::HybridGrid<R>& hybrid_grid, 
                                                    time_type evolution_time, 
                                                    size_type maximum_number_of_events) const
{
  ARIADNE_LOG(2,"ConstraintBasedHybridEvolver::upper_reach(ConstraintBasedHybridAutomaton, HybridGridCell, HybridGrid, Time, Integer)\n");
  HybridGridCellListSet<R> final_set(hybrid_grid);
  
  HybridListSet<BS> list_set(automaton.locations());
  list_set.adjoin(initial_set);
  
  list_set=this->upper_reach(automaton,list_set,evolution_time,maximum_number_of_events);
  
  final_set.adjoin_outer_approximation(list_set);

  return final_set;
}



template<class R>
Geometry::HybridGridMaskSet<R>
Evaluation::ConstraintBasedHybridEvolver<R>::upper_evolve(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
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
  size_type integration_steps = int(ceil(time_type(evolution_time/lock_to_grid_time)));
  lock_to_grid_time = evolution_time/integration_steps;

  for(uint i=0; i!=integration_steps; ++i) {
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
Evaluation::ConstraintBasedHybridEvolver<R>::upper_reach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                    time_type evolution_time, 
                                                    size_type maximum_number_of_events) const
{
  ARIADNE_LOG(2,"ConstraintBasedHybridEvolver::upper_reach(ConstraintBasedHybridAutomaton, HybridGridMaskSet, Time, Integer)\n");
  Geometry::HybridGrid<R> hybrid_grid=initial_set.grid();
  Geometry::HybridGridMaskSet<R> current_set=initial_set;
  Geometry::HybridGridMaskSet<R> next_set(current_set);
  next_set.clear();

  Geometry::HybridGridMultiMap<R> map(initial_set.grid(),initial_set.grid());

  size_type integration_steps = int(ceil(time_type(evolution_time/this->lock_to_grid_time())));
  time_type lock_to_grid_time = evolution_time/integration_steps;

  for(typename HybridGridMaskSet<R>::const_iterator cell_iter=current_set.begin();
      cell_iter!=current_set.end(); ++cell_iter)
  {
    const HybridGridCell<R>& cell=*cell_iter;
    next_set.adjoin(this->upper_reach(automaton,cell,hybrid_grid,lock_to_grid_time,maximum_number_of_events));
  }
  current_set=next_set;
  next_set.clear();

  for(uint i=0; i!=integration_steps-1; ++i) {
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
Evaluation::ConstraintBasedHybridEvolver<R>::viable(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                     const Geometry::HybridGridMaskSet<R>& bounding_set) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
Geometry::HybridGridMaskSet<R>
Evaluation::ConstraintBasedHybridEvolver<R>::chainreach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
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
