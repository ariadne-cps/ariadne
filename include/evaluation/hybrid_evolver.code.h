/***************************************************************************
 *            hybrid_evolver.code.h
 *
 *  Copyright  2004-7  Alberto Casagrande,  Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
 
#include "../geometry/rectangle_expression.h"
#include "../geometry/set_interface.h"
#include "../geometry/hybrid_set.h"
#include "../system/hybrid_automaton.h"
#include "../evaluation/applicator.h"
#include "../evaluation/integrator.h"

#include "../output/logging.h"

#include "hybrid_evolver.h"

namespace Ariadne {
  
template<class SetInterface> 
class HybridTimedSet 
{
 public:
  HybridTimedSet(const time_type& t, const id_type& id, const SetInterface& s)
    : _time(t), _discrete_state(id), _continuous_state_set(s) { }
  const time_type& time() const { return _time; }
  const id_type& discrete_state() const { return _discrete_state; }
  const SetInterface& continuous_state_set() const { return _continuous_state_set; } 
  
  bool operator==(const HybridTimedSet& other) const { 
    return this->_time == other._time 
      && this->_discrete_state==other._discrete_state 
      && this->_continuous_state_set==other._continuous_state_set; }
  bool operator!=(const HybridTimedSet& other) const { return !(*this==other); }
  bool operator<=(const HybridTimedSet& other) const { return this->_time <= other._time; }
 private:
  time_type _time;
  id_type _discrete_state;
  SetInterface _continuous_state_set;
};

}

template<class R>
Ariadne::Evaluation::HybridEvolver<R>::~HybridEvolver()
{
}

template<class R>
Ariadne::Evaluation::HybridEvolver<R>::HybridEvolver(Applicator<R>& a, Integrator<R>& i)
  : _applicator(&a), _integrator(&i) 
{
}

template<class R>
Ariadne::Geometry::HybridGridMaskSet<R> 
Ariadne::Evaluation::HybridEvolver<R>::lower_reach(const System::HybridAutomaton<R>&, 
                                                   const Geometry::HybridGridMaskSet<R>&, 
                                                   time_type t1, 
                                                   time_type t2, 
                                                   size_type nmax)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R>
Ariadne::Geometry::HybridGridMaskSet<R> 
Ariadne::Evaluation::HybridEvolver<R>::upper_reach(const System::HybridAutomaton<R>& hybrid_automaton, 
                                                   const Geometry::HybridGridMaskSet<R>& intial_set, 
                                                   time_type t1, 
                                                   time_type t2, 
                                                   size_type nmax)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R>
Ariadne::Geometry::HybridGridMaskSet<R> 
Ariadne::Evaluation::HybridEvolver<R>::_continuous_chainreach(const System::HybridAutomaton<R>& hybrid_automaton, 
                                                              const Geometry::HybridGridMaskSet<R>& initial_set,
                                                              const Geometry::HybridGridMaskSet<R>& invariant_set,
                                                              const Geometry::HybridGridMaskSet<R>& domain_set)
{
  Geometry::HybridGridMaskSet<R> result_set(initial_set);
  
  typedef typename System::HybridAutomaton<R>::discrete_mode_iterator discrete_mode_iterator;
  result_set.clear();
  
  for(discrete_mode_iterator dm_iter=hybrid_automaton.modes().begin();
      dm_iter!=hybrid_automaton.modes().end(); ++dm_iter)
  {
    const System::DiscreteMode<R>& dm = *dm_iter;
    const Geometry::GridMaskSet<R>& domain=domain_set[dm.id()];
    Geometry::GridMaskSet<R> initial=regular_intersection(initial_set[dm.id()],domain);
    result_set[dm.id()].adjoin(this->_integrator->chainreach(dm.dynamic(),initial,domain));
    result_set[dm.id()].restrict(invariant_set[dm.id()]);
  }
  return result_set;
}


template<class R>
Ariadne::Geometry::HybridGridCellListSet<R> 
Ariadne::Evaluation::HybridEvolver<R>::discrete_step(const System::HybridAutomaton<R>& hybrid_automaton, 
                                                     const Geometry::HybridGridCellListSet<R>& initial_set)
{
  if(hybrid_automaton.locations()!=initial_set.locations()) {
    throw std::runtime_error("HybridEvolver::discrete_step(...): initial_set locations do not match hybrid_automaton modes");
  }

  Geometry::HybridGridCellListSet<R> result_set(initial_set);
  result_set.clear();
  
  typedef typename System::HybridAutomaton<R>::discrete_transition_iterator discrete_transition_iterator;
  
  for(discrete_transition_iterator dt_iter=hybrid_automaton.transitions().begin();
      dt_iter!=hybrid_automaton.transitions().end(); ++dt_iter)
  {
    const System::DiscreteTransition<R>& dt = *dt_iter;
    const Geometry::GridCellListSet<R>& source_set=initial_set[dt.source().id()];
    Geometry::GridCellListSet<R>& destination_set=result_set[dt.destination().id()];
    const Geometry::GridMaskSet<R> activation=Geometry::over_approximation(dt.activation(),source_set.grid());
    destination_set.adjoin(this->_applicator->image(dt.reset(),regular_intersection(source_set,activation),destination_set.grid()));
  }
  
  return result_set;
}


template<class R>
Ariadne::Geometry::HybridGridMaskSet<R> 
Ariadne::Evaluation::HybridEvolver<R>::continuous_chainreach(const System::HybridAutomaton<R>& hybrid_automaton, 
                                                             const Geometry::HybridGridMaskSet<R>& initial_set,
                                                             const Geometry::HybridGridMaskSet<R>& bounding_set)
{
  if(hybrid_automaton.locations()!=initial_set.locations()) {
    throw std::runtime_error("HybridEvolver::continuous_chainreach(...): initial_set locations do not match hybrid_automaton modes");
  }
  if(hybrid_automaton.locations()!=bounding_set.locations()) {
    throw std::runtime_error("HybridEvolver::continuous_chainreach(...): bounding_set locations do not match hybrid_automaton modes");
  }
  
  for(typename Geometry::HybridGridMaskSet<R>::const_iterator bs_iter=bounding_set.begin();
      bs_iter!=bounding_set.end(); ++bs_iter)
  {
    if(!bs_iter->second.bounded()) {
      throw Geometry::UnboundedSet("HybridEvolver::continuous_chainreach(...): bounding_set is not compact");
    }
  }
  
  typedef typename System::HybridAutomaton<R>::discrete_mode_iterator discrete_mode_iterator;
  
  Geometry::HybridGridMaskSet<R> invariant_set(bounding_set);
  invariant_set.clear();
  Geometry::HybridGridMaskSet<R> domain_set(bounding_set);
  for(discrete_mode_iterator dm_iter=hybrid_automaton.modes().begin();
      dm_iter!=hybrid_automaton.modes().end(); ++dm_iter)
  {
    invariant_set[dm_iter->id()].adjoin_over_approximation(dm_iter->invariant());
    domain_set[dm_iter->id()].restrict(invariant_set[dm_iter->id()]);
  }
  
  return _continuous_chainreach(hybrid_automaton,initial_set,invariant_set,domain_set);
}



template<class R>
Ariadne::Geometry::HybridGridMaskSet<R> 
Ariadne::Evaluation::HybridEvolver<R>::chainreach(const System::HybridAutomaton<R>& hybrid_automaton, 
                                                  const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                  const Geometry::HybridGridMaskSet<R>& bounding_set)
{
  if(hybrid_automaton.locations()!=initial_set.locations()) {
    throw std::runtime_error("HybridEvolver::chainreach(...): initial_set locations do not match hybrid_automaton modes");
  }
  if(hybrid_automaton.locations()!=bounding_set.locations()) {
    throw std::runtime_error("HybridEvolver::chainreach(...): bounding_set locations do not match hybrid_automaton modes");
  }
  for(typename Geometry::HybridGridMaskSet<R>::const_iterator bs_iter=bounding_set.begin();
      bs_iter!=bounding_set.end(); ++bs_iter)
  {
    if(!bs_iter->second.bounded()) {
      throw Geometry::UnboundedSet("HybridEvolver::chainreach(...): bounding_set is not compact");
    }
  }
  
  
  typedef typename System::HybridAutomaton<R>::discrete_transition_iterator discrete_transition_iterator;
  typedef typename System::HybridAutomaton<R>::discrete_mode_iterator discrete_mode_iterator;
  
  // Compute activation regions as GridCellListSets
  //FIXME: translate grids properly
  Geometry::HybridGridCellListSet<R> activations(initial_set);
  activations.clear();
  for(discrete_transition_iterator dt_iter=hybrid_automaton.transitions().begin();
      dt_iter!=hybrid_automaton.transitions().end(); ++dt_iter)
  {
    const System::DiscreteTransition<R>& dt=*dt_iter;
    Geometry::GridCellListSet<R>& activation=activations[dt.source().id()];
    activation.adjoin_over_approximation(dt.activation());
  }
  
  // Compute restricted invariant domains as GridMaskSets
  Geometry::HybridGridMaskSet<R> invariant_set(bounding_set);
  invariant_set.clear();
  Geometry::HybridGridMaskSet<R> domain_set(bounding_set);
  for(discrete_mode_iterator dm_iter=hybrid_automaton.modes().begin();
      dm_iter!=hybrid_automaton.modes().end(); ++dm_iter)
  {
    const System::DiscreteMode<R>& dm=*dm_iter;
    Geometry::GridMaskSet<R>& invariant=invariant_set[dm.id()];
    Geometry::GridMaskSet<R>& domain=domain_set[dm.id()];
    invariant.adjoin_over_approximation(dm.invariant());
    domain.restrict_over_approximation(invariant);
  }

  Geometry::HybridGridMaskSet<R> already_activated=domain_set;
  already_activated.clear();

  Geometry::HybridGridMaskSet<R> intermediate_set=this->_continuous_chainreach(hybrid_automaton,initial_set,invariant_set,domain_set);
  Geometry::HybridGridMaskSet<R> result_set=intermediate_set;
  
  Geometry::HybridGridCellListSet<R> new_activated=regular_intersection(intermediate_set,activations);
  while(!Geometry::subset(new_activated,already_activated)) {
    if(verbosity > 5) { std::clog << new_activated.size() << " activated cells, " << std::flush; }
    new_activated.unique_sort();
    if(verbosity > 5) { std::clog << "of which " << new_activated.size() << " are not duplicates, " << std::flush; }
    new_activated=Geometry::difference(new_activated,already_activated);
    if(verbosity > 5) { std::clog << "and " << new_activated.size() << " are new. " << std::endl; }
    already_activated.adjoin(new_activated);
    intermediate_set.clear();
    Geometry::HybridGridCellListSet<R> new_activated_image=this->discrete_step(hybrid_automaton,new_activated);
    intermediate_set.adjoin(new_activated_image);
    intermediate_set=this->_continuous_chainreach(hybrid_automaton,intermediate_set,invariant_set,domain_set);
    result_set.adjoin(intermediate_set);
    new_activated=regular_intersection(intermediate_set,activations);
  }
  return result_set;
}


template<class R>
Ariadne::Geometry::HybridGridMaskSet<R> 
Ariadne::Evaluation::HybridEvolver<R>::viable(const System::HybridAutomaton<R>& hybrid_automaton, 
                                              const Geometry::HybridGridMaskSet<R>& bounding_set)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
Ariadne::tribool
Ariadne::Evaluation::HybridEvolver<R>::verify(const System::HybridAutomaton<R>& hybrid_automaton, 
                                              const Geometry::HybridGridMaskSet<R>& initial_set, 
                                              const Geometry::HybridGridMaskSet<R>& bounding_set)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


