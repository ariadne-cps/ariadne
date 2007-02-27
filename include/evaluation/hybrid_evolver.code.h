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
 
#include "../declarations.h"
#include "../exceptions.h"
#include "../geometry/hybrid_set.h"
#include "../system/hybrid_automaton.h"
#include "../evaluation/applicator.h"
#include "../evaluation/integrator.h"
#include "hybrid_evolver.h"

namespace Ariadne {
  
template<class Set> 
class HybridTimedSet 
{
 public:
  HybridTimedSet(const time_type& t, const id_type& id, const Set& s)
    : _time(t), _discrete_state(id), _continuous_state_set(s) { }
  const time_type& time() const { return _time; }
  const id_type& discrete_state() const { return _discrete_state; }
  const Set& continuous_state_set() const { return _continuous_state_set; } 
  
  bool operator==(const HybridTimedSet& other) const { 
    return this->_time == other._time 
      && this->_discrete_state==other._discrete_state 
      && this->_continuous_state_set==other._continuous_state_set; }
  bool operator!=(const HybridTimedSet& other) const { return !(*this==other); }
  bool operator<=(const HybridTimedSet& other) const { return this->_time <= other._time; }
 private:
  time_type _time;
  id_type _discrete_state;
  Set _continuous_state_set;
};

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
Ariadne::Geometry::HybridGridCellListSet<R> 
Ariadne::Evaluation::HybridEvolver<R>::discrete_step(const System::HybridAutomaton<R>& hybrid_automaton, 
                                                     const Geometry::HybridGridCellListSet<R>& initial_set)
{

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
                                                             const Geometry::HybridGridMaskSet<R>& invariants)
{
  //if(verbosity>5) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }
  Geometry::HybridGridMaskSet<R> result_set(initial_set);
  
  typedef typename System::HybridAutomaton<R>::discrete_mode_iterator discrete_mode_iterator;
  
  for(discrete_mode_iterator dm_iter=hybrid_automaton.modes().begin();
      dm_iter!=hybrid_automaton.modes().end(); ++dm_iter)
  {
    const System::DiscreteMode<R>& dm = *dm_iter;
    const Geometry::GridMaskSet<R>& invariant=invariants[dm.id()];
    Geometry::GridMaskSet<R>  continuous_set=regular_intersection(result_set[dm.id()],invariant);
    result_set[dm.id()].adjoin(this->_integrator->chainreach(dm.dynamic(),continuous_set,invariant));
  }
  return result_set;
}


template<class R>
Ariadne::Geometry::HybridGridMaskSet<R> 
Ariadne::Evaluation::HybridEvolver<R>::chainreach(const System::HybridAutomaton<R>& hybrid_automaton, 
                                                  const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                  const Geometry::HybridGridMaskSet<R>& bounding_set)
{
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
    Geometry::GridCellListSet<R>& activation=activations[dt.id()];
    activation.adjoin(over_approximation(dt.activation(),activation.grid()));
    const Geometry::Set<R>& activation_set=dt.activation();
  }
  
  // Compute invariant regions as GridMaskSets
  Geometry::HybridGridMaskSet<R> invariants(initial_set);
  invariants.clear();
  for(discrete_mode_iterator dm_iter=hybrid_automaton.modes().begin();
      dm_iter!=hybrid_automaton.modes().end(); ++dm_iter)
  {
    const System::DiscreteMode<R>& dm=*dm_iter;
    Geometry::GridMaskSet<R>& invariant=invariants[dm.id()];
    invariant=bounding_set[dm.id()];
    invariant.restrict(over_approximation(dm.invariant(),invariant.grid()));
  }

  Geometry::HybridGridMaskSet<R> intermediate_set=this->continuous_chainreach(hybrid_automaton,initial_set,invariants);
  Geometry::HybridGridMaskSet<R> result_set=intermediate_set;
  
  Geometry::HybridGridCellListSet<R> new_activated=regular_intersection(intermediate_set,activations);
  while(!new_activated.empty()) {
    intermediate_set.clear();
    Geometry::HybridGridCellListSet<R> new_activated_image=this->discrete_step(hybrid_automaton,new_activated);
    intermediate_set.adjoin(new_activated_image);
    intermediate_set=this->continuous_chainreach(hybrid_automaton,intermediate_set,invariants);
    result_set.adjoin(intermediate_set);
    new_activated=regular_intersection(intermediate_set,activations);
  }
  return result_set;
}

