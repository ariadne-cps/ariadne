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

#include "../output/epsfstream.h"
#include "../output/logging.h"

#include "hybrid_evolver.h"

namespace Ariadne {
 
namespace Evaluation { static int& verbosity = hybrid_evolver_verbosity; }



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



template<class R>
Evaluation::HybridEvolver<R>::~HybridEvolver()
{
}

template<class R>
Evaluation::HybridEvolver<R>::HybridEvolver(Applicator<R>& a, Integrator<R>& i)
  : _applicator(&a), _integrator(&i) 
{
}





template<class R>
Geometry::HybridSet<R> 
Evaluation::HybridEvolver<R>::discrete_step(const System::HybridAutomaton<R>& automaton, 
                                            const Geometry::HybridSet<R>& initial_set)
{
  ARIADNE_LOG(2,"HybridSet HybridEvolver::discrete_step(HybridAutomaton automaton, HybridSet initial_set)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  if(automaton.locations()!=initial_set.locations()) {
    throw std::runtime_error("HybridEvolver::discrete_step(HybridAutomaton,HybridSet): initial_set locations do not match hybrid_automaton modes");
  }

  typedef Geometry::HybridSpace::const_iterator locations_iterator;

  Geometry::HybridSpace locations=automaton.locations();
  R grid_separation=this->_applicator->grid_size();
  
  Geometry::HybridGridCellListSet<R> grid_initial_set;
  for(locations_iterator loc_iter=locations.begin(); loc_iter!=locations.end(); ++loc_iter) {
    id_type id=loc_iter->id();
    dimension_type dim=loc_iter->dimension();
    grid_initial_set.new_location(id,Geometry::Grid<R>(LinearAlgebra::Vector<R>(dim,grid_separation)));
    grid_initial_set[id].adjoin_over_approximation(initial_set[id]);
  }
  Geometry::HybridGridCellListSet<R> grid_result=discrete_step(automaton,grid_initial_set);
  Geometry::HybridSetBase< Geometry::GridCellListSet<R> >& grid_base_result(grid_result);
  Geometry::HybridSet<R> result(grid_base_result);
  return result;
}


template<class R>
Geometry::HybridSet<R> 
Evaluation::HybridEvolver<R>::continuous_chainreach(const System::HybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridSet<R>& initial_set,
                                                    const Geometry::HybridSet<R>& bounding_set)
{
  ARIADNE_LOG(2,"HybridSet HybridEvolver::continuous_chainreach(HybridAutomaton automaton, HybridSet initial_set, HybridSet bounding_set)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n"<<"bounding_set="<<bounding_set<<"\n");
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
Geometry::HybridSet<R> 
Evaluation::HybridEvolver<R>::lower_reach(const System::HybridAutomaton<R>&, 
                                                   const Geometry::HybridSet<R>&, 
                                                   time_type t1, 
                                                   time_type t2, 
                                                   size_type nmax)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
Geometry::HybridSet<R> 
Evaluation::HybridEvolver<R>::upper_reach(const System::HybridAutomaton<R>& hybrid_automaton, 
                                                   const Geometry::HybridSet<R>& intial_set, 
                                                   time_type t1, 
                                                   time_type t2, 
                                                   size_type nmax)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
Geometry::HybridSet<R> 
Evaluation::HybridEvolver<R>::chainreach(const System::HybridAutomaton<R>& automaton, 
                                         const Geometry::HybridSet<R>& initial_set, 
                                         const Geometry::HybridSet<R>& bounding_set)
{
  using namespace Geometry;
  ARIADNE_LOG(1,"HybridEvolver::chainreach(HybridAutomaton,HybridSet,HybridSet)"<<std::endl);
  ARIADNE_LOG(2,"  initial_set="<<initial_set<<"\n  bounding_set="<<bounding_set<<std::endl);

  ARIADNE_CHECK_SAME_LOCATIONS(automaton,initial_set,"HybridEvolver::chainreach(HybridAutomaton,HybridSet,HybridSet)");
  ARIADNE_CHECK_SAME_LOCATIONS(automaton,bounding_set,"HybridEvolver::chainreach(HybridAutomaton,HybridSet,HybridSet)");
  ARIADNE_CHECK_BOUNDED(bounding_set,"HybridEvolver::chainreach(HybridAutomaton,HybridSet,HybridSet): bounding_set");
  ARIADNE_LOG(5,"Checked input"<<std::endl);
  
  HybridGridMaskSet<R> grid_bounding_set;
  HybridGridMaskSet<R> grid_initial_set;
  
  for(typename Geometry::HybridSet<R>::const_iterator bs_iter=bounding_set.begin();
      bs_iter!=bounding_set.end(); ++bs_iter)
  {
    id_type id=bs_iter->first;
    Grid<R> grid(bs_iter->second.dimension(),this->_applicator->grid_size());
    FiniteGrid<R> fgrid(grid,bs_iter->second.bounding_box());
    ARIADNE_LOG(5,"Made grid"<<std::endl);
    grid_bounding_set.new_location(id,fgrid);
    grid_initial_set.new_location(id,fgrid);
    grid_bounding_set[id].adjoin_over_approximation(bs_iter->second);
    grid_bounding_set[id].restrict_over_approximation(automaton.mode(id).invariant());
    grid_initial_set[id].adjoin_over_approximation(bs_iter->second);
    grid_initial_set[id].restrict_over_approximation(initial_set[id]);
  }
  ARIADNE_LOG(5,"Made cells"<<std::endl);


  ARIADNE_LOG(2,"  grid_initial_set="<<grid_initial_set<<"\n  grid_bounding_set="<<grid_bounding_set<<std::endl);
  HybridGridMaskSet<R> grid_chainreach_set=this->chainreach(automaton,grid_initial_set,grid_bounding_set);
  ARIADNE_LOG(2,"  grid_chainreach_set="<<grid_chainreach_set<<std::endl);
  HybridSet<R> chainreach_set(grid_chainreach_set);
  ARIADNE_LOG(2,"  chainreach_set="<<chainreach_set<<std::endl);
  return chainreach_set;
}


template<class R>
Geometry::HybridSet<R> 
Evaluation::HybridEvolver<R>::viable(const System::HybridAutomaton<R>& automaton, 
                                              const Geometry::HybridSet<R>& bounding_set)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}
     

template<class R>
tribool 
Evaluation::HybridEvolver<R>::verify(const System::HybridAutomaton<R>& automaton, 
                                              const Geometry::HybridSet<R>& initial_set, 
                                              const Geometry::HybridSet<R>& safe_set)

{
  throw NotImplemented(__PRETTY_FUNCTION__);
}






template<class R>
Geometry::HybridGridMaskSet<R> 
Evaluation::HybridEvolver<R>::lower_reach(const System::HybridAutomaton<R>&, 
                                                   const Geometry::HybridGridMaskSet<R>&, 
                                                   time_type t1, 
                                                   time_type t2, 
                                                   size_type nmax)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R>
Geometry::HybridGridMaskSet<R> 
Evaluation::HybridEvolver<R>::upper_reach(const System::HybridAutomaton<R>& hybrid_automaton, 
                                                   const Geometry::HybridGridMaskSet<R>& intial_set, 
                                                   time_type t1, 
                                                   time_type t2, 
                                                   size_type nmax)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}




template<class R>
Geometry::HybridGridMaskSet<R> 
Evaluation::HybridEvolver<R>::_discrete_step(const System::HybridAutomaton<R>& hybrid_automaton, 
                                             const Geometry::HybridGridMaskSet<R>& initial_set, 
                                             const Geometry::HybridGridMaskSet<R>& domain_set)
{
  Geometry::HybridGridMaskSet<R> result_set(initial_set);
  result_set.clear();
  
  typedef typename System::HybridAutomaton<R>::discrete_transition_iterator discrete_transition_iterator;
  
  for(discrete_transition_iterator dt_iter=hybrid_automaton.transitions().begin();
      dt_iter!=hybrid_automaton.transitions().end(); ++dt_iter)
  {
    const System::DiscreteTransition<R>& dt = *dt_iter;
    const Geometry::GridCellListSet<R>& source_set=initial_set[dt.source().id()];
    Geometry::GridMaskSet<R>& destination_set=result_set[dt.destination().id()];
    const Geometry::GridMaskSet<R> activation=Geometry::over_approximation(dt.activation(),source_set.grid());
    const Geometry::GridMaskSet<R> active_cells=Geometry::regular_intersection(source_set,activation);
    ARIADNE_LOG(4,"discrete_step of transition "<<dt.id()<<" from mode "<<dt.source().id()<<" to mode "<<dt.destination().id()<<":\n")
      ARIADNE_LOG(4,"  "<<active_cells.size()<<" activated cells");
    const Geometry::GridMaskSet<R> image_cells=this->_applicator->image(dt.reset(),active_cells,destination_set.grid());
    ARIADNE_LOG(4,", "<<image_cells.size()<<" cells in image\n");
    destination_set.adjoin(image_cells);
  }
  
  return result_set;
}


template<class R>
Geometry::HybridGridMaskSet<R> 
Evaluation::HybridEvolver<R>::_continuous_chainreach(const System::HybridAutomaton<R>& hybrid_automaton, 
                                                     const Geometry::HybridGridMaskSet<R>& initial_set,
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
    ARIADNE_LOG(4,"continuous_chainreach in mode "<<dm.id()<<":\n  "<<initial.size()<<" initial cells,");
    result_set[dm.id()].adjoin(this->_integrator->chainreach(dm.dynamic(),initial,domain));
    result_set[dm.id()].restrict(domain);
    ARIADNE_LOG(4," reached "<<result_set[dm.id()].size()<<" cells\n");
  }
  return result_set;
}


template<class R>
Geometry::HybridGridCellListSet<R> 
Evaluation::HybridEvolver<R>::discrete_step(const System::HybridAutomaton<R>& hybrid_automaton, 
                                            const Geometry::HybridGridCellListSet<R>& initial_set)
{
  ARIADNE_LOG(2,"HybridGridCellListSet HybridEvolver::discrete_step(HybridAutomaton automaton, HybridGridCellListSet initial_set)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  ARIADNE_CHECK_SAME_LOCATIONS(hybrid_automaton,initial_set,"HybridEvolver::discrete_step(HybridAutomaton,HybridGridCellListSet)");

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
    const Geometry::GridMaskSet<R> active_cells=Geometry::regular_intersection(source_set,activation);
    ARIADNE_LOG(4,"discrete_step of transition "<<dt.id()<<" from mode "<<dt.source().id()<<" to mode "<<dt.destination().id()<<":\n")
    ARIADNE_LOG(4,"  "<<active_cells.size()<<" activated cells");
    const Geometry::GridMaskSet<R> image_cells=this->_applicator->image(dt.reset(),active_cells,destination_set.grid());
    ARIADNE_LOG(4,", "<<image_cells.size()<<" cells in image\n");
    destination_set.adjoin(image_cells);

  }
  
  return result_set;
}


template<class R>
Geometry::HybridGridMaskSet<R> 
Evaluation::HybridEvolver<R>::continuous_chainreach(const System::HybridAutomaton<R>& hybrid_automaton, 
                                                    const Geometry::HybridGridMaskSet<R>& initial_set,
                                                    const Geometry::HybridGridMaskSet<R>& bounding_set)
{
  ARIADNE_LOG(2,"HybridGridMaskSet HybridEvolver::continuous_chainreach(HybridAutomaton automaton, HybridGridMaskSet initial_set, HybridGridMaskSet bounding_set)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n"<<"bounding_set="<<bounding_set<<"\n");
  ARIADNE_LOG(3,"checking input... ");
  ARIADNE_CHECK_SAME_LOCATIONS(hybrid_automaton,initial_set,"HybridEvolver::continuous_chainreach(HybridAutomaton,HybridGridMaskSet,HybridGridMaskSet)");
  ARIADNE_CHECK_SAME_LOCATIONS(hybrid_automaton,bounding_set,"HybridEvolver::continuous_chainreach(HybridAutomaton,HybridGridMaskSet,HybridGridMaskSet)");
  
  ARIADNE_CHECK_BOUNDED(bounding_set,"HybridEvolver::continuous_chainreach(HybridAutomaton,HybridGridMaskSet,HybridGridMaskSet): bounding_set");
  
  ARIADNE_LOG(3,"successful\n");

  typedef typename System::HybridAutomaton<R>::discrete_mode_iterator discrete_mode_iterator;
  
  Geometry::HybridGridMaskSet<R> domain_set(bounding_set);
  for(discrete_mode_iterator dm_iter=hybrid_automaton.modes().begin();
      dm_iter!=hybrid_automaton.modes().end(); ++dm_iter)
  {
    ARIADNE_LOG(3,"computing domain set in mode "<<dm_iter->id()<<"... invariant="<<dm_iter->invariant()<<"\n");
    domain_set[dm_iter->id()].restrict_over_approximation(dm_iter->invariant());
  }
  ARIADNE_LOG(3,"domain_set="<<domain_set<<"\n");

  return _continuous_chainreach(hybrid_automaton,initial_set,domain_set);
}



template<class R>
Geometry::HybridGridMaskSet<R> 
Evaluation::HybridEvolver<R>::chainreach(const System::HybridAutomaton<R>& hybrid_automaton, 
                                         const Geometry::HybridGridMaskSet<R>& initial_set, 
                                         const Geometry::HybridGridMaskSet<R>& bounding_set)
{
  using namespace Geometry;
  ARIADNE_LOG(2,"HybridGridMaskSet HybridEvolver::chainreach(HybridAutomaton automaton, HybridGridMaskSet initial_set, HybridGridMaskSet bounding_set)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n"<<"bounding_set="<<bounding_set<<"\n");
  ARIADNE_CHECK_SAME_LOCATIONS(hybrid_automaton,initial_set,"HybridEvolver::chainreach(HybridAutomaton,HybridGridMaskSet,HybridGridMaskSet)");
  ARIADNE_CHECK_SAME_LOCATIONS(hybrid_automaton,bounding_set,"HybridEvolver::chainreach(HybridAutomaton,HybridGridMaskSet,HybridGridMaskSet)");
  
  ARIADNE_CHECK_BOUNDED(bounding_set,"HybridEvolver::continuous_chainreach(HybridAutomaton,HybridGridMaskSet,HybridGridMaskSet): bounding_set");

  typedef typename System::HybridAutomaton<R>::discrete_transition_iterator discrete_transition_iterator;
  typedef typename System::HybridAutomaton<R>::discrete_mode_iterator discrete_mode_iterator;
  
  // Compute restricted invariant domains as GridMaskSets
  Geometry::HybridGridMaskSet<R> domain_set(bounding_set);
  for(discrete_mode_iterator dm_iter=hybrid_automaton.modes().begin();
      dm_iter!=hybrid_automaton.modes().end(); ++dm_iter)
  {
    domain_set[dm_iter->id()].restrict_over_approximation(dm_iter->invariant());
  }

  HybridGridMaskSet<R> result_set=domain_set;
  HybridGridMaskSet<R> integrated_set=domain_set;
  integrated_set.clear();
  result_set.restrict(initial_set);
  ARIADNE_LOG(4,result_set.size()<<" cells in initial set\n");
  HybridGridMaskSet<R> found_set=result_set;
  uint step=0;
  while(!found_set.empty()) {
    integrated_set=this->_continuous_chainreach(hybrid_automaton,found_set,domain_set);
    ARIADNE_LOG(4,"\nchainreach found "<<found_set.size()<<" cell by continuous evolution, \n"<<std::endl);
    found_set=this->_discrete_step(hybrid_automaton,integrated_set,domain_set);
    ARIADNE_LOG(4,"\nchainreach found "<<found_set.size()<<" cells by discrete step");
    found_set.remove(result_set);
    ARIADNE_LOG(4," of which "<<found_set.size()<<" are new; ");
    result_set.adjoin(integrated_set);
    result_set.adjoin(found_set);
    ARIADNE_LOG(4,"reached "<<result_set.size()<< " cells in total.\n" << std::endl);
    ++step;
    if(verbosity>=4) {
      std::stringstream filename;
      filename << "hybrid_chainreach-"<<step<<".eps";
      const GridMaskSet<R>& dom=domain_set.begin()->second;
      const GridMaskSet<R>& res=result_set.begin()->second;
      const GridMaskSet<R>& fnd=found_set.begin()->second;
      Output::epsfstream eps;
      eps.open(filename.str().c_str(),dom.bounding_box());
      eps.set_fill_colour("cyan"); eps<<dom.extent();
      eps.set_fill_colour("green"); eps<<res;
      eps.set_fill_colour("red"); eps<<fnd;
      eps.close();
    }
  }
  return result_set;
}


template<class R>
Geometry::HybridGridMaskSet<R> 
Evaluation::HybridEvolver<R>::viable(const System::HybridAutomaton<R>& hybrid_automaton, 
                                              const Geometry::HybridGridMaskSet<R>& bounding_set)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
tribool
Evaluation::HybridEvolver<R>::verify(const System::HybridAutomaton<R>& hybrid_automaton, 
                                              const Geometry::HybridGridMaskSet<R>& initial_set, 
                                              const Geometry::HybridGridMaskSet<R>& bounding_set)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

}
