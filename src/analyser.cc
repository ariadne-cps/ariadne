/***************************************************************************
 *            analyser.cc
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
 
#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "exceptions.h"

#include "numeric.h"

#include "vector.h"
#include "matrix.h"

#include "box.h"
#include "list_set.h"
#include "grid_set.h"

#include "orbit.h"

#include "hybrid_automaton.h"

#include "evolution_parameters.h"
#include "evolver_interface.h"
#include "discretiser_interface.h"

#include "discretiser.h"
#include "analyser.h"
#include "logging.h"


namespace Ariadne {
  
 
static const double DEFAULT_MAXIMUM_ENCLOSURE_RADIUS=0.25;
static const double DEFAULT_GRID_LENGTH=0.125;

static int verbosity=global_verbosity;


HybridAnalyser::
HybridAnalyser(const EvolutionParameters& parameters, 
               const EvolverInterface<HybridAutomaton,DefaultHybridEnclosureType>& evolver)
  : _parameters(new EvolutionParameters(parameters))
  , _discretiser(new HybridDiscretiser<DefaultEnclosureType>(evolver))
{
}






// Helper functions for operators on lists of sets.
HybridGridTreeSet
HybridAnalyser::_upper_reach(const HybridAutomaton& sys, 
                             const HybridGridTreeSet& set, 
                             const HybridTime& time, 
                             const int accuracy) const 
{
  HybridGridTreeSet result(set.grid()); 
  HybridGridTreeSet cells=set; 
  cells.mince(accuracy);    
  for(HybridGridTreeSet::locations_const_iterator loc_iter=cells.locations_begin(); loc_iter!=cells.locations_end(); ++loc_iter) {
    DiscreteState const& q=loc_iter->first;
    GridTreeSet const& loc_cells=loc_iter->second;
    for(GridTreeSet::const_iterator cell_iter=loc_cells.begin(); cell_iter!=loc_cells.end(); ++cell_iter) {
      HybridGridCell hybrid_cell(q,*cell_iter);
      Orbit<HybridGridCell> orbit=this->_discretiser->upper_evolution(sys,hybrid_cell,time,accuracy);
      HybridGridCellListSet reach=orbit.reach();
      result.adjoin(reach); 
    }
  }
  return result; 
}


HybridGridTreeSet
HybridAnalyser::_upper_evolve(const HybridAutomaton& sys, 
                              const HybridGridTreeSet& set, 
                              const HybridTime& time, 
                              const int accuracy) const 
{
  std::cerr<<"  upper_evolve"<<std::endl;
  GTS result(set.grid()); GTS cells=set; cells.mince(accuracy); 
  for(HybridGridTreeSet::locations_const_iterator loc_iter=cells.locations_begin(); loc_iter!=cells.locations_end(); ++loc_iter) {
    for(GridTreeSet::const_iterator cell_iter=loc_iter->second.begin(); cell_iter!=loc_iter->second.end(); ++cell_iter) {
      result.adjoin(this->_discretiser->upper_evolution(sys,HybridGridCell(loc_iter->first,*cell_iter),time,accuracy).final()); 
    } 
  }
  return result; 
}


std::pair<HybridGridTreeSet,HybridGridTreeSet>
HybridAnalyser::_upper_reach_evolve(const HybridAutomaton& sys, 
                                    const HybridGridTreeSet& set, 
                                    const HybridTime& time, 
                                    const int accuracy) const 
{
  std::cerr<<"  upper_reach_evolve"<<std::endl;
  std::pair<GTS,GTS> result(GTS(set.grid()),GTS(set.grid()));
  GTS& reach=result.first; GTS& evolve=result.second;
  GTS cells=set; cells.mince(accuracy); 
  for(HybridGridTreeSet::locations_const_iterator loc_iter=cells.locations_begin(); loc_iter!=cells.locations_end(); ++loc_iter) {
    for(GridTreeSet::const_iterator cell_iter=loc_iter->second.begin(); cell_iter!=loc_iter->second.end(); ++cell_iter) {
      Orbit<HybridGridCell> evolution=this->_discretiser->upper_evolution(sys,HybridGridCell(loc_iter->first,*cell_iter),time,accuracy);
      reach.adjoin(evolution.reach()); 
      evolve.adjoin(evolution.final()); 
    } 
  }
  return result; 
}




HybridAnalyser::ConcreteSetType*
HybridAnalyser::
lower_evolve(const SystemType& system, 
             const OvertSetType& initial_set,
             const TimeType& time) const
{
  int grid_depth = this->_parameters->grid_depth;
  GTS initial; GTS& final=*new GTS();

  // Improve accuracy of initial set for lower computations
  initial.adjoin_lower_approximation(initial_set,grid_depth+4);

  for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
    final.adjoin(this->_discretiser->lower_evolution(system,*bs_iter,time,grid_depth).final());
  }
  return &final;
}



HybridAnalyser::ConcreteSetType*
HybridAnalyser::
lower_reach(const SystemType& system, 
            const OvertSetType& initial_set,
            const TimeType& time) const
{
  int grid_depth = this->_parameters->grid_depth+4;
  GTS initial; GTS& reach=*new GTS();
  
  // Improve accuracy of initial set for lower computations
  initial.adjoin_lower_approximation(initial_set,grid_depth+4);
 
  for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
    reach.adjoin(this->_discretiser->lower_evolution(system,*bs_iter,time,grid_depth).intermediate());
  }
  return &reach;
}



std::pair<HybridAnalyser::ConcreteSetType*,HybridAnalyser::ConcreteSetType*>
HybridAnalyser::
lower_reach_evolve(const SystemType& system, 
                   const OvertSetType& initial_set,
                   const TimeType& time) const
{
  int grid_depth = this->_parameters->grid_depth;
  GTS initial; 
  
  GTS& reach=*new GTS; GTS& evolve=*new GTS;

  // Improve accuracy of initial set for lower computations
  initial.adjoin_lower_approximation(initial_set,grid_depth+4);

  for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
    Orbit<GC> orbit = this->_discretiser->lower_evolution(system,*bs_iter,time,grid_depth);
    reach.adjoin(orbit.reach());
    evolve.adjoin(orbit.final());
  }
  return make_pair(&reach,&evolve);
}


HybridAnalyser::ConcreteSetType*
HybridAnalyser::
upper_evolve(const SystemType& system, 
             const CompactSetType& initial_set,
             const TimeType& time) const
{
  verbosity=0;
  ARIADNE_LOG(2,"HybridAnalyser::upper_evolve(...)\n");
  GTS initial;
  int grid_depth = this->_parameters->grid_depth;
  initial.adjoin_outer_approximation(initial_set,grid_depth);
  ARIADNE_LOG(4,"initial="<<initial<<"\n");
  GTS& evolve=*new GTS(initial);
  ARIADNE_LOG(4,"initial_evolve="<<evolve<<"\n");
  Float real_time=time.first;
  Float lock_to_grid_time=this->_parameters->lock_to_grid_time;
  uint time_steps=uint(real_time/lock_to_grid_time);
  Float remainder_time=real_time-time_steps*lock_to_grid_time;
  HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,time_steps);
  HybridTime hybrid_remainder_time(remainder_time,time_steps);
  for(uint i=0; i!=time_steps; ++i) {
    evolve=this->_upper_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
  }
  evolve=this->_upper_evolve(system,evolve,hybrid_remainder_time,grid_depth);
  evolve.recombine();
  ARIADNE_LOG(4,"final_evolve="<<evolve<<"\n");
  return &evolve;
}



HybridAnalyser::ConcreteSetType*
HybridAnalyser::
upper_reach(const SystemType& system, 
            const CompactSetType& initial_set,
            const TimeType& time) const
{
  verbosity=6;
  ARIADNE_LOG(2,"HybridAnalyser::upper_reach(system,set,time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  GTS initial;
  int grid_depth = this->_parameters->grid_depth;
  ARIADNE_LOG(3,"grid_depth="<<grid_depth<<"\n");
  initial.adjoin_outer_approximation(initial_set,grid_depth);
  ARIADNE_LOG(4,"initial"<<initial<<"\n");
  GTS evolve=initial;
  ARIADNE_LOG(4,"evolve="<<evolve<<"\n");
  GTS& reach=*new GTS(evolve);
  ARIADNE_LOG(4,"reach="<<reach<<"\n");
  Float real_time=time.first;
  Float lock_to_grid_time = this->_parameters->lock_to_grid_time;
  uint time_steps=uint(real_time/lock_to_grid_time);
  Float remainder_time=real_time-time_steps*lock_to_grid_time;
  HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,time_steps);
  HybridTime hybrid_remainder_time(remainder_time,time_steps);
  ARIADNE_LOG(5,"evolve="<<evolve<<"\n");
  for(uint i=0; i!=time_steps; ++i) {
    reach.adjoin(this->_upper_reach(system,evolve,hybrid_lock_to_grid_time,grid_depth));
    ARIADNE_LOG(5,"reach="<<reach<<"\n");
    reach.recombine();
    ARIADNE_LOG(6,"reach.recombine()="<<reach<<"\n");
    evolve=this->_upper_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
    ARIADNE_LOG(5,"evolve="<<evolve<<"\n");
  }
  reach.adjoin(this->_upper_reach(system,evolve,hybrid_remainder_time,grid_depth));
  reach.recombine();
  ARIADNE_LOG(4,"final_reach="<<reach<<"\n");
  return &reach;
}


std::pair<HybridAnalyser::ConcreteSetType*,HybridAnalyser::ConcreteSetType*>
HybridAnalyser::
upper_reach_evolve(const SystemType& system, 
                   const CompactSetType& initial_set,
                   const TimeType& time) const
{
  verbosity=6;
  ARIADNE_LOG(2,"HybridAnalyser::upper_reach(system,set,time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  GTS initial;
  int grid_depth = this->_parameters->grid_depth;
  ARIADNE_LOG(3,"grid_depth="<<grid_depth<<"\n");
  initial.adjoin_outer_approximation(initial_set,grid_depth);
  ARIADNE_LOG(4,"initial"<<initial<<"\n");
  GTS& evolve=*new GTS(initial);
  ARIADNE_LOG(4,"evolve="<<evolve<<"\n");
  GTS& reach=*new GTS(initial);
  ARIADNE_LOG(4,"reach="<<reach<<"\n");
  Float real_time=time.first;
  Float lock_to_grid_time = this->_parameters->lock_to_grid_time;
  uint time_steps=uint(real_time/lock_to_grid_time);
  Float remainder_time=real_time-time_steps*lock_to_grid_time;
  HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,time_steps);
  HybridTime hybrid_remainder_time(remainder_time,time_steps);
  ARIADNE_LOG(5,"evolve="<<evolve<<"\n");
  for(uint i=0; i!=time_steps; ++i) {
    reach.adjoin(this->_upper_reach(system,evolve,hybrid_lock_to_grid_time,grid_depth));
    ARIADNE_LOG(5,"reach="<<reach<<"\n");
    reach.recombine();
    ARIADNE_LOG(6,"reach.recombine()="<<reach<<"\n");
    evolve=this->_upper_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
    ARIADNE_LOG(5,"evolve="<<evolve<<"\n");
  }
  reach.adjoin(this->_upper_reach(system,evolve,hybrid_remainder_time,grid_depth));
  reach.recombine();
  ARIADNE_LOG(4,"reach="<<reach<<"\n");
  evolve.recombine();
  ARIADNE_LOG(4,"evolve="<<evolve<<"\n");
  return std::make_pair(&reach,&evolve);
}






HybridAnalyser::ConcreteSetType*
HybridAnalyser::
chain_reach(const SystemType& system, 
            const CompactSetType& initial_set,
            const BoundingSetType& bounding_domain) const
{
  // FIXME: Use tree sets throughout
  uint verbosity=0;

  HybridTime hybrid_transient_time(this->_parameters->transient_time,this->_parameters->transient_steps);
  HybridTime hybrid_lock_to_grid_time(this->_parameters->lock_to_grid_time,this->_parameters->lock_to_grid_steps);
  int grid_depth = this->_parameters->grid_depth;

  ARIADNE_LOG(5,"bounding_domain="<<bounding_domain<<"\n");
  ARIADNE_LOG(5,"initial_set="<<initial_set<<"\n");

  //GTS bounding; bounding.adjoin_inner_approximation(bounding_domain,grid_depth);
  GTS bounding; bounding.adjoin_outer_approximation(bounding_domain,grid_depth);
  ARIADNE_LOG(5,"bounding_size="<<bounding.size()<<"\n");

  GTS initial; initial.adjoin_outer_approximation(initial_set,grid_depth);
  ARIADNE_LOG(5,"initial_size="<<initial.size()<<"\n");

  GTS transient=this->_upper_reach(system,initial,hybrid_transient_time,grid_depth);
  ARIADNE_LOG(5,"transient_size="<<transient.size()<<"\n");
  GTS evolve= this->_upper_evolve(system,initial,hybrid_transient_time,grid_depth);
  ARIADNE_LOG(5,"evolve_size="<<evolve.size()<<"\n");
  GTS reach=evolve;
  GTS found=reach;
  while(!found.empty()) {
    ARIADNE_LOG(5,"found.empty()="<<found.empty()<<"\n");
    reach.adjoin(found);
    ARIADNE_LOG(5,"result_size="<<found.size()<<"\n");
    found=this->_upper_evolve(system,found,hybrid_lock_to_grid_time,grid_depth);
    ARIADNE_LOG(5,"found_size="<<found.size()<<"\n");
    found.remove(reach);
    ARIADNE_LOG(5,"new_size="<<found.size()<<"\n");
    found.restrict(bounding);
    ARIADNE_LOG(5,"bounded_new_size="<<found.size()<<"\n");
  }
  GTS* result=new GTS(reach);
  return result;
}





HybridAnalyser::ConcreteSetType*
HybridAnalyser::
viable(const SystemType& system, 
       const CompactSetType& bounding_set) const
{
  ARIADNE_NOT_IMPLEMENTED;
}



tribool
HybridAnalyser::
verify(const SystemType& system, 
       const LocatedSetType& initial_set, 
       const RegularSetType& safe_set) const
{
  ARIADNE_NOT_IMPLEMENTED;
}












} // namespace Ariadne
