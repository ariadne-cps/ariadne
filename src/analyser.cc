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

HybridAnalyser::ConcreteSetType*
HybridAnalyser::
lower_evolve(const SystemType& system, 
             const OvertSetType& initial_set,
             const TimeType& time) const
{
  GTS initial; GCLS final;
  initial.adjoin_lower_approximation(initial_set);
  int grid_depth = this->_parameters->grid_depth;
  for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
    final.adjoin(this->_discretiser->lower_evolve(system,*bs_iter,time,grid_depth).final());
  }
  return new GTS(final);
}



HybridAnalyser::ConcreteSetType*
HybridAnalyser::
lower_reach(const SystemType& system, 
            const OvertSetType& initial_set,
            const TimeType& time) const
{
  GTS initial; GCLS reach;
  initial.adjoin_lower_approximation(initial_set);
  int grid_depth = this->_parameters->grid_depth;
  for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
    reach.adjoin(this->_discretiser->lower_evolve(system,*bs_iter,time,grid_depth).intermediate());
  }
  return new GTS(reach);
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
  GCLS evolve=initial;
  ARIADNE_LOG(4,"initial_evolve="<<evolve<<"\n");
  Float real_time=time.first;
  Float lock_to_grid_time=this->_parameters->lock_to_grid_time;
  uint time_steps=real_time/lock_to_grid_time;
  Float remainder_time=real_time-time_steps*lock_to_grid_time;
  HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,time_steps);
  for(uint i=0; i!=time_steps; ++i) {
    //FIXME!! evolve=this->_discretiser->upper_evolve(system,evolve,hybrid_lock_to_grid_time).final();
  }
  //FIXME!! evolve=this->_upper_evolve(system,evolve,remainder_time).final();
  ARIADNE_LOG(4,"final_evolve="<<evolve<<"\n");
  ConcreteSetType* final_ptr = new GTS(evolve);
  ARIADNE_LOG(4,"final="<<*final_ptr<<"\n");
  return final_ptr;
}



HybridAnalyser::ConcreteSetType*
HybridAnalyser::
upper_reach(const SystemType& system, 
            const CompactSetType& initial_set,
            const TimeType& time) const
{
  verbosity=0;
  ARIADNE_LOG(2,"HybridAnalyser::upper_reach(system,set,time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  GTS initial;
  int grid_depth = this->_parameters->grid_depth;
  ARIADNE_LOG(3,"grid_depth="<<grid_depth<<"\n");
  initial.adjoin_outer_approximation(initial_set,grid_depth);
  ARIADNE_LOG(4,"initial"<<initial<<"\n");
  GCLS evolve=initial;
  ARIADNE_LOG(4,"evolve="<<evolve<<"\n");
  GCLS reach=evolve;
  ARIADNE_LOG(4,"reach="<<reach<<"\n");
  Float real_time=time.first;
  Float lock_to_grid_time = this->_parameters->lock_to_grid_time;
  uint time_steps=real_time/lock_to_grid_time;
  Float remainder_time=real_time-time_steps*lock_to_grid_time;
  HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,time_steps);
  HybridTime hybrid_remainder_time(remainder_time,time_steps);
  for(uint i=0; i!=time_steps; ++i) {
    reach.adjoin(this->_upper_reach(system,evolve,hybrid_lock_to_grid_time,grid_depth));
    evolve=this->_upper_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
  }
  reach.adjoin(this->_upper_reach(system,evolve,hybrid_remainder_time,grid_depth));
  ARIADNE_LOG(4,"reach="<<reach<<"\n");
  GTS* final_ptr = new GTS(reach);
  ARIADNE_LOG(4,"final="<<*final_ptr<<"\n");
  return final_ptr;
}





HybridAnalyser::ConcreteSetType*
HybridAnalyser::
chain_reach(const SystemType& system, 
            const CompactSetType& initial_set,
            const BoundingSetType& bounding_domain) const
{
  // FIXME: Use tree sets throughout
  uint verbosity=0;
  ARIADNE_LOG(5,"bounding_domain="<<bounding_domain<<"\n");
  Float lock_to_grid_time = this->_parameters->lock_to_grid_time;
  uint lock_to_grid_steps = this->_parameters->lock_to_grid_steps;
  HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,lock_to_grid_steps);
  int grid_depth = this->_parameters->grid_depth;
  ARIADNE_LOG(5,"initial_set="<<initial_set<<"\n");
  GTS initial_tree;
  initial_tree.adjoin_outer_approximation(initial_set,grid_depth);
  GTS bounding_tree;
  //FIXME: Need inner approximation
  initial_tree.adjoin_outer_approximation(bounding_domain,grid_depth);
  //FIXME: Need this for multiple cells
  initial_tree.locations_begin()->second.mince(grid_depth);
  ARIADNE_LOG(5,"initial_tree_size="<<initial_tree.size()<<"\n");
  GCLS initial(initial_tree);
  GCLS reach=initial;
  GCLS found=initial;
  ARIADNE_LOG(5,"initial_size="<<initial.size()<<"\n");
  found=this->_upper_reach(system,found,hybrid_lock_to_grid_time,grid_depth);
  ARIADNE_LOG(5,"reach_size="<<found.size()<<"\n");
  while(!found.empty()) {
    ARIADNE_LOG(5,"found.empty()="<<found.empty()<<"\n");
    reach.adjoin(found);
    ARIADNE_LOG(5,"result_size="<<found.size()<<"\n");
    found=this->_upper_evolve(system,found,hybrid_lock_to_grid_time,grid_depth);
    ARIADNE_LOG(5,"found_size="<<found.size()<<"\n");
    found.unique_sort();
    ARIADNE_LOG(5,"unique_size="<<found.size()<<"\n");
    found.remove(reach);
    ARIADNE_LOG(5,"new_size="<<found.size()<<"\n");
    found.restrict(bounding_domain);
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
