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
  for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
    final.adjoin(this->_discretiser->lower_evolve(system,*bs_iter,time).final());
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
  for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
    reach.adjoin(this->_discretiser->lower_evolve(system,*bs_iter,time).intermediate());
  }
  return new GTS(reach);
}



HybridAnalyser::ConcreteSetType*
HybridAnalyser::
upper_evolve(const SystemType& system, 
             const CompactSetType& initial_set,
             const TimeType& time) const
{
  GTS initial;
  initial.adjoin_outer_approximation(initial_set);
  GCLS evolve=initial;
  Float real_time=time.first;
  Float lock_to_grid_time=this->_parameters->lock_to_grid_time;
  uint time_steps=real_time/lock_to_grid_time;
  Float remainder_time=real_time-time_steps*lock_to_grid_time;
  HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,time_steps);
  for(uint i=0; i!=time_steps; ++i) {
    //FIXME!! evolve=this->_discretiser->upper_evolve(system,evolve,hybrid_lock_to_grid_time).final();
  }
  //FIXME!! evolve=this->_upper_evolve(system,evolve,remainder_time).final();
  return new GTS(evolve);
}



HybridAnalyser::ConcreteSetType*
HybridAnalyser::
upper_reach(const SystemType& system, 
            const CompactSetType& initial_set,
            const TimeType& time) const
{
  ARIADNE_LOG(2,"HybridAnalyser::upper_reach(system,set,time)\n");
  ARIADNE_LOG(3,"intial_set="<<initial_set<<"\n");
  GTS initial;
  initial.adjoin_outer_approximation(initial_set);
  GCLS evolve=initial;
  GCLS reach=evolve;
  Float real_time=time.first;
  Float lock_to_grid_time = this->_parameters->lock_to_grid_time;
  uint time_steps=real_time/lock_to_grid_time;
  Float remainder_time=real_time-time_steps*lock_to_grid_time;
  HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,time_steps);
  for(uint i=0; i!=time_steps; ++i) {
    //FIXME!! evolve=this->_upper_evolve(system,evolve,lock_to_grid_time);
    //FIXME!! reach.adjoin(this->_upper_reach(system,evolve,lock_to_grid_time));
  }
  //FIXME!! reach.adjoin(this->_upper_reach(system,evolve,remainder_time));
  return new GTS(reach);
}





HybridAnalyser::ConcreteSetType*
HybridAnalyser::
chain_reach(const SystemType& system, 
            const CompactSetType& initial_set,
            const BoundingSetType& bounding_domain) const
{
  uint verbosity=0;
  ARIADNE_LOG(5,"bounding_domain="<<bounding_domain);
  Float time=this->_parameters->lock_to_grid_time;
  GTS initial;
  initial.adjoin_outer_approximation(initial_set);
  GTS* result=new GTS;
  ARIADNE_LOG(5,"initial_set="<<initial_set);
  GCLS found=initial;
  ARIADNE_LOG(5,"initial_size="<<found.size());
  //FIXME!  found=this->_upper_reach(system,found,time);
  ARIADNE_LOG(5,"reach_size="<<found.size());
  while(!found.empty()) {
    result->adjoin(found);
    ARIADNE_LOG(5,"result_size="<<found.size());
    //FIXME!  found=this->_upper_evolve(system,found,time);
    ARIADNE_LOG(5,"found_size="<<found.size());
    found.unique_sort();
    ARIADNE_LOG(5,"unique_size="<<found.size());
    found.remove(*result);
    ARIADNE_LOG(5,"new_size="<<found.size());
    found.restrict(bounding_domain);
    ARIADNE_LOG(5,"bounded_new_size="<<found.size());
  }
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
