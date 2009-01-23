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

#include "discretiser.h"
#include "reachability_analyser.h"
#include "logging.h"


namespace Ariadne {
  
 
static const double DEFAULT_MAXIMUM_ENCLOSURE_RADIUS=0.25;
static const double DEFAULT_GRID_LENGTH=0.125;

static int verbosity=global_verbosity;

HybridReachabilityAnalyser::
~HybridReachabilityAnalyser()
{
}


HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver)
    : _parameters(new EvolutionParametersType())
    , _discretiser(new HybridDiscretiser<EnclosureType>(evolver))
{
}


HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const EvolutionParametersType& parameters, 
                           const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver)
    : _parameters(new EvolutionParametersType(parameters))
    , _discretiser(new HybridDiscretiser<EnclosureType>(evolver))
{
}






// Helper functions for operators on lists of sets.
HybridGridTreeSet
HybridReachabilityAnalyser::_upper_reach(const HybridAutomaton& sys, 
                                         const HybridGridTreeSet& set, 
                                         const HybridTime& time, 
                                         const int accuracy) const 
{
    HybridGridTreeSet result(set.grid()); 
    HybridGridTreeSet cells=set; 
    cells.mince(accuracy);    
    for(HybridGridTreeSet::const_iterator iter=cells.begin(); iter!=cells.end(); ++iter) {
        result.adjoin(this->_discretiser->upper_evolution(sys,*iter,time,accuracy).reach());
    }
    return result; 
}


HybridGridTreeSet
HybridReachabilityAnalyser::_upper_evolve(const HybridAutomaton& sys, 
                                          const HybridGridTreeSet& set, 
                                          const HybridTime& time, 
                                          const int accuracy) const 
{
    GTS result(set.grid()); GTS cells=set; cells.mince(accuracy); 
    for(HybridGridTreeSet::const_iterator iter=cells.begin(); iter!=cells.end(); ++iter) {
        result.adjoin(this->_discretiser->upper_evolution(sys,*iter,time,accuracy).final()); 
    }
    return result; 
}


std::pair<HybridGridTreeSet,HybridGridTreeSet>
HybridReachabilityAnalyser::_upper_reach_evolve(const HybridAutomaton& sys, 
                                                const HybridGridTreeSet& set, 
                                                const HybridTime& time, 
                                                const int accuracy) const 
{
    std::pair<GTS,GTS> result(GTS(set.grid()),GTS(set.grid()));
    GTS& reach=result.first; GTS& evolve=result.second;
    GTS cells=set; cells.mince(accuracy); 

    /*
      for(HybridGridTreeSet::locations_const_iterator loc_iter=cells.locations_begin(); loc_iter!=cells.locations_end(); ++loc_iter) {
      for(GridTreeSet::const_iterator cell_iter=loc_iter->second.begin(); cell_iter!=loc_iter->second.end(); ++cell_iter) {
      Orbit<HybridGridCell> evolution=this->_discretiser->upper_evolution(sys,HybridGridCell(loc_iter->first,*cell_iter),time,accuracy);
      reach.adjoin(evolution.reach()); 
      evolve.adjoin(evolution.final()); 
      } 
      }
    */

    for(HybridGridTreeSet::const_iterator iter=cells.begin(); iter!=cells.end(); ++iter) {
        Orbit<HybridGridCell> evolution=this->_discretiser->upper_evolution(sys,*iter,time,accuracy);
        reach.adjoin(evolution.reach()); 
        evolve.adjoin(evolution.final()); 
    }
    return result; 
}




HybridGridTreeSet*
HybridReachabilityAnalyser::
lower_evolve(const SystemType& system, 
             const OvertSetInterfaceType& initial_set,
             const TimeType& time) const
{
    int grid_depth = this->_parameters->maximum_grid_depth;
    int grid_height = this->_parameters->maximum_grid_height;
    GTS initial; GTS& final=*new GTS();

    // Improve accuracy of initial set for lower computations
    initial.adjoin_lower_approximation(initial_set,grid_height,grid_depth+4);

    for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
        final.adjoin(this->_discretiser->lower_evolution(system,*bs_iter,time,grid_depth).final());
    }
    return &final;
}



HybridReachabilityAnalyser::SetApproximationType*
HybridReachabilityAnalyser::
lower_reach(const SystemType& system, 
            const OvertSetInterfaceType& initial_set,
            const TimeType& time) const
{
    int initial_grid_depth = this->_parameters->initial_grid_depth;
    int maximum_grid_depth = this->_parameters->maximum_grid_depth;
    int maximum_grid_height = this->_parameters->maximum_grid_height;
    GTS initial; GTS& reach=*new GTS();
  
    // Improve accuracy of initial set for lower computations
    initial.adjoin_lower_approximation(initial_set,maximum_grid_height,initial_grid_depth);
 
    for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
        reach.adjoin(this->_discretiser->lower_evolution(system,*bs_iter,time,maximum_grid_depth).reach());
    }
    return &reach;
}



std::pair<HybridReachabilityAnalyser::SetApproximationType*,HybridReachabilityAnalyser::SetApproximationType*>
HybridReachabilityAnalyser::
lower_reach_evolve(const SystemType& system, 
                   const OvertSetInterfaceType& initial_set,
                   const TimeType& time) const
{
    int initial_grid_depth = this->_parameters->initial_grid_depth;
    int maximum_grid_depth = this->_parameters->maximum_grid_depth;
    int maximum_grid_height = this->_parameters->maximum_grid_height;
    GTS initial; 
  
    GTS& reach=*new GTS; GTS& evolve=*new GTS;

    // Improve accuracy of initial set for lower computations
    initial.adjoin_lower_approximation(initial_set,maximum_grid_height,initial_grid_depth);

    for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
        Orbit<GC> orbit = this->_discretiser->lower_evolution(system,*bs_iter,time,maximum_grid_depth);
        reach.adjoin(orbit.reach());
        evolve.adjoin(orbit.final());
    }
    return make_pair(&reach,&evolve);
}


HybridReachabilityAnalyser::SetApproximationType*
HybridReachabilityAnalyser::
upper_evolve(const SystemType& system, 
             const CompactSetInterfaceType& initial_set,
             const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_evolve(...)\n");
    GTS initial;
    int grid_depth = this->_parameters->maximum_grid_depth;
    initial.adjoin_outer_approximation(initial_set,grid_depth);
    ARIADNE_LOG(4,"initial="<<initial<<"\n");
    GTS& evolve=*new GTS(initial);
    ARIADNE_LOG(4,"initial_evolve="<<evolve<<"\n");
    Float real_time=time.continuous_time;
    uint discrete_steps=time.discrete_time;
    Float lock_to_grid_time=this->_parameters->lock_to_grid_time;
    uint time_steps=uint(real_time/lock_to_grid_time);
    Float remainder_time=real_time-time_steps*lock_to_grid_time;
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);
    for(uint i=0; i!=time_steps; ++i) {
        evolve=this->_upper_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
    }
    evolve=this->_upper_evolve(system,evolve,hybrid_remainder_time,grid_depth);
    evolve.recombine();
    ARIADNE_LOG(4,"final_evolve="<<evolve<<"\n");
    return &evolve;
}



HybridReachabilityAnalyser::SetApproximationType*
HybridReachabilityAnalyser::
upper_reach(const SystemType& system, 
            const CompactSetInterfaceType& initial_set,
            const TimeType& time) const
{
    verbosity=global_verbosity;
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach(system,set,time)\n");
    ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
    GTS initial;
    int grid_depth = this->_parameters->maximum_grid_depth;
    ARIADNE_LOG(3,"grid_depth="<<grid_depth<<"\n");
    initial.adjoin_outer_approximation(initial_set,grid_depth);
    ARIADNE_LOG(4,"initial"<<initial<<"\n");
    GTS evolve=initial;
    ARIADNE_LOG(4,"evolve="<<evolve<<"\n");
    GTS& reach=*new GTS(evolve);
    ARIADNE_LOG(4,"reach="<<reach<<"\n");
    Float real_time=time.continuous_time;
    uint discrete_steps=time.discrete_time;
    Float lock_to_grid_time = this->_parameters->lock_to_grid_time;
    uint time_steps=uint(real_time/lock_to_grid_time);
    Float remainder_time=real_time-time_steps*lock_to_grid_time;
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);
    ARIADNE_LOG(5,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(5,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");
    for(uint i=0; i!=time_steps; ++i) {
        ARIADNE_LOG(5,"computing "<<i<<"-th reachability step...\n");
        reach.adjoin(this->_upper_reach(system,evolve,hybrid_lock_to_grid_time,grid_depth));
        ARIADNE_LOG(5,"reach="<<reach<<"\n");
        reach.recombine();
        ARIADNE_LOG(6,"reach.recombine()="<<reach<<"\n");
        evolve=this->_upper_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
        ARIADNE_LOG(5,"evolve="<<evolve<<"\n");
    }
    ARIADNE_LOG(5,"remainder_time="<<remainder_time<<"\n");
    reach.adjoin(this->_upper_reach(system,evolve,hybrid_remainder_time,grid_depth));
    reach.recombine();
    ARIADNE_LOG(4,"final_reach="<<reach<<"\n");
    return &reach;
}


std::pair<HybridReachabilityAnalyser::SetApproximationType*,HybridReachabilityAnalyser::SetApproximationType*>
HybridReachabilityAnalyser::
upper_reach_evolve(const SystemType& system, 
                   const CompactSetInterfaceType& initial_set,
                   const TimeType& time) const
{
    verbosity=global_verbosity;
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach(system,set,time)\n");
    ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
    GTS initial;
    int grid_depth = this->_parameters->maximum_grid_depth;
    ARIADNE_LOG(3,"grid_depth="<<grid_depth<<"\n");
    initial.adjoin_outer_approximation(initial_set,grid_depth);
    ARIADNE_LOG(4,"initial"<<initial<<"\n");
    GTS& evolve=*new GTS(initial);
    ARIADNE_LOG(4,"evolve="<<evolve<<"\n");
    GTS& reach=*new GTS(initial);
    ARIADNE_LOG(4,"reach="<<reach<<"\n");
    Float real_time=time.continuous_time;
    uint discrete_steps=time.discrete_time;
    Float lock_to_grid_time = this->_parameters->lock_to_grid_time;
    uint time_steps=uint(real_time/lock_to_grid_time);
    Float remainder_time=real_time-time_steps*lock_to_grid_time;
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);
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






HybridReachabilityAnalyser::SetApproximationType*
HybridReachabilityAnalyser::
chain_reach(const SystemType& system, 
            const CompactSetInterfaceType& initial_set) const
{
    HybridTime hybrid_transient_time(this->_parameters->transient_time,this->_parameters->transient_steps);
    HybridTime hybrid_lock_to_grid_time(this->_parameters->lock_to_grid_time,this->_parameters->lock_to_grid_steps);
    int maximum_grid_depth = this->_parameters->maximum_grid_depth;
    int maximum_grid_height = this->_parameters->maximum_grid_height;

    ARIADNE_LOG(5,"initial_set="<<initial_set<<"\n");

    GTS initial; initial.adjoin_outer_approximation(initial_set,maximum_grid_depth);
    ARIADNE_LOG(5,"initial_size="<<initial.size()<<"\n");

    GTS transient,evolve;
    make_lpair(transient,evolve)=this->_upper_reach_evolve(system,initial,hybrid_transient_time,maximum_grid_depth);
    ARIADNE_LOG(5,"transient_size="<<transient.size()<<"\n");
    ARIADNE_LOG(5,"evolve_size="<<evolve.size()<<"\n");
    GTS reach=evolve;
    GTS found=this->_upper_reach(system,reach,hybrid_transient_time,maximum_grid_depth);
    while(!found.empty()) {
        ARIADNE_LOG(5,"found.empty()="<<found.empty()<<"\n");
        reach.adjoin(found);
        ARIADNE_LOG(5,"result_size="<<found.size()<<"\n");
        found=this->_upper_evolve(system,found,hybrid_lock_to_grid_time,maximum_grid_depth);
        ARIADNE_LOG(5,"found_size="<<found.size()<<"\n");
        found.remove(reach);
        ARIADNE_LOG(5,"new_size="<<found.size()<<"\n");
        found.restrict_to_height(maximum_grid_height);
        ARIADNE_LOG(5,"bounded_new_size="<<found.size()<<"\n");
    }
    reach.adjoin(transient);
    reach.recombine();
    GTS* result=new GTS(reach);
    return result;
}


HybridReachabilityAnalyser::SetApproximationType*
HybridReachabilityAnalyser::
chain_reach(const SystemType& system, 
            const CompactSetInterfaceType& initial_set,
            const BoundingSetType& bounding_set) const
{
    // FIXME: Use tree sets throughout

    HybridBoxes bounding_domain=bounding_set;
    HybridTime hybrid_transient_time(this->_parameters->transient_time,this->_parameters->transient_steps);
    HybridTime hybrid_lock_to_grid_time(this->_parameters->lock_to_grid_time,this->_parameters->lock_to_grid_steps);
    int maximum_grid_depth = this->_parameters->maximum_grid_depth;

    ARIADNE_LOG(5,"bounding_domain="<<bounding_domain<<"\n");
    ARIADNE_LOG(5,"initial_set="<<initial_set<<"\n");

    //GTS bounding; bounding.adjoin_inner_approximation(bounding_domain,maximum_grid_depth); 
    GTS bounding; bounding.adjoin_outer_approximation(bounding_domain,maximum_grid_depth); bounding.recombine();
    ARIADNE_LOG(5,"bounding_size="<<bounding.size()<<"\n");

    GTS initial; initial.adjoin_outer_approximation(initial_set,maximum_grid_depth);
    ARIADNE_LOG(5,"initial_size="<<initial.size()<<"\n");

    GTS transient,evolve;
    make_lpair(transient,evolve)=this->_upper_reach_evolve(system,initial,hybrid_transient_time,maximum_grid_depth);
    ARIADNE_LOG(5,"transient_size="<<transient.size()<<"\n");
    ARIADNE_LOG(5,"evolve_size="<<evolve.size()<<"\n");
    GTS reach=evolve;
    GTS found=this->_upper_reach(system,reach,hybrid_transient_time,maximum_grid_depth);
    while(!found.empty()) {
        ARIADNE_LOG(5,"found.empty()="<<found.empty()<<"\n");
        reach.adjoin(found);
        ARIADNE_LOG(5,"result_size="<<found.size()<<"\n");
        found=this->_upper_evolve(system,found,hybrid_lock_to_grid_time,maximum_grid_depth);
        ARIADNE_LOG(5,"found_size="<<found.size()<<"\n");
        found.remove(reach);
        ARIADNE_LOG(5,"new_size="<<found.size()<<"\n");
        found.restrict(bounding);
        ARIADNE_LOG(5,"bounded_new_size="<<found.size()<<"\n");
    }
    reach.adjoin(transient);
    reach.recombine();
    GTS* result=new GTS(reach);
    return result;
}





HybridReachabilityAnalyser::SetApproximationType*
HybridReachabilityAnalyser::
viable(const SystemType& system, 
       const CompactSetInterfaceType& bounding_set) const
{
    ARIADNE_NOT_IMPLEMENTED;
}



tribool
HybridReachabilityAnalyser::
verify(const SystemType& system, 
       const LocatedSetInterfaceType& initial_set, 
       const RegularSetInterfaceType& safe_set) const
{
    ARIADNE_NOT_IMPLEMENTED;
}












} // namespace Ariadne
