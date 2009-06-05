/***************************************************************************
 *            reachability_analyser.cc
 *
 *  Copyright  2006-9  Alberto Casagrande, Pieter Collins, Davide Bresolin
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

#include "hybrid_time.h"
#include "hybrid_automaton.h"

#include "evolution_parameters.h"
#include "evolver_interface.h"

#include "discretiser.h"
#include "reachability_analyser.h"
#include "logging.h"

#include "graphics.h"


namespace Ariadne {

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
    HybridGrid grid=sys.grid();
    HybridGridTreeSet result(grid);
    HybridGridTreeSet cells=set;
    cells.mince(accuracy);
    for(HybridGridTreeSet::const_iterator iter=cells.begin(); iter!=cells.end(); ++iter) {
        result.adjoin(this->_discretiser->reach(sys,*iter,time,grid,accuracy,UPPER_SEMANTICS));
    }
    return result;
}


HybridGridTreeSet
HybridReachabilityAnalyser::_upper_evolve(const HybridAutomaton& sys,
                                          const HybridGridTreeSet& set,
                                          const HybridTime& time,
                                          const int accuracy) const
{
    HybridGrid grid=sys.grid();
    GTS result(grid); GTS cells=set; cells.mince(accuracy);
    for(HybridGridTreeSet::const_iterator iter=cells.begin(); iter!=cells.end(); ++iter) {
        ARIADNE_LOG(5,"Evolving cell = "<<*iter<<"\n");
        result.adjoin(this->_discretiser->evolve(sys,*iter,time,grid,accuracy,UPPER_SEMANTICS));
    }
    ARIADNE_LOG(4,"_upper_evolve result size = "<<result.size()<<"\n");
    return result;
}


std::pair<HybridGridTreeSet,HybridGridTreeSet>
HybridReachabilityAnalyser::_upper_reach_evolve(const HybridAutomaton& sys,
                                                const HybridGridTreeSet& set,
                                                const HybridTime& time,
                                                const int accuracy) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::_upper_reach_evolve(...)\n");
    HybridGrid grid=sys.grid();
    std::pair<GTS,GTS> result=make_pair(GTS(grid),GTS(grid));
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
        Orbit<HybridGridCell> evolution=this->_discretiser->evolution(sys,*iter,time,grid,accuracy,UPPER_SEMANTICS);
        ARIADNE_LOG(4,"  evolution reach size= "<<evolution.reach().size()<<"\n");
        ARIADNE_LOG(4,"  evolution final size= "<<evolution.final().size()<<"\n");
        reach.adjoin(evolution.reach());
        evolve.adjoin(evolution.final());
    }
    ARIADNE_LOG(3,"  final reach size = "<<reach.size()<<"\n");
    ARIADNE_LOG(3,"  final evolve size = "<<evolve.size()<<"\n");
    ARIADNE_LOG(2,"Done.\n");
    return result;
}




HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_evolve(const SystemType& system,
             const OvertSetInterfaceType& initial_set,
             const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_evolve(...)\n");
    int grid_depth = this->_parameters->maximum_grid_depth;
    int grid_height = this->_parameters->maximum_grid_height;
    HybridGrid grid=system.grid();
    GTS initial; GTS final;

    // Improve accuracy of initial set for lower computations
    initial.adjoin_lower_approximation(initial_set,grid_height,grid_depth+4);
    ARIADNE_LOG(3,"initial.size()="<<initial.size()<<"\n");
    ARIADNE_LOG(3,"computing lower evolution.");
    for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
        ARIADNE_LOG(3,".");
        GC cell=*bs_iter;
        GTS cell_final=this->_discretiser->evolve(system,cell,time,grid,grid_depth,LOWER_SEMANTICS);
        final.adjoin(cell_final);
    }
    ARIADNE_LOG(3,"\n");
    return final;
}

HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_evolve(const SystemType& system,
             const HybridEnclosureType& initial_set,
             const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_evolve(...)\n");
    int grid_depth = this->_parameters->maximum_grid_depth;
    HybridGrid grid=system.grid();

    ARIADNE_LOG(3,"computing lower evolution...");
    GTS reach=this->_discretiser->evolve(system,initial_set,time,grid,grid_depth,LOWER_SEMANTICS);
    ARIADNE_LOG(3,"\n");
    return reach;
}



HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_reach(const SystemType& system,
            const OvertSetInterfaceType& initial_set,
            const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_reach(...)\n");
    int grid_depth = this->_parameters->maximum_grid_depth;
    int grid_height = this->_parameters->maximum_grid_height;
    HybridGrid grid=system.grid();
    GTS initial; GTS reach;

    ARIADNE_LOG(3,"Adjoining initial set to the grid...\n");
    // Improve accuracy of initial set for lower computations
    initial.adjoin_lower_approximation(initial_set,grid_height,grid_depth+4);
    ARIADNE_LOG(3,"initial.size()="<<initial.size()<<"\n");
    ARIADNE_LOG(3,"Computing lower reach set...");
    for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
        ARIADNE_LOG(3,".");
        GC cell=*bs_iter;
        GTS cell_reach=this->_discretiser->reach(system,cell,time,grid,grid_depth,LOWER_SEMANTICS);
        reach.adjoin(cell_reach);
    }
    ARIADNE_LOG(3,"\n");
    return reach;
}

HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_reach(const SystemType& system,
            const HybridEnclosureType& initial_set,
            const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_reach(...)\n");
    int grid_depth = this->_parameters->maximum_grid_depth;
    HybridGrid grid=system.grid();

    ARIADNE_LOG(3,"Computing lower reach set...");
    GTS reach=this->_discretiser->reach(system,initial_set,time,grid,grid_depth,LOWER_SEMANTICS);
    ARIADNE_LOG(3,"\n");
    return reach;
}


std::pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
lower_reach_evolve(const SystemType& system,
                   const OvertSetInterfaceType& initial_set,
                   const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_reach_evolve(...)\n");
    int grid_depth = this->_parameters->maximum_grid_depth;
    int grid_height = this->_parameters->maximum_grid_height;

    HybridGrid grid=system.grid();

    GTS initial;

    GTS reach=(grid); GTS evolve(grid);

    // Improve accuracy of initial set for lower computations
    initial.adjoin_lower_approximation(initial_set,grid_height,grid_depth+4);
    ARIADNE_LOG(3,"initial.size()="<<initial.size()<<"\n");
    ARIADNE_LOG(3,"computing lower evolution.");
    for(GTS::const_iterator bs_iter=initial.begin(); bs_iter!=initial.end(); ++bs_iter) {
        ARIADNE_LOG(3,".");
        Orbit<GC> orbit = this->_discretiser->evolution(system,*bs_iter,time,grid,grid_depth,LOWER_SEMANTICS);
        reach.adjoin(orbit.reach());
        evolve.adjoin(orbit.final());
    }
    ARIADNE_LOG(3,"\n");
    return make_pair(reach,evolve);
}

std::pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
lower_reach_evolve(const SystemType& system,
                   const HybridEnclosureType& initial_set,
                   const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_reach_evolve(...)\n");
    int grid_depth = this->_parameters->maximum_grid_depth;

    HybridGrid grid=system.grid();

    ARIADNE_LOG(3,"computing lower evolution.");
    Orbit<GC> orbit = this->_discretiser->evolution(system,initial_set,time,grid,grid_depth,LOWER_SEMANTICS);
    GTS reach(orbit.reach());
    GTS evolve(orbit.final());
    return make_pair(reach,evolve);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_evolve(const SystemType& system,
             const CompactSetInterfaceType& initial_set,
             const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_evolve(...)\n");
    HybridGrid grid=system.grid();
    GTS evolve = *new GTS(grid);
    int grid_depth = this->_parameters->maximum_grid_depth;
    evolve.adjoin_outer_approximation(initial_set,grid_depth);
    ARIADNE_LOG(4,"initial_evolve.size()="<<evolve.size()<<"\n");
    Float real_time=time.continuous_time;
    uint discrete_steps=time.discrete_time;
    Float lock_to_grid_time=this->_parameters->lock_to_grid_time;
    uint time_steps=uint(real_time/lock_to_grid_time);
    Float remainder_time=real_time-time_steps*lock_to_grid_time;
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);
    ARIADNE_LOG(3,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(3,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");
    for(uint i=0; i!=time_steps; ++i) {
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...\n");
        evolve=this->_upper_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
    }
    ARIADNE_LOG(3,"remainder_time="<<remainder_time<<"\n");
    if(!evolve.empty() && remainder_time > 0) {
        ARIADNE_LOG(3,"computing evolution for remainder time...\n");
        evolve=this->_upper_evolve(system,evolve,hybrid_remainder_time,grid_depth);
    }
    evolve.recombine();
    ARIADNE_LOG(4,"final_evolve.size()="<<evolve.size()<<"\n");
    return evolve;
}

HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_evolve(const SystemType& system,
             const HybridEnclosureType& initial_set,
             const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_evolve(...)\n");
    HybridGrid grid=system.grid();
    GTS evolve = *new GTS(grid);
    int grid_depth = this->_parameters->maximum_grid_depth;
    Float real_time=time.continuous_time;
    uint discrete_steps=time.discrete_time;
    Float lock_to_grid_time=this->_parameters->lock_to_grid_time;
    uint time_steps=uint(real_time/lock_to_grid_time);
    Float remainder_time=real_time-time_steps*lock_to_grid_time;
    if(time_steps <= 0) {       // evolution time is smaller than lock_to_grid_time
        lock_to_grid_time = remainder_time;
        remainder_time = 0.0;
    }
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);
    ARIADNE_LOG(3,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(3,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");
    ARIADNE_LOG(3,"computing first reachability step...\n");
    evolve = this->_discretiser->evolve(system,initial_set,hybrid_lock_to_grid_time,grid,grid_depth,UPPER_SEMANTICS);
    for(uint i=1; i<time_steps; ++i) {
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...\n");
        evolve=this->_upper_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
    }
    ARIADNE_LOG(3,"remainder_time="<<remainder_time<<"\n");
    if(!evolve.empty() && remainder_time > 0.0) {
        ARIADNE_LOG(3,"computing evolution for remainder time...\n");
        evolve=this->_upper_evolve(system,evolve,hybrid_remainder_time,grid_depth);
    }
    evolve.recombine();
    ARIADNE_LOG(4,"final_evolve.size()="<<evolve.size()<<"\n");
    return evolve;
}



HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_reach(const SystemType& system,
            const CompactSetInterfaceType& initial_set,
            const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach(system,set,time)\n");
    ARIADNE_LOG(4,"initial_set="<<initial_set<<"\n");
    HybridGrid grid=system.grid();
    GTS evolve(grid);
    int grid_depth = this->_parameters->maximum_grid_depth;
    ARIADNE_LOG(4,"grid_depth="<<grid_depth<<"\n");
    evolve.adjoin_outer_approximation(initial_set,grid_depth);
    ARIADNE_LOG(4,"initial size = "<<evolve.size()<<"\n");
    GTS reach(evolve);
    ARIADNE_LOG(4,"reach size ="<<reach.size()<<"\n");
    Float real_time=time.continuous_time;
    uint discrete_steps=time.discrete_time;
    Float lock_to_grid_time = this->_parameters->lock_to_grid_time;
    uint time_steps=uint(real_time/lock_to_grid_time);
    Float remainder_time=real_time-time_steps*lock_to_grid_time;
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);
    ARIADNE_LOG(3,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(3,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");
    ARIADNE_LOG(3,"discrete_steps="<<discrete_steps<<"\n");
    GTS found;
    for(uint i=0; i!=time_steps; ++i) {
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...\n");
        make_lpair(found,evolve)=this->_upper_reach_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
        ARIADNE_LOG(5,"found.size()="<<found.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve.size()<<"\n");
        evolve.remove(reach);
        reach.adjoin(found);
        ARIADNE_LOG(3,"  found "<<found.size()<<" cells, of which "<<evolve.size()<<" are new.\n");
        if(evolve.empty()) break;
        ARIADNE_LOG(6,"evolve="<<evolve<<"\n");
    }
    ARIADNE_LOG(3,"remainder_time="<<remainder_time<<"\n");
    if(!evolve.empty() && remainder_time > 0) {
        ARIADNE_LOG(3,"computing evolution for remainder time...\n");
        reach.adjoin(this->_upper_reach(system,evolve,hybrid_remainder_time,grid_depth));
    }
    reach.recombine();
    ARIADNE_LOG(4,"final_reach size = "<<reach.size()<<"\n");
    return reach;
}

HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_reach(const SystemType& system,
            const HybridEnclosureType& initial_set,
            const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach(system,set,time)\n");
    ARIADNE_LOG(4,"initial_set="<<initial_set<<"\n");
    HybridGrid grid=system.grid();
    int grid_depth = this->_parameters->maximum_grid_depth;
    Float real_time=time.continuous_time;
    uint discrete_steps=time.discrete_time;
    Float lock_to_grid_time = this->_parameters->lock_to_grid_time;
    uint time_steps=uint(real_time/lock_to_grid_time);
    Float remainder_time=real_time-time_steps*lock_to_grid_time;
    if(time_steps <= 0) {       // evolution time is smaller than lock_to_grid_time
        lock_to_grid_time = remainder_time;
        remainder_time = 0.0;
    }
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);
    ARIADNE_LOG(3,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(3,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");
    ARIADNE_LOG(3,"discrete_steps="<<discrete_steps<<"\n");
    ARIADNE_LOG(3,"computing first reachability step...\n");
    Orbit<GC> orbit = this->_discretiser->evolution(system,initial_set,time,grid,grid_depth,UPPER_SEMANTICS);
    GTS reach(orbit.reach());
    GTS evolve(orbit.final());
    GTS found;
    for(uint i=1; i<time_steps; ++i) {
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...\n");
        make_lpair(found,evolve)=this->_upper_reach_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
        ARIADNE_LOG(5,"found.size()="<<found.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve.size()<<"\n");
        evolve.remove(reach);
        reach.adjoin(found);
        ARIADNE_LOG(3,"  found "<<found.size()<<" cells, of which "<<evolve.size()<<" are new.\n");
        if(evolve.empty()) break;
        ARIADNE_LOG(6,"evolve="<<evolve<<"\n");
    }
    ARIADNE_LOG(3,"remainder_time="<<remainder_time<<"\n");
    if(!evolve.empty() && remainder_time > 0) {
        ARIADNE_LOG(3,"computing evolution for remainder time...\n");
        reach.adjoin(this->_upper_reach(system,evolve,hybrid_remainder_time,grid_depth));
    }
    reach.recombine();
    ARIADNE_LOG(4,"final_reach size = "<<reach.size()<<"\n");
    return reach;
}




std::pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
upper_reach_evolve(const SystemType& system,
                   const CompactSetInterfaceType& initial_set,
                   const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach_evolve(system,set,time)\n");
    ARIADNE_LOG(4,"initial_set="<<initial_set<<"\n");
    HybridGrid grid=system.grid();
    GTS evolve(grid);
    int grid_depth = this->_parameters->maximum_grid_depth;
    ARIADNE_LOG(4,"grid_depth="<<grid_depth<<"\n");
    evolve.adjoin_outer_approximation(initial_set,grid_depth);
    ARIADNE_LOG(4,"initial_evolve"<<evolve<<"\n");
    GTS reach(evolve);
    ARIADNE_LOG(4,"reach="<<reach<<"\n");
    Float real_time=time.continuous_time;
    uint discrete_steps=time.discrete_time;
    Float lock_to_grid_time = this->_parameters->lock_to_grid_time;
    uint time_steps=uint(real_time/lock_to_grid_time);
    Float remainder_time=real_time-time_steps*lock_to_grid_time;
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);
    ARIADNE_LOG(3,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(3,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");
    GTS found(grid);
    for(uint i=0; i!=time_steps; ++i) {
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...\n");
        make_lpair(found,evolve) = this->_upper_reach_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
        ARIADNE_LOG(5,"found.size()="<<found.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve.size()<<"\n");
        reach.adjoin(found);
        ARIADNE_LOG(3,"  found "<<found.size()<<" cells.\n");
    }
    ARIADNE_LOG(3,"remainder_time="<<remainder_time<<"\n");
    if(!evolve.empty() && remainder_time > 0) {
        ARIADNE_LOG(3,"computing evolution for remainder time...\n");
        make_lpair(found,evolve) = this->_upper_reach_evolve(system,evolve,hybrid_remainder_time,grid_depth);
        reach.adjoin(found);
    }
    reach.recombine();
    ARIADNE_LOG(4,"reach="<<reach<<"\n");
    evolve.recombine();
    ARIADNE_LOG(4,"evolve="<<evolve<<"\n");
    return std::make_pair(reach,evolve);
}

std::pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
upper_reach_evolve(const SystemType& system,
            const HybridEnclosureType& initial_set,
            const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach(system,set,time)\n");
    ARIADNE_LOG(4,"initial_set="<<initial_set<<"\n");
    HybridGrid grid=system.grid();
    int grid_depth = this->_parameters->maximum_grid_depth;
    Float real_time=time.continuous_time;
    uint discrete_steps=time.discrete_time;
    Float lock_to_grid_time = this->_parameters->lock_to_grid_time;
    uint time_steps=uint(real_time/lock_to_grid_time);
    Float remainder_time=real_time-time_steps*lock_to_grid_time;
    if(time_steps <= 0) {       // evolution time is smaller than lock_to_grid_time
        lock_to_grid_time = remainder_time;
        remainder_time = 0.0;
    }
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);
    ARIADNE_LOG(3,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(3,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");
    ARIADNE_LOG(3,"discrete_steps="<<discrete_steps<<"\n");
    ARIADNE_LOG(3,"computing first reachability step...\n");
    Orbit<GC> orbit = this->_discretiser->evolution(system,initial_set,time,grid,grid_depth,UPPER_SEMANTICS);
    GTS reach(orbit.reach());
    GTS evolve(orbit.final());
    GTS found;
    for(uint i=1; i<time_steps; ++i) {
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...\n");
        make_lpair(found,evolve)=this->_upper_reach_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
        ARIADNE_LOG(5,"found.size()="<<found.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve.size()<<"\n");
        evolve.remove(reach);
        reach.adjoin(found);
        ARIADNE_LOG(3,"  found "<<found.size()<<" cells, of which "<<evolve.size()<<" are new.\n");
        if(evolve.empty()) break;
        ARIADNE_LOG(6,"evolve="<<evolve<<"\n");
    }
    ARIADNE_LOG(3,"remainder_time="<<remainder_time<<"\n");
    if(!evolve.empty() && remainder_time > 0) {
        ARIADNE_LOG(3,"computing evolution for remainder time...\n");
        make_lpair(found,evolve) = this->_upper_reach_evolve(system,evolve,hybrid_remainder_time,grid_depth);
        reach.adjoin(found);
    }
    reach.recombine();
    ARIADNE_LOG(4,"reach="<<reach<<"\n");
    evolve.recombine();
    ARIADNE_LOG(4,"evolve="<<evolve<<"\n");
    return std::make_pair(reach,evolve);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
chain_reach(const SystemType& system,
            const CompactSetInterfaceType& initial_set) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::chain_reach(system,set)\n");
    Float transient_time = this->_parameters->transient_time;
    int transient_steps = this->_parameters->transient_steps;
    HybridTime hybrid_transient_time(transient_time, transient_steps);
    Float lock_to_grid_time=this->_parameters->lock_to_grid_time;
    int lock_to_grid_steps=this->_parameters->lock_to_grid_steps;
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,lock_to_grid_steps);
    int maximum_grid_depth = this->_parameters->maximum_grid_depth;
    int maximum_grid_height = this->_parameters->maximum_grid_height;
    ARIADNE_LOG(3,"transient_time=("<<transient_time<<","<<transient_steps<<")\n");
    ARIADNE_LOG(3,"lock_to_grid_time=("<<lock_to_grid_time<<","<<lock_to_grid_steps<<")\n");

    ARIADNE_LOG(5,"initial_set="<<initial_set<<"\n");

    HybridGrid grid=system.grid();

    GTS evolve(grid);
    evolve.adjoin_outer_approximation(initial_set,maximum_grid_depth);
    ARIADNE_LOG(5,"initial_size="<<evolve.size()<<"\n");

    GTS reach = GTS(grid);
    if(transient_time > 0.0 || transient_steps > 0) {
        ARIADNE_LOG(3,"Computing transient evolution...\n");
        make_lpair(reach,evolve)=this->_upper_reach_evolve(system,evolve,hybrid_transient_time,maximum_grid_depth);
        ARIADNE_LOG(5,"transient_reach_size="<<reach.size()<<"\n");
        ARIADNE_LOG(5,"evolve_size="<<evolve.size()<<"\n");
        ARIADNE_LOG(3,"  found "<<reach.size()<<" cells.\n");
    }
    GTS found(grid);
    ARIADNE_LOG(3,"Computing recurrent evolution...\n");
    while(!evolve.empty()) {
        make_lpair(found,evolve)=this->_upper_reach_evolve(system,evolve,hybrid_lock_to_grid_time,maximum_grid_depth);
        ARIADNE_LOG(5,"found.size()="<<found.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve.size()<<"\n");
        evolve.remove(reach);
        found.restrict_to_height(maximum_grid_height);
        reach.adjoin(found);
        ARIADNE_LOG(3,"  found "<<found.size()<<" cells, of which "<<evolve.size()<<" are new.\n");
        // ARIADNE_LOG(5,"bounded_new_size="<<found.size()<<"\n");
    }
    reach.recombine();
    return reach;
}

HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
chain_reach(const SystemType& system,
            const HybridEnclosureType& initial_set) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::chain_reach(system,set)\n");
    Float transient_time = this->_parameters->transient_time;
    int transient_steps = this->_parameters->transient_steps;
    HybridTime hybrid_transient_time(transient_time, transient_steps);
    Float lock_to_grid_time=this->_parameters->lock_to_grid_time;
    int lock_to_grid_steps=this->_parameters->lock_to_grid_steps;
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,lock_to_grid_steps);
    int maximum_grid_depth = this->_parameters->maximum_grid_depth;
    int maximum_grid_height = this->_parameters->maximum_grid_height;
    ARIADNE_LOG(3,"transient_time=("<<transient_time<<","<<transient_steps<<")\n");
    ARIADNE_LOG(3,"lock_to_grid_time=("<<lock_to_grid_time<<","<<lock_to_grid_steps<<")\n");

    ARIADNE_LOG(5,"initial_set="<<initial_set<<"\n");

    HybridGrid grid=system.grid();

    GTS evolve(grid);
    GTS reach(grid);
    
    if(transient_time > 0.0 || transient_steps > 0) {
        ARIADNE_LOG(3,"Computing transient evolution...\n");
        Orbit<GC> orbit = this->_discretiser->evolution(system,initial_set,hybrid_transient_time,grid,maximum_grid_depth,UPPER_SEMANTICS);
        reach.adjoin(orbit.reach());
        evolve.adjoin(orbit.final());
        ARIADNE_LOG(5,"transient_reach_size="<<reach.size()<<"\n");
        ARIADNE_LOG(5,"evolve_size="<<evolve.size()<<"\n");
        ARIADNE_LOG(3,"  found "<<reach.size()<<" cells.\n");
    } else {
        ARIADNE_LOG(3,"Computing first step of recurrent evolution...\n");
        Orbit<GC> orbit = this->_discretiser->evolution(system,initial_set,hybrid_lock_to_grid_time,grid,maximum_grid_depth,UPPER_SEMANTICS);
        reach.adjoin(orbit.reach());
        evolve.adjoin(orbit.final());
        ARIADNE_LOG(5,"transient_reach_size="<<reach.size()<<"\n");
        ARIADNE_LOG(5,"evolve_size="<<evolve.size()<<"\n");
        ARIADNE_LOG(3,"  found "<<reach.size()<<" cells.\n");
    }    
    GTS found(grid);
    ARIADNE_LOG(3,"Computing recurrent evolution...\n");
    while(!evolve.empty()) {
        make_lpair(found,evolve)=this->_upper_reach_evolve(system,evolve,hybrid_lock_to_grid_time,maximum_grid_depth);
        ARIADNE_LOG(5,"found.size()="<<found.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve.size()<<"\n");
        evolve.remove(reach);
        found.restrict_to_height(maximum_grid_height);
        reach.adjoin(found);
        ARIADNE_LOG(3,"  found "<<found.size()<<" cells, of which "<<evolve.size()<<" are new.\n");
        // ARIADNE_LOG(5,"bounded_new_size="<<found.size()<<"\n");
    }
    reach.recombine();
    return reach;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
chain_reach(const SystemType& system,
            const CompactSetInterfaceType& initial_set,
            const BoundingSetType& bounding_set) const
{
    // FIXME: Use tree sets throughout

    HybridBoxes bounding_domain=bounding_set;
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::chain_reach(system,set,bounding_set)\n");
    Float transient_time = this->_parameters->transient_time;
    int transient_steps = this->_parameters->transient_steps;
    HybridTime hybrid_transient_time(transient_time, transient_steps);
    Float lock_to_grid_time=this->_parameters->lock_to_grid_time;
    int lock_to_grid_steps=this->_parameters->lock_to_grid_steps;
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,lock_to_grid_steps);
    int maximum_grid_depth = this->_parameters->maximum_grid_depth;
    ARIADNE_LOG(3,"transient_time=("<<transient_time<<","<<transient_steps<<")\n");
    ARIADNE_LOG(3,"lock_to_grid_time=("<<lock_to_grid_time<<","<<lock_to_grid_steps<<")\n");

    ARIADNE_LOG(5,"bounding_domain="<<bounding_domain<<"\n");
    ARIADNE_LOG(5,"initial_set="<<initial_set<<"\n");

    HybridGrid grid=system.grid();
    //GTS bounding; bounding.adjoin_inner_approximation(bounding_domain,maximum_grid_depth);
    GTS bounding(grid);
    bounding.adjoin_outer_approximation(bounding_domain,maximum_grid_depth); bounding.recombine();
    ARIADNE_LOG(5,"bounding_size="<<bounding.size()<<"\n");

    GTS evolve(grid);
    evolve.adjoin_outer_approximation(initial_set,maximum_grid_depth);
    evolve.restrict(bounding);
    ARIADNE_LOG(5,"initial_size="<<evolve.size()<<"\n");

    GTS reach(grid);
    if(transient_time > 0.0 || transient_steps > 0) {
        ARIADNE_LOG(3,"Computing transient evolution...\n");
        make_lpair(reach,evolve)=this->_upper_reach_evolve(system,evolve,hybrid_transient_time,maximum_grid_depth);
        evolve.restrict(bounding);
        ARIADNE_LOG(5,"transient_reach_size="<<reach.size()<<"\n");
        ARIADNE_LOG(5,"evolve_size="<<evolve.size()<<"\n");
        ARIADNE_LOG(3,"  found "<<reach.size()<<" cells.\n");
    }
    GTS found(grid);
    ARIADNE_LOG(3,"Computing recurrent evolution...\n");

    while(!evolve.empty()) {
        make_lpair(found,evolve)=this->_upper_reach_evolve(system,evolve,hybrid_lock_to_grid_time,maximum_grid_depth);
        ARIADNE_LOG(5,"found.size()="<<found.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve.size()<<"\n");
        evolve.remove(reach);
        evolve.restrict(bounding);
        reach.adjoin(found);
        ARIADNE_LOG(3,"  found "<<found.size()<<" cells, of which "<<evolve.size()<<" are new.\n");
        // ARIADNE_LOG(5,"bounded_new_size="<<found.size()<<"\n");
    }
    reach.recombine();
    reach.restrict(bounding);
    return reach;
}





HybridReachabilityAnalyser::SetApproximationType
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
