/***************************************************************************
 *            hybrid_reachability_analyser.cc
 *
 *  Copyright  2006-11  Alberto Casagrande, Pieter Collins, Davide Bresolin
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

#include "config.h"

#include "exceptions.h"

#include "numeric.h"

#include "vector.h"
#include "matrix.h"

#include "box.h"
#include "list_set.h"
#include "grid_set.h"

#include "evolution_parameters.h"

#include "hybrid_automaton_interface.h"

#include "hybrid_time.h"
#include "hybrid_space.h"
#include "hybrid_orbit.h"
#include "hybrid_set.h"
#include "hybrid_evolver.h"
#include "hybrid_reachability_analyser.h"

#include "logging.h"
#include "graphics.h"
#include <include/linear_programming.h>


namespace Ariadne {

static const double DEFAULT_MAXIMUM_ENCLOSURE_RADIUS=0.25;
static const double DEFAULT_GRID_LENGTH=0.125;

HybridReachabilityAnalyser::
~HybridReachabilityAnalyser()
{
}


HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const HybridEvolverInterface& evolver)
    : _parameters(new EvolutionParametersType())
    , _evolver(evolver.clone())
    , _scaling(new HybridScaling())
{
}

HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const DiscreteEvolutionParameters& parameters,
                           const HybridEvolverInterface& evolver)
    : _parameters(new DiscreteEvolutionParameters(parameters))
    , _evolver(evolver.clone())
    , _scaling(new HybridScaling())
{
}

HybridReachabilityAnalyser*
HybridReachabilityAnalyser::
clone() const
{
    return new HybridReachabilityAnalyser(*this);
}



void
HybridReachabilityAnalyser::
set_scaling(const HybridScalingInterface& scaling)
{
    this->_scaling=shared_ptr<HybridScalingInterface>(scaling.clone());
}






// Helper functions for operators on lists of sets.
HybridGridTreeSet
HybridReachabilityAnalyser::_upper_reach(const HybridAutomatonInterface& sys,
                                         const HybridGridTreeSet& set,
                                         const HybridTime& time,
                                         const int accuracy) const
{
    HybridGrid grid=set.grid();
    HybridGridTreeSet result(grid);
    HybridGridTreeSet cells=set;
    cells.mince(accuracy);
    for(HybridGridTreeSet::const_iterator cell_iter=cells.begin(); cell_iter!=cells.end(); ++cell_iter) {
        //HybridEnclosure initial_enclosure(cell_iter->box(),this->_evolver->function_factory());
        HybridEnclosure initial_enclosure=this->_evolver->enclosure(cell_iter->box());
        ListSet<HybridEnclosure> reach = this->_evolver->reach(sys,initial_enclosure,time,UPPER_SEMANTICS);
        for(ListSet<HybridEnclosure>::const_iterator enclosure_iter=reach.begin(); enclosure_iter!=reach.end(); ++enclosure_iter) {
            enclosure_iter->adjoin_outer_approximation_to(result,accuracy);
            //result.adjoin_outer_approximation(*enclosure_iter,accuracy);
        }
    }
    return result;
}


HybridGridTreeSet
HybridReachabilityAnalyser::_upper_evolve(const HybridAutomatonInterface& sys,
                                          const HybridGridTreeSet& set,
                                          const HybridTime& time,
                                          const int accuracy) const
{
    HybridGrid grid=set.grid();
    HybridGridTreeSet result(grid); HybridGridTreeSet cells=set; cells.mince(accuracy);
    for(HybridGridTreeSet::const_iterator cell_iter=cells.begin(); cell_iter!=cells.end(); ++cell_iter) {
        ARIADNE_LOG(5,"Evolving cell = "<<*cell_iter<<"\n");
        HybridEnclosure initial_enclosure=this->_evolver->enclosure(cell_iter->box());
        ListSet<HybridEnclosure> final = this->_evolver->evolve(sys,initial_enclosure,time,UPPER_SEMANTICS);
        for(ListSet<HybridEnclosure>::const_iterator enclosure_iter=final.begin(); enclosure_iter!=final.end(); ++enclosure_iter) {
            enclosure_iter->adjoin_outer_approximation_to(result,accuracy);
            //result.adjoin_outer_approximation(*enclosure_iter,accuracy);
        }
    }
    ARIADNE_LOG(4,"_upper_evolve result size = "<<result.size()<<"\n");
    return result;
}


Pair<HybridGridTreeSet,HybridGridTreeSet>
HybridReachabilityAnalyser::_upper_reach_evolve(const HybridAutomatonInterface& sys,
                                                const HybridGridTreeSet& set,
                                                const HybridTime& time,
                                                const int accuracy) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::_upper_reach_evolve(...)\n");
    HybridGrid grid=set.grid();
    std::pair<HybridGridTreeSet,HybridGridTreeSet> result=make_pair(HybridGridTreeSet(grid),HybridGridTreeSet(grid));
    HybridGridTreeSet& reach_cells=result.first; HybridGridTreeSet& evolve_cells=result.second;
    HybridGridTreeSet cells=set; cells.mince(accuracy);

    /*
      for(HybridGridTreeSet::locations_const_iterator loc_iter=cells.locations_begin(); loc_iter!=cells.locations_end(); ++loc_iter) {
      for(GridTreeSet::const_iterator cell_iter=loc_iter->second.begin(); cell_iter!=loc_iter->second.end(); ++cell_iter) {
      Orbit<HybridGridCell> evolution=this->_discretiser->upper_evolution(sys,HybridGridCell(loc_iter->first,*cell_iter),time,accuracy);
      reach.adjoin(evolution.reach());
      evolve.adjoin(evolution.final());
      }
      }
    */

    for(HybridGridTreeSet::const_iterator cell_iter=cells.begin(); cell_iter!=cells.end(); ++cell_iter) {
        ARIADNE_LOG(5,"Evolving cell = "<<*cell_iter<<"\n");
        HybridEnclosure initial_enclosure=this->_evolver->enclosure(cell_iter->box());
        Orbit<HybridEnclosure> orbit = this->_evolver->orbit(sys,initial_enclosure,time,UPPER_SEMANTICS);
        ListSet<HybridEnclosure> const& reach_enclosures=orbit.reach();
        ListSet<HybridEnclosure> const& final_enclosures=orbit.final();
        for(ListSet<HybridEnclosure>::const_iterator enclosure_iter=reach_enclosures.begin(); enclosure_iter!=reach_enclosures.end(); ++enclosure_iter) {
            enclosure_iter->adjoin_outer_approximation_to(reach_cells,accuracy);
        }
        for(ListSet<HybridEnclosure>::const_iterator enclosure_iter=final_enclosures.begin(); enclosure_iter!=final_enclosures.end(); ++enclosure_iter) {
            enclosure_iter->adjoin_outer_approximation_to(evolve_cells,accuracy);
            //result.adjoin_outer_approximation(*enclosure_iter,accuracy);
        }
    }
    ARIADNE_LOG(3,"  final reach size = "<<reach_cells.size()<<"\n");
    ARIADNE_LOG(3,"  final evolve size = "<<evolve_cells.size()<<"\n");
    ARIADNE_LOG(2,"Done.\n");
    return result;
}

HybridGrid
HybridReachabilityAnalyser::
_hybrid_grid(const Sys& sys) const
{
    return HybridGrid(sys.state_space(),HybridScaling());
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
    HybridGrid grid=this->_hybrid_grid(system);
    HybridGridTreeSet initial_cells(grid); HybridGridTreeSet final_cells(grid);

    // Improve accuracy of initial set for lower computations
    initial_cells.adjoin_lower_approximation(initial_set,grid_height,grid_depth+4);
    ARIADNE_LOG(3,"initial_cells.size()="<<initial_cells.size()<<"\n");
    ARIADNE_LOG(3,"computing lower evolution.");
    for(HybridGridTreeSet::const_iterator cell_iter=initial_cells.begin(); cell_iter!=initial_cells.end(); ++cell_iter) {
        ARIADNE_LOG(3,".");
        HybridGridCell cell=*cell_iter;
        HybridEnclosure initial_enclosure=this->_evolver->enclosure(cell_iter->box());
        ListSet<HybridEnclosure> final_enclosures=this->_evolver->evolve(system,initial_enclosure,time,LOWER_SEMANTICS);
        for(ListSet<HybridEnclosure>::const_iterator enclosure_iter=final_enclosures.begin(); enclosure_iter!=final_enclosures.end(); ++enclosure_iter) {
            enclosure_iter->adjoin_outer_approximation_to(final_cells,grid_depth);
        }
    }
    ARIADNE_LOG(3,"\n");
    return final_cells;
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
    HybridGrid grid=this->_hybrid_grid(system);
    HybridGridTreeSet initial_cells(grid); HybridGridTreeSet reach_cells(grid);

    ARIADNE_LOG(3,"Adjoining initial set to the grid...\n");
    // Improve accuracy of initial set for lower computations
    initial_cells.adjoin_lower_approximation(initial_set,grid_height,grid_depth+4);
    ARIADNE_LOG(3,"initial_cells.size()="<<initial_cells.size()<<"\n");
    ARIADNE_LOG(3,"Computing lower reach set...");
    for(HybridGridTreeSet::const_iterator cell_iter=initial_cells.begin(); cell_iter!=initial_cells.end(); ++cell_iter) {
        ARIADNE_LOG(3,".");
        HybridGridCell cell=*cell_iter;
        HybridEnclosure initial_enclosure=this->_evolver->enclosure(cell_iter->box());
        ListSet<HybridEnclosure> reach_enclosures=this->_evolver->reach(system,initial_enclosure,time,LOWER_SEMANTICS);
        for(ListSet<HybridEnclosure>::const_iterator enclosure_iter=reach_enclosures.begin(); enclosure_iter!=reach_enclosures.end(); ++enclosure_iter) {
            enclosure_iter->adjoin_outer_approximation_to(reach_cells,grid_depth);
        }
    }
    ARIADNE_LOG(3,"\n");
    return reach_cells;
}



Pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
lower_reach_evolve(const SystemType& system,
                   const OvertSetInterfaceType& initial_set,
                   const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_reach_evolve(...)\n");
    int grid_depth = this->_parameters->maximum_grid_depth;
    int grid_height = this->_parameters->maximum_grid_height;

    HybridGrid grid=this->_hybrid_grid(system);

    HybridGridTreeSet initial_cells(grid);

    HybridGridTreeSet reach=(grid); HybridGridTreeSet evolve_cells(grid);

    // Improve accuracy of initial set for lower computations
    initial_cells.adjoin_lower_approximation(initial_set,grid_height,grid_depth+4);
    ARIADNE_LOG(3,"initial_cells.size()="<<initial_cells.size()<<"\n");
    ARIADNE_LOG(3,"computing lower evolution.");
    for(HybridGridTreeSet::const_iterator cell_iter=initial_cells.begin(); cell_iter!=initial_cells.end(); ++cell_iter) {
        ARIADNE_LOG(3,".");
        HybridEnclosure initial_enclosure=this->_evolver->enclosure(cell_iter->box());
        ListSet<HybridEnclosure> reach_enclosures;
        ListSet<HybridEnclosure> final_enclosures;
        make_lpair(reach_enclosures,final_enclosures) = this->_evolver->reach_evolve(system,initial_enclosure,time,LOWER_SEMANTICS);
        for(ListSet<HybridEnclosure>::const_iterator enclosure_iter=reach_enclosures.begin(); enclosure_iter!=reach_enclosures.end(); ++enclosure_iter) {
            enclosure_iter->adjoin_outer_approximation_to(reach,grid_depth);
        }
        for(ListSet<HybridEnclosure>::const_iterator enclosure_iter=final_enclosures.begin(); enclosure_iter!=final_enclosures.end(); ++enclosure_iter) {
            enclosure_iter->adjoin_outer_approximation_to(evolve_cells,grid_depth);
        }
    }
    ARIADNE_LOG(3,"\n");
    return make_pair(reach,evolve_cells);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_evolve(const SystemType& system,
             const CompactSetInterfaceType& initial_set,
             const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_evolve(...)\n");
    HybridGrid grid=this->_hybrid_grid(system);
    HybridGridTreeSet evolve_cells = *new HybridGridTreeSet(grid);
    int grid_depth = this->_parameters->maximum_grid_depth;
    evolve_cells.adjoin_outer_approximation(initial_set,grid_depth);
    ARIADNE_LOG(4,"initial_evolve.size()="<<evolve_cells.size()<<"\n");
    Float real_time=time.continuous_time();
    uint discrete_steps=time.discrete_time();
    Float lock_to_grid_time=this->_parameters->lock_to_grid_time;
    uint time_steps=integer_cast<uint>(real_time/lock_to_grid_time);
    Float remainder_time=real_time-time_steps*lock_to_grid_time;
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);
    ARIADNE_LOG(3,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(3,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");
    for(uint i=0; i!=time_steps; ++i) {
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...\n");
        evolve_cells=this->_upper_evolve(system,evolve_cells,hybrid_lock_to_grid_time,grid_depth);
    }
    ARIADNE_LOG(3,"remainder_time="<<remainder_time<<"\n");
    if(!evolve_cells.empty() && remainder_time > 0) {
        ARIADNE_LOG(3,"computing evolution for remainder time...\n");
        evolve_cells=this->_upper_evolve(system,evolve_cells,hybrid_remainder_time,grid_depth);
    }
    evolve_cells.recombine();
    ARIADNE_LOG(4,"final_evolve.size()="<<evolve_cells.size()<<"\n");
    return evolve_cells;
}



HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_reach(const SystemType& system,
            const CompactSetInterfaceType& initial_set,
            const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach(system,set,time)\n");
    ARIADNE_LOG(4,"initial_set="<<initial_set<<"\n");
    HybridGrid grid=this->_hybrid_grid(system);
    HybridGridTreeSet evolve_cells(grid);
    int grid_depth = this->_parameters->maximum_grid_depth;
    ARIADNE_LOG(4,"grid_depth="<<grid_depth<<"\n");
    evolve_cells.adjoin_outer_approximation(initial_set,grid_depth);
    ARIADNE_LOG(4,"initial size = "<<evolve_cells.size()<<"\n");
    HybridGridTreeSet reach_cells(evolve_cells);
    ARIADNE_LOG(4,"reach size ="<<reach_cells.size()<<"\n");
    Float real_time=time.continuous_time();
    uint discrete_steps=time.discrete_time();
    Float lock_to_grid_time = this->_parameters->lock_to_grid_time;
    uint time_steps=integer_cast<uint>(real_time/lock_to_grid_time);
    Float remainder_time=real_time-time_steps*lock_to_grid_time;
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);
    ARIADNE_LOG(3,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(3,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");
    ARIADNE_LOG(3,"discrete_steps="<<discrete_steps<<"\n");
    HybridGridTreeSet found_cells(grid), accumulated_evolve_cells(grid);
    for(uint i=0; i!=time_steps; ++i) {
        accumulated_evolve_cells.adjoin(evolve_cells);
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...\n");
        make_lpair(found_cells,evolve_cells)=this->_upper_reach_evolve(system,evolve_cells,hybrid_lock_to_grid_time,grid_depth);
        ARIADNE_LOG(5,"found.size()="<<found_cells.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve_cells.size()<<"\n");
        evolve_cells.remove(accumulated_evolve_cells);
        reach_cells.adjoin(found_cells);
        ARIADNE_LOG(3,"  found "<<found_cells.size()<<" cells, and "<<evolve_cells.size()<<" are new intermediate.\n");
        if(evolve_cells.empty()) break;
        ARIADNE_LOG(6,"evolve_cells="<<evolve_cells<<"\n");
    }
    ARIADNE_LOG(3,"remainder_time="<<remainder_time<<"\n");
    if(!evolve_cells.empty() && remainder_time > 0) {
        ARIADNE_LOG(3,"computing evolution for remainder time...\n");
        make_lpair(found_cells,evolve_cells)=this->_upper_reach_evolve(system,evolve_cells,hybrid_remainder_time,grid_depth);
        reach_cells.adjoin(found_cells);
    }
    // This last step is necessary to add the final set to the result.
    reach_cells.adjoin(evolve_cells);
    reach_cells.recombine();
    ARIADNE_LOG(4,"final_reach size = "<<reach_cells.size()<<"\n");
    return reach_cells;
}



Pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
upper_reach_evolve(const SystemType& system,
                   const CompactSetInterfaceType& initial_set,
                   const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach_evolve(system,set,time)\n");
    ARIADNE_LOG(4,"initial_set="<<initial_set<<"\n");
    HybridGrid grid=this->_hybrid_grid(system);
    HybridGridTreeSet evolve_cells(grid);
    int grid_depth = this->_parameters->maximum_grid_depth;
    ARIADNE_LOG(4,"grid_depth="<<grid_depth<<"\n");
    evolve_cells.adjoin_outer_approximation(initial_set,grid_depth);
    ARIADNE_LOG(4,"initial_evolve"<<evolve_cells<<"\n");
    HybridGridTreeSet reach_cells(evolve_cells);
    ARIADNE_LOG(4,"reach="<<reach_cells<<"\n");
    Float real_time=time.continuous_time();
    uint discrete_steps=time.discrete_time();
    Float lock_to_grid_time = this->_parameters->lock_to_grid_time;
    uint time_steps=integer_cast<uint>(real_time/lock_to_grid_time);
    Float remainder_time=real_time-time_steps*lock_to_grid_time;
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);
    ARIADNE_LOG(3,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(3,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");
    HybridGridTreeSet found_cells(grid);
    for(uint i=0; i!=time_steps; ++i) {
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...\n");
        make_lpair(found_cells,evolve_cells) = this->_upper_reach_evolve(system,evolve_cells,hybrid_lock_to_grid_time,grid_depth);
        ARIADNE_LOG(5,"found.size()="<<found_cells.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve_cells.size()<<"\n");
        reach_cells.adjoin(found_cells);
        ARIADNE_LOG(3,"  found "<<found_cells.size()<<" cells.\n");
    }
    ARIADNE_LOG(3,"remainder_time="<<remainder_time<<"\n");
    if(!evolve_cells.empty() && remainder_time > 0) {
        ARIADNE_LOG(3,"computing evolution for remainder time...\n");
        make_lpair(found_cells,evolve_cells) = this->_upper_reach_evolve(system,evolve_cells,hybrid_remainder_time,grid_depth);
        reach_cells.adjoin(found_cells);
    }
    reach_cells.recombine();
    ARIADNE_LOG(4,"reach="<<reach_cells<<"\n");
    evolve_cells.recombine();
    ARIADNE_LOG(4,"evolve="<<evolve_cells<<"\n");
    return std::make_pair(reach_cells,evolve_cells);
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

    HybridGrid grid=this->_hybrid_grid(system);
    HybridGridTreeSet evolve_cells(grid), accumulated_evolve_cells(grid);
    evolve_cells.adjoin_outer_approximation(initial_set,maximum_grid_depth);
    ARIADNE_LOG(5,"initial_size="<<evolve_cells.size()<<"\n");
    accumulated_evolve_cells.adjoin(evolve_cells);

    HybridGridTreeSet reach_cells = HybridGridTreeSet(grid);
    if(transient_time > 0.0 || transient_steps > 0) {
        ARIADNE_LOG(3,"Computing transient evolution...\n");
        make_lpair(reach_cells,evolve_cells)=this->_upper_reach_evolve(system,evolve_cells,hybrid_transient_time,maximum_grid_depth);
        ARIADNE_LOG(5,"transient_reach_size="<<reach_cells.size()<<"\n");
        ARIADNE_LOG(5,"evolve_size="<<evolve_cells.size()<<"\n");
        ARIADNE_LOG(3,"  found "<<reach_cells.size()<<" cells.\n");
    }
    HybridGridTreeSet found_cells(grid);
    ARIADNE_LOG(3,"Computing recurrent evolution...\n");
    while(!evolve_cells.empty()) {
        accumulated_evolve_cells.adjoin(evolve_cells);
        make_lpair(found_cells,evolve_cells)=this->_upper_reach_evolve(system,evolve_cells,hybrid_lock_to_grid_time,maximum_grid_depth);
        ARIADNE_LOG(5,"found.size()="<<found_cells.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve_cells.size()<<"\n");
        evolve_cells.remove(accumulated_evolve_cells);
        evolve_cells.restrict_to_height(maximum_grid_height);
        found_cells.restrict_to_height(maximum_grid_height);
        reach_cells.adjoin(found_cells);
        ARIADNE_LOG(3,"  found "<<found_cells.size()<<" cells, of which "<<evolve_cells.size()<<" are new.\n");
        // ARIADNE_LOG(5,"bounded_new_size="<<found_cells.size()<<"\n");
    }
    reach_cells.recombine();
    return reach_cells;
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

    HybridGrid grid=this->_hybrid_grid(system);
    //HybridGridTreeSet bounding; bounding.adjoin_inner_approximation(bounding_domain,maximum_grid_depth);
    HybridGridTreeSet bounding(grid);
    bounding.adjoin_outer_approximation(bounding_domain,maximum_grid_depth); bounding.recombine();
    ARIADNE_LOG(5,"bounding_size="<<bounding.size()<<"\n");

    HybridGridTreeSet evolve_cells(grid), accumulated_evolve_cells(grid);
    evolve_cells.adjoin_outer_approximation(initial_set,maximum_grid_depth);
    evolve_cells.restrict(bounding);
    accumulated_evolve_cells.adjoin(evolve_cells);
    ARIADNE_LOG(5,"initial_size="<<evolve_cells.size()<<"\n");

    HybridGridTreeSet reach_cells(grid);
    if(transient_time > 0.0 || transient_steps > 0) {
        ARIADNE_LOG(3,"Computing transient evolution...\n");
        make_lpair(reach_cells,evolve_cells)=this->_upper_reach_evolve(system,evolve_cells,hybrid_transient_time,maximum_grid_depth);
        evolve_cells.restrict(bounding);
        ARIADNE_LOG(5,"transient_reach_size="<<reach_cells.size()<<"\n");
        ARIADNE_LOG(5,"evolve_size="<<evolve_cells.size()<<"\n");
        ARIADNE_LOG(3,"  found "<<reach_cells.size()<<" cells.\n");
    }
    HybridGridTreeSet found_cells(grid);
    ARIADNE_LOG(3,"Computing recurrent evolution...\n");

    while(!evolve_cells.empty()) {
        accumulated_evolve_cells.adjoin(evolve_cells);
        make_lpair(found_cells,evolve_cells)=this->_upper_reach_evolve(system,evolve_cells,hybrid_lock_to_grid_time,maximum_grid_depth);
        ARIADNE_LOG(5,"found.size()="<<found_cells.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve_cells.size()<<"\n");
        evolve_cells.remove(accumulated_evolve_cells);
        evolve_cells.restrict(bounding);
        reach_cells.adjoin(found_cells);
        ARIADNE_LOG(3,"  found "<<found_cells.size()<<" cells, of which "<<evolve_cells.size()<<" are new.\n");
        // ARIADNE_LOG(5,"bounded_new_size="<<found_cells.size()<<"\n");
    }
    reach_cells.recombine();
    reach_cells.restrict(bounding);
    return reach_cells;
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
