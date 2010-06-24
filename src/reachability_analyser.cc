/*****************************************************************************************
 *            reachability_analyser.cc
 *
 *  Copyright  2006-10  Alberto Casagrande, Pieter Collins, Davide Bresolin, Luca Geretti
 *
 *****************************************************************************************/

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

#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
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
#include "evolution_statistics.h"
#include "evolver_interface.h"
#include "taylor_calculus.h"

#include "discretiser.h"
#include "reachability_analyser.h"
#include "logging.h"

#include "graphics.h"

#include "workers.h"


namespace Ariadne {

HybridReachabilityAnalyser::
~HybridReachabilityAnalyser()
{
}


HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const HybridDiscretiser<HybridEvolver::ContinuousEnclosureType>& discretiser)
    : _parameters(new EvolutionParametersType())
    , _discretiser(discretiser.clone())
{
	plot_verify_results = false;
}

list<HybridBasicSet<TaylorSet> >
HybridReachabilityAnalyser::
_split_initial_set(const HybridImageSet initial_set,
				   const HybridGrid grid) const
{
	// The initial enclosures
	list<EnclosureType> initial_enclosures;

	// For each location, split the original enclosure into a list of enclosures
	// having widths lesser than the grid lengths at maximum depth
	for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin();
		 loc_iter!=initial_set.locations_end(); ++loc_iter)
		{
		// Get the grid lengths at the maximum depth (i.e. the minimum cell widths)
		Vector<Float> mincell_widths = grid[loc_iter->first].lengths()/(1 << _parameters->maximum_grid_depth);

		// Keep a list of the continuous enclosures to split
		list<ContinuousEnclosureType> continuous_enclosures;
		// Insert the original continuous enclosure
		continuous_enclosures.push_front(ContinuousEnclosureType(loc_iter->second));

		/* While there are continuous enclosures:
		 * a) pick one, pop it and check it against the mincell_lengths
		 * 	 i) if larger, split it and put the couple into the continuous_enclosures
		 *   ii) if not, put the enclosure (paired with the location) into the initial_enclosures
		 */
		while (!continuous_enclosures.empty()) {
			// Pick and pop from the back (interpreted as the least recent)
			ContinuousEnclosureType enclosure = continuous_enclosures.back();
			continuous_enclosures.pop_back();
			// Keep a flag to identify if any splitting occurred (defaults to false)
			bool hasBeenSplit = false;
			// For each dimension
			for (uint i=0; i < enclosure.dimension(); i++) {
				// If the width of the range for the model on the i-th dimension is greater than the minimum cell length
				if (enclosure[i].range().width() > mincell_widths[i]) {
					// Split on the dimension and put the resulting enclosures in the continuous enclosures
					pair<ContinuousEnclosureType,ContinuousEnclosureType> split_result = enclosure.split(i);
					continuous_enclosures.push_front(split_result.first);
					continuous_enclosures.push_front(split_result.second);
					// Set the flag as true and skip the remaining checks
					hasBeenSplit = true;
					break;
				}
			}
			// If it has not been split, put the hybrid enclosure into the initial enclosures
			if (!hasBeenSplit)
				initial_enclosures.push_back(EnclosureType(loc_iter->first,enclosure));
		}
	}

	ARIADNE_LOG(6,"Split initial enclosures size = " << initial_enclosures.size() << "\n");

	// Returns
	return initial_enclosures;
}

// Helper functions for operators on lists of sets.
HybridGridTreeSet
HybridReachabilityAnalyser::_upper_reach(const HybridAutomaton& sys,
                                         const HybridGridTreeSet& set,
                                         const HybridTime& time,
                                         const int accuracy) const
{
    Gr grid=sys.grid();
    HybridGridTreeSet result(grid);
    HybridGridTreeSet cells=set;
    cells.mince(accuracy);
    for(HybridGridTreeSet::const_iterator iter=cells.begin(); iter!=cells.end(); ++iter) {
        EnclosureType enclosure=_discretiser->enclosure(*iter);
        result.adjoin(_discretiser->reach(sys,enclosure,time,accuracy,UPPER_SEMANTICS));
    }
    return result;
}


HybridGridTreeSet
HybridReachabilityAnalyser::_lower_reach(const HybridAutomaton& system,
							 			 std::list<EnclosureType>& initial_enclosures,
										 const HybridTime& time,
										 const HybridTime& lock_time) const
{
	typedef std::list<EnclosureType> EL;
	typedef std::map<DiscreteState,uint> HUM;

	// Get the actual concurrency
	const uint concurrency = boost::thread::hardware_concurrency() - free_cores;
	// Check that the concurrency is greater than zero
	ARIADNE_ASSERT_MSG(concurrency>0,"Error: concurrency must be at least 1.");

	// The final enclosures
	EL final_enclosures;
	HybridGrid grid = system.grid();

	// The result
	GTS reach(grid);

	// For each grid lock
	for (uint i=0;i<((uint)time.discrete_time()/lock_time.discrete_time());i++)
	{
		// The sizes of the adjoined (the cells) or superposed (the enclosures) evolve
		std::map<DiscreteState,uint> adjoined_evolve_sizes;
		std::map<DiscreteState,uint> superposed_evolve_sizes;
		// The evolve
		GTS evolve;
		// The current reach
		GTS current_reach;

		// Create the worker
		LowerReachWorker worker(_discretiser,initial_enclosures,system,lock_time,
							    _parameters->bounding_domain,_parameters->maximum_grid_depth,concurrency);

		// Compute and get the result
		make_ltuple<HUM,HUM,EL,GTS>(adjoined_evolve_sizes,superposed_evolve_sizes,final_enclosures,current_reach) = worker.get_result();

		// Adjoin the current reach
		reach.adjoin(current_reach);

		// Pruning of the dump of the final enclosures into the initial enclosures
		while (!final_enclosures.empty())
		{
			// Pop the current enclosure
			EnclosureType encl = final_enclosures.front();
			final_enclosures.pop_front();

			// If pruning is to be performed
			if (_parameters->enable_lower_pruning)
			{
				// Get the ratio between the adjoined evolve size and the superposed evolve size
				Float ratio = (Float)adjoined_evolve_sizes[encl.location()]/(Float)superposed_evolve_sizes[encl.location()];

				// At least 2 enclosures are inserted, then the enclosures are pruned as long as the number of enclosures is at least twice the number of evolve cells
				if (initial_enclosures.size() <= 2 || rand() < 2*ratio*RAND_MAX)
					initial_enclosures.push_back(encl);
			}
			else
				initial_enclosures.push_back(encl);
		}
	}

    return reach;
}


HybridGridTreeSet
HybridReachabilityAnalyser::_upper_evolve(const HybridAutomaton& sys,
                                          const HybridGridTreeSet& set,
                                          const HybridTime& time,
                                          const int accuracy) const
{
    Gr grid=sys.grid();
    GTS result(grid); GTS cells=set; cells.mince(accuracy);
    for(HybridGridTreeSet::const_iterator iter=cells.begin(); iter!=cells.end(); ++iter) {
        ARIADNE_LOG(6,"Evolving cell = "<<*iter<<"\n");
        EnclosureType enclosure=_discretiser->enclosure(*iter);
        result.adjoin(_discretiser->evolve(sys,enclosure,time,accuracy,UPPER_SEMANTICS));
    }
    ARIADNE_LOG(5,"_upper_evolve result size = "<<result.size()<<"\n");
    return result;
}


std::pair<HybridGridTreeSet,HybridGridTreeSet>
HybridReachabilityAnalyser::_upper_reach_evolve(const HybridAutomaton& sys,
                                                const HybridGridTreeSet& set,
                                                const HybridTime& time,
                                                const int accuracy) const
{
    ARIADNE_LOG(6,"HybridReachabilityAnalyser::_upper_reach_evolve(...)\n");

	// Copy the set into the initial set and mince it to the desired precision
	GTS initial_set = set;
	initial_set.mince(accuracy);

	// Get the actual concurrency
	const uint concurrency = boost::thread::hardware_concurrency() - free_cores;

	// Check that the concurrency is greater than zero
	ARIADNE_ASSERT_MSG(concurrency>0,"Error: concurrency must be at least 1.");

	// Create the worker
	UpperReachEvolveWorker worker(_discretiser,sys,set,time,accuracy,concurrency);

	// Compute and get the result
	std::pair<GTS,GTS> result = worker.get_result();

    ARIADNE_LOG(6,"Reach size = "<<result.first.size()<<"\n");
    ARIADNE_LOG(6,"Final size = "<<result.second.size()<<"\n");
    return result;
}

std::pair<HybridGridTreeSet,HybridGridTreeSet>
HybridReachabilityAnalyser::_upper_reach_evolve_continuous(const HybridAutomaton& sys,
                                                		   const list<EnclosureType>& initial_enclosures,
                                                		   const HybridTime& time,
                                                		   const int accuracy) const
{
    ARIADNE_LOG(6,"HybridReachabilityAnalyser::_upper_reach_evolve_continuous(...)\n");

	// Get the actual concurrency
	const uint concurrency = boost::thread::hardware_concurrency() - free_cores;

	// Check that the concurrency is greater than zero
	ARIADNE_ASSERT_MSG(concurrency>0,"Error: concurrency must be at least 1.");

	// Create the worker
	UpperReachEvolveContinuousWorker worker(_discretiser,sys,initial_enclosures,_parameters->bounding_domain,time,accuracy,concurrency);

	ARIADNE_LOG(6,"Evolving and discretising...\n");
	// Compute and get the result
	std::pair<GTS,GTS> result = worker.get_result();

    ARIADNE_LOG(6,"Reach size = "<<result.first.size()<<"\n");
    ARIADNE_LOG(6,"Final size = "<<result.second.size()<<"\n");
    return result;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_evolve(const SystemType& system,
             const HybridImageSet& initial_set,
             const TimeType& time) const
{
    ARIADNE_LOG(4,"HybridReachabilityAnalyser::lower_evolve(...)\n");

	_discretiser->reset_lower_statistics(); // Resets the continuous statistics

    GTS evolve;

    // Split the initial set into enclosures that are smaller than the minimum cell, given the grid
    list<EnclosureType> initial_enclosures = _split_initial_set(initial_set,system.grid());

	ARIADNE_LOG(5,"Computing evolution...\n");
    // For each initial enclosure, perform evolution and adjoin to the reach and evolve regions
    for (list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); encl_it++) {
        GTS cell_final=_discretiser->evolve(system,*encl_it,time,_parameters->maximum_grid_depth,LOWER_SEMANTICS);
        evolve.adjoin(cell_final);
    }

    return evolve;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_reach(const SystemType& system,
            const HybridImageSet& initial_set,
            const TimeType& time) const
{
    ARIADNE_LOG(4,"HybridReachabilityAnalyser::lower_reach(...)\n");

	_statistics->lower().reset(); // Reset the discrete statistics
	_discretiser->reset_lower_statistics(); // Reset the continuous statistics

	// Get the actual concurrency
	const uint concurrency = boost::thread::hardware_concurrency() - free_cores;
	// Check that the concurrency is greater than zero
	ARIADNE_ASSERT_MSG(concurrency>0,"Error: concurrency must be at least 1.");

	// Get the lock time, in the case pruning is to be performed
	TimeType lock_time(_parameters->lock_to_grid_time,_parameters->lock_to_grid_steps);

    // Split the initial set into enclosures that are smaller than the minimum cell, given the grid
    list<EnclosureType> initial_enclosures = _split_initial_set(initial_set,system.grid());

    ARIADNE_LOG(6,"Evolving and discretising...\n");
    // Calculates the reached region
    GTS reach = _lower_reach(system,initial_enclosures,time,lock_time);

	// Copies the reached region into the statistics
	_statistics->lower().reach = reach;

	return reach;
}


std::pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
lower_reach_evolve(const SystemType& system,
                   const HybridImageSet& initial_set,
                   const TimeType& time) const
{
    ARIADNE_LOG(4,"HybridReachabilityAnalyser::lower_reach_evolve(...)\n");

	_statistics->lower().reset(); // Resets the discrete statistics
	_discretiser->reset_lower_statistics(); // Resets the continuous statistics

    int grid_depth = _parameters->maximum_grid_depth;

    Gr grid=system.grid();
    GTS reach; GTS evolve;

    // Split the initial set into enclosures that are smaller than the minimum cell, given the grid
    list<EnclosureType> initial_enclosures = _split_initial_set(initial_set,grid);

	ARIADNE_LOG(5,"Computing evolution...\n");
    // For each initial enclosure, perform evolution and adjoin to the reach and evolve regions
    for (list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); encl_it++) {
        GTS cell_reach,cell_final;
        make_lpair(cell_reach,cell_final)=_discretiser->evolution(system,*encl_it,time,grid_depth,LOWER_SEMANTICS);
        reach.adjoin(cell_reach);
        evolve.adjoin(cell_final);
    }

	// Copies the reached region
	_statistics->lower().reach = reach;

    return make_pair(reach,evolve);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_evolve(const SystemType& system,
             const HybridImageSet& initial_set,
             const TimeType& time) const
{
    ARIADNE_LOG(4,"HybridReachabilityAnalyser::upper_evolve(...)\n");
 
	_statistics->upper().reset(); // Resets the discrete statistics
	_discretiser->reset_upper_statistics(); // Resets the continuous statistics

    Gr grid=system.grid();
    GTS evolve(grid);
    int grid_depth = _parameters->maximum_grid_depth;
    Float real_time=time.continuous_time();
    uint discrete_steps=time.discrete_time();
    Float lock_to_grid_time=_parameters->lock_to_grid_time;
    uint time_steps=uint(real_time/lock_to_grid_time);
    Float remaining_time=real_time-time_steps*lock_to_grid_time;
    if(time_steps == 0) {
        time_steps=1;
        remaining_time=0.0;
        lock_to_grid_time=real_time;
    }
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remaining_time(remaining_time,discrete_steps);
    ARIADNE_LOG(5,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(5,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");

    // Split the initial set into enclosures that are smaller than the minimum cell, given the grid
    list<EnclosureType> initial_enclosures = _split_initial_set(initial_set,grid);

	ARIADNE_LOG(5,"Computing initial evolution...\n");
    // For each initial enclosure, perform evolution and adjoin to the reach and evolve regions
    for (list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); encl_it++) {
        GTS cell_final=_discretiser->evolve(system,*encl_it,hybrid_lock_to_grid_time,grid_depth,UPPER_SEMANTICS);
        evolve.adjoin(cell_final);
    }

    // time steps evolution loop
    for(uint i=1; i<time_steps; ++i) {
        ARIADNE_LOG(5,"Computing "<<i+1<<"-th reachability step...\n");
        evolve=_upper_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
    }

    ARIADNE_LOG(5,"remaining_time="<<remaining_time<<"\n");
    if(!evolve.empty() && remaining_time > 0) {
        ARIADNE_LOG(5,"Computing evolution for remaining time...\n");
        evolve=_upper_evolve(system,evolve,hybrid_remaining_time,grid_depth);
    }
    evolve.recombine();
    ARIADNE_LOG(6,"final_evolve.size()="<<evolve.size()<<"\n");
    return evolve;
}



HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_reach(const SystemType& system,
            const HybridImageSet& initial_set,
            const TimeType& time) const
{
    ARIADNE_LOG(4,"HybridReachabilityAnalyser::upper_reach(system,set,time)\n");

	GTS reach, evolve;
	make_lpair(reach,evolve) = upper_reach_evolve(system, initial_set, time); // Runs the upper_reach_evolve routine on its behalf

    return reach; // Returns the reached region only
}



std::pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
upper_reach_evolve(const SystemType& system,
                   const HybridImageSet& initial_set,
                   const TimeType& time) const
{
    ARIADNE_LOG(4,"HybridReachabilityAnalyser::upper_reach_evolve(system,set,time)\n");
    ARIADNE_LOG(5,"initial_set="<<initial_set<<"\n");

	_statistics->upper().reset(); // Reset the discrete statistics
	_discretiser->reset_upper_statistics(); // Reset the continuous statistics

    Gr grid=system.grid();
    GTS found(grid),evolve(grid),reach(grid);
    int grid_depth = _parameters->maximum_grid_depth;
    Float real_time=time.continuous_time();
    uint discrete_steps=time.discrete_time();
    Float lock_to_grid_time=_parameters->lock_to_grid_time;
    uint time_steps=uint(real_time/lock_to_grid_time);
    Float remaining_time=real_time-time_steps*lock_to_grid_time;
    if(time_steps == 0) {
        time_steps=1;
        remaining_time=0.0;
        lock_to_grid_time=real_time;
    }
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remaining_time(remaining_time,discrete_steps);
    ARIADNE_LOG(5,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(5,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");

    // Split the initial set into enclosures that are smaller than the minimum cell, given the grid
    list<EnclosureType> initial_enclosures = _split_initial_set(initial_set,grid);

	ARIADNE_LOG(5,"Computing initial evolution...\n");
    // For each initial enclosure, perform evolution and adjoin to the reach and evolve regions
    for (list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); encl_it++) {
        GTS cell_reach,cell_final;
        make_lpair(cell_reach,cell_final)=_discretiser->evolution(system,*encl_it,hybrid_lock_to_grid_time,grid_depth,UPPER_SEMANTICS);
        reach.adjoin(cell_reach);
        evolve.adjoin(cell_final);
    }

    // time steps evolution loop        
    for(uint i=1; i<time_steps; ++i) {
        ARIADNE_LOG(5,"computing "<<i+1<<"-th reachability step...\n");
        make_lpair(found,evolve) = _upper_reach_evolve(system,evolve,hybrid_lock_to_grid_time,grid_depth);
        ARIADNE_LOG(6,"found.size()="<<found.size()<<"\n");
        ARIADNE_LOG(6,"evolve.size()="<<evolve.size()<<"\n");
        reach.adjoin(found);
        ARIADNE_LOG(5,"found "<<found.size()<<" cells.\n");
    }

    ARIADNE_LOG(5,"remaining_time="<<remaining_time<<"\n");
    if(!evolve.empty() && remaining_time > 0) {
        ARIADNE_LOG(5,"computing evolution for the remaining time...\n");
        make_lpair(found,evolve) = _upper_reach_evolve(system,evolve,hybrid_remaining_time,grid_depth);
        reach.adjoin(found);
    }

    reach.recombine();
	_statistics->upper().reach = reach;
    ARIADNE_LOG(5,"reach="<<reach<<"\n");
    evolve.recombine();
    ARIADNE_LOG(5,"evolve="<<evolve<<"\n");
    return std::make_pair(reach,evolve);
}




HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
chain_reach(const SystemType& system,
            const HybridImageSet& initial_set) const
{
    ARIADNE_LOG(4,"HybridReachabilityAnalyser::chain_reach(system,initial_set)\n");

	// Reset statistics
	_statistics->upper().reset();
	_discretiser->reset_upper_statistics();

	// Assign local variables
    HybridBoxes bounding_domain = _parameters->bounding_domain;
    Float transient_time = _parameters->transient_time;
    int transient_steps = _parameters->transient_steps;
    Float lock_to_grid_time = _parameters->lock_to_grid_time;
    int lock_to_grid_steps = _parameters->lock_to_grid_steps;
    int maximum_grid_depth = _parameters->maximum_grid_depth;

    ARIADNE_LOG(6,"transient_time=("<<transient_time<<","<<transient_steps<<")\n");
    ARIADNE_LOG(6,"lock_to_grid_time=("<<lock_to_grid_time<<","<<lock_to_grid_steps<<")\n");

	// Checks consistency of the bounding domain in respect to the state space
	HybridSpace hspace = system.state_space();
	// If the DiscreteState was not found or otherwise if the continuous space sizes mismatch, throws an error
	for (HybridSpace::locations_const_iterator hs_it = hspace.locations_begin(); hs_it != hspace.locations_end(); ++hs_it) {
		if (bounding_domain.find(hs_it->first) == bounding_domain.end()) {
			ARIADNE_FAIL_MSG("Error: the system state space and the bounding domain space do not match on the discrete space."); }		
		else if (hs_it->second != bounding_domain[hs_it->first].size()) {
			ARIADNE_FAIL_MSG("Error: the system state space and the bounding domain space do not match on the continuous space."); }}

    Gr grid=system.grid();
    GTS bounding(grid), evolve(grid), reach(grid), found(grid), intermediate(grid);

    bounding.adjoin_outer_approximation(bounding_domain,0u); 
    ARIADNE_LOG(6,"bounding_size(pre recombine)="<<bounding.size()<<"\n");
	bounding.recombine();
	HybridBoxes bounding_box = bounding.bounding_box(); // Used for the restriction check
    ARIADNE_LOG(6,"bounding_size(post recombine)="<<bounding.size()<<"\n");

    if(transient_time <= 0.0 || transient_steps <= 0) {
        transient_time = lock_to_grid_time;
        transient_steps = lock_to_grid_steps; }

    HybridTime hybrid_transient_time(transient_time, transient_steps);
    ARIADNE_LOG(5,"Computing first evolution step...\n");

    // Split the initial set into enclosures that are smaller than the minimum cell, given the grid
    list<EnclosureType> initial_enclosures = _split_initial_set(initial_set,grid);

	ARIADNE_LOG(5,"Computing transient evolution...\n");
    // For each initial enclosure, perform evolution and adjoin to the reach and evolve regions
    for (list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); encl_it++) {
        GTS cell_reach,cell_final;
        make_lpair(cell_reach,cell_final)=_discretiser->evolution(system,*encl_it,hybrid_transient_time,maximum_grid_depth,UPPER_SEMANTICS);
        reach.adjoin(cell_reach);
        evolve.adjoin(cell_final);
    }
	 
    // If the evolve region is not a subset of the bounding region, the region will be restricted (NOTE: for efficiency, only performed if the region is currently considered unrestricted)
	if (!_statistics->upper().has_restriction_occurred)
		if (!evolve.subset(bounding_box)) _statistics->upper().has_restriction_occurred = true;

    evolve.restrict(bounding);

    ARIADNE_LOG(5,"found "<<reach.size()<<" cells, of which "<<evolve.size()<<" are new.\n");
    
    ARIADNE_LOG(5,"Computing recurrent evolution...\n");
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,lock_to_grid_steps);

	// While the final set has new cells in respect to the previous reach set, process them and increase the number of locks (starting from 1 due to the initial transient phase)
    for (uint i = 0; !evolve.empty(); i++) {

    	ARIADNE_LOG(5,"Iteration "<<i<<"\n");
		// Add the evolve region to the intermediate region
		intermediate.adjoin(evolve);  

        make_lpair(found,evolve)=_upper_reach_evolve(system,evolve,hybrid_lock_to_grid_time,maximum_grid_depth);
        ARIADNE_LOG(6,"found.size()="<<found.size()<<"\n");
        ARIADNE_LOG(6,"evolve.size()="<<evolve.size()<<"\n");

		// If the evolve region is not a subset of the bounding region, the region will be restricted (NOTE: for efficiency, only performed if the region is currently considered unrestricted)
		if (!_statistics->upper().has_restriction_occurred)
			if (!evolve.subset(bounding_box)) _statistics->upper().has_restriction_occurred = true;

        evolve.remove(intermediate); // Remove only the cells that are intermediate
        evolve.restrict(bounding);
        reach.adjoin(found);

        ARIADNE_LOG(6,"found "<<found.size()<<" cells, of which "<<evolve.size()<<" are new.\n");
    }
    reach.recombine();

    // If the evolve region is not a subset of the bounding region, the region will be restricted (NOTE: for efficiency, only performed if the region is currently considered unrestricted)
	if (!_statistics->upper().has_restriction_occurred)
		if (!evolve.subset(bounding_box)) _statistics->upper().has_restriction_occurred = true;

    reach.restrict(bounding);
	_statistics->upper().reach = reach;

    return reach;
}

HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
chain_reach(const SystemType& system,
            const HybridImageSet& initial_set,
            const HybridBoxes& bounding_set) const
{
	// Assigns the input bounding_set to the bounding domain
	_parameters->bounding_domain = bounding_set;

	// Returns the result of chain_reach with implicit bounding set
	return chain_reach(system,initial_set);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_chain_reach(const SystemType& system,
            const HybridImageSet& initial_set) const
{
    ARIADNE_LOG(4,"HybridReachabilityAnalyser::upper_chain_reach(system,initial_set)\n");

	// Reset statistics
	_statistics->upper().reset();
	_discretiser->reset_upper_statistics();

	// Helper class for operations on Taylor sets
	TaylorCalculus tc;

	// Assign local variables
    Float lock_to_grid_time = _parameters->lock_to_grid_time;
    int lock_to_grid_steps = _parameters->lock_to_grid_steps;
    int maximum_grid_depth = _parameters->maximum_grid_depth;

    Gr grid=system.grid();
    GTS new_final(grid), new_reach(grid), reach(grid), intermediate(grid);

    // The information about previous and current locations in use (only for dumping the chain reach)
	std::map<DiscreteState,bool> previousLocationsInUse;
	std::map<DiscreteState,bool> currentLocationsInUse;

    // If dumping must be performed
    if (chain_reach_dumping) {
    	// Initialize (with true values) the information about locations in use on the previous iteration and
    	// initialize (with false values) the information about locations in use on the current iteration
		for (HybridGrid::locations_const_iterator grid_it = grid.locations_begin(); grid_it != grid.locations_end(); grid_it++) {
			previousLocationsInUse[grid_it->first] = true;
			currentLocationsInUse[grid_it->first] = false;
		}
		// For each location, set to true the corresponding currentLocationsInUse value
		for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin(); loc_iter!=initial_set.locations_end(); ++loc_iter)
			currentLocationsInUse[loc_iter->first] = true;
    }

    // Split the initial set into enclosures that are smaller than the minimum cell, given the grid
    list<EnclosureType> initial_enclosures = _split_initial_set(initial_set,grid);

    ARIADNE_LOG(5,"Computing recurrent evolution...\n");
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,lock_to_grid_steps);

	// While the final set has new cells in respect to the previous reach set, process them and increase the number of locks
    for (uint i=0; !initial_enclosures.empty(); i++)
	{
    	ARIADNE_LOG(5,"Iteration " << i << "\n");

	    // If dumping must be performed
	    if (chain_reach_dumping) {
	    	// Prepare the initial string
	    	string locations = "Working in locations: ";
			// Dump into files those reached regions associated with newly-found not-in-use locations, and
	    	// load from files those reached regions associated with newly-found in-use locations
			for (std::map<DiscreteState,bool>::const_iterator loc_it = currentLocationsInUse.begin(); loc_it != currentLocationsInUse.end(); loc_it++) {
				// Add the location in use
				if (loc_it->second)
					locations += loc_it->first.name() + " ";

				// Get the file name string as [systemname].[locationname].dump
				string filename = system.name() + "." + loc_it->first.name() + ".dump";
				const char * filename_ch = filename.c_str();
				// If it is not currently in use but was in use before, dump it
				// If it is currently in use but was not in use before, load it
				if (!loc_it->second && previousLocationsInUse[loc_it->first]) {
					ARIADNE_LOG(6,"Dumping reached region for location " << loc_it->first.name() << "\n");
					reach[loc_it->first].export_to_file(filename_ch);
				}
				else if (loc_it->second && !previousLocationsInUse[loc_it->first]) {
					ARIADNE_LOG(6,"Restoring reached region for location " << loc_it->first.name() << "\n");
					reach[loc_it->first].import_from_file(filename_ch);
				}
			}
			ARIADNE_LOG(6,locations << "\n");
			// Copy the current locations in use into the previous ones
			previousLocationsInUse = currentLocationsInUse;
	    }

        // Reset (with false values) the information about locations in use on the current iteration
        std::map<DiscreteState,bool> nextLocationsInUse;
        for (HybridGrid::locations_const_iterator grid_it = grid.locations_begin(); grid_it != grid.locations_end(); grid_it++)
        	currentLocationsInUse[grid_it->first] = false;

		// Evolve the initial enclosures
        make_lpair(new_reach,new_final)=_upper_reach_evolve_continuous(system,initial_enclosures,hybrid_lock_to_grid_time,maximum_grid_depth);
		// Clear the initial enclosures
		initial_enclosures.clear();

		// Remove from the result those cells that have already been evolved from or checked for activations
        new_final.remove(intermediate);
		new_reach.remove(reach);

	    // If dumping must be performed
	    if (chain_reach_dumping) {
			// Check the new final region and set as in-use all the locations whose reached region is not empty
			for (GTS::locations_iterator loc_it = new_final.locations_begin(); loc_it != new_final.locations_end(); loc_it++)
				currentLocationsInUse[loc_it->first] = !loc_it->second.empty();
	    }

        ARIADNE_LOG(6,"Reach size after removal = "<<new_reach.size()<<"\n");
        ARIADNE_LOG(6,"Final size after removal = "<<new_final.size()<<"\n");

		// Mince the final cells, then for each of them add the enclosure into the new initial enclosures
		new_final.mince(maximum_grid_depth);
		for (GTS::const_iterator cell_it = new_final.begin(); cell_it != new_final.end(); cell_it++)
			initial_enclosures.push_back(_discretiser->enclosure(*cell_it));

		// Mince the reach cells, then for each of them and for each transition, check the activation:
		// If it is possibly active, create the enclosure, apply the reset and put the result into the new initial enclosures
		new_reach.mince(maximum_grid_depth);
		ARIADNE_LOG(6,"Reach size after mincing = "<<new_reach.size()<<"\n");
		ARIADNE_LOG(6,"Final size after mincing = "<<new_final.size()<<"\n");
		for (GTS::const_iterator cell_it = new_reach.begin(); cell_it != new_reach.end(); cell_it++)
		{
			// Get the transitions for the corresponding location
			list<DiscreteTransition> transitions = system.transitions(cell_it->first);
			// Get the enclosure of the cell
			EnclosureType encl = _discretiser->enclosure(*cell_it);

			ARIADNE_LOG(7,"Enclosure = "<<encl<<"\n");

			// For each of them
			for (list<DiscreteTransition>::const_iterator trans_it = transitions.begin(); trans_it != transitions.end(); trans_it++)
			{
				ARIADNE_LOG(8,"Transition = "<<*trans_it<<"\n");
				// Get the time model of the guard
    			TaylorModel guard_time_model = apply(trans_it->activation(),encl.continuous_state_set())[0];

    			ARIADNE_LOG(8,"Guard time model = "<<guard_time_model<<" (range: " << guard_time_model.range() << ")\n");

				// If possibly active, get the set model after the reset, enclose it and add it to the initial enclosures
            	if (!(guard_time_model.range().upper()<0)) {
                	TaylorSet jump_set_model= tc.reset_step(trans_it->reset(),encl.continuous_state_set());
                	initial_enclosures.push_back(EnclosureType(trans_it->target(),jump_set_model));
                	currentLocationsInUse[trans_it->target()] = true; // Set the in-use value for the target location as true
            	}	
			}
		}

		ARIADNE_LOG(6,"Enclosures size after transition checking = "<<initial_enclosures.size()<<"\n");

		// Adjoin to the accumulated reach and intermediate cells
        reach.adjoin(new_reach);
		intermediate.adjoin(new_final);

		// Recombine both
		reach.recombine();
		intermediate.recombine();
    }

    // Restore those reached regions dumped in the previous iteration
    // (please note that currentLocationsInUse must feature false values only)
    if (chain_reach_dumping) {
		for (std::map<DiscreteState,bool>::const_iterator loc_it = previousLocationsInUse.begin(); loc_it != previousLocationsInUse.end(); loc_it++) {
			// Get the file name string as [systemname].[locationname].dump
			string filename = system.name() + "." + loc_it->first.name() + ".dump";
			const char * filename_ch = filename.c_str();
			// If the location was not used
			if (!loc_it->second) {
				ARIADNE_LOG(6,"Restoring reached region for location " << loc_it->first.name() << "\n");
				reach[loc_it->first].import_from_file(filename_ch);
			}
		}
    }

    // Save the reached region
	_statistics->upper().reach = reach;

	ARIADNE_LOG(5,"Found a total of " << reach.size() << " reached cells.\n");

    return reach;
}

HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_chain_reach(const SystemType& system,
				  const HybridImageSet& initial_set) const
{
	typedef std::list<EnclosureType> EL;
	typedef std::map<DiscreteState,uint> HUM;

    ARIADNE_LOG(4,"HybridReachabilityAnalyser::lower_chain_reach(system,initial_set)\n");

	_statistics->lower().reset(); // Reset the discrete statistics
	_discretiser->reset_lower_statistics(); // Reset the continuous statistics

	// Get the actual concurrency
	const uint concurrency = boost::thread::hardware_concurrency() - free_cores;
	// Check that the concurrency is greater than zero
	ARIADNE_ASSERT_MSG(concurrency>0,"Error: concurrency must be at least 1.");

	// The global reach region
    GTS reach(system.grid());

	// Get the lock time, in the case pruning is to be performed
	TimeType lock_time(_parameters->lock_to_grid_time,_parameters->lock_to_grid_steps);

    // Split the initial set into enclosures that are smaller than the minimum cell, given the grid
    EL initial_enclosures = _split_initial_set(initial_set,system.grid());

	// Perform loop
	for (uint i=0;;i++)
	{
		ARIADNE_LOG(5,"Iteration " << i << "\n");

		// The final enclosures
		EL final_enclosures;

		// The sizes of the adjoined (the cells) or superposed (the enclosures) evolve
		std::map<DiscreteState,uint> adjoined_evolve_sizes;
		std::map<DiscreteState,uint> superposed_evolve_sizes;
		// The evolve
		GTS evolve;
		// The new reach
		GTS new_reach;

		ARIADNE_LOG(6,"Initial enclosures size = " << initial_enclosures.size() << "\n");

		// Create the worker
		LowerReachWorker worker(_discretiser,initial_enclosures,system,lock_time,
								_parameters->bounding_domain,_parameters->maximum_grid_depth,concurrency);

		ARIADNE_LOG(6,"Evolving and discretising...\n");

		// Compute and get the result
		make_ltuple<HUM,HUM,EL,GTS>(adjoined_evolve_sizes,superposed_evolve_sizes,final_enclosures,new_reach) = worker.get_result();

		ARIADNE_LOG(6,"Reach size before removal = " << new_reach.size() << "\n");

		// Remove from the partial reached region the global reached region
		new_reach.remove(reach);

		ARIADNE_LOG(6,"Reach size after removal  = " << new_reach.size() << "\n");

		// Check for inclusion: if no new cells are found, terminate
		if (new_reach.empty())
			break;

		// Otherwise proceed and adjoin the new reach
		reach.adjoin(new_reach);

		// Pruning of the dump of the final enclosures into the initial enclosures
		while (!final_enclosures.empty())
		{
			// Pop the current enclosure
			EnclosureType encl = final_enclosures.front();
			final_enclosures.pop_front();

			/* If pruning is to be performed, push just a fraction of the final_enclosures into the initial_enclosures;
			 * otherwise, push indiscriminately.
			 */
			if (_parameters->enable_lower_pruning)
			{
				// Get the ratio between the adjoined evolve size and the superposed evolve size
				Float ratio = (Float)adjoined_evolve_sizes[encl.location()]/(Float)superposed_evolve_sizes[encl.location()];

				// At least 2 enclosures are inserted, then the enclosures are pruned as long as the number of enclosures is at least twice the number of evolve cells
				if (initial_enclosures.size() <= 2 || rand() < 2*ratio*RAND_MAX)
					initial_enclosures.push_back(encl);
			}
			else
				initial_enclosures.push_back(encl);
		}
	}

	// Copies the reached region
	_statistics->lower().reach = reach;

	return reach;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
viable(const SystemType& system,
       const HybridImageSet& bounding_set) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


bool 
HybridReachabilityAnalyser::
_safe(const SystemType& system, 
	  const HybridImageSet& initial_set, 
	  const HybridBoxes& safe_box)
{
	ARIADNE_LOG(4,"Safety analysis...\n");

	HybridGridTreeSet reach = upper_chain_reach(system,initial_set); // Perform the chain reachability analysis

	// Get the bounds of the reach set
	HybridBoxes bounds = reach.bounding_box();

	// Check for each location that the bounds are included into the domain
	bool is_inside_domain = true; // Initial value for the reached set being inside the domain
	for (HybridBoxes::const_iterator loc_it = bounds.begin(); loc_it != bounds.end(); loc_it++) {
		if (!(loc_it->second.subset(_parameters->bounding_domain[loc_it->first]))) {
			is_inside_domain = false;
			break;
		}
	}

	// If the reached region is inside the domain
	if (is_inside_domain)
	{
		// If the reached region is definitely inside the hybrid safe box, the result is safe 
		bool result = definitely(reach.subset(safe_box));

		ARIADNE_LOG(4, (result ? "Safe.\n" : "Not safe.\n") );

		return result;
	}
	// Otherwise notify and return false
	else
	{
		ARIADNE_LOG(4,"Not checked due to domain bounds.\n");
		return false;
	}
}


bool 
HybridReachabilityAnalyser::
_unsafe(const SystemType& system, 
		const HybridImageSet& initial_set, 
		const HybridBoxes& safe_box)
{
	ARIADNE_LOG(4,"Unsafety analysis...\n");

	// Get the size of the continuous space (NOTE: assumed equal for all locations)
	const uint css = system.state_space().locations_begin()->second; 

	// Get the lower reach
	HybridGridTreeSet lowerreach = lower_chain_reach(system,initial_set);

	ARIADNE_LOG(5, "Largest enclosure cell: " << _discretiser->statistics().lower().largest_enclosure_cell << "\n");

	if (lowerreach.size() > 0)
	{
		// Create a safe set enlarged for the falsification
		HybridBoxes safe_box_enl = safe_box;
		for (HybridBoxes::iterator hbx_it = safe_box_enl.begin(); hbx_it != safe_box_enl.end(); hbx_it++)
			hbx_it->second.widen();

		// Get a copy of the grid
		HybridGrid hg = system.grid();
		// Get the safe box enlarged by a grid cell and half a maximum enclosure cell
		for (HybridGrid::locations_const_iterator hg_it = hg.locations_begin(); hg_it != hg.locations_end(); hg_it++)
			for (uint i=0;i<css;i++)
				safe_box_enl[hg_it->first][i] += Interval(-_discretiser->statistics().lower().largest_enclosure_cell[i]/2,_discretiser->statistics().lower().largest_enclosure_cell[i]/2)+Interval(-hg_it->second.lengths()[i],hg_it->second.lengths()[i])/(1<<_parameters->maximum_grid_depth);

		// If the reach is definitely not included in the enlarged safe box, then it is unsafe
		bool result = definitely(!lowerreach.subset(safe_box_enl));

		ARIADNE_LOG(4, (result ? "Unsafe.\n" : "Not unsafe.\n") );

		return result;
	}
	else
		return false;
}


tribool
HybridReachabilityAnalyser::
verify(const SystemType& system,
       const HybridImageSet& initial_set,
       const HybridBoxes& safe_box)
{
		ARIADNE_LOG(3, "Verification...\n");

		// Perform the safety analysis
		if (_safe(system,initial_set,safe_box))
		{
			ARIADNE_LOG(3, "Safe.\n");
			return true;
		}

		// Perform the unsafety analysis
		if (_unsafe(system,initial_set,safe_box))
		{
			ARIADNE_LOG(3, "Unsafe.\n");
			return false;
		}

		ARIADNE_LOG(3, "Indeterminate.\n");

		// Return indeterminate if both failed
		return indeterminate;
}


HybridFloatVector 
HybridReachabilityAnalyser::
_getHybridMaximumAbsoluteDerivatives(const HybridAutomaton& system) const
{
	// Gets the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = system.state_space().locations_begin()->second;

	// The variable for the bounding box of the derivatives
	Vector<Interval> der;
    // The variable for the result
	HybridFloatVector hmad;

	// Get the HybridGridTreeSet
	HybridGridTreeSet& hybridreach = _statistics->upper().reach;

	// For each mode
	for (list<DiscreteMode>::const_iterator modes_it = system.modes().begin(); modes_it != system.modes().end(); modes_it++)
	{
		// Get the location
		const DiscreteState& loc = modes_it->location();

		// Insert the corresponding hmad pair, initialized with zero maximum absolute derivatives
		hmad.insert(pair<DiscreteState,Vector<Float> >(loc,Vector<Float>(css)));

		// If the reached region for the location exists and is not empty, check its cells, otherwise use the whole domain
		if (hybridreach.has_location(loc) && !hybridreach[loc].empty()) {
			// Get the GridTreeSet
			GridTreeSet& reach = hybridreach[loc];
			// For each of its hybrid cells
			for (GridTreeSet::const_iterator cells_it = reach.begin(); cells_it != reach.end(); cells_it++)
			{
				// Gets the derivative bounds
				der = modes_it->dynamic()(cells_it->box());

				// For each variable, sets the maximum value
				for (uint i=0;i<css;i++)
					hmad[loc][i] = max(hmad[loc][i], abs(der[i]).upper());
			}
		} else {
			// Gets the first order derivatives in respect to the dynamic of the mode, applied to the domain of the corresponding location
			der = modes_it->dynamic()(_parameters->bounding_domain.find(loc)->second);

			// Gets the maximum absolute derivatives
			for (uint i=0;i<css;i++)
				hmad[loc][i] = abs(der[i]).upper();
		}
	}

	// Returns
	return hmad;
}


void
HybridReachabilityAnalyser::
_setMaximumEnclosureCell(const HybridGrid& hgrid)
{
	// Gets the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = hgrid.locations_begin()->second.lengths().size();

	// Initializes the result
	Vector<Float> mec(css);

	// For each location and dimension of the space
	for (HybridGrid::locations_const_iterator hg_it = hgrid.locations_begin(); hg_it != hgrid.locations_end(); hg_it++)
		for (uint i=0;i<css;i++)
			if (hg_it->second.lengths()[i] > mec[i])
				mec[i] = hg_it->second.lengths()[i];

	// Scales the cell in respect to the maximum grid depth
	for (uint i=0;i<css;i++)
		mec[i] /= (1<<_parameters->maximum_grid_depth);

	// Assigns (NOTE: it is preferable to have the multiplier slightly lesser than an integer multiple of the grid cell)
	_discretiser->parameters().maximum_enclosure_cell = 1.9*mec;
}


void 
HybridReachabilityAnalyser::
_setEqualizedHybridMaximumStepSize(const HybridFloatVector& hmad, const HybridGrid& hgrid)
{
	// Gets the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = hmad.begin()->second.size();

	// Initialize the hybrid maximum step size
	std::map<DiscreteState,Float> hmss;

	// Initializes the maximum step size
	Float mss = 0.0;
	// For each couple DiscreteState,Vector<Float>
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) 
	{
		// For each dimension of the space, if the derivative is not zero,
		// evaluates the ratio between the minimum cell length and the derivative itself
		for (uint i=0;i<css;i++)
			if (hfv_it->second[i] > 0)
				mss = max(mss,hgrid[hfv_it->first].lengths()[i]/(1 << _parameters->maximum_grid_depth)/hfv_it->second[i]);
	}

	// Populates the map (twice the values since the maximum enclosure is set as ~2 times the grid cell)
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
		hmss.insert(std::pair<DiscreteState,Float>(hfv_it->first,2*mss));

	// Assigns
	_discretiser->parameters().hybrid_maximum_step_size = hmss;
}


void
HybridReachabilityAnalyser::
_setHybridMaximumStepSize(const HybridFloatVector& hmad, const HybridGrid& hgrid)
{
	// Gets the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = hmad.begin()->second.size();

	// Initialize the hybrid maximum step size
	std::map<DiscreteState,Float> hmss;

	// For each couple DiscreteState,Vector<Float>
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
	{
		// Initializes the maximum step size
		Float mss = 0.0;
		// For each dimension of the space, if the derivative is not zero,
		// evaluates the ratio between the minimum cell length and the derivative itself
		for (uint i=0;i<css;i++)
			if (hfv_it->second[i] > 0)
				mss = max(mss,hgrid[hfv_it->first].lengths()[i]/(1 << _parameters->maximum_grid_depth)/hfv_it->second[i]);

		// Inserts the value (twice the value since the maximum enclosure is set as ~2 the grid cell)
		hmss.insert(std::pair<DiscreteState,Float>(hfv_it->first,2*mss));
	}

	// Assigns
	_discretiser->parameters().hybrid_maximum_step_size = hmss;
}


HybridGrid 
HybridReachabilityAnalyser::
_getHybridGrid(const HybridFloatVector& hmad) const
{
	// Get the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = hmad.begin()->second.size();

	// Initialize the hybrid grid
	HybridGrid hg;	

	// For each couple DiscreteState,Vector<Float>
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
	{
		// Keep the location
		const DiscreteState& loc = hfv_it->first;

		// Initialize the maximum derivative/domain_width ratio
		Float maxratio = 0.0;
		// Initialize the gridlengths
		Vector<Float> gridlengths(css);
		// Initialize the minimum length to use for zero-derivative variables (will be the minimum among both grid lengths and domain widths)
		Float minlength = std::numeric_limits<double>::infinity();

		// For each dimension of the continuous space
		for (uint i=0;i<css;i++)
		{	
			maxratio = max(maxratio,hfv_it->second[i]/_parameters->bounding_domain.find(loc)->second[i].width()); // Get the largest ratio
			minlength = min(minlength,_parameters->bounding_domain.find(loc)->second[i].width()); // Get the smallest domain width
		}

		// Assign the lengths and check the minimum length to be assigned
		for (uint i=0;i<css;i++)
		{
			if (hfv_it->second[i] > 0) // If the derivative is greater than zero
			{
				gridlengths[i] = hfv_it->second[i]/maxratio; // Assign the length
				minlength = min(minlength,gridlengths[i]); // Reduce the minimum length by comparing it to the current grid length
			}
		}

		// Assign the minimum length for zero-derivative variables
		for (uint i=0;i<css;i++)
			if (gridlengths[i] == 0) // If it has zero grid length, assign minlength
				gridlengths[i] = minlength;

		/* Choose the center for the grid: the center of the reached region, if a non-empty reached region exists,
		 * otherwise the center of the domain */
		Vector<Float> centre(css);
		if (_statistics->upper().reach.has_location(loc) && !_statistics->upper().reach[loc].empty()) {
			centre = _statistics->upper().reach[loc].bounding_box().centre();
		} else {
			centre = _parameters->bounding_domain.find(loc)->second.centre();
		}
		hg[loc] = Grid(centre,gridlengths);
	}

	// Return
	return hg;
}


HybridGrid 
HybridReachabilityAnalyser::
_getEqualizedHybridGrid(const HybridFloatVector& hmad) const
{
	// Get the size of the space (NOTE: taken as equal for all locations)
	const uint css = hmad.begin()->second.size();

	// Initialize the maximum derivative/domain_width ratio
	Float maxratio = 0.0;
	// Initialize the minimum length to use for zero-derivative variables
	Float minlength = std::numeric_limits<double>::infinity();
	// Initialize the gridlengths
	Vector<Float> gridlengths(css);		
	// Initialize the maximum absolute derivatives
	Vector<Float> mad(css);

	// Get the maximum absolute derivatives
	for (uint i=0;i<css;i++)
		// For each couple DiscreteState,Vector<Float>
		for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
			mad[i] = max(mad[i],hfv_it->second[i]);

	// For each dimension
	for (uint i=0;i<css;i++)
		// For each couple DiscreteState,Vector<Float>
		for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
		{
			maxratio = max(maxratio,mad[i]/_parameters->bounding_domain.find(hfv_it->first)->second[i].width()); // Check for maximum ratio
			minlength = min(minlength,_parameters->bounding_domain.find(hfv_it->first)->second[i].width()); // Get the minimum domain length
		}

	// Assign the lengths and checks the minimum length to be assigned
	for (uint i=0;i<css;i++)
	{
		// If the derivative is greater than zero
		if (mad[i] > 0)
		{	
			gridlengths[i] = mad[i]/maxratio; // Assign the length
			minlength = min(minlength,gridlengths[i]); // Get the minimum between the current minimum and the current length
		}		
	}

	// Correct the lengths that are zero
	for (uint i=0;i<css;i++)
		// If it has zero grid length, assign minlength
		if (gridlengths[i] == 0)
			gridlengths[i] = minlength;

	// Initialize the hybrid grid
	HybridGrid hg;
	// Populate it, centered on the centre of the domain
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
		hg[hfv_it->first] = Grid(_parameters->bounding_domain.find(hfv_it->first)->second.centre(),gridlengths);
	
	// Return
	return hg;
} 


void 
HybridReachabilityAnalyser::
_setInitialParameters(SystemType& system, const HybridBoxes& domain)
{
	// Set the domain
	_parameters->bounding_domain = domain;
}


void 
HybridReachabilityAnalyser::
_tuneParameters(SystemType& system)
{
	// Evaluate the maximum absolute derivatives
	HybridFloatVector hmad = _getHybridMaximumAbsoluteDerivatives(system);
	ARIADNE_LOG(3, "Hybrid maximum absolute derivatives: " << hmad << "\n");
	// Grid
	system.set_grid(_getHybridGrid(hmad));
	ARIADNE_LOG(3, "Hybrid grid: " << system.grid() << "\n");
	// Maximum enclosure cell
	_setMaximumEnclosureCell(system.grid());
	ARIADNE_LOG(3, "Maximum enclosure cell: " << _discretiser->parameters().maximum_enclosure_cell << "\n");
	// Maximum step size
	_setHybridMaximumStepSize(hmad,system.grid());
	ARIADNE_LOG(3, "Hybrid maximum step size: " << _discretiser->parameters().hybrid_maximum_step_size << "\n");
}


tribool 
HybridReachabilityAnalyser::
verify_iterative(SystemType& system, 
				 const HybridImageSet& initial_set, 
				 const HybridBoxes& safe_box, 
				 const HybridBoxes& domain)
{
	ARIADNE_LOG(2,"\nIterative verification...\n");

	time_t mytime;
	time(&mytime);
	string foldername = system.name()+"-png";
	char mgd_char[10];
	string filename;

	// Save the folder name as a function of the automaton name and of the current timestamp, then create the main folder and the verification run folder
	if (plot_verify_results)
	{
		mkdir(foldername.c_str(),0777);
		string timestring = asctime(localtime(&mytime));
		timestring.erase(std::remove(timestring.begin(), timestring.end(), '\n'), timestring.end());
		foldername = foldername+"/"+timestring;
		mkdir(foldername.c_str(),0777);
	}

	// Set the initial parameters
	_setInitialParameters(system, domain);

    for(_parameters->maximum_grid_depth = _parameters->lowest_maximum_grid_depth;
    	_parameters->maximum_grid_depth <= _parameters->highest_maximum_grid_depth;
    	_parameters->maximum_grid_depth++)
	{ 
    	sprintf(mgd_char,"%i", _parameters->maximum_grid_depth);
		ARIADNE_LOG(2, "DEPTH " << _parameters->maximum_grid_depth << "\n");

		// Tune the parameters for the current iteration
		_tuneParameters(system);

		// Perform the verification
		tribool result = verify(system,initial_set,safe_box);

		// Plot the results, if desired
		if (plot_verify_results)
		{
			// Plot the upper region
			filename = "upper-";
			plot(foldername,filename + mgd_char, _statistics->upper().reach);
			// Plot the lower region, if the result is not true (otherwise the lower region has not been computed on this iteration)
			filename = "lower-";
			if (!definitely(result))
				plot(foldername,filename + mgd_char, _statistics->lower().reach);
		}

		// Return the result, if it is not indeterminate
		if (!indeterminate(result)) return result;

    }

	// Return indeterminate
	return indeterminate;
}

Interval
HybridReachabilityAnalyser::
safety_parametric(SystemType& system, 
				  const HybridImageSet& initial_set, 
				  const HybridBoxes& safe_box, 
				  const HybridBoxes& domain,
				  const RealConstant& parameter,
				  const Interval& parameter_interval,
				  const Float& tolerance)	
{
	// Copy the parameter and interval for local operations
	RealConstant param = parameter;
	Interval param_int = parameter_interval;

	// Check the lower bound
	ARIADNE_LOG(1,"\nChecking lower interval bound... ");

	// Set the parameter
	param.set_value(parameter_interval.lower());
	// Substitute the value
	system.substitute(param);
	// Perform the verification
	tribool lower_result = verify_iterative(system,initial_set,safe_box,domain);

	if (definitely(lower_result)) { ARIADNE_LOG(1,"Safe.\n"); }
	else { ARIADNE_LOG(1,"Not safe.\n"); }

	// Check the upper bound
	ARIADNE_LOG(1,"Checking upper interval bound... ");

	// Set the parameter
	param.set_value(parameter_interval.upper());
	// Substitute the value
	system.substitute(param);
	// Perform the verification
	tribool upper_result = verify_iterative(system,initial_set,safe_box,domain);

	if (definitely(upper_result)) { ARIADNE_LOG(1,"Safe.\n"); }
	else { ARIADNE_LOG(1,"Not safe.\n"); }

	// Analyse the results

	// The source of update
	bool updateFromBottom;

	// Create an empty interval
	Interval empty_int;
	empty_int.make_empty();

	// If both extremes are definitely safe, no more verification is involved
	if (definitely(lower_result) && definitely(upper_result)) {
		return parameter_interval;
	}
	// If both extremes are not definitely safe, no more verification is involved
	else if (!definitely(lower_result) && !definitely(upper_result)) {
		return empty_int;
	}
	// Otherwise it updates from the bottom or the top depending on the lower_result being safe or not
	else updateFromBottom = definitely(lower_result);		

	// While the tolerance bound has not been hit
	while (param_int.width() > tolerance)
	{
		// Set the parameter as the midpoint of the interval
		param.set_value(param_int.midpoint());
		// Substitute the value
		system.substitute(param);

		ARIADNE_LOG(1,"Checking " << param_int << " (midpoint: " << param_int.midpoint() << ", width: " << param_int.width() << ") ... ");

		// Perform the verification
		tribool result = verify_iterative(system,initial_set,safe_box,domain);

		if (definitely(result)) {
			if (updateFromBottom) {
				ARIADNE_LOG(1,"Safe, refining upwards.\n");
				param_int.set_lower(param.value()); }
			else {
				ARIADNE_LOG(1,"Safe, refining downwards.\n");
				param_int.set_upper(param.value()); }}
		else {
			if (updateFromBottom) {
				ARIADNE_LOG(1,"Not safe, refining downwards.\n");
				param_int.set_upper(param.value()); }
			else {
				ARIADNE_LOG(1,"Not safe, refining upwards.\n");
				param_int.set_lower(param.value()); }}}

	if (updateFromBottom)
		return Interval(parameter_interval.lower(),param_int.lower());
	else
		return Interval(param_int.upper(),parameter_interval.upper());
}

std::pair<Interval,Interval>
HybridReachabilityAnalyser::
safety_unsafety_parametric(SystemType& system, 
						   const HybridImageSet& initial_set, 
						   const HybridBoxes& safe_box, 
						   const HybridBoxes& domain,
						   const RealConstant& parameter,
						   const Interval& parameter_interval,
						   const Float& tolerance)	
{
	// Copy the parameter for local operations
	RealConstant param = parameter;
	// Create the safety and unsafety intervals: they represent the search intervals, not the intervals where the system is proved safe or unsafe
	Interval safety_int = parameter_interval;
	Interval unsafety_int = parameter_interval;

	// Check the lower bound
	ARIADNE_LOG(1,"\nChecking lower interval bound... ");

	// Set the parameter
	param.set_value(parameter_interval.lower());
	// Substitute the value
	system.substitute(param);
	// Perform the verification
	tribool lower_result = verify_iterative(system,initial_set,safe_box,domain);

	if (definitely(lower_result)) { ARIADNE_LOG(1,"Safe.\n"); }
	else if (!possibly(lower_result)) { ARIADNE_LOG(1,"Unsafe.\n"); }
	else ARIADNE_LOG(1,"Indeterminate.\n");

	// Check the upper bound
	ARIADNE_LOG(1,"Checking upper interval bound... ");

	// Set the parameter
	param.set_value(parameter_interval.upper());
	// Substitute the value
	system.substitute(param);
	// Perform the verification
	tribool upper_result = verify_iterative(system,initial_set,safe_box,domain);

	if (definitely(upper_result)) { ARIADNE_LOG(1,"Safe.\n"); }
	else if (!possibly(upper_result)) { ARIADNE_LOG(1,"Unsafe.\n"); }
	else ARIADNE_LOG(1,"Indeterminate.\n");

	// Analyze the results

	// Where the safe value is found
	bool safeOnBottom;

	// Create an empty interval
	Interval empty_int;
	empty_int.make_empty();

	// If both extremes are safe, no more verification is involved
	if (definitely(lower_result) && definitely(upper_result)) {
		return make_pair<Interval,Interval>(parameter_interval,empty_int); }
	// If both extremes are unsafe, no more verification is involved
	else if (!possibly(lower_result) && !possibly(upper_result)) {
		return make_pair<Interval,Interval>(empty_int,parameter_interval); }
	// If both extremes are indeterminate, no verification is possible
	else if (indeterminate(lower_result) && indeterminate(upper_result)) {
		return make_pair<Interval,Interval>(empty_int,empty_int); }
	// If the lower extreme is safe or the upper extreme is unsafe, the safe values are on the bottom
	else if (definitely(lower_result) || !possibly(upper_result)) {
		safeOnBottom = true;		
		// If there are indeterminate values, reset the corresponding intervals as empty
		if (indeterminate(lower_result)) safety_int = empty_int;
		if (indeterminate(upper_result)) unsafety_int = empty_int; }		
	// If the upper extreme is safe or the lower extreme is unsafe, the safe values are on the top
	else {
		safeOnBottom = false;		
		// If there are indeterminate values, reset the corresponding intervals as empty
		if (indeterminate(lower_result)) unsafety_int = empty_int;
		if (indeterminate(upper_result)) safety_int = empty_int; }

	// Verification loop
	while (true) 
	{
		// The verification result
		tribool result;

		// Safety interval check
		if (!safety_int.empty()) 
		{
			// Set the parameter as the midpoint of the interval
			param.set_value(safety_int.midpoint());
			// Substitute the value
			system.substitute(param);

			ARIADNE_LOG(1,"Checking safety interval " << safety_int << " (midpoint: " << safety_int.midpoint() << ", width: " << safety_int.width() << ") ... ");

			// Perform the verification
			result = verify_iterative(system,initial_set,safe_box,domain);

			// If safe
			if (definitely(result)) {
				if (safeOnBottom) {
					// If the unsafety interval is the same as the safety interval, update it too
					if (equal(unsafety_int,safety_int)) unsafety_int.set_lower(param.value());

					ARIADNE_LOG(1,"Safe, refining upwards.\n");
					safety_int.set_lower(param.value()); }
				else {
					// If the unsafety interval is the same as the safety interval, update it too
					if (equal(unsafety_int,safety_int)) unsafety_int.set_upper(param.value());

					ARIADNE_LOG(1,"Safe, refining downwards.\n");
					safety_int.set_upper(param.value()); }}
			// If unsafe
			else if (!possibly(result)) {
				if (safeOnBottom) {
					ARIADNE_LOG(1,"Unsafe, refining downwards and resetting the unsafety.\n");
					safety_int.set_upper(param.value()); }
				else {
					ARIADNE_LOG(1,"Unsafe, refining upwards and resetting the unsafety.\n");
					safety_int.set_lower(param.value()); }

				// The unsafety interval now becomes the same as the safety interval
				unsafety_int = safety_int; }
			// If indeterminate
			else {
				if (safeOnBottom) {
					// If the unsafety interval is the same as the safety interval, update it too
					if (equal(unsafety_int,safety_int)) unsafety_int.set_lower(param.value());

					ARIADNE_LOG(1,"Indeterminate, refining downwards.\n");
					safety_int.set_upper(param.value()); }
				else {
					// If the unsafety interval is the same as the safety interval, update it too
					if (equal(unsafety_int,safety_int))	unsafety_int.set_upper(param.value());

					ARIADNE_LOG(1,"Indeterminate, refining upwards.\n");
					safety_int.set_lower(param.value()); }}

			/* Break if the safety interval is lesser than the tolerance or, if the unsafety interval is not empty,
			 if the minimum distance between the safe and unsafe values is lesser than the tolerance (which is a more relaxed condition) */
			if ((safety_int.width() <= tolerance) || (!unsafety_int.empty() &&
				((safeOnBottom && (unsafety_int.upper() - safety_int.lower() <= tolerance)) ||
				(!safeOnBottom && (safety_int.upper() - unsafety_int.lower() <= tolerance)))))
				break;
		}

		// Unsafety interval check
		if (!unsafety_int.empty())
		{
			// Set the parameter as the midpoint of the interval
			param.set_value(unsafety_int.midpoint());
			// Substitute the value
			system.substitute(param);

			ARIADNE_LOG(1,"Checking unsafety interval " << unsafety_int << " (midpoint: " << unsafety_int.midpoint() << ", width: " << unsafety_int.width() << ") ... ");

			// Perform the verification
			result = verify_iterative(system,initial_set,safe_box,domain);

			// If safe
			if (definitely(result)) {
				if (safeOnBottom) {
					ARIADNE_LOG(1,"Safe, refining upwards and resetting the safety.\n");
					unsafety_int.set_lower(param.value()); }
				else {
					ARIADNE_LOG(1,"Safe, refining downwards and resetting the safety.\n");
					unsafety_int.set_upper(param.value()); }

				// The safety interval now becomes the same as the unsafety interval
				safety_int = unsafety_int; }
			// If unsafe
			else if (!possibly(result)) {
				if (safeOnBottom) {
					// If the unsafety interval is the same as the safety interval, update it too
					if (equal(unsafety_int,safety_int)) safety_int.set_upper(param.value());

					ARIADNE_LOG(1,"Unsafe, refining downwards.\n");
					unsafety_int.set_upper(param.value()); }
				else {
					// If the unsafety interval is the same as the safety interval, update it too
					if (equal(unsafety_int,safety_int)) safety_int.set_lower(param.value());

					ARIADNE_LOG(1,"Unsafe, refining upwards.\n");
					unsafety_int.set_lower(param.value()); }}
			// If indeterminate
			else {
				if (safeOnBottom) {
					// If the unsafety interval is the same as the safety interval, update it too
					if (equal(unsafety_int,safety_int)) safety_int.set_upper(param.value());

					ARIADNE_LOG(1,"Indeterminate, refining upwards.\n");
					unsafety_int.set_lower(param.value()); }
				else {
					// If the unsafety interval is the same as the safety interval, update it too
					if (equal(unsafety_int,safety_int)) safety_int.set_lower(param.value());

					ARIADNE_LOG(1,"Indeterminate, refining downwards.\n");
					unsafety_int.set_upper(param.value()); }}

			/* Break if the safety interval is lesser than the tolerance or, if the unsafety interval is not empty,
			 if the minimum distance between the safe and unsafe values is lesser than the tolerance (which is a more relaxed condition) */
			if ((unsafety_int.width() <= tolerance) || (!safety_int.empty() &&
				((safeOnBottom && (unsafety_int.upper() - safety_int.lower() <= tolerance)) ||
				(!safeOnBottom && (safety_int.upper() - unsafety_int.lower() <= tolerance)))))
				break;
		}
	}

	// The result intervals for safe and unsafe values
	Interval safe_result, unsafe_result;

	// Get the safe and unsafe intervals
	if (safeOnBottom) {
		safe_result = (safety_int.empty() ? safety_int : Interval(parameter_interval.lower(),safety_int.lower()));
		unsafe_result = (unsafety_int.empty() ? unsafety_int : Interval(unsafety_int.upper(),parameter_interval.upper())); }	
	else {
		safe_result = (safety_int.empty() ? safety_int : Interval(safety_int.upper(),parameter_interval.upper()));	
		unsafe_result = (unsafety_int.empty() ? unsafety_int : Interval(parameter_interval.lower(),unsafety_int.lower())); }

	return make_pair<Interval,Interval>(safe_result,unsafe_result);
}


} // namespace Ariadne
