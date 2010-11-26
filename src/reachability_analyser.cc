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


void
HybridReachabilityAnalyser::
_splitAndCreateTargetEnclosures(bool& isValid,
								std::list<EnclosureType>& initial_enclosures,
								const ContinuousEnclosureType& encl,
								const Vector<Float>& minTargetCellWidths,
								const Box& target_bounding,
								const DiscreteTransition& trans,
								const TaylorCalculus& tc) const
{
	// Get the target enclosure
	ContinuousEnclosureType target_encl = tc.reset_step(trans.reset(),encl);
	// Get its box
	const Box& target_encl_box = target_encl.bounding_box();

	// If the cell box lies outside the bounding domain (i.e. not inside and disjoint)
	// of the target location, ignore it and notify
	if (!target_encl_box.inside(target_bounding) &&
		target_encl_box.disjoint(target_bounding)) {
		ARIADNE_LOG(7,"Discarding target enclosure " << target_encl_box <<
					  " being outside the domain in location " << trans.target().name() << ".\n");
		isValid = false;
		return;
	}

	// Otherwise try splitting the enclosure
	bool hasSplit = false;
	for (uint i=0; i < minTargetCellWidths.size(); i++) {
		// If the enclosure has width larger than that of the minimum cell, split on that dimension
		if (encl.bounding_box()[i].width() > minTargetCellWidths[i]) {
			hasSplit = true;
			std::pair<ContinuousEnclosureType,ContinuousEnclosureType> split_sets = encl.split(i);
			// Call recursively on the two enclosures
			_splitAndCreateTargetEnclosures(isValid,initial_enclosures,split_sets.first,minTargetCellWidths,target_bounding,trans,tc);
			_splitAndCreateTargetEnclosures(isValid,initial_enclosures,split_sets.second,minTargetCellWidths,target_bounding,trans,tc);
		}
	}

	// If we could not split, we put the hybrid enclosure into the new initial_enclosures
	if (!hasSplit)
		initial_enclosures.push_back(EnclosureType(trans.target(),target_encl));
}

void
HybridReachabilityAnalyser::
_splitTargetEnclosures(bool& isValid,
					   std::list<EnclosureType>& initial_enclosures,
					   const DiscreteState& target_loc,
					   const ContinuousEnclosureType& target_encl,
					   const Vector<Float>& minTargetCellWidths,
					   const Box& target_bounding) const
{
	// Get the target enclosure box
	const Box& target_encl_box = target_encl.bounding_box();

	// If the cell box lies outside the bounding domain (i.e. not inside and disjoint)
	// of the target location, ignore it and notify
	if (!target_encl_box.inside(target_bounding) &&
		target_encl_box.disjoint(target_bounding)) {
		ARIADNE_LOG(7,"Discarding target enclosure " << target_encl_box <<
					  " being outside the domain in location " << target_loc.name() << ".\n");
		isValid = false;
		return;
	}

	// Otherwise try splitting the enclosure
	bool hasSplit = false;
	for (uint i=0; i < minTargetCellWidths.size(); i++) {
		// If the enclosure has width larger than that of the minimum cell, split on that dimension
		if (target_encl_box[i].width() > minTargetCellWidths[i]) {
			hasSplit = true;
			std::pair<ContinuousEnclosureType,ContinuousEnclosureType> split_sets = target_encl.split(i);
			// Call recursively on the two enclosures
			_splitTargetEnclosures(isValid,initial_enclosures,target_loc,split_sets.first,minTargetCellWidths,target_bounding);
			_splitTargetEnclosures(isValid,initial_enclosures,target_loc,split_sets.second,minTargetCellWidths,target_bounding);
		}
	}

	// If we could not split, we put the hybrid enclosure into the new initial_enclosures
	if (!hasSplit)
		initial_enclosures.push_back(EnclosureType(target_loc,target_encl));
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
	// Get the actual concurrency
	const uint concurrency = boost::thread::hardware_concurrency() - free_cores;

	// Check that the concurrency is greater than zero
	ARIADNE_ASSERT_MSG(concurrency>0,"Error: concurrency must be at least 1.");

	// Create the worker
	UpperReachEvolveContinuousWorker worker(_discretiser,sys,initial_enclosures,time,accuracy,concurrency);

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

    // The reached region
    GTS reach(system.grid());

    // Split the initial set into enclosures that are smaller than the minimum cell, given the grid
    list<EnclosureType> initial_enclosures = _split_initial_set(initial_set,system.grid());

    ARIADNE_LOG(6,"Evolving and discretising...\n");

    // For each initial enclosure, perform evolution and adjoin to the reach and evolve regions
    for (list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); encl_it++) {
        reach.adjoin(_discretiser->reach(system,*encl_it,time,_parameters->maximum_grid_depth,LOWER_SEMANTICS));
    }

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

	// The reached and evolved regions
    GTS reach; GTS evolve;

    // Split the initial set into enclosures that are smaller than the minimum cell, given the grid
    list<EnclosureType> initial_enclosures = _split_initial_set(initial_set,system.grid());

	ARIADNE_LOG(5,"Computing evolution...\n");
    // For each initial enclosure, perform evolution and adjoin to the reach and evolve regions
    for (list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); encl_it++) {
        GTS cell_reach,cell_final;
        make_lpair(cell_reach,cell_final)=_discretiser->evolution(system,*encl_it,time,_parameters->maximum_grid_depth,LOWER_SEMANTICS);
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
    ARIADNE_LOG(4,"HybridReachabilityAnalyser::upper_evolve(system,set,time)\n");
 
	GTS reach, evolve;
	make_lpair(reach,evolve) = upper_reach_evolve(system, initial_set, time); // Runs the upper_reach_evolve routine on its behalf

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

        evolve.remove(intermediate); // Remove only the cells that are intermediate
        evolve.restrict(bounding);
        reach.adjoin(found);

        ARIADNE_LOG(6,"found "<<found.size()<<" cells, of which "<<evolve.size()<<" are new.\n");
    }
    reach.recombine();

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


std::pair<HybridReachabilityAnalyser::SetApproximationType,bool>
HybridReachabilityAnalyser::
upper_chain_reach(const SystemType& system,
            const HybridImageSet& initial_set) const
{
	/* General procedure:
		1) Build the working sets from the initial enclosures
		2) While a working set exists, evolve the current working set for a step, where the following exceptions are handled:
			a) If the evolution time limit is reached, adjoin both the reached set and the final set
			b) If the enclosure of the initial set is larger than the maximum enclosure allowed, skip the evolution and adjoin the final set
			c) If a blocking event is definitely initially active, or definitely finally active due to flow, adjoin the reached set only and discard the final set
			d) If the bounding box of the final set is not a subset of the domain box, adjoin the reached set only and discard the final set (would be plausible to discard the reached set also)
		3) Discretise the reached and final sets into cells
		4) Remove the previous intermediate cells from the final cells and remove the previous reached cells from the reached cells
		5) Mince the resulting final cells and put their enclosures into the new initial enclosures
		6) Check activations of transitions for each cell (NOTE: the maximum allowed depth is related to the cells of the target location):
			a) If a transition is definitely active for a cell: mince the cell, then for each resulting cell
				i) enclose it
				ii) apply the transition
				iii) put the resulting enclosure into the new initial enclosures
			b) If a transition is possibly active at the maximum allowed depth:
				i) enclose it
				ii) apply the transition
				iii) put the resulting enclosure into the new initial enclosures
			c) If a transition is possibly active at a depth lesser than the maximum allowed: split the cell and repeat a) for each of the two cells
		7) Adjoin the reached cells into the previous reached cells, adjoin the final cells into the intermediate cells, then recombine both the resulting cells sets
		8) If new initial enclosures exist, restart from 1), otherwise terminate.
	*/

	// Initialize the validity flag (will be invalidated if any restriction of the reached region is performed)
	bool isValid = true;

	// Helper class for operations on Taylor sets
	TaylorCalculus tc;

	// Assign local variables
    Float lock_to_grid_time = _parameters->lock_to_grid_time;
    int lock_to_grid_steps = _parameters->lock_to_grid_steps;
    int maximum_grid_depth = _parameters->maximum_grid_depth;

    // The number of divisions to a cell
    long numCellDivisions = (1<<maximum_grid_depth);

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

        ARIADNE_LOG(6,"Initial enclosures size = " << initial_enclosures.size() << "\n");

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

		/* Check the reached cells for transitions: we should check the largest enclosures possible, but also provide target
		 * enclosures that are smaller than the smallest cell of the target location.
		 * (NOTE ON SPLITTING: it is consequently done in respect to the smallest cell of the target location:
		 *  we allow splitting as long as there is any dimension having width larger than the width of the smallest cell)
		 * Procedure:
		 * a) Recombine the reach cells, in order to start from the largest enclosures;
		 * b) For each cell:
		 * 	  i) Create the hybrid enclosure from it and push-front its continuous component into a split_list (initially empty);
		 * 	  ii) For each transition of the involved location:
		 * 		   *) While the split_list is not empty;
		 * 				1) Pop-back an enclosure from the split_list;
		 * 		   		2) If the transition is not active, drop the enclosure;
		 * 		   		3) If the transition is definitely active,
		 * 					x) split the enclosures as finely as possible;
		 * 					xx) on each enclosure of minimum size, apply the reset and create the target enclosure;
		 * 		   		4) If the transition is possibly active, try splitting the enclosure:
		 * 			  		x) if a splitting is possible, push-front the split enclosures into the split_list;
		 * 		   			xx) if no splitting is possible, apply the reset and create the target enclosure.
		 */

		// Recombine to have the largest enclosures available at the original location
		new_reach.recombine();
		ARIADNE_LOG(6,"Reach size after recombining = "<<new_reach.size()<<"\n");
		for (GTS::const_iterator cell_it = new_reach.begin(); cell_it != new_reach.end(); cell_it++)
		{
			// Get the location
			const DiscreteState& loc = cell_it->first;
			// Get the box
			const Box& bx = cell_it->second.box();

			// If the cell box lies outside the bounding domain (i.e. not inside and disjoint)
			// of the target location, ignore it and notify
			if (!bx.inside(_parameters->bounding_domain[loc]) &&
				bx.disjoint(_parameters->bounding_domain[loc])) {
				ARIADNE_LOG(7,"Discarding reach enclosure " << bx <<
							  " being outside the domain in location " << loc.name() << ".\n");
				isValid = false;

				// Return false if we skip when the system is unprovable
				if (_parameters->skip_if_unprovable)
					return make_pair<GTS,bool>(HybridGridTreeSet(),false);
			}

			// Initialize the split list
			std::list<ContinuousEnclosureType> split_list;
			// Push the continuous enclosure of the cell into the list
			split_list.push_front(ContinuousEnclosureType(bx));

			// Get the transitions for the corresponding location
			list<DiscreteTransition> transitions = system.transitions(loc);

			ARIADNE_LOG(7,"Checking transitions for box "<< bx <<" in location " << loc.name() << "\n");

			// For each transition
			for (list<DiscreteTransition>::const_iterator trans_it = transitions.begin(); trans_it != transitions.end(); trans_it++)
			{
				ARIADNE_LOG(8,"Transition = "<<*trans_it<<"\n");

				// Get the minimum cell widths in the target location
				Vector<Float> minTargetCellWidths = grid[trans_it->target()].lengths()/numCellDivisions;
				// Get the target bounding domain
				const Box& target_bounding = _parameters->bounding_domain[trans_it->target()];

				// While there is an enclosure into the split list
				while (!split_list.empty()) {

					// Get an enclosure and remove it from the list
					ContinuousEnclosureType encl = split_list.back();
					split_list.pop_back();

					// Get the time model of the guard
					TaylorModel guard_time_model = apply(trans_it->activation(),encl)[0];

					ARIADNE_LOG(9,"Guard time model = "<<guard_time_model<<" (range: " << guard_time_model.range() << ")\n");

					// If definitely active, add the target enclosure on the splittings of the enclosure
					// If possibly active, split the enclosure if possible, otherwise create the target enclosure
					// (if otherwise the transition is not active, do nothing)

					if (guard_time_model.range().lower()>0) {

						time_t t1,t2;
						time(&t1);

						// Populate the initial enclosures with the hybrid enclosures derived from the splittings of the original enclosure
						//_splitAndCreateTargetEnclosures(isValid,initial_enclosures,encl,minTargetCellWidths,target_bounding,*trans_it,tc);

						// Get the target continuous enclosure
						ContinuousEnclosureType target_encl = tc.reset_step(trans_it->reset(),encl);
						// Populate the initial enclosures with the hybrid enclosures derived from the splittings of the target enclosure
						_splitTargetEnclosures(isValid,initial_enclosures,trans_it->target(),target_encl,minTargetCellWidths,target_bounding);

						// Return false if we skip when the system is unprovable
						if (!isValid && _parameters->skip_if_unprovable)
							return make_pair<GTS,bool>(HybridGridTreeSet(),false);

						time(&t2);
						ARIADNE_LOG(7,"Time required for splitting and creating: " << difftime(t2,t1) << ".\n");
						// Set the in-use value for the target location as true
						currentLocationsInUse[trans_it->target()] = true;

					} else if (!(guard_time_model.range().upper()<0)) {

						// Try splitting the enclosure
						bool hasSplit = false;
						for (uint i=0; i < minTargetCellWidths.size(); i++) {
							// If the enclosure has width larger than that of the minimum cell, split on that dimension
							if (encl.bounding_box()[i].width() > minTargetCellWidths[i]) {
								hasSplit = true;
								std::pair<ContinuousEnclosureType,ContinuousEnclosureType> split_sets = encl.split(i);
								split_list.push_front(split_sets.first);
								split_list.push_front(split_sets.second);
								break;
							}
						}
						// If we could not split
						if (!hasSplit) {

							// Get the target continuous enclosure
							ContinuousEnclosureType target_encl = tc.reset_step(trans_it->reset(),encl);
							// Get its bounding box
							const Box& target_encl_box = target_encl.bounding_box();

							// If the cell box lies outside the bounding domain (i.e. not inside and disjoint)
							// of the target location, ignore it and notify
							if (!target_encl_box.inside(target_bounding) &&
								target_encl_box.disjoint(target_bounding)) {
								ARIADNE_LOG(7,"Discarding target enclosure " << target_encl_box <<
											  " being outside the domain in location " << trans_it->target().name() << ".\n");
								isValid = false;

								// Return false if we skip when the system is unprovable
								if (_parameters->skip_if_unprovable)
									return make_pair<GTS,bool>(HybridGridTreeSet(),false);
							} else {
								// We apply the reset and put the hybrid enclosure into the new initial_enclosures
								initial_enclosures.push_back(EnclosureType(trans_it->target(),target_encl));
								// Set the in-use value for the target location as true
								currentLocationsInUse[trans_it->target()] = true;
							}
						}
					}
				}
			}
		}

		/* Mince the final cells, then for each of them add the enclosure into the new initial enclosures
		 * if their bounding box is not outside (i.e. inside, or not inside but not disjoint) the bounding domain */
		new_final.mince(maximum_grid_depth);
		ARIADNE_LOG(6,"Final size after mincing = "<<new_final.size()<<"\n");
		for (GTS::const_iterator cell_it = new_final.begin(); cell_it != new_final.end(); cell_it++) {
			// Get the location
			const DiscreteState& loc = cell_it->first;
			// If inside the domain
			if (cell_it->second.box().inside(_parameters->bounding_domain[loc])) {
				initial_enclosures.push_back(_discretiser->enclosure(*cell_it));
			} else if (!cell_it->second.box().disjoint(_parameters->bounding_domain[loc])) {
				initial_enclosures.push_back(_discretiser->enclosure(*cell_it));
			} else {
				ARIADNE_LOG(7,"Discarding enclosure " << cell_it->second.box() <<
							" from final cell outside the domain in location " << loc.name() <<".\n");
				isValid = false;

				// Return false if we skip when the system is unprovable
				if (_parameters->skip_if_unprovable)
					return make_pair<GTS,bool>(HybridGridTreeSet(),false);
			}

		}

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

    return make_pair<GTS,bool>(reach,isValid);
}

std::pair<HybridReachabilityAnalyser::SetApproximationType,bool>
HybridReachabilityAnalyser::
lower_chain_reach(const SystemType& system,
				  const HybridImageSet& initial_set) const
{
	typedef std::list<EnclosureType> EL;
	typedef std::map<DiscreteState,uint> HUM;

	// Get the actual concurrency
	const uint concurrency = boost::thread::hardware_concurrency() - free_cores;
	// Check that the concurrency is greater than zero
	ARIADNE_ASSERT_MSG(concurrency>0,"Error: concurrency must be at least 1.");

	// The global reach region
    GTS reach(system.grid());

	// Get the lock time, in the case pruning is to be performed
	TimeType lock_time(_parameters->lock_to_grid_time,_parameters->lock_to_grid_steps);

	// Get the widened safe hybrid box
	HybridBoxes disprove_bounds;
	for (HybridBoxes::const_iterator loc_it = _parameters->safe_region.begin(); loc_it != _parameters->safe_region.end(); loc_it++) {
		// Copy the related box and widen it
		Box bx = loc_it->second;
		bx.widen();
		// Insert the couple with widened box
		disprove_bounds.insert(make_pair<DiscreteState,Box>(loc_it->first,bx));
	}

    // The disproving result
    bool isDisproved = false;

    // Split the initial set into enclosures that are smaller than the minimum cell, given the grid
    EL initial_enclosures = _split_initial_set(initial_set,system.grid());

    ARIADNE_LOG(5,"Computing recurrent evolution...\n");

	// Perform an infinite loop (exception taken if no initial enclosures exist)
	for (uint i=0;!initial_enclosures.empty();i++)
	{
		ARIADNE_LOG(5,"Iteration " << i << "\n");

		// The final enclosures
		EL final_enclosures;

		// The evolve sizes resulting from evolution
		std::pair<HUM,HUM> evolve_sizes;

		// The sizes of the adjoined (the cells) or superposed (the enclosures) evolve
		HUM& adjoined_evolve_sizes = evolve_sizes.first;
		HUM& superposed_evolve_sizes = evolve_sizes.second;

		// The evolve
		GTS evolve;
		// The new reach
		GTS new_reach;
		// The new disproving result
		bool new_isDisproved = false;

		ARIADNE_LOG(6,"Initial enclosures size = " << initial_enclosures.size() << "\n");

		// Create the worker
		LowerChainReachWorker worker(_discretiser,initial_enclosures,system,lock_time,
									 _parameters->maximum_grid_depth,concurrency, disprove_bounds, _parameters->skip_if_disproved);

		ARIADNE_LOG(6,"Evolving and discretising...\n");

		// Compute and get the result
		make_ltuple<std::pair<HUM,HUM>,EL,GTS,bool>(evolve_sizes,final_enclosures,
																new_reach,new_isDisproved) = worker.get_result();

		// Update the disproving result flag
		isDisproved = isDisproved || new_isDisproved;

		ARIADNE_LOG(6,"Reach size before removal = " << new_reach.size() << "\n");

		// Remove from the partial reached region the global reached region
		new_reach.remove(reach);

		ARIADNE_LOG(6,"Reach size after removal  = " << new_reach.size() << "\n");

		// Check for inclusion: if no new cells are found, terminate
		if (new_reach.empty())
			break;

		// Otherwise proceed and adjoin the new reach
		reach.adjoin(new_reach);

		ARIADNE_LOG(6,"Final enclosures size = " << final_enclosures.size() << "\n");

		// Pruning of the dump of the final enclosures into the initial enclosures
		while (!final_enclosures.empty())
		{
			// Pop the current enclosure
			EnclosureType encl = final_enclosures.front();
			final_enclosures.pop_front();

			// Get the location
			const DiscreteState& loc = encl.location();

			/* If the enclosure lies outside the bounding domain (i.e. not inside and disjoint), then the
			 * domain is not proper and an error should be thrown. */
			if (!encl.continuous_state_set().bounding_box().inside(_parameters->bounding_domain[loc]) &&
				encl.continuous_state_set().bounding_box().disjoint(_parameters->bounding_domain[loc])) {
				ARIADNE_FAIL_MSG("Found an enclosure in location " << loc.name() << " with bounding box " << encl.continuous_state_set().bounding_box() <<
								 " lying outside the domain in lower semantics: the domain is incorrect.\n");
			}

			/* If pruning is to be performed, push just a fraction of the final_enclosures into the initial_enclosures;
			 * otherwise, push indiscriminately.
			 */
			if (_parameters->enable_lower_pruning)
			{
				// Get the ratio between the adjoined evolve size and the superposed evolve size
				Float ratio = (Float)adjoined_evolve_sizes[loc]/(Float)superposed_evolve_sizes[loc];

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

	return make_pair<GTS,bool>(reach,isDisproved);
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
_prove(const SystemType& system,
	   const HybridImageSet& initial_set)
{
	ARIADNE_LOG(4,"Proving...\n");

	// The reachable set
	HybridGridTreeSet reach;
	// The flag that informs whether the region is valid for proving
	bool isValid;

	make_lpair<HybridGridTreeSet,bool>(reach,isValid) = upper_chain_reach(system,initial_set); // Perform the chain reachability analysis

	// If the reached region is valid (i.e. it has not been restricted), perform checking, otherwise return false
	if (isValid) {
		// If the reached region is definitely inside the hybrid safe region, the safety property is proved
		bool result = definitely(reach.subset(_parameters->safe_region));
		ARIADNE_LOG(4, (result ? "Proved.\n" : "Not proved.\n") );
		return result;
	} else {
		ARIADNE_LOG(4,"Not checked due to restriction inside the domain.\n");
		return false;
	}
}


bool 
HybridReachabilityAnalyser::
_disprove(const SystemType& system,
		  const HybridImageSet& initial_set)
{
	ARIADNE_LOG(4,"Disproving...\n");

	// The reach
	HybridGridTreeSet reach;
	// The flag which notifies if the system has been disproved
	bool isDisproved;

	// Get the result
	make_lpair<HybridGridTreeSet,bool>(reach,isDisproved) = lower_chain_reach(system,initial_set);

	// Notify
	ARIADNE_LOG(4, (isDisproved ? "Disproved.\n" : "Not disproved.\n") );

	return isDisproved;
}


tribool
HybridReachabilityAnalyser::
verify(const SystemType& system,
       const HybridImageSet& initial_set)
{
		ARIADNE_LOG(3, "Verification...\n");

		// Perform the proving
		if (_prove(system,initial_set))
		{
			ARIADNE_LOG(3, "Safe.\n");
			return true;
		}

		// Perform the disproving
		if (_disprove(system,initial_set))
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
_getLooselyConstrainedHybridGrid(const HybridFloatVector& hmad) const
{
	// Get the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = hmad.begin()->second.size();

	// Initialize the hybrid grid
	HybridGrid hg;	

	// Initialize the hybrid grid lengths
	std::map<DiscreteState,Vector<Float> > hybridgridlengths;

	// Keep the minimum non-zero length for each variable
	Vector<Float> minNonZeroLengths(css);
	// Initialize the values to infinity
	for (uint i=0;i<css;i++)
		minNonZeroLengths[i] = std::numeric_limits<double>::infinity();

	// Get the grid lengths from the derivatives
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
	{
		// Keep the location
		const DiscreteState& loc = hfv_it->first;
		// Keep the domain box
		const Box& domain_box = _parameters->bounding_domain.find(loc)->second;

		// Initialize the lengths
		Vector<Float> gridlengths(css);

		// Initialize the maximum derivative/domain_width ratio among variables
		Float maxratio = 0.0;

		// For each dimension of the continuous space
		for (uint i=0;i<css;i++)
			maxratio = max(maxratio,hfv_it->second[i]/domain_box[i].width()); // Get the largest ratio

		// Assign the lengths and check the minimum length to be assigned
		for (uint i=0;i<css;i++)
		{
			// Update the minimum in respect to the domain width
			minNonZeroLengths[i] = min(minNonZeroLengths[i],domain_box[i].width());
			// If the derivative is greater than zero, get the grid length and update the minimum in respect to the grid length
			if (hfv_it->second[i] > 0) {
				gridlengths[i] = hfv_it->second[i]/maxratio;
				minNonZeroLengths[i] = min(minNonZeroLengths[i],hfv_it->second[i]/maxratio);
			}
		}

		// Add the pair to the hybrid lengths
		hybridgridlengths.insert(make_pair<DiscreteState,Vector<Float> >(loc,gridlengths));
	}

	// Now that the minimum lengths have been evaluated, we can assign the lengths in those cases where the variables
	// have zero derivative, then complete the construction of the hybrid grid
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) {

		// Keep the location
		const DiscreteState& loc = hfv_it->first;

		/* Choose the center for the grid: the center of the reached region, if a non-empty reached region exists,
		 * otherwise the center of the domain box */
		Vector<Float> centre(css);
		if (_statistics->upper().reach.has_location(loc) && !_statistics->upper().reach[loc].empty()) {
			centre = _statistics->upper().reach[loc].bounding_box().centre();
		} else {
			centre = _parameters->bounding_domain.find(loc)->second.centre();
		}

		// For each variable having zero length, assign the minimum value
		for (uint i=0;i<css;i++)
			if (hybridgridlengths[loc][i] == 0)
				hybridgridlengths[loc][i] = minNonZeroLengths[i];

		// Add the related grid
		hg[loc] = Grid(centre,hybridgridlengths[loc]);
	}

	// Return
	return hg;
}


HybridGrid 
HybridReachabilityAnalyser::
_getHybridGrid(const HybridFloatVector& hmad) const
{
	// Get the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = hmad.begin()->second.size();

	// Initialize the hybrid grid
	HybridGrid hg;

	// Initialize the hybrid grid lengths
	std::map<DiscreteState,Vector<Float> > hybridgridlengths;

	// Get the minimum domain length for each variable
	Vector<Float> minDomainLengths(css);
	// Initialize the values to infinity
	for (uint i=0;i<css;i++) {
		minDomainLengths[i] = std::numeric_limits<double>::infinity();
	}
	// Get the actual values
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) {
		for (uint i=0;i<css;i++) {
			minDomainLengths[i] = min(minDomainLengths[i], _parameters->bounding_domain.find(hfv_it->first)->second[i].width());
		}
	}

	// Initialize the minimum non-zero length for each variable as the minimum domain lengths
	Vector<Float> minNonZeroLengths = minDomainLengths;

	// Get the maximum derivative/domainlength ratio among variables and locations
	Float maxratio = 0.0;
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) {
		for (uint i=0;i<css;i++) {
			maxratio = max(maxratio,hfv_it->second[i]/minDomainLengths[i]);
		}
	}

	// Update the minimum non zero lengths
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) {
		for (uint i=0;i<css;i++) {
			if (hfv_it->second[i] > 0) {
				minNonZeroLengths[i] = min(minNonZeroLengths[i],hfv_it->second[i]/maxratio);
			}
		}
	}

	// Get the grid lengths from the derivatives
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
	{
		// Initialize the lengths
		Vector<Float> gridlengths(css);

		// Assign the lengths
		for (uint i=0;i<css;i++) {
			gridlengths[i] = (hfv_it->second[i] > 0) ? hfv_it->second[i]/maxratio : minNonZeroLengths[i];
		}

		// Add the pair to the hybrid lengths
		hybridgridlengths.insert(make_pair<DiscreteState,Vector<Float> >(hfv_it->first,gridlengths));
	}

	// We complete the construction of the hybrid grid by defining the grid center
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) {

		// Keep the location
		const DiscreteState& loc = hfv_it->first;

		/* Choose the center for the grid: the center of the reached region, if a non-empty reached region exists,
		 * otherwise the center of the domain box */
		Vector<Float> centre(css);
		if (_statistics->upper().reach.has_location(loc) && !_statistics->upper().reach[loc].empty()) {
			centre = _statistics->upper().reach[loc].bounding_box().centre();
		} else {
			centre = _parameters->bounding_domain.find(loc)->second.centre();
		}

		// Add the related grid
		hg[loc] = Grid(centre,hybridgridlengths[loc]);
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
	Vector<Float> minlengths(css,std::numeric_limits<double>::infinity());
	// Initialize the gridlengths
	Vector<Float> gridlengths(css);		
	// Initialize the maximum absolute derivatives
	Vector<Float> mad(css);

	// Get the maximum absolute derivatives
	for (uint i=0;i<css;i++)
		for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
			mad[i] = max(mad[i],hfv_it->second[i]);

	// For each dimension
	for (uint i=0;i<css;i++)
		// For each couple DiscreteState,Vector<Float>
		for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
		{
			maxratio = max(maxratio,mad[i]/_parameters->bounding_domain.find(hfv_it->first)->second[i].width()); // Check for maximum ratio
			minlengths[i] = min(minlengths[i],_parameters->bounding_domain.find(hfv_it->first)->second[i].width()); // Get the minimum domain length
		}

	// Assign the lengths
	for (uint i=0;i<css;i++)
		gridlengths[i] = (mad[i] > 0) ? mad[i]/maxratio : minlengths[i];

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
_setInitialParameters(SystemType& system, const HybridBoxes& domain, const HybridBoxes& safe_region)
{
	// Set the domain
	_parameters->bounding_domain = domain;
	// Set the safe region
	_parameters->safe_region = safe_region;
}


void 
HybridReachabilityAnalyser::
_tuneParameters(SystemType& system)
{
	/* Tune the parameters:
	 * a) evaluate the derivatives from the domain or the latest upper reached region;
	 * b) get the grid with lengths proportional to the derivatives,
	 * 	  and scale such that the grid cell can be included into the domain;
	 * c) get the maximum enclosure cell as a multiple of the minimum cell;
	 * d) get the maximum step size as the step size which makes a variable evolve a distance equal
	 * 	  to twice the minimum cell length (for each dimension),
	 *    under the assumption of moving at speed equal to the maximum derivative.
	 */

	// Evaluate the maximum absolute derivatives
	HybridFloatVector hmad = _getHybridMaximumAbsoluteDerivatives(system);
	ARIADNE_LOG(3, "Maximum absolute derivatives: " << hmad << "\n");
	// Grid
	system.set_grid(_getHybridGrid(hmad));
	ARIADNE_LOG(3, "Grid: " << system.grid() << "\n");
	// Maximum enclosure cell
	_setMaximumEnclosureCell(system.grid());
	ARIADNE_LOG(3, "Maximum enclosure cell: " << _discretiser->parameters().maximum_enclosure_cell << "\n");
	// Maximum step size
	_setHybridMaximumStepSize(hmad,system.grid());
	ARIADNE_LOG(3, "Maximum step size: " << _discretiser->parameters().hybrid_maximum_step_size << "\n");
}


tribool 
HybridReachabilityAnalyser::
verify_iterative(SystemType& system, 
				 const HybridImageSet& initial_set, 
				 const HybridBoxes& safe,
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
	_setInitialParameters(system, domain, safe);

    for(_parameters->maximum_grid_depth = _parameters->lowest_maximum_grid_depth;
    	_parameters->maximum_grid_depth <= _parameters->highest_maximum_grid_depth;
    	_parameters->maximum_grid_depth++)
	{ 
    	sprintf(mgd_char,"%i", _parameters->maximum_grid_depth);
		ARIADNE_LOG(2, "DEPTH " << _parameters->maximum_grid_depth << "\n");

		// Tune the parameters for the current iteration
		_tuneParameters(system);

		// Perform the verification
		tribool result = verify(system,initial_set);

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
				  const HybridBoxes& safe,
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
	tribool lower_result = verify_iterative(system,initial_set,safe,domain);

	if (definitely(lower_result)) { ARIADNE_LOG(1,"Safe.\n"); }
	else { ARIADNE_LOG(1,"Not safe.\n"); }

	// Check the upper bound
	ARIADNE_LOG(1,"Checking upper interval bound... ");

	// Set the parameter
	param.set_value(parameter_interval.upper());
	// Substitute the value
	system.substitute(param);
	// Perform the verification
	tribool upper_result = verify_iterative(system,initial_set,safe,domain);

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
	while (param_int.width() > tolerance*parameter_interval.width())
	{
		// Set the parameter as the midpoint of the interval
		param.set_value(param_int.midpoint());
		// Substitute the value
		system.substitute(param);

		ARIADNE_LOG(1,"Checking " << param_int << " (midpoint: " << param_int.midpoint() << ", width: " << param_int.width() << ") ... ");

		// Perform the verification
		tribool result = verify_iterative(system,initial_set,safe,domain);

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
				param_int.set_lower(param.value()); 
				 }
		}
	}

	// Returns the system to its original parameter value
	system.substitute(parameter);

	if (updateFromBottom)
		return Interval(parameter_interval.lower(),param_int.lower());
	else
		return Interval(param_int.upper(),parameter_interval.upper());
}

std::pair<Interval,Interval>
HybridReachabilityAnalyser::
safety_unsafety_parametric(SystemType& system, 
						   const HybridImageSet& initial_set, 
						   const HybridBoxes& safe,
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
	tribool lower_result = verify_iterative(system,initial_set,safe,domain);

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
	tribool upper_result = verify_iterative(system,initial_set,safe,domain);

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
			result = verify_iterative(system,initial_set,safe,domain);

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
					safety_int.set_lower(param.value()); }
			}

			// If both intervals are not empty
			if (!safety_int.empty() && !unsafety_int.empty())
			{
				// Breaks if the minimum distance between the safe and unsafe values is lesser than the tolerance
				// (which is a more relaxed condition than the ones on the interval width)
				if ((safeOnBottom && (unsafety_int.upper() - safety_int.lower() <= tolerance*parameter_interval.width())) ||
					(!safeOnBottom && (safety_int.upper() - unsafety_int.lower() <= tolerance*parameter_interval.width())))
						break;
			}
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
			result = verify_iterative(system,initial_set,safe,domain);

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
					unsafety_int.set_upper(param.value()); }
			}

			// If both intervals are not empty
			if (!safety_int.empty() && !unsafety_int.empty())
			{
				// Breaks if the minimum distance between the safe and unsafe values is lesser than the tolerance
				// (which is a more relaxed condition than the ones on the interval width)
				if ((safeOnBottom && (unsafety_int.upper() - safety_int.lower() <= tolerance*parameter_interval.width())) ||
					(!safeOnBottom && (safety_int.upper() - unsafety_int.lower() <= tolerance*parameter_interval.width())))
						break;
			}
		}

		// Breaks if the safety interval is not empty and the safety interval width is lesser than the tolerance
		if (!safety_int.empty() && safety_int.width() <= tolerance*parameter_interval.width())
			break;
		// Breaks if the unsafety interval is not empty and the unsafety interval width is lesser than the tolerance
		if (!unsafety_int.empty() && unsafety_int.width() <= tolerance*parameter_interval.width())
			break;
	}

	// Returns the system to its original parameter value
	system.substitute(parameter);

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

void HybridReachabilityAnalyser::parametric_2d(SystemType& system,
											   const HybridImageSet& initial_set,
											   const HybridBoxes& safe,
											   const HybridBoxes& domain,
											   const RealConstant& xParam,
											   const RealConstant& yParam,
											   const Interval& xBounds,
											   const Interval& yBounds,
											   const unsigned& numPointsPerAxis,
											   const Float& tolerance)
{
	// FIXME: in the end, the system will be modified with the upper values of the
	// XY bounds. We must be able to recover (or pass) the initial values in order to restore them.

	// Generates the file name
	std::string filename = system.name();
	filename = filename + "-" + xParam.name() + "-" + yParam.name();

	// Copies the parameters into constants
	RealConstant xConstant = xParam;
	RealConstant yConstant = yParam;

	// Initializes the results
	Parametric2DAnalysisResults results(filename.c_str(),xBounds,yBounds,numPointsPerAxis);

	ARIADNE_LOG(1,"\nSweeping on " << xParam.name() << " in " << xBounds << ":\n");

	// Sweeps in the X direction
	for (unsigned i=0; i<numPointsPerAxis; i++)
	{
		// Changes the value of X
		xConstant.set_value(xBounds.lower()+i*xBounds.width()/(numPointsPerAxis-1));
		// Modifies the system with the new X
		system.substitute(xConstant);

		ARIADNE_LOG(1,"\nAnalysis with " << xParam.name() << " = " << xConstant.value().lower() << "...\n");

		// Performs the analysis on Y and adds to the results of X
		std::pair<Interval,Interval> result = safety_unsafety_parametric(system,initial_set,safe,domain,yParam,yBounds,tolerance);
		results.insertXValue(result);

		ARIADNE_LOG(1,"Obtained safety in " << result.first << " and unsafety in " << result.second << ".\n");
	}

	ARIADNE_LOG(1,"\nSweeping on " << yParam.name() << " in " << yBounds << ":\n");

	// Sweeps in the Y direction
	for (unsigned i=0; i<numPointsPerAxis; i++)
	{
		// Changes the value of Y
		yConstant.set_value(yBounds.lower()+i*yBounds.width()/(numPointsPerAxis-1));
		// Modifies the system with the new Y
		system.substitute(yConstant);

		ARIADNE_LOG(1,"\nAnalysis with " << yParam.name() << " = " << yConstant.value().lower() << "...\n");

		// Performs the analysis on X and adds to the results of Y
		std::pair<Interval,Interval> result = safety_unsafety_parametric(system,initial_set,safe,domain,xParam,xBounds,tolerance);
		results.insertYValue(result);

		ARIADNE_LOG(1,"Obtained safety in " << result.first << " and unsafety in " << result.second << ".\n");
	}

	// Modifies the system to return to its original parameter values
	system.substitute(xParam);
	system.substitute(yParam);

	// Dumps the results into a file
	results.dump();
	// Draws the result
	results.draw();
}

} // namespace Ariadne
