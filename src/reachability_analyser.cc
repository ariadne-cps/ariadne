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

SystemVerificationInfo::SystemVerificationInfo(HybridAutomaton& system,
						   HybridImageSet& initial_set,
						   HybridBoxes& domain,
						   HybridBoxes& safe_region) :
							_system(system), _initial_set(initial_set),
							_domain(domain), _safe_region(safe_region),
							_projection(std::vector<uint>(0))
{
	_check_fields();
}

SystemVerificationInfo::SystemVerificationInfo(HybridAutomaton& system,
						   HybridImageSet& initial_set,
						   HybridBoxes& domain,
						   std::vector<uint>& projection) :
							_system(system), _initial_set(initial_set),
							 _domain(domain), _safe_region(unbounded_hybrid_box(system.state_space())),
							 _projection(projection)
{
	_check_fields();
}

SystemVerificationInfo::SystemVerificationInfo(HybridAutomaton& system,
						   HybridImageSet& initial_set,
						   HybridBoxes& domain,
						   HybridBoxes& safe_region,
						   std::vector<uint>& projection) :
							_system(system), _initial_set(initial_set),
							_domain(domain), _safe_region(safe_region),
							_projection(projection)
{
	_check_fields();
}

void SystemVerificationInfo::_check_fields() const
{
	HybridSpace hspace = _system.state_space();
	for (HybridImageSet::const_iterator it = _initial_set.begin(); it != _initial_set.end(); ++it) {
		HybridSpace::const_iterator hspace_it = hspace.find(it->first);
		ARIADNE_ASSERT_MSG(hspace_it != hspace.end(),
						   "The location " << it->first.name() << "is not present into the system.");
		ARIADNE_ASSERT_MSG(hspace_it->second == it->second.dimension(),
						   "The dimension of the continuous space in the initial set for location " << it->first.name() << " does not match the system space");
	}
	for (HybridSpace::const_iterator hspace_it = hspace.begin(); hspace_it != hspace.end(); ++hspace_it) {
		HybridBoxes::const_iterator domain_it = _domain.find(hspace_it->first);
		ARIADNE_ASSERT_MSG(domain_it != _domain.end(),
						   "The location " << hspace_it->first.name() << "is not present into the domain.");
		ARIADNE_ASSERT_MSG(hspace_it->second == domain_it->second.dimension(),
						   "The dimension of the continuous space in the domain for location " << hspace_it->first.name() << " does not match the system space");
		HybridBoxes::const_iterator safe_it = _safe_region.find(hspace_it->first);
		ARIADNE_ASSERT_MSG(safe_it != _safe_region.end(),
						   "The location " << hspace_it->first.name() << "is not present into the safe region.");
		ARIADNE_ASSERT_MSG(hspace_it->second == safe_it->second.dimension(),
						   "The dimension of the continuous space in the safe region for location " << hspace_it->first.name() << " does not match the system space");
	}
}


std::ostream&
SystemVerificationInfo::write(std::ostream& os) const
{
	os << "(System: " << _system << "; Initial set: " << _initial_set << "; Domain: " <<
			_domain << "; Safe region: " << _safe_region << "; Projection: " << _projection << ")";
	return os;
}

HybridReachabilityAnalyser::
~HybridReachabilityAnalyser()
{
}


HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const HybridDiscretiser<HybridEvolver::ContinuousEnclosureType>& discretiser)
    : _parameters(new EvolutionParametersType())
    , _discretiser(discretiser.clone())
	, plot_verify_results(false)
{
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
splitTargetEnclosures(bool& isValid,
					   std::list<EnclosureType>& initial_enclosures,
					   const DiscreteState& target_loc,
					   const ContinuousEnclosureType& target_encl,
					   const Vector<Float>& minTargetCellWidths,
					   const Box& target_bounding)
{
	// Get the target enclosure box
	const Box& target_encl_box = target_encl.bounding_box();

	// If the cell box lies outside the bounding domain (i.e. not inside and disjoint)
	// of the target location, ignore it and notify
	if (!target_encl_box.inside(target_bounding) &&
		target_encl_box.disjoint(target_bounding)) {
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
			splitTargetEnclosures(isValid,initial_enclosures,target_loc,split_sets.first,minTargetCellWidths,target_bounding);
			splitTargetEnclosures(isValid,initial_enclosures,target_loc,split_sets.second,minTargetCellWidths,target_bounding);
		}
	}

	// If we could not split, we put the hybrid enclosure into the new initial_enclosures
	if (!hasSplit)
		initial_enclosures.push_back(EnclosureType(target_loc,target_encl));
}

RealConstantIntMap
getSplitFactorsOfConstants(SystemType& system,
						   const RealConstantSet& locked_constants,
						   const Float& targetRatioPerc,
						   const HybridBoxes& bounding_domain)
{
	/*! Procedure:
	 * 1) Get the derivatives bounds for the system with constants set as their intervals
	 * 2) Get the derivatives bounds for the system with constants set as their midpoints
	 * 3) Get the maximum ratio from this bounds sets and set the target ratio as a fraction of this value
	 * 4) While the ratio is greater than the target ratio:
 	 *	 a) For each non-singleton accessible constant
 	 *	   i) Get the derivatives bounds for the system where the range of the chosen constant has been halved
 	 *	   ii) Update the best ratio
 	 *	 b) Increase the factor of the constant having the best ratio and update the system accordingly by halving the range of the constant
	 */

	RealConstantIntMap result;

	RealConstantSet working_constants = system.nonsingleton_accessible_constants();

	// Remove the locked constants
	for (RealConstantSet::const_iterator locked_constant_it = locked_constants.begin();
										 locked_constant_it != locked_constants.end();
										 ++locked_constant_it) {
		RealConstantSet::iterator original_constant_it = working_constants.find(*locked_constant_it);

		if (original_constant_it != working_constants.end())
			working_constants.erase(*original_constant_it);
	}

	if (working_constants.empty())
		return result;

	// Initializes the result and sets the system with the midpoints of the corresponding intervals
	for (RealConstantSet::const_iterator constant_it = working_constants.begin();
												 constant_it != working_constants.end();
												 ++constant_it) {
			result.insert(std::pair<RealConstant,int>(*constant_it,1));
			system.substitute(*constant_it,constant_it->value().midpoint());
	}

	// Gets the derivative widths corresponding to all accessible constants having midpoint value
	HybridFloatVector mid_der_widths = getDerivativeWidths(system, bounding_domain);
	// Restores the system to the original values
	system.substitute(working_constants);

	// While the ratio is sufficiently high, gets the best constant and substitutes half its interval into the system
	// If the maximum ratio is zero, then no constant affects the derivatives and splitting them would neither be necessary nor correct
	// given the current unfair implementation of _getBestConstantToSplit
	Float maxRatio = getMaxDerivativeWidthRatio(system, mid_der_widths,bounding_domain);
	if (maxRatio > 0) {
		Float ratio = maxRatio;
		while (ratio > targetRatioPerc*maxRatio) {
			RealConstant bestConstant = getBestConstantToSplit(system, working_constants, mid_der_widths, bounding_domain);
			Interval originalInterval = bestConstant.value();
			Float quarterIntervalWidth = originalInterval.width()/4;
			Interval halvedInterval = Interval(originalInterval.midpoint()-quarterIntervalWidth,
											   originalInterval.midpoint()+quarterIntervalWidth);
			system.substitute(bestConstant,Real(halvedInterval));
			result[bestConstant]++;
			ratio = getMaxDerivativeWidthRatio(system, mid_der_widths, bounding_domain);
		}

		system.substitute(working_constants);
	}

	return result;
}

RealConstant
getBestConstantToSplit(SystemType& system,
					   const RealConstantSet& working_constants,
					   const HybridFloatVector& referenceWidths,
					   const HybridBoxes& bounding_domain)
{
	RealConstant bestConstant = *working_constants.begin();
	Float bestLocalRatio = std::numeric_limits<Float>::infinity();

	for (RealConstantSet::const_iterator constant_it = working_constants.begin();
												 constant_it != working_constants.end();
												 ++constant_it) {
		// Modifies the system in order to have the range of the original given constant halved
		Real originalValue = system.accessible_constant_value(constant_it->name());
		Float quarterIntervalWidth = originalValue.width()/4;
		Interval halvedInterval = Interval(constant_it->value().midpoint()-quarterIntervalWidth,
										   constant_it->value().midpoint()+quarterIntervalWidth);
		system.substitute(*constant_it,Real(halvedInterval));

		Float localRatio = getMaxDerivativeWidthRatio(system,referenceWidths,bounding_domain);
		if (localRatio < bestLocalRatio) {
			bestLocalRatio = localRatio;
			bestConstant = RealConstant(constant_it->name(),originalValue);
		}

		// Restores the related constant to its original value
		system.substitute(*constant_it,originalValue);
	}

	return bestConstant;
}

HybridFloatVector
getDerivativeWidths(const HybridAutomaton& system,
				    const HybridBoxes& bounding_domain)
{
	ARIADNE_ASSERT_MSG(bounding_domain.size() == system.state_space().size(), "The bounding domain must be defined");

	HybridFloatVector result;

	// Gets the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = system.state_space().locations_begin()->second;

	for (list<DiscreteMode>::const_iterator modes_it = system.modes().begin(); modes_it != system.modes().end(); modes_it++) {
		const DiscreteState& loc = modes_it->location();

		// Gets the first order derivatives in respect to the dynamic of the mode, applied to the domain of the corresponding location
		Vector<Interval> der = modes_it->dynamic()(bounding_domain.find(loc)->second);

		Vector<Float> der_widths(css);
		for (uint i=0;i<css;i++)
			der_widths[i] = der[i].width();

		result.insert(pair<DiscreteState,Vector<Float> >(loc,der_widths));
	}

	return result;
}

Float
getMaxDerivativeWidthRatio(const HybridAutomaton& system,
						   const HybridFloatVector& referenceWidths,
						   const HybridBoxes& bounding_domain)
{
	ARIADNE_ASSERT_MSG(bounding_domain.size() == system.state_space().size(), "The bounding domain must be defined.");

	Float result = 0;

	// Gets the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = system.state_space().locations_begin()->second;

	// For each location and dimension, updates the result with the (w - wm)/wm derivative width ratio, excluding the undefined wm = 0 case
	for (list<DiscreteMode>::const_iterator modes_it = system.modes().begin(); modes_it != system.modes().end(); modes_it++) {
		const DiscreteState& loc = modes_it->location();

		Vector<Interval> der = modes_it->dynamic()(bounding_domain.find(loc)->second);

		for (uint i=0; i<css; ++i) {
			Float referenceWidth = referenceWidths.find(loc)->second[i];
			if (referenceWidth != 0)
				result = max(result,(der[i].width()-referenceWidth)/referenceWidth);
		}

	}

	return result;
}

/*! \brief Splits a RealConstant \a con into \a numParts parts.
 * \details Orders the subintervals by putting the second leftmost subinterval up to the rightmost, followed by the leftmost. */
std::vector<RealConstant>
split(const RealConstant& con, uint numParts)
{
	Interval bounds;
	Float lower, upper;

	std::vector<RealConstant> result(numParts,con);

	String name = con.name();
	Float intervalWidth = con.value().width();

	// Puts the first element
	lower = con.value().lower();
	upper = con.value().lower() + intervalWidth/numParts;
	result[numParts-1] = RealConstant(name,Interval(lower,upper));
	// Puts the last to the second element, in inverse order
	for (uint i=numParts;i>1;--i) {
		lower = con.value().lower() + intervalWidth*(i-1)/numParts;
		upper = con.value().lower() + intervalWidth*i/numParts;
		result[i-2] = RealConstant(name,Interval(lower,upper));
	}

	return result;
}

void _fillSplitSet(const std::vector<std::vector<RealConstant> >& src,
				   std::vector<std::vector<RealConstant> >::iterator col_it,
				   std::vector<RealConstant>::iterator row_it,
				   RealConstantSet s,
				   std::list<RealConstantSet>& dest)
{
	if (col_it != src.end() && row_it != col_it->end()) {
		_fillSplitSet(src,col_it,row_it+1,s,dest);
		s.insert(*row_it);
		col_it++;
		row_it = col_it->begin();

		if (col_it != src.end())
			_fillSplitSet(src,col_it,row_it,s,dest);
		else
			dest.push_back(s);
	}
}

std::list<RealConstantSet>
HybridReachabilityAnalyser::_getSplitConstantsSet() const
{
	std::list<RealConstantSet> result;

	if (_parameters->split_factors.empty())
		return result;

	// Creates a vector for all the interval splits (i.e. a jagged matrix)
	std::vector<std::vector<RealConstant> > split_intervals_set(_parameters->split_factors.size());
	uint i=0;
	for (RealConstantIntMap::const_iterator factor_it = _parameters->split_factors.begin();
											factor_it != _parameters->split_factors.end();
											++factor_it)
		split_intervals_set[i++] = split(factor_it->first, factor_it->second);

	// Generates all the possible split combinations
	RealConstantSet initial_combination;
	std::vector<std::vector<RealConstant> >::iterator initial_col_it = split_intervals_set.begin();
	std::vector<RealConstant>::iterator initial_row_it = initial_col_it->begin();
	_fillSplitSet(split_intervals_set,initial_col_it,initial_row_it,initial_combination,result);

	return result;
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
	_parameters->bounding_domain = bounding_set;

	return chain_reach(system,initial_set);
}


std::pair<HybridReachabilityAnalyser::SetApproximationType,bool>
HybridReachabilityAnalyser::
_upper_chain_reach_forward_domainCheck(const SystemType& system,
            const HybridImageSet& initial_set) const
{
	/* Complete procedure:
		1) Build the working sets from the initial enclosures
		2) While a working set exists, evolve the current working set, with the following termination conditions:
			a) If the evolution time limit is reached, adjoin both the reached set and the final set
			b) If the enclosure of the initial set is larger than the maximum enclosure allowed, skip the evolution and adjoin the final set
			c) If a blocking event is definitely initially active, or definitely finally active due to flow, adjoin the reached set only and discard the final set
		3) Discretise the reached and final sets into cells
		4) Remove the previous intermediate cells from the final cells and remove the previous reached cells from the reached cells
		5) Put the enclosures of the transitioned reach cells into the new initial enclosures, removing those outside the domain in either the source or destination location
		6) Put the enclosures of the final cells into the new initial enclosures, removing those outside the domain
		7) If new initial enclosures exist, restart from 1), otherwise terminate.
	*/

	// Initialize the validity flag (will be invalidated if any restriction of the reached region is performed)
	bool isValid = true;

    const Float& lock_to_grid_time = _parameters->lock_to_grid_time;
    const int& lock_to_grid_steps = _parameters->lock_to_grid_steps;
    const int& maximum_grid_depth = _parameters->maximum_grid_depth;

    HybridGrid grid=system.grid();
    HybridGridTreeSet new_final(grid), new_reach(grid), reach(grid), intermediate(grid);

    // Split the initial set into enclosures that are smaller than the minimum cell, given the grid
    list<EnclosureType> working_enclosures = _split_initial_set(initial_set,grid);

    ARIADNE_LOG(5,"Computing recurrent evolution...\n");
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,lock_to_grid_steps);

	// While the final set has new cells in respect to the previous reach set, process them
    for (uint i=0; !working_enclosures.empty(); i++)
	{
    	ARIADNE_LOG(5,"Iteration " << i << "\n");

        ARIADNE_LOG(6,"Initial enclosures size = " << working_enclosures.size() << "\n");

		// Evolve the initial enclosures
        make_lpair(new_reach,new_final)=_upper_reach_evolve_continuous(system,working_enclosures,hybrid_lock_to_grid_time,maximum_grid_depth);
        working_enclosures.clear();

		// Remove from the result those cells that have already been evolved from or checked for activations
        new_final.remove(intermediate);
		new_reach.remove(reach);

	    ARIADNE_LOG(6,"Reach size after removal = "<<new_reach.size()<<"\n");
	    ARIADNE_LOG(6,"Final size after removal = "<<new_final.size()<<"\n");

	    // Recombines in order to start checking activations on the largest possible enclosures
		new_reach.recombine();
		ARIADNE_LOG(6,"Reach size after recombining = "<<new_reach.size()<<"\n");
		bool reachValid = _upper_chain_reach_forward_pushTargetCells(new_reach,system,working_enclosures);
		if (_upper_chain_reach_break(reachValid,isValid))
			break;

		// Minces in order to identify possible enclosures lying outside the domain
		new_final.mince(maximum_grid_depth);
		ARIADNE_LOG(6,"Final size after mincing = "<<new_final.size()<<"\n");
		bool finalValid = _upper_chain_reach_pushLocalFinalCells(new_final,working_enclosures);
		if (_upper_chain_reach_break(finalValid,isValid))
			break;

		ARIADNE_LOG(6,"Final enclosures size = " << working_enclosures.size() << "\n");

        reach.adjoin(new_reach);
		intermediate.adjoin(new_final);

		reach.recombine();
		intermediate.recombine();
    }

	ARIADNE_LOG(5,"Found a total of " << reach.size() << " reached cells.\n");

    return make_pair<GTS,bool>(reach,isValid);
}

bool
HybridReachabilityAnalyser::_upper_chain_reach_break(const bool& flagToCheck, bool& flagToUpdate) const
{
	bool hasToBreak = false;

	if (_parameters->skip_if_unprovable && !flagToCheck) {
		flagToUpdate = false;
		hasToBreak = true;
	}

	return hasToBreak;
}

bool
HybridReachabilityAnalyser::_upper_chain_reach_forward_pushTargetCells(const HybridGridTreeSet& reachCells,
															  const SystemType& system,
															  std::list<EnclosureType>& destination) const
{
	bool reachIsValid = true;

	HybridGrid grid = system.grid();

	for (GTS::const_iterator cell_it = reachCells.begin(); cell_it != reachCells.end(); cell_it++)
	{
		// Get the location
		const DiscreteState& loc = cell_it->first;
		// Get the box
		const Box& bx = cell_it->second.box();

		// If the cell box lies outside the bounding domain (i.e. not inside and disjoint)
		// of the target location, notify
		if (!bx.inside(_parameters->bounding_domain[loc]) &&
			bx.disjoint(_parameters->bounding_domain[loc])) {
			ARIADNE_LOG(7,"Discarding reach enclosure " << bx <<
						  " being outside the domain in location " << loc.name() << ".\n");
			reachIsValid = false;

			// Return false if we skip when the system is unprovable
			if (_parameters->skip_if_unprovable)
				break;
		}

		ARIADNE_LOG(7,"Checking transitions for box "<< bx <<" in location " << loc.name() << "\n");
		bool transitionsEnclosuresAreValid = _upper_chain_reach_pushTargetEnclosures(system.transitions(loc),ContinuousEnclosureType(bx),grid,destination);
		if (_upper_chain_reach_break(transitionsEnclosuresAreValid,reachIsValid))
			break;
	}

	return reachIsValid;
}

bool
HybridReachabilityAnalyser::_upper_chain_reach_pushTargetEnclosures(const std::list<DiscreteTransition>& transitions,
																	const ContinuousEnclosureType& source,
																	const HybridGrid& grid,
																	std::list<EnclosureType>& destination) const
{
	bool enclosuresAreValid = true;

    long numCellDivisions = (1<<_parameters->maximum_grid_depth);

	// For each transition
	for (list<DiscreteTransition>::const_iterator trans_it = transitions.begin(); trans_it != transitions.end(); trans_it++)
	{
		ARIADNE_LOG(8,"Target: "<<trans_it->target()<<"\n");

		Vector<Float> minTargetCellWidths = grid[trans_it->target()].lengths()/numCellDivisions;

		bool currentTransitionEnclosuresAreValid = _upper_chain_reach_pushTargetEnclosuresOfTransition(*trans_it,source,minTargetCellWidths,destination);
		if (_upper_chain_reach_break(currentTransitionEnclosuresAreValid,enclosuresAreValid))
			break;
	}

	return enclosuresAreValid;
}

bool
HybridReachabilityAnalyser::_upper_chain_reach_pushTargetEnclosuresOfTransition(const DiscreteTransition& trans,
																				const ContinuousEnclosureType& source,
																				const Vector<Float>& minTargetCellWidths,
																				std::list<EnclosureType>& destination) const
{
	bool enclosuresAreValid = true;

	const Box& target_bounding = _parameters->bounding_domain[trans.target()];

	// Helper class for operations on Taylor sets
	TaylorCalculus tc;

	std::list<ContinuousEnclosureType> split_list;
	split_list.push_front(source);

	while (!split_list.empty()) {

		ContinuousEnclosureType current_encl = split_list.back();
		split_list.pop_back();

		TaylorModel guard_time_model = apply(trans.activation(),current_encl)[0];
		ARIADNE_LOG(9,"Guard time model = "<<guard_time_model<<" (range: " << guard_time_model.range() << ")\n");

		// If definitely active, add the target enclosure on the splittings of the enclosure
		// If possibly active, split the enclosure if possible, otherwise create the target enclosure
		// (if otherwise the transition is not active, do nothing)

		if (guard_time_model.range().lower()>0) {
			ARIADNE_LOG(10,"Definitely active, adding the splittings of the target enclosure to the destination list.\n");
			// Get the target continuous enclosure
			ContinuousEnclosureType target_encl = tc.reset_step(trans.reset(),current_encl);
			// Populate the initial enclosures with the hybrid enclosures derived from the splittings of the target enclosure
			bool currentEnclosuresAreValid = true;
			splitTargetEnclosures(currentEnclosuresAreValid,destination,trans.target(),target_encl,minTargetCellWidths,target_bounding);

			if (_upper_chain_reach_break(currentEnclosuresAreValid,enclosuresAreValid))
				break;
		} else if (!(guard_time_model.range().upper()<0)) {
			// Try splitting the enclosure
			bool hasSplit = false;
			for (uint i=0; i < minTargetCellWidths.size(); i++) {
				// If the enclosure has width larger than that of the minimum cell, split on that dimension
				if (current_encl.bounding_box()[i].width() > minTargetCellWidths[i]) {
					hasSplit = true;
					std::pair<ContinuousEnclosureType,ContinuousEnclosureType> split_sets = current_encl.split(i);
					split_list.push_front(split_sets.first);
					split_list.push_front(split_sets.second);
					break;
				}
			}
			// If we could not split
			if (!hasSplit) {
				ARIADNE_LOG(10,"Possibly active but minimum enclosure size reached, adding the splittings of the target enclosure to the destination list.\n");
				bool transitionedEnclosureIsValid = _upper_chain_reach_pushTransitioningEnclosure(trans,current_encl,tc,destination);
				if (_upper_chain_reach_break(transitionedEnclosureIsValid,enclosuresAreValid))
					break;
			}
			else {
				ARIADNE_LOG(10,"Possibly active, adding the splittings of the current enclosure to the evaluation list.\n");
			}
		} else {
			ARIADNE_LOG(10,"Inactive, discarding.\n");
		}
	}

	return enclosuresAreValid;
}

bool
HybridReachabilityAnalyser::_upper_chain_reach_pushTransitioningEnclosure(const DiscreteTransition& trans,
																	      const ContinuousEnclosureType& encl,
																		  const TaylorCalculus& tc,
																		  std::list<EnclosureType>& destination) const
{
	bool isValid = true;

	ContinuousEnclosureType target_encl = tc.reset_step(trans.reset(),encl);
	const Box& target_encl_box = target_encl.bounding_box();
	const Box& target_bounding = _parameters->bounding_domain[trans.target()];

	// If the cell box lies outside the bounding domain (i.e. not inside and disjoint)
	// of the target location, notify, otherwise add it
	if (!target_encl_box.inside(target_bounding) && target_encl_box.disjoint(target_bounding)) {
		ARIADNE_LOG(10,"Discarding target enclosure " << target_encl_box <<
					  " being outside the domain in location " << trans.target().name() << ".\n");
		isValid = false;
	} else {
		destination.push_back(EnclosureType(trans.target(),target_encl));
	}

	return isValid;
}


bool
HybridReachabilityAnalyser::_upper_chain_reach_pushLocalFinalCells(const HybridGridTreeSet& finalCells,
																		    std::list<EnclosureType>& destination) const
{
	bool isValid = true;

	for (GTS::const_iterator cell_it = finalCells.begin(); cell_it != finalCells.end(); ++cell_it) {
		const DiscreteState& loc = cell_it->first;
		const Box& bounding = _parameters->bounding_domain[loc];

		if (!cell_it->second.box().inside(bounding) && cell_it->second.box().disjoint(bounding)) {
			ARIADNE_LOG(7,"Discarding enclosure " << cell_it->second.box() <<
						" from final cell outside the domain in location " << loc.name() <<".\n");
			isValid = false;

			if (_parameters->skip_if_unprovable)
				break;
		} else {
			destination.push_back(_discretiser->enclosure(*cell_it));
		}
	}

	return isValid;
}


void
HybridReachabilityAnalyser::_upper_chain_reach_pushFinalCells_noDomainCheck(const HybridGridTreeSet& finalCells,
																		    std::list<EnclosureType>& destination) const
{
	for (GTS::const_iterator cell_it = finalCells.begin(); cell_it != finalCells.end(); ++cell_it)
			destination.push_back(_discretiser->enclosure(*cell_it));
}


std::pair<HybridReachabilityAnalyser::SetApproximationType,bool>
HybridReachabilityAnalyser::
upper_chain_reach_forward(SystemType& system,
				  const HybridImageSet& initial_set) const
{
	UpperChainReachFuncPtr func = &Ariadne::HybridReachabilityAnalyser::_upper_chain_reach_forward_domainCheck;
	return _upper_chain_reach(system,initial_set,func);
}

std::pair<HybridReachabilityAnalyser::SetApproximationType,bool>
HybridReachabilityAnalyser::
_upper_chain_reach(SystemType& system,
				   const HybridImageSet& initial_set,UpperChainReachFuncPtr func) const
{
	// The reachable set
	HybridGridTreeSet reach;
	// The flag that informs whether the region is valid for proving
	bool isValid;

	std::list<RealConstantSet> split_set = _getSplitConstantsSet();

	// If no split set exists, performs the reachability analysis on the system, otherwise checks each element in the set
	if (split_set.empty()) {
		make_lpair<HybridGridTreeSet,bool>(reach,isValid) = _upper_chain_reach_forward_domainCheck(system,initial_set);
	} else {
		isValid = true;
		RealConstantSet original_constants = system.nonsingleton_accessible_constants();

		uint i = 0;
		// Progressively adds the results for each subsystem
		for (std::list<RealConstantSet>::const_iterator set_it = split_set.begin(); set_it != split_set.end(); ++set_it) {
			ARIADNE_LOG(5,"<Constants set #" << i++ << " : " << *set_it << " >\n");

			system.substitute(*set_it);

			std::pair<HybridGridTreeSet,bool> reachAndValidity = (this->*func)(system,initial_set);
			const HybridGridTreeSet& local_reach = reachAndValidity.first;
			const bool& local_isValid = reachAndValidity.second;

			reach.adjoin(local_reach);

			isValid = isValid && local_isValid;

			// We skip if allowed and if either the partial result is not valid or outside the safe region
			if (_parameters->skip_if_unprovable) {
				if (!isValid || (isValid && definitely(!local_reach.subset(_parameters->safe_region)))) {
					isValid = false;
					break;
				}
			}
		}

		system.substitute(original_constants);
	}

	return std::pair<HybridReachabilityAnalyser::SetApproximationType,bool>(reach,isValid);
}

std::pair<HybridReachabilityAnalyser::SetApproximationType,DisproveData>
HybridReachabilityAnalyser::
_lower_chain_reach(const SystemType& system,
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

    DisproveData globalFalsInfo(system.state_space());

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
		// The falsification info
		DisproveData newFalsInfo(system.state_space());

		ARIADNE_LOG(6,"Initial enclosures size = " << initial_enclosures.size() << "\n");

		// Create the worker
		LowerChainReachWorker worker(_discretiser,initial_enclosures,system,lock_time,
									 _parameters->maximum_grid_depth,concurrency, disprove_bounds, _parameters->skip_if_disproved);

		ARIADNE_LOG(6,"Evolving and discretising...\n");

		// Compute and get the result
		make_ltuple<std::pair<HUM,HUM>,EL,GTS,DisproveData>(evolve_sizes,final_enclosures,
																new_reach,newFalsInfo) = worker.get_result();


		ARIADNE_LOG(6,"Reach size before removal = " << new_reach.size() << "\n");

		// Update the falsification info
		globalFalsInfo.updateWith(newFalsInfo);

		new_reach.remove(reach);
		ARIADNE_LOG(6,"Reach size after removal  = " << new_reach.size() << "\n");
		if (new_reach.empty())
			break;

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

	return make_pair<GTS,DisproveData>(reach,globalFalsInfo);
}

std::pair<HybridReachabilityAnalyser::SetApproximationType,DisproveData>
HybridReachabilityAnalyser::
lower_chain_reach(SystemType& system,
				  const HybridImageSet& initial_set) const
{
	HybridGridTreeSet reach;
	DisproveData disproveData(system.state_space());

	std::list<RealConstantSet> split_set = _getSplitConstantsSet();

	// If no split set exists, performs the reachability analysis on the system, otherwise checks each element in the set
	if (split_set.empty()) {
		make_lpair<HybridGridTreeSet,DisproveData>(reach,disproveData) = _lower_chain_reach(system,initial_set);
	} else {
		RealConstantSet original_constants = system.nonsingleton_accessible_constants();

		uint i = 0;
		// Progressively adds the results for each subsystem
		for (std::list<RealConstantSet>::const_iterator set_it = split_set.begin(); set_it != split_set.end(); ++set_it)
		{
			ARIADNE_LOG(5,"<Constants set #" << i++ << " : " << *set_it << " >\n");

			system.substitute(*set_it);

			std::pair<HybridGridTreeSet,DisproveData> reachAndDisproveData = _lower_chain_reach(system,initial_set);

			reach.adjoin(reachAndDisproveData.first);
			disproveData.updateWith(reachAndDisproveData.second);

			if (_parameters->skip_if_disproved && reachAndDisproveData.second.getIsDisproved())
				break;
		}
		system.substitute(original_constants);
	}

	_statistics->lower().reach = reach;

	return std::pair<HybridGridTreeSet,DisproveData>(reach,disproveData);
}

bool 
HybridReachabilityAnalyser::
_prove(SystemType& system,
	   const HybridImageSet& initial_set)
{
	bool result;

	ARIADNE_LOG(4,"Proving... " << (verbosity == 4 ? "" : "\n"));

	std::pair<HybridGridTreeSet,bool> reachAndIsValid = upper_chain_reach_forward(system,initial_set);

	// Proved iff the reached region is valid (i.e. it has not been restricted) and is inside the safe region
	result = (reachAndIsValid.second && definitely(reachAndIsValid.first.subset(_parameters->safe_region)));

	ARIADNE_LOG((verbosity == 4 ? 1 : 4), (result ? "Proved.\n" : "Not proved.\n") );

	return result;
}

bool 
HybridReachabilityAnalyser::
_disprove(SystemType& system,
		  const HybridImageSet& initial_set)
{
	ARIADNE_LOG(4,"Disproving... " << (verbosity == 4 ? "" : "\n"));

	std::pair<HybridGridTreeSet,DisproveData> reachAndDisproveData = lower_chain_reach(system,initial_set);
	const bool& isDisproved = reachAndDisproveData.second.getIsDisproved();

	ARIADNE_LOG((verbosity == 4 ? 1 : 4), (isDisproved ? "Disproved.\n" : "Not disproved.\n") );

	return isDisproved;
}


tribool
HybridReachabilityAnalyser::
verify(SystemType& system,
       const HybridImageSet& initial_set)
{
		ARIADNE_LOG(3, "Verification...\n");

		// Return true if proved
		if (_prove(system,initial_set)) {
			ARIADNE_LOG(3, "Safe.\n");
			return true;
		}

		// Return false if disproved
		if (_disprove(system,initial_set)) {
			ARIADNE_LOG(3, "Unsafe.\n");
			return false;
		}

		// Return indeterminate otherwise
		ARIADNE_LOG(3, "Indeterminate.\n");
		return indeterminate;
}


HybridFloatVector 
HybridReachabilityAnalyser::
_getHybridMaximumAbsoluteDerivatives(const HybridAutomaton& system) const
{
	HybridFloatVector result;

	// Get the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = system.state_space().locations_begin()->second;

	// The variable for the bounding box of the derivatives
	Vector<Interval> der;

	HybridGridTreeSet& hybridreach = _statistics->upper().reach;

	// For each mode
	for (list<DiscreteMode>::const_iterator modes_it = system.modes().begin(); modes_it != system.modes().end(); modes_it++) {

		const DiscreteState& loc = modes_it->location();

		// Insert the corresponding hmad pair, initialized with zero maximum absolute derivatives
		result.insert(pair<DiscreteState,Vector<Float> >(loc,Vector<Float>(css)));

		// If the reached region for the location exists and is not empty, check its cells, otherwise use the whole domain
		if (hybridreach.has_location(loc) && !hybridreach[loc].empty()) {
			// Get the GridTreeSet
			GridTreeSet& reach = hybridreach[loc];
			// For each of its hybrid cells
			for (GridTreeSet::const_iterator cells_it = reach.begin(); cells_it != reach.end(); cells_it++) {
				// Gets the derivative bounds
				der = modes_it->dynamic()(cells_it->box());

				// For each variable, sets the maximum value
				for (uint i=0;i<css;i++)
					result[loc][i] = max(result[loc][i], abs(der[i]).upper());
			}
		} else {
			// Gets the first order derivatives in respect to the dynamic of the mode, applied to the domain of the corresponding location
			der = modes_it->dynamic()(_parameters->bounding_domain.find(loc)->second);

			// Gets the maximum absolute derivatives
			for (uint i=0;i<css;i++)
				result[loc][i] = abs(der[i]).upper();
		}
	}

	// Returns
	return result;
}

void
HybridReachabilityAnalyser::
_setLockToGridTime(const HybridAutomaton& system) const
{
	// Get the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = system.state_space().locations_begin()->second;

	// The variable for the bounding box of the derivatives
	Vector<Interval> der;
    // The variable for the result
	Float result = 0;

	// For each mode
	for (list<DiscreteMode>::const_iterator modes_it = system.modes().begin(); modes_it != system.modes().end(); modes_it++)
	{
		// Gets the location
		const DiscreteState& loc = modes_it->location();

		// Gets the domain for this mode
		const Box& domain = _parameters->bounding_domain.find(loc)->second;

		// Gets the first order derivatives in respect to the dynamic of the mode, applied to the domain of the corresponding location
		der = modes_it->dynamic()(domain);

		// Updates the lock time
		for (uint i=0;i<css;i++)
		{
			Float maxAbsDer = abs(der[i]).upper();
			if (maxAbsDer > 0)
				result = max(result,domain[i].width()/maxAbsDer);
		}
	}

	_parameters->lock_to_grid_time = result;
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

	// Assigns the cell
	// NOTE: it is preferable to have the multiplier slightly lesser than an integer multiple of the grid cell, so that
	// the overapproximation error due to discretization is minimized
	_discretiser->parameters().maximum_enclosure_cell = 1.9*mec;
}


std::map<DiscreteState,Float>
HybridReachabilityAnalyser::
_getHybridMaximumStepSize(const HybridFloatVector& hmad, const HybridGrid& hgrid)
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

	return hmss;
}

HybridGrid 
HybridReachabilityAnalyser::
_getHybridGrid(const HybridFloatVector& hmad) const
{
	// Get the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = hmad.begin()->second.size();

	HybridGrid hg;

	// The lengths of the grid cell
	std::map<DiscreteState,Vector<Float> > hybridgridlengths;

	// Get the minimum domain length for each variable
	Vector<Float> minDomainLengths(css);
	for (uint i=0;i<css;i++) {
		minDomainLengths[i] = std::numeric_limits<double>::infinity();
	}
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

	// Populate the grid, centered on the centre of the domain
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) {
		const DiscreteState& loc = hfv_it->first;
		hg[loc] = Grid(_parameters->bounding_domain.find(loc)->second.centre(),hybridgridlengths[loc]);
	}

	return hg;
}

void
HybridReachabilityAnalyser::
_setInitialParameters(SystemType& system,
					  const HybridBoxes& domain,
					  const HybridBoxes& safe_region,
					  const RealConstantSet& locked_constants)
{
	// Set the domain
	_parameters->bounding_domain = domain;
	ARIADNE_LOG(3, "Domain: " << _parameters->bounding_domain << "\n");
	// Set the safe region
	_parameters->safe_region = safe_region;
	ARIADNE_LOG(3, "Safe region: " << _parameters->safe_region << "\n");
	// Set the lock to grid time
	_setLockToGridTime(system);
	ARIADNE_LOG(3, "Lock to grid time: " << _parameters->lock_to_grid_time << "\n");
	// Set the split factors
	_parameters->split_factors = getSplitFactorsOfConstants(system,locked_constants,0.01,_parameters->bounding_domain);
	ARIADNE_LOG(3, "Split factors: " << _parameters->split_factors << "\n");
}


void 
HybridReachabilityAnalyser::
_tuneIterativeStepParameters(SystemType& system)
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
	_discretiser->parameters().hybrid_maximum_step_size = _getHybridMaximumStepSize(hmad,system.grid());
	ARIADNE_LOG(3, "Maximum step size: " << _discretiser->parameters().hybrid_maximum_step_size << "\n");
}

void
HybridReachabilityAnalyser::
_setDominanceParameters(SystemVerificationInfo& verInfo, const RealConstantSet& lockedConstants)
{
	// Clears the previously reached regions
	this->_statistics->upper().reach = HybridGridTreeSet();
	this->_statistics->lower().reach = HybridGridTreeSet();
	// Domain
	_parameters->bounding_domain = verInfo.getDomain();
	// Dummy safe region
	_parameters->safe_region = verInfo.getSafeRegion();
	// General parameters
	_tuneIterativeStepParameters(verInfo.getSystem());
	// Lock to grid time
	_setLockToGridTime(verInfo.getSystem());
	// Split factors
	_parameters->split_factors = getSplitFactorsOfConstants(verInfo.getSystem(),lockedConstants,0.01,_parameters->bounding_domain);
}

tribool
HybridReachabilityAnalyser::
verify_iterative(SystemVerificationInfo& verInfo)
{
	RealConstantSet locked_constants;

	return _verify_iterative(verInfo,locked_constants);
}

tribool
HybridReachabilityAnalyser::
verify_iterative(SystemVerificationInfo& verInfo, const RealConstant& parameter)
{
	HybridAutomaton& system = verInfo.getSystem();

	Real original_value = system.accessible_constant_value(parameter.name());

	system.substitute(parameter);
	tribool result = verify_iterative(verInfo);
	system.substitute(parameter,original_value);

	return result;
}

tribool
HybridReachabilityAnalyser::
verify_iterative(SystemVerificationInfo& verInfo, const RealConstant& parameter, const Float& value)
{
	const RealConstant modifiedParameter(parameter.name(),Interval(value));

	return verify_iterative(verInfo, modifiedParameter);
}

tribool 
HybridReachabilityAnalyser::
_verify_iterative(SystemVerificationInfo& verInfo,
				 const RealConstantSet& locked_constants)
{
	ARIADNE_LOG(2,"\nIterative verification...\n");

	HybridAutomaton& system = verInfo.getSystem();

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

	// Reset the reach statistics
	_statistics->upper().reach = HybridGridTreeSet();
	_statistics->lower().reach = HybridGridTreeSet();

	// Set the initial parameters
	_setInitialParameters(system, verInfo.getDomain(), verInfo.getSafeRegion(), locked_constants);

	int& depth = _parameters->maximum_grid_depth;
    for(depth = _parameters->lowest_maximum_grid_depth; depth <= _parameters->highest_maximum_grid_depth; ++depth)
	{ 
    	sprintf(mgd_char,"%i", depth);
		ARIADNE_LOG(2, "DEPTH " << depth << "\n");

		// Tune the parameters for the current iteration
		_tuneIterativeStepParameters(system);

		// Perform the verification
		tribool result = verify(system,verInfo.getInitialSet());

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


std::pair<Interval,Interval>
HybridReachabilityAnalyser::
parametric_verification_1d_bisection(SystemVerificationInfo& verInfo,
						   	   	   	 const RealConstant& parameter,
						   	   	   	 const Float& tolerance)
{
	ARIADNE_ASSERT(tolerance > 0 && tolerance < 1);

	HybridAutomaton& system = verInfo.getSystem();

	// Get the original value of the parameter and the related range
	Real original_value = system.accessible_constant_value(parameter.name());
	Interval parameter_range(parameter.value());

	ARIADNE_ASSERT(parameter_range.width() > 0);

	// Create the safety and unsafety intervals: they represent the search intervals,
	// NOT the intervals where the system is proved safe or unsafe
	Interval safety_int = parameter_range;
	Interval unsafety_int = parameter_range;

	ARIADNE_LOG(1,"\nChecking parameter " << parameter.name() << " in " << parameter_range << " with " << tolerance*100 << "% tolerance\n");

	// Check the lower bound
	ARIADNE_LOG(1,"\nChecking lower interval bound... ");
	tribool lower_result = verify_iterative(verInfo,parameter,parameter_range.lower());
	_log_verification_result(lower_result);

	// Check the upper bound
	ARIADNE_LOG(1,"Checking upper interval bound... ");
	tribool upper_result = verify_iterative(verInfo,parameter,parameter_range.upper());
	_log_verification_result(upper_result);

	// If we must proceed with bisection refining
	bool proceed;
	// If the safe value is found on the lower extreme of the parameter interval
	bool safeOnBottom;

	// Updates the initial values for the safety/unsafety intervals, and decides whether to proceed and the eventual ordering
	make_lpair<bool,bool>(proceed,safeOnBottom) = this->_process_initial_bisection_results(safety_int,unsafety_int,parameter_range,
																						 lower_result,upper_result);

	// Verification loop
	bool proceedSafety = proceed;
	bool proceedUnsafety = proceed;
	while (proceedSafety || proceedUnsafety)
	{
		// Safety interval check
		if (!safety_int.empty() && safety_int.width() > tolerance*parameter_range.width())
		{
			ARIADNE_LOG(1,"Checking safety (positive) interval " << safety_int << " (midpoint: " << safety_int.midpoint() <<
					      ", width ratio: " << safety_int.width()/parameter_range.width()*100 << "%) ... ");

			const Float current_value = safety_int.midpoint();
			tribool result = verify_iterative(verInfo,parameter,current_value);

			_process_positive_bisection_result(result,safety_int,unsafety_int,current_value,safeOnBottom);
		}
		else
			proceedSafety = false;

		// Unsafety interval check
		if (!unsafety_int.empty() && unsafety_int.width() > tolerance*parameter_range.width())
		{
			ARIADNE_LOG(1,"Checking unsafety (negative) interval " << unsafety_int << " (midpoint: " << unsafety_int.midpoint() <<
						  ", width ratio: " << unsafety_int.width()/parameter_range.width()*100 << "%) ... ");

			const Float current_value = unsafety_int.midpoint();
			tribool result = verify_iterative(verInfo,parameter,current_value);

			_process_negative_bisection_result(result,safety_int,unsafety_int,current_value,safeOnBottom);
		}
		else
			proceedUnsafety = false;
	}

	system.substitute(parameter,original_value);

	return _pos_neg_bounds_from_search_intervals(safety_int,unsafety_int,parameter_range,safeOnBottom);
}


std::pair<Interval,Interval>
HybridReachabilityAnalyser::
parametric_dominance_1d_bisection(SystemVerificationInfo& dominating,
								  SystemVerificationInfo& dominated,
						   	   	  const RealConstant& parameter,
						   	   	  const Float& tolerance)
{
	ARIADNE_ASSERT(tolerance > 0 && tolerance < 1);

	HybridAutomaton& system = dominating.getSystem();

	// Get the original value of the parameter and the related range
	Real original_value = system.accessible_constant_value(parameter.name());
	Interval parameter_range(parameter.value());

	ARIADNE_ASSERT(parameter_range.width() > 0);

	// Create the dominating and nondominating intervals: they represent the search intervals,
	// NOT the intervals where the system is proved dominating or nondominating
	Interval dominating_int = parameter_range;
	Interval nondominating_int = parameter_range;

	ARIADNE_LOG(1,"\nChecking parameter " << parameter.name() << " in " << parameter_range << " with " << tolerance*100 << "% tolerance\n");

	// Check the lower bound
	ARIADNE_LOG(1,"\nChecking lower interval bound... ");
	tribool lower_result = dominance(dominating,dominated,parameter,parameter_range.lower());
	_log_verification_result(lower_result);

	// Check the upper bound
	ARIADNE_LOG(1,"Checking upper interval bound... ");
	tribool upper_result = dominance(dominating,dominated,parameter,parameter_range.upper());
	_log_verification_result(upper_result);

	// If we must proceed with bisection refining
	bool proceed;
	// If the dominating value is found on the lower extreme of the parameter interval
	bool dominatingOnBottom;

	// Updates the initial values for the dominating/nondominating intervals, and decides whether to proceed and the eventual ordering
	make_lpair<bool,bool>(proceed,dominatingOnBottom) = this->_process_initial_bisection_results(dominating_int,nondominating_int,parameter_range,
																						 lower_result,upper_result);

	// Verification loop
	bool proceedDominating = proceed;
	bool proceedNondominating = proceed;
	while (proceedDominating || proceedNondominating)
	{
		// Dominating interval check
		if (!dominating_int.empty() && dominating_int.width() > tolerance*parameter_range.width())
		{
			ARIADNE_LOG(1,"Checking dominating (positive) interval " << dominating_int << " (midpoint: " << dominating_int.midpoint() <<
					      ", width ratio: " << dominating_int.width()/parameter_range.width()*100 << "%) ... ");

			const Float current_value = dominating_int.midpoint();
			tribool result = dominance(dominating,dominated,parameter,current_value);

			_process_positive_bisection_result(result,dominating_int,nondominating_int,current_value,dominatingOnBottom);
		}
		else
			proceedDominating = false;

		// Nondominating interval check
		if (!nondominating_int.empty() && nondominating_int.width() > tolerance*parameter_range.width())
		{
			ARIADNE_LOG(1,"Checking nondominating (negative) interval " << nondominating_int << " (midpoint: " << nondominating_int.midpoint() <<
						  ", width ratio: " << nondominating_int.width()/parameter_range.width()*100 << "%) ... ");

			const Float current_value = nondominating_int.midpoint();
			tribool result = dominance(dominating,dominated,parameter,current_value);

			_process_negative_bisection_result(result,dominating_int,nondominating_int,current_value,dominatingOnBottom);
		}
		else
			proceedNondominating = false;
	}

	system.substitute(parameter,original_value);

	return _pos_neg_bounds_from_search_intervals(dominating_int,nondominating_int,parameter_range,dominatingOnBottom);
}

Parametric2DBisectionResults
HybridReachabilityAnalyser::parametric_dominance_2d_bisection(SystemVerificationInfo& dominating,
															  SystemVerificationInfo& dominated,
											   const RealConstantSet& params,
											   const Float& tolerance,
											   const unsigned& numPointsPerAxis)
{
	ARIADNE_ASSERT_MSG(params.size() == 2,"Provide exactly two parameters.");

	RealConstantSet::const_iterator param_it = params.begin();

	const RealConstant& xParam = *param_it;
	const RealConstant& yParam = *(++param_it);

	return parametric_dominance_2d_bisection(dominating,dominated,xParam,yParam,tolerance,numPointsPerAxis);
}

Parametric2DBisectionResults
HybridReachabilityAnalyser::parametric_dominance_2d_bisection(SystemVerificationInfo& dominating,
															  SystemVerificationInfo& dominated,
											   const RealConstant& xParam,
											   const RealConstant& yParam,
											   const Float& tolerance,
											   const unsigned& numPointsPerAxis)
{
	// Generates the file name
	std::string filename = dominating.getSystem().name() + "&" + dominated.getSystem().name();
	filename = filename + "[" + xParam.name() + "," + yParam.name() + "]";

	// Initializes the results
	Parametric2DBisectionResults results(filename,xParam.value(),yParam.value(),numPointsPerAxis);

	// Sweeps on each axis
	_parametric_dominance_2d_bisection_sweep(results, dominating, dominated, xParam, yParam, tolerance, numPointsPerAxis, true);
	_parametric_dominance_2d_bisection_sweep(results, dominating, dominated, xParam, yParam, tolerance, numPointsPerAxis, false);

	return results;
}

void HybridReachabilityAnalyser::_parametric_dominance_2d_bisection_sweep(Parametric2DBisectionResults& results,
					  	  	  	    							SystemVerificationInfo& dominating,
					  	  	  	    							SystemVerificationInfo& dominated,
					  	  	  	    							RealConstant xParam,
					  	  	  	    							RealConstant yParam,
					  	  	  	    							const Float& tolerance,
					  	  	  	    							const uint& numPointsPerAxis,
					  	  	  	    							bool sweepOnX)
{
	HybridAutomaton& system = dominating.getSystem();

	RealConstant& sweepParam = (sweepOnX ? xParam : yParam);
	RealConstant& otherParam = (sweepOnX ? yParam : xParam);
	Real originalSweepValue = system.accessible_constant_value(sweepParam.name());
	Interval sweepBounds = sweepParam.value();

	ARIADNE_LOG(1,"\nSweeping on " << sweepParam.name() << " in " << sweepBounds << ":\n");

	// Sweeps in the given direction
	for (uint i=0; i<numPointsPerAxis; ++i)
	{
		// Changes the value
		sweepParam.set_value(sweepBounds.lower()+i*sweepBounds.width()/(numPointsPerAxis-1));
		// Modifies the system with the new value
		system.substitute(sweepParam);

		ARIADNE_LOG(1,"\nAnalysis with " << sweepParam.name() << " = " << sweepParam.value().lower() << "...\n");

		// Performs the analysis on the other axis and adds to the results of the sweep variable
		std::pair<Interval,Interval> result = parametric_dominance_1d_bisection(dominating,dominated,otherParam,tolerance);
		if (sweepOnX)
			results.insertXValue(result);
		else
			results.insertYValue(result);

		ARIADNE_LOG(1,"Obtained dominance in " << result.first << " and nondominance in " << result.second << ".\n");
	}

	// Restores the original value
	system.substitute(sweepParam,originalSweepValue);
}

void
HybridReachabilityAnalyser::_log_verification_result(const tribool& result) const
{
	if (definitely(result)) { ARIADNE_LOG(1,"True.\n"); }
	else if (!possibly(result)) { ARIADNE_LOG(1,"False.\n"); }
	else ARIADNE_LOG(1,"Indeterminate.\n");
}

std::pair<bool,bool>
HybridReachabilityAnalyser::_process_initial_bisection_results(Interval& positive_int,
															 Interval& negative_int,
															 const Interval& parameter_range,
															 const tribool& lower_result,
															 const tribool& upper_result) const
{
	bool proceed;
	bool positiveOnBottom;

	// If both extremes are safe, no more verification is involved
	if (definitely(lower_result) && definitely(upper_result)) {
		proceed = false;
		positiveOnBottom = false;
		positive_int = Interval(parameter_range.lower());
		negative_int.make_empty();
	}
	// If both extremes are unsafe, no more verification is involved
	else if (!possibly(lower_result) && !possibly(upper_result)) {
		proceed = false;
		positiveOnBottom = false;
		positive_int.make_empty();
		negative_int = Interval(parameter_range.upper());
	}
	// If both extremes are indeterminate, no verification is possible
	else if (indeterminate(lower_result) && indeterminate(upper_result)) {
		proceed = false;
		positiveOnBottom = false;
		positive_int.make_empty();
		negative_int.make_empty();
	}
	// If the lower extreme is safe or the upper extreme is unsafe, the safe values are on the bottom
	else if (definitely(lower_result) || !possibly(upper_result)) {
		proceed = true;
		positiveOnBottom = true;
		// If there are indeterminate values, reset the corresponding intervals as empty
		if (indeterminate(lower_result)) positive_int.make_empty();
		if (indeterminate(upper_result)) negative_int.make_empty();
	}
	// If the upper extreme is safe or the lower extreme is unsafe, the safe values are on the top
	else {
		proceed = true;
		positiveOnBottom = false;
		// If there are indeterminate values, reset the corresponding intervals as empty
		if (indeterminate(lower_result)) negative_int.make_empty();
		if (indeterminate(upper_result)) positive_int.make_empty();
	}

	return std::pair<bool,bool>(proceed,positiveOnBottom);
}

void
HybridReachabilityAnalyser::_process_positive_bisection_result(const tribool& result,
															   Interval& positive_int,
															   Interval& negative_int,
															   const Float& current_value,
															   const bool& positiveOnBottom) const
{
	if (definitely(result)) {
		if (positiveOnBottom) {
			// If the negative interval is the same as the positive one, update it too
			if (equal(negative_int,positive_int)) negative_int.set_lower(current_value);

			ARIADNE_LOG(1,"True, refining upwards.\n");
			positive_int.set_lower(current_value);
		} else {
			// If the negative interval is the same as the positive one, update it too
			if (equal(negative_int,positive_int)) negative_int.set_upper(current_value);

			ARIADNE_LOG(1,"True, refining downwards.\n");
			positive_int.set_upper(current_value);
		}
	}
	else if (!possibly(result)) {
		if (positiveOnBottom) {
			ARIADNE_LOG(1,"False, refining downwards and resetting the unsafety.\n");
			positive_int.set_upper(current_value);
		}
		else {
			ARIADNE_LOG(1,"False, refining upwards and resetting the unsafety.\n");
			positive_int.set_lower(current_value);
		}

		negative_int = positive_int;
	}
	else {
		if (positiveOnBottom) {
			// If the negative interval is the same as the positive interval, update it too
			if (equal(negative_int,positive_int)) negative_int.set_lower(current_value);

			ARIADNE_LOG(1,"Indeterminate, refining downwards.\n");
			positive_int.set_upper(current_value); }
		else {
			// If the negative interval is the same as the positive interval, update it too
			if (equal(negative_int,positive_int)) negative_int.set_upper(current_value);

			ARIADNE_LOG(1,"Indeterminate, refining upwards.\n");
			positive_int.set_lower(current_value); }
	}
}

void
HybridReachabilityAnalyser::_process_negative_bisection_result(const tribool& result,
															   Interval& positive_int,
															   Interval& negative_int,
															   const Float& current_value,
															   const bool& positiveOnBottom) const
{
	if (definitely(result)) {
		if (positiveOnBottom) {
			ARIADNE_LOG(1,"True, refining upwards and resetting the safety.\n");
			negative_int.set_lower(current_value);
		} else {
			ARIADNE_LOG(1,"True, refining downwards and resetting the safety.\n");
			negative_int.set_upper(current_value);
		}

		positive_int = negative_int;
	}

	else if (!possibly(result)) {
		if (positiveOnBottom) {
			// If the negative interval is the same as the positive interval, update it too
			if (equal(negative_int,positive_int)) positive_int.set_upper(current_value);

			ARIADNE_LOG(1,"False, refining downwards.\n");
			negative_int.set_upper(current_value);
		} else {
			// If the negative interval is the same as the positive interval, update it too
			if (equal(negative_int,positive_int)) positive_int.set_lower(current_value);

			ARIADNE_LOG(1,"False, refining upwards.\n");
			negative_int.set_lower(current_value);
		}
	} else {
		if (positiveOnBottom) {
			// If the negative interval is the same as the positive interval, update it too
			if (equal(negative_int,positive_int)) positive_int.set_upper(current_value);

			ARIADNE_LOG(1,"Indeterminate, refining upwards.\n");
			negative_int.set_lower(current_value);
		} else {
			// If the negative interval is the same as the positive interval, update it too
			if (equal(negative_int,positive_int)) positive_int.set_lower(current_value);

			ARIADNE_LOG(1,"Indeterminate, refining downwards.\n");
			negative_int.set_upper(current_value);
		}
	}
}


std::pair<Interval,Interval>
HybridReachabilityAnalyser::_pos_neg_bounds_from_search_intervals(const Interval& positive_int,
																	  const Interval& negative_int,
																	  const Interval& parameter_range,
																	  const bool& positiveOnBottom) const
{
	Interval positive_result, negative_result;

	if (positiveOnBottom) {
		positive_result = (positive_int.empty() ? positive_int : Interval(parameter_range.lower(),positive_int.lower()));
		negative_result = (negative_int.empty() ? negative_int : Interval(negative_int.upper(),parameter_range.upper())); }
	else {
		positive_result = (positive_int.empty() ? positive_int : Interval(positive_int.upper(),parameter_range.upper()));
		negative_result = (negative_int.empty() ? negative_int : Interval(parameter_range.lower(),negative_int.lower())); }

	return make_pair<Interval,Interval>(positive_result,negative_result);
}

Parametric2DBisectionResults
HybridReachabilityAnalyser::parametric_verification_2d_bisection(SystemVerificationInfo& verInfo,
											   const RealConstantSet& params,
											   const Float& tolerance,
											   const unsigned& numPointsPerAxis)
{
	ARIADNE_ASSERT_MSG(params.size() == 2,"Provide exactly two parameters.");

	RealConstantSet::const_iterator param_it = params.begin();

	const RealConstant& xParam = *param_it;
	const RealConstant& yParam = *(++param_it);

	return parametric_verification_2d_bisection(verInfo,xParam,yParam,tolerance,numPointsPerAxis);
}

Parametric2DBisectionResults
HybridReachabilityAnalyser::parametric_verification_2d_bisection(SystemVerificationInfo& verInfo,
											   const RealConstant& xParam,
											   const RealConstant& yParam,
											   const Float& tolerance,
											   const unsigned& numPointsPerAxis)
{
	// Generates the file name
	std::string filename = verInfo.getSystem().name();
	filename = filename + "[" + xParam.name() + "," + yParam.name() + "]";

	// Initializes the results
	Parametric2DBisectionResults results(filename,xParam.value(),yParam.value(),numPointsPerAxis);

	// Sweeps on each axis
	_parametric_verification_2d_bisection_sweep(results, verInfo, xParam, yParam, tolerance, numPointsPerAxis, true);
	_parametric_verification_2d_bisection_sweep(results, verInfo, xParam, yParam, tolerance, numPointsPerAxis, false);

	return results;
}

void HybridReachabilityAnalyser::_parametric_verification_2d_bisection_sweep(Parametric2DBisectionResults& results,
					  	  	  	    							SystemVerificationInfo& verInfo,
					  	  	  	    							RealConstant xParam,
					  	  	  	    							RealConstant yParam,
					  	  	  	    							const Float& tolerance,
					  	  	  	    							const uint& numPointsPerAxis,
					  	  	  	    							bool sweepOnX)
{
	HybridAutomaton& system = verInfo.getSystem();

	RealConstant& sweepParam = (sweepOnX ? xParam : yParam);
	RealConstant& otherParam = (sweepOnX ? yParam : xParam);
	Real originalSweepValue = system.accessible_constant_value(sweepParam.name());
	Interval sweepBounds = sweepParam.value();

	ARIADNE_LOG(1,"\nSweeping on " << sweepParam.name() << " in " << sweepBounds << ":\n");

	// Sweeps in the given direction
	for (uint i=0; i<numPointsPerAxis; ++i)
	{
		// Changes the value
		sweepParam.set_value(sweepBounds.lower()+i*sweepBounds.width()/(numPointsPerAxis-1));
		// Modifies the system with the new value
		system.substitute(sweepParam);

		ARIADNE_LOG(1,"\nAnalysis with " << sweepParam.name() << " = " << sweepParam.value().lower() << "...\n");

		// Performs the analysis on the other axis and adds to the results of the sweep variable
		std::pair<Interval,Interval> result = parametric_verification_1d_bisection(verInfo,otherParam,tolerance);
		if (sweepOnX)
			results.insertXValue(result);
		else
			results.insertYValue(result);

		ARIADNE_LOG(1,"Obtained safety in " << result.first << " and unsafety in " << result.second << ".\n");
	}

	// Restores the original value
	system.substitute(sweepParam,originalSweepValue);
}

ParametricPartitioningOutcomeList
HybridReachabilityAnalyser::parametric_verification_partitioning(SystemVerificationInfo& verInfo,
											  const RealConstantSet& params,
											  const Float& minPartitioningRatio)
{
	ARIADNE_ASSERT_MSG(params.size() > 0, "Provide at least one parameter.");
	ARIADNE_ASSERT(minPartitioningRatio > 0 && minPartitioningRatio < 1);

	ParametricPartitioningOutcomeList result(verInfo.getSystem().name(),params);

	RealConstantSet original_constants = verInfo.getSystem().accessible_constants();

	std::list<RealConstantSet> splittings = maximally_split_parameters_onTolerance(params,minPartitioningRatio);
	uint i=0;
	for (std::list<RealConstantSet>::const_iterator splitting_it = splittings.begin();
													 splitting_it != splittings.end();
													 ++splitting_it)
	{
		ARIADNE_LOG(1,"<Split parameters set #" << ++i << "/" << splittings.size() << ">\n");
		RealConstantSet current_params = *splitting_it;
		ARIADNE_LOG(1,"Parameter values: " << current_params << " ");
		verInfo.getSystem().substitute(current_params);
		tribool outcome = _verify_iterative(verInfo,params);
		ARIADNE_LOG(1,"Outcome: " << outcome << "\n");
		result.push_back(ParametricPartitioningOutcome(current_params,outcome));
	}

	verInfo.getSystem().substitute(original_constants);

	return result;
}

tribool
HybridReachabilityAnalyser::dominance(SystemVerificationInfo& dominating,
		  	  	  	  	  	  	  	  SystemVerificationInfo& dominated)
{
	const RealConstantSet dominatingLockedConstants;
	return _dominance(dominating,dominated,dominatingLockedConstants);
}

tribool
HybridReachabilityAnalyser::dominance(SystemVerificationInfo& dominating,
		  	  	  	  	  	  	  	  SystemVerificationInfo& dominated,
		  	  	  	  	  	  	  	  const RealConstant& parameter)
{
	HybridAutomaton& system = dominating.getSystem();

	Real original_value = system.accessible_constant_value(parameter.name());

	system.substitute(parameter);
	tribool result = dominance(dominating,dominated);
	system.substitute(parameter,original_value);

	return result;
}

tribool
HybridReachabilityAnalyser::dominance(SystemVerificationInfo& dominating,
		  	  	  	  	  	  	  	  SystemVerificationInfo& dominated,
		  	  	  	  	  	  	  	  const RealConstant& parameter,
		  	  	  	  	  	  	  	  const Float& value)
{
	const RealConstant modifiedParameter(parameter.name(),Interval(value));

	return dominance(dominating,dominated,modifiedParameter);
}

tribool
HybridReachabilityAnalyser::_dominance(SystemVerificationInfo& dominating,
									  SystemVerificationInfo& dominated,
									  const RealConstantSet& dominatingLockedConstants)
{
	ARIADNE_ASSERT(dominating.getProjection().size() == dominated.getProjection().size());

	ARIADNE_LOG(1, "Dominance checking...\n");

	// We are not allowed to skip if disproved, since we need as much reached region as possible
	// We are however allowed to skip if unprovable, since we could not determine dominance anyway
	_parameters->skip_if_disproved = false;
	_parameters->skip_if_unprovable = true;

	int& depth = _parameters->maximum_grid_depth;
    for(depth = _parameters->lowest_maximum_grid_depth; depth <= _parameters->highest_maximum_grid_depth; ++depth)
	{
		ARIADNE_LOG(2, "DEPTH " << depth << "\n");

		if (_dominance_positive(dominating, dominated, dominatingLockedConstants)) {
			ARIADNE_LOG(3, "Dominates.\n");
			return true;
		}

		if (_dominance_negative(dominating, dominated, dominatingLockedConstants)) {
			ARIADNE_LOG(3, "Does not dominate.\n");
			return false;
		}
    }

	// Return indeterminate otherwise
	ARIADNE_LOG(3, "Indeterminate.\n");
	return indeterminate;
}

bool
HybridReachabilityAnalyser::_dominance_positive(SystemVerificationInfo& dominating,
		  	  	  	  	  	  	  	  	  	  	SystemVerificationInfo& dominated,
		  	  	  	  	  	  	  	  	  	  	const RealConstantSet& dominatingLockedConstants)
{
	ARIADNE_LOG(3,"Looking for a positive answer...\n");

	bool result;

	HybridGridTreeSet dominating_reach, dominated_reach;
	bool isValid;
	DisproveData disproveData(dominated.getSystem().state_space());

	ARIADNE_LOG(3,"Getting the outer approximation of the dominating system...\n");

	_setDominanceParameters(dominating,dominatingLockedConstants);
	make_lpair<HybridGridTreeSet,bool>(dominating_reach,isValid) = upper_chain_reach_forward(dominating.getSystem(),dominating.getInitialSet());

	if (!isValid) {
		ARIADNE_LOG(3,"The reached region is unbounded.\n");
		result = false;
	}
	else
	{
		ARIADNE_LOG(3,"Getting the lower approximation of the dominated system...\n");

		RealConstantSet emptyLockedConstants;
		_setDominanceParameters(dominated,emptyLockedConstants);
		make_lpair<HybridGridTreeSet,DisproveData>(dominated_reach,disproveData) = lower_chain_reach(dominated.getSystem(),dominated.getInitialSet());

		// We must shrink the lower approximation of the the dominated system, but underapproximating in terms of rounding
		HybridBoxes shrinked_dominated_bounds = Ariadne::shrink_in(disproveData.getReachBounds(),disproveData.getEpsilon());

		Box projected_dominating_bounds = Ariadne::project(dominating_reach.bounding_box(),dominating.getProjection());
		Box projected_shrinked_dominated_bounds = Ariadne::project(shrinked_dominated_bounds,dominated.getProjection());

		ARIADNE_LOG(4,"Epsilon: " << disproveData.getEpsilon() << "\n");
		ARIADNE_LOG(4,"Projected dominating bounds: " << projected_dominating_bounds << "\n");
		ARIADNE_LOG(4,"Projected shrinked dominated bounds: " << projected_shrinked_dominated_bounds << "\n");

		result = inside(projected_dominating_bounds,projected_shrinked_dominated_bounds);
	}

	return result;
}

bool
HybridReachabilityAnalyser::_dominance_negative(SystemVerificationInfo& dominating,
	  	  	  	  								SystemVerificationInfo& dominated,
	  	  	  	  								const RealConstantSet& dominatingLockedConstants)
{
	ARIADNE_LOG(3,"Looking for a negative answer...\n");

	bool result;

	HybridGridTreeSet dominating_reach, dominated_reach;
	bool isValid;
	DisproveData disproveData(dominating.getSystem().state_space());

	ARIADNE_LOG(3,"Getting the outer approximation of the dominated system...\n");

	RealConstantSet emptyLockedConstants;
	_setDominanceParameters(dominated,emptyLockedConstants);
	make_lpair<HybridGridTreeSet,bool>(dominated_reach,isValid) = upper_chain_reach_forward(dominated.getSystem(),dominated.getInitialSet());

	if (!isValid) {
		ARIADNE_LOG(3,"The reached region is unbounded.\n");
		result = false;
	}
	else
	{
		ARIADNE_LOG(3,"Getting the lower approximation of the dominating system...\n");

		_setDominanceParameters(dominating,dominatingLockedConstants);

		make_lpair<HybridGridTreeSet,DisproveData>(dominating_reach,disproveData) = lower_chain_reach(dominating.getSystem(),dominating.getInitialSet());

		// We must shrink the lower approximation of the the dominating system, but overapproximating in terms of rounding
		HybridBoxes shrinked_dominating_bounds = Ariadne::shrink_out(disproveData.getReachBounds(),disproveData.getEpsilon());

		Box projected_shrinked_dominating_bounds = Ariadne::project(shrinked_dominating_bounds,dominating.getProjection());
		Box projected_dominated_bounds = Ariadne::project(dominated_reach.bounding_box(),dominated.getProjection());

		ARIADNE_LOG(4,"Epsilon: " << disproveData.getEpsilon() << "\n");
		ARIADNE_LOG(4,"Projected shrinked dominating bounds: " << projected_shrinked_dominating_bounds << "\n");
		ARIADNE_LOG(4,"Projected dominated bounds: " << projected_dominated_bounds << "\n");

		result = !inside(projected_shrinked_dominating_bounds,projected_dominated_bounds);
	}

	return result;
}

ParametricPartitioningOutcomeList
HybridReachabilityAnalyser::parametric_dominance_partitioning(SystemVerificationInfo& dominating,
												 SystemVerificationInfo& dominated,
											  	 const RealConstantSet& dominating_params,
											  	 const Float& minPartitioningRatio)
{
	ARIADNE_ASSERT_MSG(dominating_params.size() > 0, "Provide at least one parameter.");
	ARIADNE_ASSERT(minPartitioningRatio > 0 && minPartitioningRatio < 1);

	std::string name = dominating.getSystem().name() + "&" + dominated.getSystem().name();
	ParametricPartitioningOutcomeList result(name,dominating_params);

	RealConstantSet original_constants = dominating.getSystem().accessible_constants();

	std::list<RealConstantSet> splittings = maximally_split_parameters_onTolerance(dominating_params,minPartitioningRatio);
	uint i=0;
	for (std::list<RealConstantSet>::const_iterator splitting_it = splittings.begin();
													 splitting_it != splittings.end();
													 ++splitting_it)
	{
		ARIADNE_LOG(1,"<Split parameters set #" << ++i << "/" << splittings.size() << ">\n");
		RealConstantSet current_params = *splitting_it;
		ARIADNE_LOG(1,"Parameter values: " << current_params << " ");
		dominating.getSystem().substitute(current_params);
		tribool outcome = _dominance(dominating,dominated,dominating_params);
		ARIADNE_LOG(1,"Outcome: " << outcome << "\n");
		result.push_back(ParametricPartitioningOutcome(current_params,outcome));
	}

	dominating.getSystem().substitute(original_constants);

	return result;
}

std::list<RealConstantSet>
maximally_split_parameters_onTolerance(const RealConstantSet& params,
									   const Float& tolerance)
{
	std::list<RealConstantSet> result;

	std::list<RealConstantSet> working;
	working.push_back(params);

	while(!working.empty()) {
		RealConstantSet currentParams = working.back();
		working.pop_back();

		bool hasSplit = false;
		for (RealConstantSet::const_iterator param_it = currentParams.begin();
											 param_it != currentParams.end();
											 ++param_it)
		{
			Float currentWidth = param_it->value().width();
			Float initialWidth = params.find(*param_it)->value().width();

			if (currentWidth/initialWidth > tolerance)
			{
				hasSplit = true;
				currentParams.erase(*param_it);
				RealConstantSet newConfigurationLeft = currentParams;
				RealConstantSet newConfigurationRight = currentParams;

				const Real& currentInterval = param_it->value();
				newConfigurationLeft.insert(RealConstant(param_it->name(),
													  Interval(currentInterval.lower(),currentInterval.midpoint())));
				newConfigurationRight.insert(RealConstant(param_it->name(),
													  Interval(currentInterval.midpoint(),currentInterval.upper())));

				working.push_back(newConfigurationLeft);
				working.push_back(newConfigurationRight);
				break;
			}
		}

		if (!hasSplit)
			result.push_back(currentParams);
	}

	return result;
}

} // namespace Ariadne
