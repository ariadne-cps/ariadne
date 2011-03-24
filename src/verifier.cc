/***************************************************************************
 *            verifier.cc
 *
 *  Copyright 2011  Luca Geretti
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

#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <string>

#include "verifier.h"

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
							 _domain(domain), _safe_region(unbounded_hybrid_boxes(system.state_space())),
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

Verifier::Verifier(const HybridReachabilityAnalyser& outer_analyser,
				   const HybridReachabilityAnalyser& lower_analyser) :
				   _outer_analyser(outer_analyser.clone()),
				   _lower_analyser(lower_analyser.clone()),
				   plot_results(false),
				   maximum_parameter_depth(3)
{
}

Verifier::Verifier(const HybridReachabilityAnalyser& analyser) :
				   _outer_analyser(analyser.clone()),
				   _lower_analyser(analyser.clone()),
				   plot_results(false),
				   maximum_parameter_depth(3)
{
}

Verifier::~Verifier()
{
}


bool
Verifier::
_safety_positive_once(SystemType& system,
	  const HybridImageSet& initial_set,
	  const HybridBoxes& safe_region) const
{
	bool result;

	ARIADNE_LOG(4,"Proving...\n");

	DiscreteEvolutionParameters& params = _outer_analyser->parameters();

	if (!_is_grid_depth_within_bounds(UPPER_SEMANTICS))
		return false;

	_tuneIterativeStepParameters(system,_outer_analyser->statistics().upper().reach,UPPER_SEMANTICS);

	HybridGridTreeSet reach;

	try
	{
		reach = _outer_analyser->outer_chain_reach_quick_proving(system,initial_set,safe_region);

		if (params.domain_enforcing_policy == OFFLINE && definitely(!reach.subset(params.bounding_domain)))
			result = false;
		else
			result = definitely(reach.subset(safe_region));
	}
	catch (ReachOutOfDomainException ex)
	{
		ARIADNE_LOG(5, "The reached set could not be bounded (" << ex.what() << ").\n");
		result = false;
	}

	_outer_analyser->statistics().upper().reach = reach;

	if (plot_results)
		_plot(reach,UPPER_SEMANTICS);

	ARIADNE_LOG(4, (result ? "Proved.\n" : "Not proved.\n") );

	return result;
}


bool
Verifier::
_safety_negative_once(SystemType& system,
		 const HybridImageSet& initial_set,
		 const HybridBoxes& safe_region) const
{
	ARIADNE_LOG(4,"Disproving...\n");

	if (!_is_grid_depth_within_bounds(LOWER_SEMANTICS))
		return false;

	_tuneIterativeStepParameters(system,_outer_analyser->statistics().upper().reach,LOWER_SEMANTICS);

	std::pair<HybridGridTreeSet,DisproveData> reachAndDisproveData =
			_lower_analyser->lower_chain_reach(system,initial_set,safe_region);
	const HybridGridTreeSet& reach = reachAndDisproveData.first;
	const DisproveData& disproveData = reachAndDisproveData.second;
	const bool& isDisproved = disproveData.getIsDisproved();

	if (plot_results)
		_plot(reach,LOWER_SEMANTICS);

	ARIADNE_LOG(5,"Disprove data: " << disproveData << "\n");

	ARIADNE_LOG(4, (isDisproved ? "Disproved.\n" : "Not disproved.\n") );

	return isDisproved;
}


tribool
Verifier::
_safety_once(SystemType& system,
			 const HybridImageSet& initial_set,
			 const HybridBoxes& safe_region) const
{
		ARIADNE_LOG(3, "Verification...\n");

		if (_safety_positive_once(system,initial_set,safe_region)) {
			ARIADNE_LOG(3, "Safe.\n");
			return true;
		}

		if (_safety_negative_once(system,initial_set,safe_region)) {
			ARIADNE_LOG(3, "Unsafe.\n");
			return false;
		}

		ARIADNE_LOG(3, "Indeterminate.\n");
		return indeterminate;
}

tribool
Verifier::
safety(SystemVerificationInfo& verInfo) const
{
	RealConstantSet locked_constants;

	return _safety_nosplitting(verInfo,locked_constants);
}


tribool
Verifier::
_safety(SystemVerificationInfo& verInfo, const RealConstant& parameter) const
{
	HybridAutomaton& system = verInfo.getSystem();

	Real originalParameterValue = system.accessible_constant_value(parameter.name());

	system.substitute(parameter);
	tribool result = safety(verInfo);
	system.substitute(parameter,originalParameterValue);

	return result;
}

tribool
Verifier::
_safety(SystemVerificationInfo& verInfo, const RealConstant& parameter, const Float& value) const
{
	const RealConstant modifiedParameter(parameter.name(),Interval(value));

	return _safety(verInfo, modifiedParameter);
}


tribool
Verifier::
_safety_nosplitting(SystemVerificationInfo& verInfo,
				    const RealConstantSet& locked_constants) const
{
	ARIADNE_LOG(2,"\nIterative verification...\n");

	HybridAutomaton& system = verInfo.getSystem();

	if (plot_results)
		_plot_dirpath_init(system);

	_outer_analyser->resetStatistics();

	// Set the initial parameters
	_setInitialParameters(system,verInfo.getDomain(),verInfo.getSafeRegion(),locked_constants,UPPER_SEMANTICS);
	_setInitialParameters(system,verInfo.getDomain(),verInfo.getSafeRegion(),locked_constants,LOWER_SEMANTICS);

	int initial_depth = min(_outer_analyser->parameters().lowest_maximum_grid_depth,
							_lower_analyser->parameters().lowest_maximum_grid_depth);
	int final_depth = max(_outer_analyser->parameters().highest_maximum_grid_depth,
						  _lower_analyser->parameters().highest_maximum_grid_depth);
    for (int depth = initial_depth; depth <= final_depth; ++depth)
	{
    	_outer_analyser->parameters().maximum_grid_depth = depth;
    	_lower_analyser->parameters().maximum_grid_depth = depth;

		ARIADNE_LOG(2, "DEPTH " << depth << "\n");

		_tuneIterativeStepParameters(system,_outer_analyser->statistics().upper().reach,UPPER_SEMANTICS);
		_tuneIterativeStepParameters(system,_outer_analyser->statistics().upper().reach,LOWER_SEMANTICS);

		tribool result = _safety_once(system,verInfo.getInitialSet(),verInfo.getSafeRegion());

		if (!indeterminate(result))
			return result;
    }

	// Return indeterminate
	return indeterminate;
}

bool
Verifier::
_is_grid_depth_within_bounds(Semantics semantics) const
{
	const DiscreteEvolutionParameters& parameters = (semantics == UPPER_SEMANTICS ? _outer_analyser->parameters() : _lower_analyser->parameters());

	if (parameters.maximum_grid_depth < parameters.lowest_maximum_grid_depth) {
		ARIADNE_LOG(4,"Skipped verification since the depth is lower than the lowest allowed.\n");
		return false;
	}
	if (parameters.maximum_grid_depth > parameters.highest_maximum_grid_depth) {
		ARIADNE_LOG(4,"Skipped verification since the depth is higher than the highest allowed.\n");
		return false;
	}

	return true;
}

std::pair<Interval,Interval>
Verifier::
parametric_safety_1d_bisection(SystemVerificationInfo& verInfo,
						   	   	   	 const RealConstant& parameter) const
{
	float tolerance = 1.0/(1 << this->maximum_parameter_depth);

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
	tribool lower_result = _safety(verInfo,parameter,parameter_range.lower());
	ARIADNE_LOG(1,pretty_print(lower_result) << ".\n");

	// Check the upper bound
	ARIADNE_LOG(1,"Checking upper interval bound... ");
	tribool upper_result = _safety(verInfo,parameter,parameter_range.upper());
	ARIADNE_LOG(1,pretty_print(upper_result) << ".\n");

	// If we must proceed with bisection refining
	bool proceed;
	// If the safe value is found on the lower extreme of the parameter interval
	bool safeOnBottom;

	// Updates the initial values for the safety/unsafety intervals, and decides whether to proceed and the eventual ordering
	make_lpair<bool,bool>(proceed,safeOnBottom) = process_initial_bisection_results(safety_int,unsafety_int,parameter_range,
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
			tribool result = _safety(verInfo,parameter,current_value);

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
			tribool result = _safety(verInfo,parameter,current_value);

			_process_negative_bisection_result(result,safety_int,unsafety_int,current_value,safeOnBottom);
		}
		else
			proceedUnsafety = false;
	}

	system.substitute(parameter,original_value);

	return pos_neg_bounds_from_search_intervals(safety_int,unsafety_int,parameter_range,safeOnBottom);
}

Parametric2DBisectionResults
Verifier::
parametric_safety_2d_bisection(SystemVerificationInfo& verInfo,
									 const RealConstantSet& params) const
{
	ARIADNE_ASSERT_MSG(params.size() == 2,"Provide exactly two parameters.");

	RealConstantSet::const_iterator param_it = params.begin();

	const RealConstant& xParam = *param_it;
	const RealConstant& yParam = *(++param_it);

	return parametric_safety_2d_bisection(verInfo,xParam,yParam);
}

Parametric2DBisectionResults
Verifier::
parametric_safety_2d_bisection(SystemVerificationInfo& verInfo,
									 const RealConstant& xParam,
									 const RealConstant& yParam) const
{
	// Generates the file name
	std::string filename = verInfo.getSystem().name();
	filename = filename + "[" + xParam.name() + "," + yParam.name() + "]";

	// Initializes the results
	uint numPointsPerAxis = 1+(1<<this->maximum_parameter_depth);
	Parametric2DBisectionResults results(filename,xParam.value(),yParam.value(),numPointsPerAxis);

	// Sweeps on each axis
	_parametric_safety_2d_bisection_sweep(results, verInfo, xParam, yParam, true);
	_parametric_safety_2d_bisection_sweep(results, verInfo, xParam, yParam, false);

	return results;
}

std::list<ParametricOutcome>
Verifier::
parametric_safety(SystemVerificationInfo& verInfo,
									 const RealConstantSet& params) const
{
	ARIADNE_ASSERT_MSG(params.size() > 0, "Provide at least one parameter.");

	std::list<ParametricOutcome> result;

	RealConstantSet original_constants = verInfo.getSystem().accessible_constants();

	std::list<RealConstantSet> splittings = maximally_split_parameters(params,this->maximum_parameter_depth);
	uint i=0;
	for (std::list<RealConstantSet>::const_iterator splitting_it = splittings.begin();
													 splitting_it != splittings.end();
													 ++splitting_it)
	{
		ARIADNE_LOG(1,"<Split parameters set #" << ++i << "/" << splittings.size() << ">\n");
		RealConstantSet current_params = *splitting_it;
		ARIADNE_LOG(1,"Parameter values: " << current_params << " ");
		verInfo.getSystem().substitute(current_params);
		tribool outcome = _safety_nosplitting(verInfo,params);
		ARIADNE_LOG(1,"Outcome: " << pretty_print(outcome) << "\n");
		result.push_back(ParametricOutcome(current_params,outcome));
	}

	verInfo.getSystem().substitute(original_constants);

	return result;
}

tribool
Verifier::
dominance(SystemVerificationInfo& dominating,
		  SystemVerificationInfo& dominated) const
{
	const RealConstantSet dominatingLockedConstants;
	return _dominance(dominating,dominated,dominatingLockedConstants);
}

tribool
Verifier::
_dominance(SystemVerificationInfo& dominating,
		  	SystemVerificationInfo& dominated,
		  	const RealConstant& parameter) const
{
	HybridAutomaton& system = dominating.getSystem();

	Real original_value = system.accessible_constant_value(parameter.name());

	system.substitute(parameter);
	tribool result = dominance(dominating,dominated);
	system.substitute(parameter,original_value);

	return result;
}

tribool
Verifier::
_dominance(SystemVerificationInfo& dominating,
		  SystemVerificationInfo& dominated,
		  const RealConstant& parameter,
		  const Float& value) const
{
	const RealConstant modifiedParameter(parameter.name(),Interval(value));

	return _dominance(dominating,dominated,modifiedParameter);
}

std::pair<Interval,Interval>
Verifier::
parametric_dominance_1d_bisection(SystemVerificationInfo& dominating,
								  SystemVerificationInfo& dominated,
						   	   	  const RealConstant& parameter) const
{
	float tolerance = 1.0/(1 << this->maximum_parameter_depth);

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
	tribool lower_result = _dominance(dominating,dominated,parameter,parameter_range.lower());
	ARIADNE_LOG(1,pretty_print(lower_result) << ".\n");

	// Check the upper bound
	ARIADNE_LOG(1,"Checking upper interval bound... ");
	tribool upper_result = _dominance(dominating,dominated,parameter,parameter_range.upper());
	ARIADNE_LOG(1,pretty_print(upper_result) << ".\n");

	// If we must proceed with bisection refining
	bool proceed;
	// If the dominating value is found on the lower extreme of the parameter interval
	bool dominatingOnBottom;

	// Updates the initial values for the dominating/nondominating intervals, and decides whether to proceed and the eventual ordering
	make_lpair<bool,bool>(proceed,dominatingOnBottom) = process_initial_bisection_results(dominating_int,nondominating_int,parameter_range,
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
			tribool result = _dominance(dominating,dominated,parameter,current_value);

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
			tribool result = _dominance(dominating,dominated,parameter,current_value);

			_process_negative_bisection_result(result,dominating_int,nondominating_int,current_value,dominatingOnBottom);
		}
		else
			proceedNondominating = false;
	}

	system.substitute(parameter,original_value);

	return pos_neg_bounds_from_search_intervals(dominating_int,nondominating_int,parameter_range,dominatingOnBottom);
}

Parametric2DBisectionResults
Verifier::
parametric_dominance_2d_bisection(SystemVerificationInfo& dominating,
								  SystemVerificationInfo& dominated,
								  const RealConstantSet& params) const
{
	ARIADNE_ASSERT_MSG(params.size() == 2,"Provide exactly two parameters.");

	RealConstantSet::const_iterator param_it = params.begin();

	const RealConstant& xParam = *param_it;
	const RealConstant& yParam = *(++param_it);

	return parametric_dominance_2d_bisection(dominating,dominated,xParam,yParam);
}

Parametric2DBisectionResults
Verifier::
parametric_dominance_2d_bisection(SystemVerificationInfo& dominating,
								  SystemVerificationInfo& dominated,
								  const RealConstant& xParam,
								  const RealConstant& yParam) const
{
	// Generates the file name
	std::string filename = dominating.getSystem().name() + "&" + dominated.getSystem().name();
	filename = filename + "[" + xParam.name() + "," + yParam.name() + "]";

	// Initializes the results
	uint numPointsPerAxis = (1 << this->maximum_parameter_depth);
	Parametric2DBisectionResults results(filename,xParam.value(),yParam.value(),numPointsPerAxis);

	// Sweeps on each axis
	_parametric_dominance_2d_bisection_sweep(results, dominating, dominated, xParam, yParam, true);
	_parametric_dominance_2d_bisection_sweep(results, dominating, dominated, xParam, yParam, false);

	return results;
}


std::list<ParametricOutcome>
Verifier::
parametric_dominance(SystemVerificationInfo& dominating,
								  SystemVerificationInfo& dominated,
								  const RealConstantSet& dominating_params) const
{
	ARIADNE_ASSERT_MSG(dominating_params.size() > 0, "Provide at least one parameter.");

	std::list<ParametricOutcome> result;

	RealConstantSet original_constants = dominating.getSystem().accessible_constants();

	std::list<RealConstantSet> splittings = maximally_split_parameters(dominating_params,this->maximum_parameter_depth);
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
		ARIADNE_LOG(1,"Outcome: " << pretty_print(outcome) << "\n");
		result.push_back(ParametricOutcome(current_params,outcome));
	}

	dominating.getSystem().substitute(original_constants);

	return result;
}

void
Verifier::
_process_positive_bisection_result(const tribool& result,
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
Verifier::_process_negative_bisection_result(const tribool& result,
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


void Verifier::_parametric_safety_2d_bisection_sweep(Parametric2DBisectionResults& results,
					  	  	  	    					   SystemVerificationInfo& verInfo,
					  	  	  	    					   RealConstant xParam,
					  	  	  	    					   RealConstant yParam,
					  	  	  	    					   bool sweepOnX) const
{
	HybridAutomaton& system = verInfo.getSystem();

	uint numPointsPerAxis = 1+(1 << this->maximum_parameter_depth);

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
		std::pair<Interval,Interval> result = parametric_safety_1d_bisection(verInfo,otherParam);
		if (sweepOnX)
			results.insertXValue(result);
		else
			results.insertYValue(result);

		ARIADNE_LOG(1,"Obtained safety in " << result.first << " and unsafety in " << result.second << ".\n");
	}

	// Restores the original value
	system.substitute(sweepParam,originalSweepValue);
}

void
Verifier::
_parametric_dominance_2d_bisection_sweep(Parametric2DBisectionResults& results,
					  	  	  	    	 SystemVerificationInfo& dominating,
					  	  	  	    	 SystemVerificationInfo& dominated,
					  	  	  	    	 RealConstant xParam,
					  	  	  	    	 RealConstant yParam,
					  	  	  	    	 bool sweepOnX) const
{
	HybridAutomaton& system = dominating.getSystem();

	uint numPointsPerAxis = (1 << this->maximum_parameter_depth);

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
		std::pair<Interval,Interval> result = parametric_dominance_1d_bisection(dominating,dominated,otherParam);
		if (sweepOnX)
			results.insertXValue(result);
		else
			results.insertYValue(result);

		ARIADNE_LOG(1,"Obtained dominance in " << result.first << " and nondominance in " << result.second << ".\n");
	}

	// Restores the original value
	system.substitute(sweepParam,originalSweepValue);
}

tribool
Verifier::_dominance(SystemVerificationInfo& dominating,
					 SystemVerificationInfo& dominated,
					 const RealConstantSet& dominatingLockedConstants) const
{
	ARIADNE_ASSERT(dominating.getProjection().size() == dominated.getProjection().size());

	ARIADNE_LOG(1, "Dominance checking...\n");

	// We are not allowed to skip as soon as disproved, since we need as much reached region as possible
	// We are however allowed to perform quick proving, since we could not determine dominance anyway
	_lower_analyser->parameters().enable_quick_disproving = false;
	_outer_analyser->parameters().enable_quick_proving = true;

	int initial_depth = max(_outer_analyser->parameters().lowest_maximum_grid_depth,
							_lower_analyser->parameters().lowest_maximum_grid_depth);
	int final_depth = min(_outer_analyser->parameters().highest_maximum_grid_depth,
						  _lower_analyser->parameters().highest_maximum_grid_depth);
    for (int depth = initial_depth; depth <= final_depth; ++depth)
	{
    	_outer_analyser->parameters().maximum_grid_depth = depth;
    	_lower_analyser->parameters().maximum_grid_depth = depth;

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
Verifier::_dominance_positive(SystemVerificationInfo& dominating,
		  	  	  	  	  	  SystemVerificationInfo& dominated,
		  	  	  	  	  	  const RealConstantSet& dominatingLockedConstants) const
{
	ARIADNE_LOG(3,"Looking for a positive answer...\n");

	bool result;

	try {
		ARIADNE_LOG(3,"Getting the outer approximation of the dominating system...\n");

		_setDominanceParameters(dominating,dominatingLockedConstants,UPPER_SEMANTICS);

		HybridGridTreeSet dominating_reach = _outer_analyser->outer_chain_reach(dominating.getSystem(),dominating.getInitialSet());
		Box projected_dominating_bounds = Ariadne::project(dominating_reach.bounding_box(),dominating.getProjection());

		ARIADNE_LOG(4,"Projected dominating bounds: " << projected_dominating_bounds << "\n");

		ARIADNE_LOG(3,"Getting the lower approximation of the dominated system...\n");

		HybridGridTreeSet dominated_reach;
		RealConstantSet emptyLockedConstants;
		DisproveData disproveData(dominated.getSystem().state_space());
		_setDominanceParameters(dominated,emptyLockedConstants,LOWER_SEMANTICS);
		make_lpair<HybridGridTreeSet,DisproveData>(dominated_reach,disproveData) =
				_lower_analyser->lower_chain_reach(dominated.getSystem(),dominated.getInitialSet());

		// We must shrink the lower approximation of the the dominated system, but underapproximating in terms of rounding
		HybridBoxes shrinked_dominated_bounds = Ariadne::shrink_in(disproveData.getReachBounds(),disproveData.getEpsilon());

		Box projected_shrinked_dominated_bounds = Ariadne::project(shrinked_dominated_bounds,dominated.getProjection());

		ARIADNE_LOG(4,"Epsilon: " << disproveData.getEpsilon() << "\n");

		ARIADNE_LOG(4,"Projected shrinked dominated bounds: " << projected_shrinked_dominated_bounds << "\n");

		result = inside(projected_dominating_bounds,projected_shrinked_dominated_bounds);
	}
	catch (ReachOutOfDomainException ex) {
		ARIADNE_LOG(3,"The outer reached region of the dominating system is unbounded.\n");
		result = false;
	}

	return result;
}

bool
Verifier::_dominance_negative(SystemVerificationInfo& dominating,
	  	  	  	  			  SystemVerificationInfo& dominated,
	  	  	  	  			  const RealConstantSet& dominatingLockedConstants) const
{
	ARIADNE_LOG(3,"Looking for a negative answer...\n");

	bool result;

	try {
		ARIADNE_LOG(3,"Getting the outer approximation of the dominated system...\n");

		RealConstantSet emptyLockedConstants;
		_setDominanceParameters(dominated,emptyLockedConstants,UPPER_SEMANTICS);
		HybridGridTreeSet dominated_reach = _outer_analyser->outer_chain_reach(dominated.getSystem(),dominated.getInitialSet());
		Box projected_dominated_bounds = Ariadne::project(dominated_reach.bounding_box(),dominated.getProjection());

		ARIADNE_LOG(4,"Projected dominated bounds: " << projected_dominated_bounds << "\n");

		ARIADNE_LOG(3,"Getting the lower approximation of the dominating system...\n");

		HybridGridTreeSet dominating_reach;
		DisproveData disproveData(dominating.getSystem().state_space());
		_setDominanceParameters(dominating,dominatingLockedConstants,LOWER_SEMANTICS);
		make_lpair<HybridGridTreeSet,DisproveData>(dominating_reach,disproveData) =
				_lower_analyser->lower_chain_reach(dominating.getSystem(),dominating.getInitialSet());

		// We must shrink the lower approximation of the the dominating system, but overapproximating in terms of rounding
		HybridBoxes shrinked_dominating_bounds = Ariadne::shrink_out(disproveData.getReachBounds(),disproveData.getEpsilon());

		Box projected_shrinked_dominating_bounds = Ariadne::project(shrinked_dominating_bounds,dominating.getProjection());


		ARIADNE_LOG(4,"Epsilon: " << disproveData.getEpsilon() << "\n");
		ARIADNE_LOG(4,"Projected shrinked dominating bounds: " << projected_shrinked_dominating_bounds << "\n");

		result = !inside(projected_shrinked_dominating_bounds,projected_dominated_bounds);
	}
	catch (ReachOutOfDomainException ex) {
		ARIADNE_LOG(3,"The outer reached region of the dominated system is unbounded.\n");
		result = false;
	}

	return result;
}

void
Verifier::
_setInitialParameters(HybridAutomaton& system,
					  const HybridBoxes& domain,
					  const HybridBoxes& safe_region,
					  const RealConstantSet& locked_constants,
					  Semantics semantics) const
{
	ARIADNE_LOG(3,"Setting initial parameters of the " << (semantics == UPPER_SEMANTICS ? "outer " : "lower ") << "analyser...\n");
	DiscreteEvolutionParameters& parameters = (semantics == UPPER_SEMANTICS ?
											   _outer_analyser->parameters() : _lower_analyser->parameters());
	parameters.bounding_domain = domain;
	ARIADNE_LOG(4, "Domain: " << domain << "\n");
	parameters.lock_to_grid_time = getLockToGridTime(system,domain);
	ARIADNE_LOG(4, "Lock to grid time: " << parameters.lock_to_grid_time << "\n");
	parameters.split_factors = getSplitFactorsOfConstants(system,locked_constants,0.1,domain);
	ARIADNE_LOG(4, "Split factors: " << parameters.split_factors << "\n");
}


void
Verifier::
_tuneIterativeStepParameters(HybridAutomaton& system,
							 const HybridGridTreeSet& bounding_reach,
							 Semantics semantics) const
{
	HybridReachabilityAnalyser& analyser = (semantics == UPPER_SEMANTICS ? *_outer_analyser : *_lower_analyser);
	HybridFloatVector hmad = getHybridMaximumAbsoluteDerivatives(system,bounding_reach,analyser.parameters().bounding_domain);
	ARIADNE_LOG(4, "Derivative bounds: " << hmad << "\n");
	system.set_grid(getHybridGrid(hmad,analyser.parameters().bounding_domain));
	ARIADNE_LOG(4, "Grid: " << system.grid() << "\n");
	analyser.tuneEvolverParameters(system,hmad,analyser.parameters().maximum_grid_depth,semantics);
}

void
Verifier::
_setDominanceParameters(SystemVerificationInfo& verInfo,
						const RealConstantSet& lockedConstants,
						Semantics semantics) const
{
	HybridReachabilityAnalyser& analyser = (semantics == UPPER_SEMANTICS ? *_outer_analyser : *_lower_analyser);

	_outer_analyser->resetStatistics();

	analyser.parameters().bounding_domain = verInfo.getDomain();
	ARIADNE_LOG(4, "Domain: " << analyser.parameters().bounding_domain << "\n");
	// TODO: exploit actual upper reach results for grid refinement
	_tuneIterativeStepParameters(verInfo.getSystem(),_outer_analyser->statistics().upper().reach,semantics);
	analyser.parameters().lock_to_grid_time = getLockToGridTime(verInfo.getSystem(),analyser.parameters().bounding_domain);
	ARIADNE_LOG(4, "Lock to grid time: " << analyser.parameters().lock_to_grid_time << "\n");
	analyser.parameters().split_factors = getSplitFactorsOfConstants(verInfo.getSystem(),lockedConstants,0.1,analyser.parameters().bounding_domain);
	ARIADNE_LOG(4, "Split factors: " << analyser.parameters().split_factors << "\n");
}

void
Verifier::
_plot_dirpath_init(const HybridAutomaton& system) const
{
	time_t mytime;
	time(&mytime);
	string foldername = system.name()+"-png";

	mkdir(foldername.c_str(),0777);
	string timestring = asctime(localtime(&mytime));
	timestring.erase(std::remove(timestring.begin(), timestring.end(), '\n'), timestring.end());
	foldername = foldername+"/"+timestring;
	mkdir(foldername.c_str(),0777);

	_plot_dirpath = foldername;
}

void
Verifier::
_plot(const HybridGridTreeSet& reach,
	  Semantics semantics) const
{
	int maximum_grid_depth = (semantics == UPPER_SEMANTICS ?
							   _outer_analyser->parameters().maximum_grid_depth : _lower_analyser->parameters().maximum_grid_depth);
	char mgd_char[10];
	sprintf(mgd_char,"%i",maximum_grid_depth);
	string filename = (semantics == UPPER_SEMANTICS ? "outer-" : "lower-");
	filename.append(mgd_char);
	plot(_plot_dirpath,filename,reach);
}

std::string
pretty_print(tribool value)
{
	if (definitely(value))
		return "True";
	if (!possibly(value))
		return "False";
	return "Indeterminate";
}

std::pair<bool,bool>
process_initial_bisection_results(Interval& positive_int,
								  Interval& negative_int,
								  const Interval& parameter_range,
		 						  const tribool& lower_result,
								  const tribool& upper_result)
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

std::pair<Interval,Interval>
pos_neg_bounds_from_search_intervals(const Interval& positive_int,
									 const Interval& negative_int,
									 const Interval& parameter_range,
									 bool positiveOnBottom)
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

std::list<RealConstantSet>
maximally_split_parameters(const RealConstantSet& params,
						   const uint& maximum_parameter_depth)
{
	std::list<RealConstantSet> source;
	std::list<RealConstantSet> destination;
	destination.push_back(params);

	for (RealConstantSet::const_iterator param_it = params.begin();
										 param_it != params.end();
										 ++param_it)
	{
		for (uint i=0; i<maximum_parameter_depth; i++) {
			source.clear();
			source.insert<std::list<RealConstantSet>::const_iterator>(source.begin(),destination.begin(),destination.end());
			destination.clear();

			while (!source.empty()) {
				RealConstantSet currentParams = source.back();
				source.pop_back();

				RealConstantSet newConfigurationLeft = currentParams;
				RealConstantSet newConfigurationRight = currentParams;
				newConfigurationLeft.erase(*param_it);
				newConfigurationRight.erase(*param_it);

				const Real& currentInterval = currentParams.find(*param_it)->value();
				newConfigurationLeft.insert(RealConstant(param_it->name(),
													  Interval(currentInterval.lower(),currentInterval.midpoint())));
				newConfigurationRight.insert(RealConstant(param_it->name(),
													  Interval(currentInterval.midpoint(),currentInterval.upper())));

				destination.push_back(newConfigurationLeft);
				destination.push_back(newConfigurationRight);
			}
		}
	}

	return destination;
}

}
