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

Verifier::Verifier(
		const HybridReachabilityAnalyser& outer_analyser,
		const HybridReachabilityAnalyser& lower_analyser) :
			_outer_analyser(outer_analyser.clone()),
			_lower_analyser(lower_analyser.clone()),
			_settings(new VerificationSettings()),
			_safety_coarse_outer_approximation(new OuterApproximationCache()),
			_dominating_coarse_outer_approximation(new OuterApproximationCache()),
			_dominated_coarse_outer_approximation(new OuterApproximationCache())
{
	_outer_analyser->verb_tab_prefix = verb_tab_prefix + verifier_max_verbosity_level_used;
	_lower_analyser->verb_tab_prefix = verb_tab_prefix + verifier_max_verbosity_level_used;
}

Verifier::Verifier(const HybridReachabilityAnalyser& analyser) :
			_outer_analyser(analyser.clone()),
			_lower_analyser(analyser.clone()),
			_settings(new VerificationSettings()),
			_safety_coarse_outer_approximation(new OuterApproximationCache()),
			_dominating_coarse_outer_approximation(new OuterApproximationCache()),
			_dominated_coarse_outer_approximation(new OuterApproximationCache())
{
	_outer_analyser->verb_tab_prefix = verb_tab_prefix + verifier_max_verbosity_level_used;
	_lower_analyser->verb_tab_prefix = verb_tab_prefix + verifier_max_verbosity_level_used;
}

Verifier::~Verifier()
{
}

bool
Verifier::
_safety_proving_once(
		SystemType& system,
		const HybridImageSet& initial_set,
		const HybridBoxes& safe_region,
		const RealConstantSet& constants) const
{
	bool result;

	ARIADNE_LOG(4,"Proving...\n");
	if (!_is_grid_depth_within_bounds(UPPER_SEMANTICS)) {
		ARIADNE_LOG(4,"Not proved.\n");
		return false;
	}

	RealConstantSet original_constants = system.accessible_constants();

	system.substitute(constants,_settings->use_param_midpoints_for_proving);

	ARIADNE_LOG(4,"Setting parameters for this proving iteration...\n");

	_tuneIterativeStepSettings(system,_safety_coarse_outer_approximation->get(),_safety_reachability_restriction,UPPER_SEMANTICS);

	ARIADNE_LOG(4,"Performing outer reachability analysis...\n");

	try {

		// We quicken the outer reachability calculation only if we already have a constraint available,
		// otherwise we would never obtain an outer approximation for tuning the grid and restricting the reachability.
		bool terminate_as_soon_as_unprovable = _settings->allow_quick_safety_proving && !_safety_reachability_restriction.empty();

		HybridGridTreeSet reach = _outer_analyser->outer_chain_reach(system,initial_set,DIRECTION_FORWARD,terminate_as_soon_as_unprovable,safe_region,NOT_INSIDE_TARGET);

		result = definitely(reach.subset(safe_region));

		ARIADNE_LOG(5, "The reachable set is " << (!result ? "not ":"") << "inside the safe region.\n");

		// We refine only if we have no result from the initial reach set and the grid has already been set using a coarse
		// outer approximation
		if (!result && _safety_coarse_outer_approximation->is_set() && _settings->enable_fb_refinement_for_safety_proving) {
			ARIADNE_LOG(5, "Performing forward-backward refinement...\n");
			if (_outer_analyser->fb_refinement_check(system,initial_set,safe_region,reach))
				result = true;
		}

		if (_safety_coarse_outer_approximation->is_set())
			_safety_reachability_restriction = reach;
		else
			_safety_coarse_outer_approximation->set(reach);

		if (_settings->plot_results)
			_plot_reach(reach,UPPER_SEMANTICS);

	} catch (ReachOutOfDomainException ex) {
		ARIADNE_LOG(5, "The outer reached region is partially out of the domain (" << ex.what() << ").\n");
		result = false;
	} catch (ReachOutOfTargetException ex) {
		ARIADNE_LOG(5, "The outer reached region is partially out of the safe region (" << ex.what() << ").\n");
		result = false;
	}

	system.substitute(original_constants);

	ARIADNE_LOG(4, (result ? "Proved.\n" : "Not proved.\n") );

	return result;
}


bool
Verifier::
_safety_disproving_once(
		SystemType& system,
		const HybridImageSet& initial_set,
		const HybridBoxes& safe_region,
		const RealConstantSet& constants) const
{
	ARIADNE_LOG(4,"Disproving...\n");
	if (!_is_grid_depth_within_bounds(LOWER_SEMANTICS)) {
		ARIADNE_LOG(4,"Not disproved.\n");
		return false;
	}

	RealConstantSet original_constants = system.accessible_constants();

	system.substitute(constants,_settings->use_param_midpoints_for_disproving);

	ARIADNE_LOG(4,"Setting parameters for this disproving iteration...\n");

	_tuneIterativeStepSettings(system,_safety_coarse_outer_approximation->get(),_safety_reachability_restriction,LOWER_SEMANTICS);

	ARIADNE_LOG(4,"Performing lower reachability analysis and getting disprove data...\n");

	// We have no benefit in getting the "whole" lower chain reachability, if we already have disproved:
	// hence we can safely quicken the termination at any depth
	bool terminate_as_soon_as_disproved = _settings->allow_quick_safety_disproving && true;

	std::pair<HybridGridTreeSet,DisproveData> reachAndDisproveData =
			_lower_analyser->lower_chain_reach(system,initial_set,safe_region,terminate_as_soon_as_disproved);
	const HybridGridTreeSet& reach = reachAndDisproveData.first;
	const DisproveData& disproveData = reachAndDisproveData.second;
	const bool& isDisproved = disproveData.getIsDisproved();

	if (_settings->plot_results)
		_plot_reach(reach,LOWER_SEMANTICS);

	ARIADNE_LOG(5,"Disprove data: " << disproveData << "\n");

	ARIADNE_LOG(4, (isDisproved ? "Disproved.\n" : "Not disproved.\n") );

	system.substitute(original_constants);

	return isDisproved;
}


tribool
Verifier::
_safety_once(
		SystemType& system,
		const HybridImageSet& initial_set,
		const HybridBoxes& safe_region,
		const RealConstantSet& constants) const
{
		ARIADNE_LOG(3, "Verification...\n");

		if (_safety_proving_once(system,initial_set,safe_region,constants)) {
			ARIADNE_LOG(3, "Safe.\n");
			return true;
		}

		// It is necessary to widen for disproving, in order to compensate for possible shrinking due to numerical rounding
		HybridBoxes widened_safe_region = widen(safe_region);

		if (_safety_disproving_once(system,initial_set,widened_safe_region,constants)) {
			ARIADNE_LOG(3, "Unsafe.\n");
			return false;
		}

		ARIADNE_LOG(3, "Indeterminate.\n");
		return indeterminate;
}

tribool
Verifier::
safety(SafetyVerificationInput& verInput) const
{
	RealConstantSet constants;

	return _safety_nosplitting(verInput,constants);
}


tribool
Verifier::
_safety(
		SafetyVerificationInput& verInput,
		const RealConstant& constant) const
{
	HybridAutomaton& system = verInput.getSystem();

	Real originalParameterValue = system.accessible_constant_value(constant.name());

	system.substitute(constant);
	tribool result = safety(verInput);
	system.substitute(constant,originalParameterValue);

	return result;
}

tribool
Verifier::
_safety(
		SafetyVerificationInput& verInput,
		const RealConstant& constant,
		const Float& value) const
{
	const RealConstant modifiedParameter(constant.name(),Interval(value));

	return _safety(verInput, modifiedParameter);
}


tribool
Verifier::
_safety_nosplitting(
		SafetyVerificationInput& verInput,
		const RealConstantSet& constants) const
{
	ARIADNE_LOG(2,"\n");
	ARIADNE_LOG(2,"Iterative verification...\n");

	HybridAutomaton& system = verInput.getSystem();

	if (_settings->plot_results)
		_plot_dirpath_init(system.name());

	_chooseInitialSafetySettings(system,verInput.getDomain(),verInput.getSafeRegion(),constants);

	int initial_depth = min(_outer_analyser->settings().lowest_maximum_grid_depth,
							_lower_analyser->settings().lowest_maximum_grid_depth);
	int final_depth = max(_outer_analyser->settings().highest_maximum_grid_depth,
						  _lower_analyser->settings().highest_maximum_grid_depth);
    for (int depth = initial_depth; depth <= final_depth; ++depth) {
    	_outer_analyser->settings().maximum_grid_depth = depth;
    	_lower_analyser->settings().maximum_grid_depth = depth;

		ARIADNE_LOG(2, "Depth " << depth << "\n");

		tribool result = _safety_once(system,verInput.getInitialSet(),verInput.getSafeRegion(),constants);

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
	const DiscreteEvolutionSettings& settings = (semantics == UPPER_SEMANTICS ? _outer_analyser->settings() : _lower_analyser->settings());

	if (settings.maximum_grid_depth < settings.lowest_maximum_grid_depth) {
		ARIADNE_LOG(4,"Skipped verification since the depth is lower than the lowest allowed.\n");
		return false;
	}
	if (settings.maximum_grid_depth > settings.highest_maximum_grid_depth) {
		ARIADNE_LOG(4,"Skipped verification since the depth is higher than the highest allowed.\n");
		return false;
	}

	return true;
}

std::pair<Interval,Interval>
Verifier::
parametric_safety_1d_bisection(
		SafetyVerificationInput& verInput,
		const RealConstant& parameter) const
{
	float tolerance = 1.0/(1 << _settings->maximum_parameter_depth);

	HybridAutomaton& system = verInput.getSystem();

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
	tribool lower_result = _safety(verInput,parameter,parameter_range.lower());
	ARIADNE_LOG(1,pretty_print(lower_result) << ".\n");

	// Check the upper bound
	ARIADNE_LOG(1,"Checking upper interval bound... ");
	tribool upper_result = _safety(verInput,parameter,parameter_range.upper());
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
	while (proceedSafety || proceedUnsafety) {
		// Safety interval check
		if (!safety_int.empty() && safety_int.width() > tolerance*parameter_range.width()) {
			ARIADNE_LOG(1,"Checking safety (positive) interval " << safety_int << " (midpoint: " << safety_int.midpoint() <<
					      ", width ratio: " << safety_int.width()/parameter_range.width()*100 << "%) ... ");

			const Float current_value = safety_int.midpoint();
			tribool result = _safety(verInput,parameter,current_value);

			_process_positive_bisection_result(result,safety_int,unsafety_int,current_value,safeOnBottom);
		} else
			proceedSafety = false;

		// Unsafety interval check
		if (!unsafety_int.empty() && unsafety_int.width() > tolerance*parameter_range.width()) {
			ARIADNE_LOG(1,"Checking unsafety (negative) interval " << unsafety_int << " (midpoint: " << unsafety_int.midpoint() <<
						  ", width ratio: " << unsafety_int.width()/parameter_range.width()*100 << "%) ... ");

			const Float current_value = unsafety_int.midpoint();
			tribool result = _safety(verInput,parameter,current_value);

			_process_negative_bisection_result(result,safety_int,unsafety_int,current_value,safeOnBottom);
		} else
			proceedUnsafety = false;
	}

	system.substitute(parameter,original_value);

	return pos_neg_bounds_from_search_intervals(safety_int,unsafety_int,parameter_range,safeOnBottom);
}

Parametric2DBisectionResults
Verifier::
parametric_safety_2d_bisection(
		SafetyVerificationInput& verInput,
		const RealConstantSet& params) const
{
	ARIADNE_ASSERT_MSG(params.size() == 2,"Provide exactly two parameters.");

	RealConstantSet::const_iterator param_it = params.begin();

	const RealConstant& xParam = *param_it;
	const RealConstant& yParam = *(++param_it);

	return parametric_safety_2d_bisection(verInput,xParam,yParam);
}

Parametric2DBisectionResults
Verifier::
parametric_safety_2d_bisection(
		SafetyVerificationInput& verInput,
		const RealConstant& xParam,
		const RealConstant& yParam) const
{
	// Generates the file name
	std::string filename = verInput.getSystem().name();
	filename = filename + "[" + xParam.name() + "," + yParam.name() + "]";

	// Initializes the results
	uint numPointsPerAxis = 1+(1<<_settings->maximum_parameter_depth);
	Parametric2DBisectionResults results(filename,xParam.value(),yParam.value(),numPointsPerAxis);

	// Sweeps on each axis
	_parametric_safety_2d_bisection_sweep(results, verInput, xParam, yParam, true);
	_parametric_safety_2d_bisection_sweep(results, verInput, xParam, yParam, false);

	return results;
}

std::list<ParametricOutcome>
Verifier::
parametric_safety(
		SafetyVerificationInput& verInput,
		const RealConstantSet& params) const
{
	ARIADNE_ASSERT_MSG(params.size() > 0, "Provide at least one parameter.");

	std::list<ParametricOutcome> result;

	std::list<RealConstantSet> splittings = maximally_split_parameters(params,_settings->maximum_parameter_depth);
	uint i=0;
	for (std::list<RealConstantSet>::const_iterator splitting_it = splittings.begin();
													 splitting_it != splittings.end();
													 ++splitting_it)
	{
		ARIADNE_LOG(1,"<Split parameters set #" << ++i << "/" << splittings.size() << ">\n");
		RealConstantSet current_params = *splitting_it;
		ARIADNE_LOG(1,"Parameter values: " << current_params << " ");
		tribool outcome = _safety_nosplitting(verInput,current_params);
		ARIADNE_LOG(1,"Outcome: " << pretty_print(outcome) << "\n");
		result.push_back(ParametricOutcome(current_params,outcome));
	}

	return result;
}

tribool
Verifier::
dominance(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated) const
{
	const RealConstantSet dominatingConstants;
	return _dominance(dominating,dominated,dominatingConstants);
}

tribool
Verifier::
_dominance(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealConstant& constant) const
{
	HybridAutomaton& system = dominating.getSystem();

	Real original_value = system.accessible_constant_value(constant.name());

	system.substitute(constant);
	tribool result = dominance(dominating,dominated);
	system.substitute(constant,original_value);

	return result;
}

tribool
Verifier::
_dominance(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealConstant& constant,
		const Float& value) const
{
	const RealConstant modifiedConstant(constant.name(),Interval(value));

	return _dominance(dominating,dominated,modifiedConstant);
}

std::pair<Interval,Interval>
Verifier::
parametric_dominance_1d_bisection(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealConstant& parameter) const
{
	float tolerance = 1.0/(1 << _settings->maximum_parameter_depth);

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
parametric_dominance_2d_bisection(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
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
parametric_dominance_2d_bisection(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealConstant& xParam,
		const RealConstant& yParam) const
{
	// Generates the file name
	std::string filename = dominating.getSystem().name() + "&" + dominated.getSystem().name();
	filename = filename + "[" + xParam.name() + "," + yParam.name() + "]";

	// Initializes the results
	uint numPointsPerAxis = (1 << _settings->maximum_parameter_depth);
	Parametric2DBisectionResults results(filename,xParam.value(),yParam.value(),numPointsPerAxis);

	// Sweeps on each axis
	_parametric_dominance_2d_bisection_sweep(results, dominating, dominated, xParam, yParam, true);
	_parametric_dominance_2d_bisection_sweep(results, dominating, dominated, xParam, yParam, false);

	return results;
}


std::list<ParametricOutcome>
Verifier::
parametric_dominance(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealConstantSet& dominating_params) const
{
	ARIADNE_ASSERT_MSG(dominating_params.size() > 0, "Provide at least one parameter.");

	std::list<ParametricOutcome> result;

	std::list<RealConstantSet> splittings = maximally_split_parameters(dominating_params,_settings->maximum_parameter_depth);
	uint i=0;
	for (std::list<RealConstantSet>::const_iterator splitting_it = splittings.begin(); splitting_it != splittings.end(); ++splitting_it) {
		ARIADNE_LOG(1,"<Split parameters set #" << ++i << "/" << splittings.size() << ">\n");
		RealConstantSet current_params = *splitting_it;
		ARIADNE_LOG(1,"Parameter values: " << current_params << " ");
		tribool outcome = _dominance(dominating,dominated,current_params);
		ARIADNE_LOG(1,"Outcome: " << pretty_print(outcome) << "\n");
		result.push_back(ParametricOutcome(current_params,outcome));
	}

	return result;
}

void
Verifier::
_process_positive_bisection_result(
		const tribool& result,
		Interval& positive_int,
		Interval& negative_int,
		const Float& current_value,
		const bool& positiveOnBottom) const
{
	if (definitely(result)) {
		if (positiveOnBottom) {
			// If the negative interval is the same as the positive one, update it too
			if (equal(negative_int,positive_int))
				negative_int.set_lower(current_value);

			ARIADNE_LOG(1,"True, refining upwards.\n");
			positive_int.set_lower(current_value);
		} else {
			// If the negative interval is the same as the positive one, update it too
			if (equal(negative_int,positive_int))
				negative_int.set_upper(current_value);

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
			if (equal(negative_int,positive_int))
				negative_int.set_lower(current_value);

			ARIADNE_LOG(1,"Indeterminate, refining downwards.\n");
			positive_int.set_upper(current_value); }
		else {
			// If the negative interval is the same as the positive interval, update it too
			if (equal(negative_int,positive_int))
				negative_int.set_upper(current_value);

			ARIADNE_LOG(1,"Indeterminate, refining upwards.\n");
			positive_int.set_lower(current_value); }
	}
}

void
Verifier::_process_negative_bisection_result(
		const tribool& result,
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
			if (equal(negative_int,positive_int))
				positive_int.set_upper(current_value);

			ARIADNE_LOG(1,"False, refining downwards.\n");
			negative_int.set_upper(current_value);
		} else {
			// If the negative interval is the same as the positive interval, update it too
			if (equal(negative_int,positive_int))
				positive_int.set_lower(current_value);

			ARIADNE_LOG(1,"False, refining upwards.\n");
			negative_int.set_lower(current_value);
		}
	} else {
		if (positiveOnBottom) {
			// If the negative interval is the same as the positive interval, update it too
			if (equal(negative_int,positive_int))
				positive_int.set_upper(current_value);

			ARIADNE_LOG(1,"Indeterminate, refining upwards.\n");
			negative_int.set_lower(current_value);
		} else {
			// If the negative interval is the same as the positive interval, update it too
			if (equal(negative_int,positive_int))
				positive_int.set_lower(current_value);

			ARIADNE_LOG(1,"Indeterminate, refining downwards.\n");
			negative_int.set_upper(current_value);
		}
	}
}


void Verifier::_parametric_safety_2d_bisection_sweep(
		Parametric2DBisectionResults& results,
		SafetyVerificationInput& verInput,
		RealConstant xParam,
		RealConstant yParam,
		bool sweepOnX) const
{
	HybridAutomaton& system = verInput.getSystem();

	uint numPointsPerAxis = 1+(1 << _settings->maximum_parameter_depth);

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
		std::pair<Interval,Interval> result = parametric_safety_1d_bisection(verInput,otherParam);
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
_parametric_dominance_2d_bisection_sweep(
		Parametric2DBisectionResults& results,
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		RealConstant xParam,
		RealConstant yParam,
		bool sweepOnX) const
{
	HybridAutomaton& system = dominating.getSystem();

	uint numPointsPerAxis = (1 << _settings->maximum_parameter_depth);

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
Verifier::_dominance(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealConstantSet& constants) const
{
	ARIADNE_ASSERT(dominating.getProjection().size() == dominated.getProjection().size());

	ARIADNE_LOG(1, "Dominance checking...\n");

	if (_settings->plot_results)
		_plot_dirpath_init(dominating.getSystem().name() + "&" + dominated.getSystem().name());

	_chooseInitialDominanceSettings();

	int initial_depth = max(_outer_analyser->settings().lowest_maximum_grid_depth,
							_lower_analyser->settings().lowest_maximum_grid_depth);
	int final_depth = min(_outer_analyser->settings().highest_maximum_grid_depth,
						  _lower_analyser->settings().highest_maximum_grid_depth);
    for (int depth = initial_depth; depth <= final_depth; ++depth)
	{
    	_outer_analyser->settings().maximum_grid_depth = depth;
    	_lower_analyser->settings().maximum_grid_depth = depth;

		ARIADNE_LOG(2, "Depth " << depth << "\n");

		if (_dominance_proving_once(dominating, dominated, constants)) {
			ARIADNE_LOG(3, "Dominates.\n");
			return true;
		}

		if (_dominance_disproving_once(dominating, dominated, constants)) {
			ARIADNE_LOG(3, "Does not dominate.\n");
			return false;
		}
    }

	ARIADNE_LOG(3, "Indeterminate.\n");
	return indeterminate;
}

bool
Verifier::_dominance_proving_once(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealConstantSet& constants) const
{
	ARIADNE_LOG(3,"Proving...\n");

	bool result;

	const RealConstantSet& original_constants = dominating.getSystem().accessible_constants();

	dominating.getSystem().substitute(constants,_settings->use_param_midpoints_for_proving);

	try {

		Box shrinked_dominated_bounds = _dominance_shrinked_lower_bounds(dominated,constants,DOMINATED_SYSTEM);

		HybridBoxes shrinked_dominated_bounds_on_dominating_space = Ariadne::project(shrinked_dominated_bounds,
				dominating.getProjection(),dominating.getSystem().state_space());

		Box dominating_bounds = _dominance_outer_bounds(
				dominating,shrinked_dominated_bounds_on_dominating_space,constants,DOMINATING_SYSTEM);

		result = inside(dominating_bounds,shrinked_dominated_bounds);

	} catch (ReachOutOfDomainException ex) {
		ARIADNE_LOG(4,"The outer reached region of the dominating system is partially out of the domain.\n");
		result = false;
	} catch (ReachOutOfTargetException ex) {
		ARIADNE_LOG(4,"The outer reached region of the dominating system is not inside " +
				"the projected shrinked lower reached region of the dominated system.\n");
		result = false;
	}

	ARIADNE_LOG(3, (result ? "Proved.\n" : "Not proved.\n") );

	dominating.getSystem().substitute(original_constants);

	return result;
}

bool
Verifier::_dominance_disproving_once(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
	  	const RealConstantSet& constants) const
{
	ARIADNE_LOG(3,"Disproving...\n");

	bool result;

	const RealConstantSet& original_constants = dominating.getSystem().accessible_constants();

	dominating.getSystem().substitute(constants,_settings->use_param_midpoints_for_disproving);

	try {

		Box shrinked_dominating_bounds = _dominance_shrinked_lower_bounds(dominating,constants,DOMINATING_SYSTEM);

		HybridBoxes shrinked_dominating_bounds_on_dominated_space = Ariadne::project(shrinked_dominating_bounds,
				dominated.getProjection(),dominated.getSystem().state_space());

		Box dominated_bounds = _dominance_outer_bounds(
				dominated,shrinked_dominating_bounds_on_dominated_space,constants,DOMINATED_SYSTEM);

		result = !inside(shrinked_dominating_bounds,dominated_bounds);

	} catch (ReachOutOfDomainException ex) {
		ARIADNE_LOG(4,"The outer reached region of the dominated system is partially out of the domain.\n");
		result = false;
	} catch (ReachEnclosesTargetException ex) {
		ARIADNE_LOG(4,"The outer reached region of the dominated system encloses " +
				"the projected shrinked lower reached region of the dominated system.\n");
		result = true;
	}

	ARIADNE_LOG(3, (result ? "Disproved.\n" : "Not disproved.\n") );

	dominating.getSystem().substitute(original_constants);

	return result;
}


Box
Verifier::
_dominance_shrinked_lower_bounds(
		DominanceVerificationInput& verInfo,
		const RealConstantSet& constants,
		DominanceSystem dominanceSystem) const
{
	string descriptor = (dominanceSystem == DOMINATING_SYSTEM ? "dominating" : "dominated");
	HybridGridTreeSet outer_approximation = (dominanceSystem == DOMINATING_SYSTEM ?
			_dominating_coarse_outer_approximation->get() : _dominated_coarse_outer_approximation->get());
	HybridGridTreeSet reachability_restriction = (dominanceSystem == DOMINATING_SYSTEM ?
			_dominating_reachability_restriction : _dominated_reachability_restriction);

	ARIADNE_LOG(4,"Choosing the settings for the lower reached region of the " << descriptor << " system...\n");

	RealConstantSet emptyLockedConstants;
	_chooseDominanceSettings(verInfo,emptyLockedConstants,outer_approximation,reachability_restriction,LOWER_SEMANTICS);

	ARIADNE_LOG(4,"Getting the lower reached region of the " << descriptor << " system...\n");

	HybridGridTreeSet reach;
	DisproveData disproveData(verInfo.getSystem().state_space());
	make_lpair<HybridGridTreeSet,DisproveData>(reach,disproveData) =
			_lower_analyser->lower_chain_reach(verInfo.getSystem(),verInfo.getInitialSet());

	// We must shrink the lower approximation of the system, but underapproximating in terms of rounding
	HybridBoxes shrinked_bounds = Ariadne::shrink_in(disproveData.getReachBounds(),disproveData.getEpsilon());

	Box projected_shrinked_bounds = Ariadne::project(shrinked_bounds,verInfo.getProjection());

	ARIADNE_LOG(5,"Epsilon: " << disproveData.getEpsilon() << "\n");
	ARIADNE_LOG(5,"Projected shrinked " << descriptor << " bounds: " << projected_shrinked_bounds << "\n");

	if (_settings->plot_results)
		_plot_dominance(reach,dominanceSystem,LOWER_SEMANTICS);

	return projected_shrinked_bounds;
}


Box
Verifier::
_dominance_outer_bounds(
		DominanceVerificationInput& verInput,
		HybridBoxes& lower_bounds_on_this_space,
		const RealConstantSet& constants,
		DominanceSystem dominanceSystem) const
{
	string descriptor = (dominanceSystem == DOMINATING_SYSTEM ? "dominating" : "dominated");
	OuterApproximationCache& outer_approximation_cache = (dominanceSystem == DOMINATING_SYSTEM ?
			*_dominating_coarse_outer_approximation : *_dominated_coarse_outer_approximation);
	HybridGridTreeSet& reachability_restriction = (dominanceSystem == DOMINATING_SYSTEM ?
			_dominating_reachability_restriction : _dominated_reachability_restriction);

	ARIADNE_LOG(4,"Choosing the settings for the outer reached region of the " << descriptor << " system...\n");

	_chooseDominanceSettings(verInput,constants,outer_approximation_cache.get(),reachability_restriction,UPPER_SEMANTICS);

	ARIADNE_LOG(4,"Getting the outer reached region of the " << descriptor << " system...\n");

	bool terminate_as_soon_as_unprovable = _settings->allow_quick_dominance_proving && !reachability_restriction.empty();

	HybridGridTreeSet reach = _outer_analyser->outer_chain_reach(verInput.getSystem(),verInput.getInitialSet(),
			DIRECTION_FORWARD,terminate_as_soon_as_unprovable,lower_bounds_on_this_space,SUPERSET_OF_TARGET);

	Box projected_bounds = Ariadne::project(reach.bounding_box(),verInput.getProjection());

	ARIADNE_LOG(5,"Projected " << descriptor << " bounds: " << projected_bounds << "\n");

	if (!outer_approximation_cache.is_set()) {
		outer_approximation_cache.set(reach);
	} else {
		reachability_restriction = reach;
	}

	if (_settings->plot_results)
		_plot_dominance(reach,dominanceSystem,UPPER_SEMANTICS);

	return projected_bounds;
}


void
Verifier::
_chooseInitialSafetySettings(
		const HybridAutomaton& system,
		const HybridBoxes& domain,
		const HybridBoxes& safe_region,
		const RealConstantSet& locked_constants) const
{
	_safety_coarse_outer_approximation->reset();

	_safety_reachability_restriction = HybridGridTreeSet();

	_chooseInitialSafetySettings(system,domain,safe_region,locked_constants,UPPER_SEMANTICS);
	_chooseInitialSafetySettings(system,domain,safe_region,locked_constants,LOWER_SEMANTICS);
}


void
Verifier::
_chooseInitialSafetySettings(
		const HybridAutomaton& system,
		const HybridBoxes& domain,
		const HybridBoxes& safe_region,
		const RealConstantSet& locked_constants,
		Semantics semantics) const
{
	ARIADNE_LOG(3,"Choosing the initial settings of the " << (semantics == UPPER_SEMANTICS ? "outer " : "lower ") << "analyser...\n");
	DiscreteEvolutionSettings& settings = (semantics == UPPER_SEMANTICS ?
											   _outer_analyser->settings() : _lower_analyser->settings());

	settings.domain_bounds = domain;
	ARIADNE_LOG(4, "Domain: " << domain << "\n");
	settings.lock_to_grid_time = getLockToGridTime(system,domain);
	ARIADNE_LOG(4, "Lock to grid time: " << settings.lock_to_grid_time << "\n");
	settings.locked_constants = locked_constants;
	ARIADNE_LOG(4, "Locked constants: " << locked_constants << "\n");
}


void
Verifier::
_tuneIterativeStepSettings(
		const HybridAutomaton& system,
		const HybridGridTreeSet& hgts_domain,
		const HybridGridTreeSet& reachability_restriction,
		Semantics semantics) const
{
	HybridReachabilityAnalyser& analyser = (semantics == UPPER_SEMANTICS ? *_outer_analyser : *_lower_analyser);

	// Passes through the correct outer approximation constraint to the analyser (used in practice for dominance,
	// where the outer analyser must deal with either the dominating or dominated system)
	analyser.settings().reachability_restriction = reachability_restriction;

	ARIADNE_LOG(5, "Derivatives evaluation policy: " << (hgts_domain.empty() ? "Domain box" : "Outer approximation") << "\n");

	HybridFloatVector hmad = getHybridMaximumAbsoluteDerivatives(system,hgts_domain,analyser.settings().domain_bounds);
	ARIADNE_LOG(5, "Derivatives bounds: " << hmad << "\n");
	analyser.settings().grid = boost::shared_ptr<HybridGrid>(
			new HybridGrid(getHybridGrid(hmad,analyser.settings().domain_bounds)));
	ARIADNE_LOG(5, "Grid lengths: " << analyser.settings().grid->lengths() << "\n");

	ARIADNE_LOG(5, "Use restriction: " << pretty_print(!reachability_restriction.empty()) << "\n");

	analyser.tuneEvolverSettings(system,hmad,analyser.settings().maximum_grid_depth,semantics);
}

void
Verifier::
_chooseInitialDominanceSettings() const
{
	_dominating_coarse_outer_approximation->reset();
	_dominated_coarse_outer_approximation->reset();

	_dominating_reachability_restriction = HybridGridTreeSet();
	_dominated_reachability_restriction = HybridGridTreeSet();
}

void
Verifier::
_chooseDominanceSettings(
		const DominanceVerificationInput& verInput,
		const RealConstantSet& locked_constants,
		const HybridGridTreeSet& outer_reach,
		const HybridGridTreeSet& outer_approx_constraint,
		Semantics semantics) const
{
	HybridReachabilityAnalyser& analyser = (semantics == UPPER_SEMANTICS ? *_outer_analyser : *_lower_analyser);

	analyser.settings().domain_bounds = verInput.getDomain();
	ARIADNE_LOG(5, "Domain: " << analyser.settings().domain_bounds << "\n");

	_tuneIterativeStepSettings(verInput.getSystem(),outer_reach,outer_approx_constraint,semantics);

	analyser.settings().lock_to_grid_time = getLockToGridTime(verInput.getSystem(),analyser.settings().domain_bounds);
	ARIADNE_LOG(5, "Lock to grid time: " << analyser.settings().lock_to_grid_time << "\n");
	analyser.settings().locked_constants = locked_constants;
	ARIADNE_LOG(5, "Locked constants: " << analyser.settings().locked_constants << "\n");
}

void
Verifier::
_plot_dirpath_init(std::string basename) const
{
	time_t mytime;
	time(&mytime);
	string foldername = basename+"-png";

	mkdir(foldername.c_str(),0777);
	string timestring = asctime(localtime(&mytime));
	timestring.erase(std::remove(timestring.begin(), timestring.end(), '\n'), timestring.end());
	foldername = foldername+"/"+timestring;
	mkdir(foldername.c_str(),0777);

	_plot_dirpath = foldername;
}


void
Verifier::
_plot_reach(
		const HybridGridTreeSet& reach,
		Semantics semantics) const
{
	int maximum_grid_depth = (semantics == UPPER_SEMANTICS ?
			_outer_analyser->settings().maximum_grid_depth :
			_lower_analyser->settings().maximum_grid_depth);

	char mgd_char[10];
	sprintf(mgd_char,"%i",maximum_grid_depth);
	string filename = (semantics == UPPER_SEMANTICS ? "outer-" : "lower-");
	filename.append(mgd_char);
	plot(_plot_dirpath,filename,reach);
}


void
Verifier::
_plot_dominance(
		const HybridGridTreeSet& reach,
		DominanceSystem dominanceSystem,
		Semantics semantics) const
{
	int maximum_grid_depth = _outer_analyser->settings().maximum_grid_depth;

	string system_descr = (dominanceSystem == DOMINATING_SYSTEM ? "dominating" : "dominated");
	string verification_descr = (semantics == UPPER_SEMANTICS ? "pos" : "neg");

	char mgd_char[10];
	sprintf(mgd_char,"%i",maximum_grid_depth);
	string filename = system_descr + "-" + verification_descr + "-";
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
process_initial_bisection_results(
		Interval& positive_int,
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
pos_neg_bounds_from_search_intervals(
		const Interval& positive_int,
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
maximally_split_parameters(
		const RealConstantSet& params,
		const uint& maximum_parameter_depth)
{
	std::list<RealConstantSet> source;
	std::list<RealConstantSet> destination;
	destination.push_back(params);

	for (RealConstantSet::const_iterator param_it = params.begin(); param_it != params.end(); ++param_it)
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
