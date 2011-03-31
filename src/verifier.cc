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

Verifier::Verifier(
		const HybridReachabilityAnalyser& outer_analyser,
		const HybridReachabilityAnalyser& lower_analyser) :
			_outer_analyser(outer_analyser.clone()),
			_lower_analyser(lower_analyser.clone()),
			_settings(new VerificationSettings()),
			_safety_hgts_domain(new OuterApproximationCache()),
			_dominating_hgts_domain(new OuterApproximationCache()),
			_dominated_hgts_domain(new OuterApproximationCache())
{
}

Verifier::Verifier(const HybridReachabilityAnalyser& analyser) :
			_outer_analyser(analyser.clone()),
			_lower_analyser(analyser.clone()),
			_settings(new VerificationSettings()),
			_safety_hgts_domain(new OuterApproximationCache()),
			_dominating_hgts_domain(new OuterApproximationCache()),
			_dominated_hgts_domain(new OuterApproximationCache())
{
}

Verifier::~Verifier()
{
}

bool
Verifier::
_safety_positive_once(
		SystemType& system,
		const HybridImageSet& initial_set,
		const HybridBoxes& safe_region,
		const RealConstantSet& constants) const
{
	bool result;
	bool obtained_outer_approximation;

	ARIADNE_LOG(4,"Proving...\n");
	if (!_is_grid_depth_within_bounds(UPPER_SEMANTICS)) {
		ARIADNE_LOG(4,"Not proved.\n");
		return false;
	}

	RealConstantSet original_constants = system.accessible_constants();

	system.substitute(constants,_settings->use_param_midpoints_for_proving);

	DiscreteEvolutionSettings& analyser_params = _outer_analyser->settings();

	ARIADNE_LOG(4,"Setting parameters for this proving iteration...\n");

	_tuneIterativeStepSettings(system,_safety_hgts_domain->get(),_outer_analyser->settings().outer_approx_constraint,UPPER_SEMANTICS);

	ARIADNE_LOG(4,"Performing outer reachability analysis...\n");

	HybridGridTreeSet reach;

	try
	{
		reach = _outer_analyser->outer_chain_reach(system,initial_set,safe_region,_settings->enable_quick_safety_proving);

		result = definitely(reach.subset(safe_region));
		obtained_outer_approximation = true;
	}
	catch (ReachOutOfDomainException ex)
	{
		ARIADNE_LOG(5, "The reached set could not be bounded (" << ex.what() << ").\n");
		result = false;
		obtained_outer_approximation = false;
	}

	system.substitute(original_constants);

	if (obtained_outer_approximation)
		_update_constraining(analyser_params,*_safety_hgts_domain,reach);

	if (_settings->plot_results)
		_plot_reach(reach,UPPER_SEMANTICS);

	ARIADNE_LOG(4, (result ? "Proved.\n" : "Not proved.\n") );

	return result;
}


bool
Verifier::
_safety_negative_once(
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

	_tuneIterativeStepSettings(system,_safety_hgts_domain->get(),_outer_analyser->settings().outer_approx_constraint,LOWER_SEMANTICS);

	ARIADNE_LOG(4,"Performing lower reachability analysis and getting disprove data...\n");

	std::pair<HybridGridTreeSet,DisproveData> reachAndDisproveData =
			_lower_analyser->lower_chain_reach(system,initial_set,safe_region,_settings->enable_quick_safety_disproving);
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

		if (_safety_positive_once(system,initial_set,safe_region,constants)) {
			ARIADNE_LOG(3, "Safe.\n");
			return true;
		}

		if (_safety_negative_once(system,initial_set,safe_region,constants)) {
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
	RealConstantSet constants;

	return _safety_nosplitting(verInfo,constants);
}


tribool
Verifier::
_safety(
		SystemVerificationInfo& verInfo,
		const RealConstant& constant) const
{
	HybridAutomaton& system = verInfo.getSystem();

	Real originalParameterValue = system.accessible_constant_value(constant.name());

	system.substitute(constant);
	tribool result = safety(verInfo);
	system.substitute(constant,originalParameterValue);

	return result;
}

tribool
Verifier::
_safety(
		SystemVerificationInfo& verInfo,
		const RealConstant& constant,
		const Float& value) const
{
	const RealConstant modifiedParameter(constant.name(),Interval(value));

	return _safety(verInfo, modifiedParameter);
}


tribool
Verifier::
_safety_nosplitting(
		SystemVerificationInfo& verInfo,
		const RealConstantSet& constants
		) const
{
	ARIADNE_LOG(2,"\nIterative verification...\n");

	HybridAutomaton& system = verInfo.getSystem();

	if (_settings->plot_results)
		_plot_dirpath_init(system.name());

	_chooseInitialSafetySettings(system,verInfo.getDomain(),verInfo.getSafeRegion(),constants);

	int initial_depth = min(_outer_analyser->settings().lowest_maximum_grid_depth,
							_lower_analyser->settings().lowest_maximum_grid_depth);
	int final_depth = max(_outer_analyser->settings().highest_maximum_grid_depth,
						  _lower_analyser->settings().highest_maximum_grid_depth);
    for (int depth = initial_depth; depth <= final_depth; ++depth)
	{
    	_outer_analyser->settings().maximum_grid_depth = depth;
    	_lower_analyser->settings().maximum_grid_depth = depth;

		ARIADNE_LOG(2, "DEPTH " << depth << "\n");

		tribool result = _safety_once(system,verInfo.getInitialSet(),verInfo.getSafeRegion(),constants);

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

void
Verifier::
_update_constraining(
		DiscreteEvolutionSettings& analyser_settings,
		OuterApproximationCache& outer_approximation_cache,
		const HybridGridTreeSet& new_outer_approximation) const
{
	/* If there is one outer approximation stored, then we have already tuned the grid based on an outer approximation.
	 * Consequently, we maintain the grid and use new_outer_approximation as a constraint.
	 * (on subsequent visits of this code path, we are guaranteed to refine outer_approx_constraint)
	 * Differently, we just save the new_outer_approximation for the subsequent tuning of the grid */
	if (outer_approximation_cache.is_set())
		analyser_settings.outer_approx_constraint = new_outer_approximation;
	else
		outer_approximation_cache.set(new_outer_approximation);
}

std::pair<Interval,Interval>
Verifier::
parametric_safety_1d_bisection(
		SystemVerificationInfo& verInfo,
		const RealConstant& parameter) const
{
	float tolerance = 1.0/(1 << _settings->maximum_parameter_depth);

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
parametric_safety_2d_bisection(
		SystemVerificationInfo& verInfo,
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
parametric_safety_2d_bisection(
		SystemVerificationInfo& verInfo,
		const RealConstant& xParam,
		const RealConstant& yParam) const
{
	// Generates the file name
	std::string filename = verInfo.getSystem().name();
	filename = filename + "[" + xParam.name() + "," + yParam.name() + "]";

	// Initializes the results
	uint numPointsPerAxis = 1+(1<<_settings->maximum_parameter_depth);
	Parametric2DBisectionResults results(filename,xParam.value(),yParam.value(),numPointsPerAxis);

	// Sweeps on each axis
	_parametric_safety_2d_bisection_sweep(results, verInfo, xParam, yParam, true);
	_parametric_safety_2d_bisection_sweep(results, verInfo, xParam, yParam, false);

	return results;
}

std::list<ParametricOutcome>
Verifier::
parametric_safety(
		SystemVerificationInfo& verInfo,
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
		tribool outcome = _safety_nosplitting(verInfo,current_params);
		ARIADNE_LOG(1,"Outcome: " << pretty_print(outcome) << "\n");
		result.push_back(ParametricOutcome(current_params,outcome));
	}

	return result;
}

tribool
Verifier::
dominance(
		SystemVerificationInfo& dominating,
		SystemVerificationInfo& dominated) const
{
	const RealConstantSet dominatingConstants;
	return _dominance(dominating,dominated,dominatingConstants);
}

tribool
Verifier::
_dominance(
		SystemVerificationInfo& dominating,
		SystemVerificationInfo& dominated,
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
		SystemVerificationInfo& dominating,
		SystemVerificationInfo& dominated,
		const RealConstant& constant,
		const Float& value) const
{
	const RealConstant modifiedConstant(constant.name(),Interval(value));

	return _dominance(dominating,dominated,modifiedConstant);
}

std::pair<Interval,Interval>
Verifier::
parametric_dominance_1d_bisection(SystemVerificationInfo& dominating,
								  SystemVerificationInfo& dominated,
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
		SystemVerificationInfo& dominating,
		SystemVerificationInfo& dominated,
		const RealConstantSet& dominating_params) const
{
	ARIADNE_ASSERT_MSG(dominating_params.size() > 0, "Provide at least one parameter.");

	std::list<ParametricOutcome> result;

	std::list<RealConstantSet> splittings = maximally_split_parameters(dominating_params,_settings->maximum_parameter_depth);
	uint i=0;
	for (std::list<RealConstantSet>::const_iterator splitting_it = splittings.begin();
													 splitting_it != splittings.end();
													 ++splitting_it)
	{
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


void Verifier::_parametric_safety_2d_bisection_sweep(
		Parametric2DBisectionResults& results,
		SystemVerificationInfo& verInfo,
		RealConstant xParam,
		RealConstant yParam,
		bool sweepOnX) const
{
	HybridAutomaton& system = verInfo.getSystem();

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
_parametric_dominance_2d_bisection_sweep(
		Parametric2DBisectionResults& results,
		SystemVerificationInfo& dominating,
		SystemVerificationInfo& dominated,
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
		SystemVerificationInfo& dominating,
		SystemVerificationInfo& dominated,
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

		ARIADNE_LOG(2, "DEPTH " << depth << "\n");

		if (_dominance_positive(dominating, dominated, constants)) {
			ARIADNE_LOG(3, "Dominates.\n");
			return true;
		}

		if (_dominance_negative(dominating, dominated, constants)) {
			ARIADNE_LOG(3, "Does not dominate.\n");
			return false;
		}
    }

	// Return indeterminate otherwise
	ARIADNE_LOG(3, "Indeterminate.\n");
	return indeterminate;
}

bool
Verifier::_dominance_positive(
		SystemVerificationInfo& dominating,
		SystemVerificationInfo& dominated,
		const RealConstantSet& constants) const
{
	ARIADNE_LOG(3,"Looking for a positive answer...\n");

	bool result;
	bool obtained_outer_approximation;

	const RealConstantSet& original_constants = dominating.getSystem().accessible_constants();

	dominating.getSystem().substitute(constants,_settings->use_param_midpoints_for_proving);

	HybridGridTreeSet dominating_reach;

	try {

		ARIADNE_LOG(4,"Choosing the settings for the outer approximation of the dominating system...\n");

		_chooseDominanceSettings(dominating,constants,_dominating_hgts_domain->get(),_dominating_constraint,UPPER_SEMANTICS);

		ARIADNE_LOG(4,"Getting the outer approximation of the dominating system...\n");

		dominating_reach = _outer_analyser->outer_chain_reach(dominating.getSystem(),dominating.getInitialSet());

		if (!_dominating_hgts_domain->is_set()) {
			_dominating_hgts_domain->set(dominating_reach);
		} else {
			_dominating_constraint = dominating_reach;
		}

		Box projected_dominating_bounds = Ariadne::project(dominating_reach.bounding_box(),dominating.getProjection());

		ARIADNE_LOG(5,"Projected dominating bounds: " << projected_dominating_bounds << "\n");

		ARIADNE_LOG(4,"Choosing the settings for the lower approximation of the dominated system...\n");

		RealConstantSet emptyLockedConstants;
		_chooseDominanceSettings(dominated,emptyLockedConstants,_dominated_hgts_domain->get(),_dominated_constraint,LOWER_SEMANTICS);

		ARIADNE_LOG(4,"Getting the lower approximation of the dominated system...\n");

		HybridGridTreeSet dominated_reach;
		DisproveData disproveData(dominated.getSystem().state_space());
		make_lpair<HybridGridTreeSet,DisproveData>(dominated_reach,disproveData) =
				_lower_analyser->lower_chain_reach(dominated.getSystem(),dominated.getInitialSet());

		// We must shrink the lower approximation of the the dominated system, but underapproximating in terms of rounding
		HybridBoxes shrinked_dominated_bounds = Ariadne::shrink_in(disproveData.getReachBounds(),disproveData.getEpsilon());

		Box projected_shrinked_dominated_bounds = Ariadne::project(shrinked_dominated_bounds,dominated.getProjection());

		if (_settings->plot_results) {
			_plot_dominance(dominating_reach,UPPER_SEMANTICS,DOMINANCE_POSITIVE);
			_plot_dominance(dominated_reach,LOWER_SEMANTICS,DOMINANCE_POSITIVE);
		}

		ARIADNE_LOG(5,"Epsilon: " << disproveData.getEpsilon() << "\n");
		ARIADNE_LOG(5,"Projected shrinked dominated bounds: " << projected_shrinked_dominated_bounds << "\n");

		result = inside(projected_dominating_bounds,projected_shrinked_dominated_bounds);
		obtained_outer_approximation = true;
	}
	catch (ReachOutOfDomainException ex) {
		ARIADNE_LOG(3,"The outer reached region of the dominating system is unbounded.\n");
		result = false;
		obtained_outer_approximation = false;
	}

	dominating.getSystem().substitute(original_constants);

	return result;
}

bool
Verifier::_dominance_negative(
		SystemVerificationInfo& dominating,
	  	SystemVerificationInfo& dominated,
	  	const RealConstantSet& constants) const
{
	ARIADNE_LOG(3,"Looking for a negative answer...\n");

	bool result;
	bool obtained_outer_approximation;

	const RealConstantSet& original_constants = dominating.getSystem().accessible_constants();

	HybridGridTreeSet dominated_reach;

	dominating.getSystem().substitute(constants,_settings->use_param_midpoints_for_disproving);

	try {

		ARIADNE_LOG(4,"Choosing the settings for the outer approximation of the dominated system...\n");

		RealConstantSet emptyLockedConstants;
		_chooseDominanceSettings(dominated,emptyLockedConstants,_dominated_hgts_domain->get(),_dominated_constraint,UPPER_SEMANTICS);

		ARIADNE_LOG(4,"Getting the outer approximation of the dominated system...\n");

		dominated_reach = _outer_analyser->outer_chain_reach(dominated.getSystem(),dominated.getInitialSet());

		if (!_dominated_hgts_domain->is_set()) {
			_dominated_hgts_domain->set(dominated_reach);
		} else {
			_dominated_constraint = dominated_reach;
		}

		Box projected_dominated_bounds = Ariadne::project(dominated_reach.bounding_box(),dominated.getProjection());

		ARIADNE_LOG(5,"Projected dominated bounds: " << projected_dominated_bounds << "\n");

		ARIADNE_LOG(4,"Choosing the settings for the lower approximation of the dominating system...\n");

		_chooseDominanceSettings(dominating,constants,_dominating_hgts_domain->get(),_dominating_constraint,LOWER_SEMANTICS);

		ARIADNE_LOG(4,"Getting the lower approximation of the dominating system...\n");

		HybridGridTreeSet dominating_reach;
		DisproveData disproveData(dominating.getSystem().state_space());
		make_lpair<HybridGridTreeSet,DisproveData>(dominating_reach,disproveData) =
				_lower_analyser->lower_chain_reach(dominating.getSystem(),dominating.getInitialSet());

		// We must shrink the lower approximation of the the dominating system, but overapproximating in terms of rounding
		HybridBoxes shrinked_dominating_bounds = Ariadne::shrink_out(disproveData.getReachBounds(),disproveData.getEpsilon());

		Box projected_shrinked_dominating_bounds = Ariadne::project(shrinked_dominating_bounds,dominating.getProjection());

		if (_settings->plot_results) {
			_plot_dominance(dominated_reach,UPPER_SEMANTICS,DOMINANCE_NEGATIVE);
			_plot_dominance(dominating_reach,LOWER_SEMANTICS,DOMINANCE_NEGATIVE);
		}

		ARIADNE_LOG(5,"Epsilon: " << disproveData.getEpsilon() << "\n");
		ARIADNE_LOG(5,"Projected shrinked dominating bounds: " << projected_shrinked_dominating_bounds << "\n");

		result = !inside(projected_shrinked_dominating_bounds,projected_dominated_bounds);
		obtained_outer_approximation = true;
	}
	catch (ReachOutOfDomainException ex) {
		ARIADNE_LOG(3,"The outer reached region of the dominated system is unbounded.\n");
		result = false;
		obtained_outer_approximation = false;
	}

	dominating.getSystem().substitute(original_constants);

	return result;
}

void
Verifier::
_chooseInitialSafetySettings(
		HybridAutomaton& system,
		const HybridBoxes& domain,
		const HybridBoxes& safe_region,
		const RealConstantSet& locked_constants) const
{
	_safety_hgts_domain->reset();

	_chooseInitialSafetySettings(system,domain,safe_region,locked_constants,UPPER_SEMANTICS);
	_chooseInitialSafetySettings(system,domain,safe_region,locked_constants,LOWER_SEMANTICS);
}


void
Verifier::
_chooseInitialSafetySettings(
		HybridAutomaton& system,
		const HybridBoxes& domain,
		const HybridBoxes& safe_region,
		const RealConstantSet& locked_constants,
		Semantics semantics) const
{
	ARIADNE_LOG(3,"Choosing the initial settings of the " << (semantics == UPPER_SEMANTICS ? "outer " : "lower ") << "analyser...\n");
	DiscreteEvolutionSettings& settings = (semantics == UPPER_SEMANTICS ?
											   _outer_analyser->settings() : _lower_analyser->settings());
	settings.domain_constraint = domain;
	ARIADNE_LOG(4, "Domain: " << domain << "\n");
	settings.lock_to_grid_time = getLockToGridTime(system,domain);
	ARIADNE_LOG(4, "Lock to grid time: " << settings.lock_to_grid_time << "\n");
	settings.locked_constants = locked_constants;
	ARIADNE_LOG(4, "Locked constants: " << locked_constants << "\n");
}


void
Verifier::
_tuneIterativeStepSettings(
		HybridAutomaton& system,
		const HybridGridTreeSet& hgts_domain,
		const HybridGridTreeSet& constraint_reach,
		Semantics semantics) const
{
	HybridReachabilityAnalyser& analyser = (semantics == UPPER_SEMANTICS ? *_outer_analyser : *_lower_analyser);

	ARIADNE_LOG(5, "Bounding policy: " << (hgts_domain.empty() ? "Domain box" : "Outer approximation") << "\n");
	analyser.settings().constraining_policy = (constraint_reach.empty() ? CONSTRAIN_NO : CONSTRAIN_YES);
	ARIADNE_LOG(5, "Constraining policy: " << (analyser.settings().constraining_policy == CONSTRAIN_NO ? "No": "Yes") << "\n");

	HybridFloatVector hmad = getHybridMaximumAbsoluteDerivatives(system,hgts_domain,analyser.settings().domain_constraint);
	ARIADNE_LOG(5, "Derivative bounds: " << hmad << "\n");
	system.set_grid(getHybridGrid(hmad,analyser.settings().domain_constraint));
	ARIADNE_LOG(5, "Grid: " << system.grid() << "\n");
	analyser.tuneEvolverParameters(system,hmad,analyser.settings().maximum_grid_depth,semantics);
}

void
Verifier::
_chooseInitialDominanceSettings() const
{
	_dominating_hgts_domain->reset();
	_dominated_hgts_domain->reset();

	_dominating_constraint = HybridGridTreeSet();
	_dominated_constraint = HybridGridTreeSet();
}

void
Verifier::
_chooseDominanceSettings(
		SystemVerificationInfo& verInfo,
		const RealConstantSet& locked_constants,
		const HybridGridTreeSet& outer_reach,
		const HybridGridTreeSet& outer_approx_constraint,
		Semantics semantics) const
{
	HybridReachabilityAnalyser& analyser = (semantics == UPPER_SEMANTICS ? *_outer_analyser : *_lower_analyser);

	analyser.settings().domain_constraint = verInfo.getDomain();
	ARIADNE_LOG(5, "Domain: " << analyser.settings().domain_constraint << "\n");

	_tuneIterativeStepSettings(verInfo.getSystem(),outer_reach,outer_approx_constraint,semantics);

	analyser.settings().lock_to_grid_time = getLockToGridTime(verInfo.getSystem(),analyser.settings().domain_constraint);
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
		Semantics semantics,
		DominanceChecking dominance_checking) const
{
	int maximum_grid_depth = _outer_analyser->settings().maximum_grid_depth;

	char mgd_char[10];
	sprintf(mgd_char,"%i",maximum_grid_depth);
	string filename;
	if (semantics == UPPER_SEMANTICS && dominance_checking == DOMINANCE_POSITIVE)
		filename = "dominating-pos-";
	else if (semantics == UPPER_SEMANTICS && dominance_checking == DOMINANCE_NEGATIVE)
		filename = "dominated-neg-";
	else if (semantics == LOWER_SEMANTICS && dominance_checking == DOMINANCE_POSITIVE)
		filename = "dominated-pos-";
	else
		filename = "dominating-neg-";

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
