/***************************************************************************
 *            verifier.h
 *
 *  Copyright  2011  Luca Geretti
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

/*! \file reachability_analyser.h
 *  \brief Methods for computing abstract reachable sets.
 */

#ifndef ARIADNE_VERIFIER_H_
#define ARIADNE_VERIFIER_H_

#include "reachability_analyser.h"
#include "logging.h"

namespace Ariadne {

class HybridReachabilityAnalyser;

enum DominanceChecking { DOMINANCE_POSITIVE, DOMINANCE_NEGATIVE };

/** Provides a convenient structure for systems and their relative info for verification. */
struct SystemVerificationInfo
{
public:
	HybridAutomaton& _system;
	HybridImageSet& _initial_set;
	HybridBoxes& _domain;
	HybridBoxes _safe_region;
	std::vector<uint> _projection;

	HybridAutomaton& getSystem() { return _system; }
	HybridImageSet& getInitialSet() { return _initial_set; }
	HybridBoxes& getDomain() { return _domain; }
	HybridBoxes& getSafeRegion() { return _safe_region; }
	std::vector<uint>& getProjection() { return _projection; }

	SystemVerificationInfo(HybridAutomaton& system,
						   HybridImageSet& initial_set,
						   HybridBoxes& domain,
						   HybridBoxes& safe_region);

	SystemVerificationInfo(HybridAutomaton& system,
						   HybridImageSet& initial_set,
						   HybridBoxes& domain,
						   std::vector<uint>& projection);

	SystemVerificationInfo(HybridAutomaton& system,
						   HybridImageSet& initial_set,
						   HybridBoxes& domain,
						   HybridBoxes& safe_region,
						   std::vector<uint>& projection);

	virtual std::ostream& write(std::ostream&) const;

private:

	/** \brief Checks fields consistency */
	void _check_fields() const;

};

inline std::ostream& operator<<(std::ostream& os, const SystemVerificationInfo& verInfo) {
    return verInfo.write(os); }


/** \brief Performs verification over reachable sets information. */
class Verifier
    : public Loggable
{
  private:

    /*! \brief Holds the cached result of an outer approximation of a system. */
    struct OuterApproximationCache
    {
      private:
    	bool _is_set;
    	HybridGridTreeSet _value;
      public:
    	OuterApproximationCache() : _is_set(false), _value(HybridGridTreeSet()) {}
    	virtual ~OuterApproximationCache() {}
    	bool is_set() const { return _is_set; }
    	HybridGridTreeSet get() const { return _value; }
    	void set(const HybridGridTreeSet& hgts) { ARIADNE_ASSERT(!_is_set); _is_set = true; _value = hgts; }
    	void reset() { _is_set = false; _value = HybridGridTreeSet(); }
    };

  private:
    mutable boost::shared_ptr< HybridReachabilityAnalyser > _outer_analyser, _lower_analyser;
    boost::shared_ptr< VerificationSettings > _settings;
    mutable std::string _plot_dirpath;

	/*! \brief Fields for caching outer approximations obtained during safety or dominance methods.
	 * \details Their presence is useful to refine the derivatives and consequently the grid on successive iterations. Set to mutable since their value is
	 * valid only transitively during one (iterative) verification method, therefore they do not add external state to the verifier. */
	mutable boost::shared_ptr< OuterApproximationCache > _safety_coarse_outer_approximation;
	mutable boost::shared_ptr< OuterApproximationCache > _dominating_coarse_outer_approximation;
	mutable boost::shared_ptr< OuterApproximationCache > _dominated_coarse_outer_approximation;

	/*! \brief Fields for holding the constraints for outer reachability analyses.
	 * \details The latter two are mandatory since the outer analyser alternatively works on the dominating and dominated systems. */
	mutable HybridGridTreeSet _safety_constraint;
	mutable HybridGridTreeSet _dominating_constraint;
	mutable HybridGridTreeSet _dominated_constraint;

  public:

    //@{
    //! \name Constructors and destructors

    /*! \brief Virtual destructor */
    virtual ~Verifier();

    /*! \brief Construct from a method for evolving basic sets.
     *  \details Explicitly allocates one analyser for outer and lower reachability. */
    Verifier(const HybridReachabilityAnalyser& outer_analyser,
			 const HybridReachabilityAnalyser& lower_analyser);

    /*! \brief Construct from one generic analyser.
     *  \details The analyser will be copied to both the outer and lower internal analysers. The outer and lower analysers
     *  are still independently modifiable afterwards. */
    Verifier(const HybridReachabilityAnalyser& analyser);

    //@}

    //@{
    //! \name Methods to set and get the settings controlling the verification

    const VerificationSettings& settings() const { return *this->_settings; }
    VerificationSettings& settings() { return *this->_settings; }
    //@}

    //@{
    //! \name Safety methods

	/*! \brief Attempt to verify that the reachable set of a system starting in an initial_set remains in a safe region inside a domain.
	 * \details This is done in an iterative way by tuning the evolution/analysis parameters. The \a verInfo contains all the information
	 * necessary for verification.
	 */
	tribool safety(SystemVerificationInfo& verInfo) const;

	/**
	 * \brief Performs a parametric verification on a set of parameters \a params, by partitioning the parameters space.
	 * \details The \a logNumIntervalsPerParam variable determines how many times any parameter is split in two.
	 * The values in \a params are substituted into the system, the latter
	 * being restored to its initial conditions by the end of the method.
	 */
	std::list<ParametricOutcome> parametric_safety(
			SystemVerificationInfo& verInfo,
			const RealConstantSet& params) const;

	/*! \brief Compute an underapproximation of the safety/unsafety intervals of \a parameter (defined as an interval) for the automaton
		\a system starting in \a initial_set, where the safe region is \a safe inside \a domain.
        \details The procedure uses the bisection method. The parameter is assumed as having separable safe and unsafe intervals in its range.
        \return The intervals of safety and unsafety. */
	std::pair<Interval,Interval> parametric_safety_1d_bisection(
			SystemVerificationInfo& verInfo,
			const RealConstant& parameter) const;

	/**
	 * \brief Performs a parametric verification on two parameters \a xParam, \a yParam.
	 * \details The procedure uses the bisection method. Saves the results in a file called "<systemName>-<xName>-<yName>" and
	 * generates a "<systemName>-<xName>-<yName>.png" plot, where <systemName> is the name of the system,
	 * <xName> is the name of xParam and <yName> is the name of yParam.
	 */
	Parametric2DBisectionResults parametric_safety_2d_bisection(
			SystemVerificationInfo& verInfo,
			const RealConstant& xParam,
			const RealConstant& yParam) const;

	/**
	 * \brief Performs a parametric verification on a set of two parameters \a params, by using bisection.
	 */
	Parametric2DBisectionResults parametric_safety_2d_bisection(
			SystemVerificationInfo& verInfo,
			const RealConstantSet& params) const;

	//@}

	//@{
	//! \name Dominance methods

	/**
	 * \brief Performs dominance checking.
	 * \details Verifies if the \a dominating system dominates the \a dominated system. Dominance is intended in terms of
	 * the outer reachability of a dominating system being inside the BOUNDING BOX of the lower reachability of a dominated
	 * system minus its epsilon.
	 */
	tribool dominance(
			SystemVerificationInfo& dominating,
			SystemVerificationInfo& dominated) const;

	/**
	 * \brief Performs a parametric dominance checking on a set of parameters \a dominating_params of the \a dominating system, by partitioning
	 * the parameters space.
	 * \details The \a logNumIntervalsPerParam variable determines how many times any parameter is split in two.
     * The values in \a dominating_params are substituted into the \a dominating
	 * system alone, the latter being restored to its initial conditions by the end of the method.
	 */
	std::list<ParametricOutcome> parametric_dominance(
			SystemVerificationInfo& dominating,
			SystemVerificationInfo& dominated,
			const RealConstantSet& dominating_params) const;

	/*! \brief Compute an underapproximation of the dominating/non-dominating intervals of \a parameter for the dominance problem.
        \details The parameter is varied on the \a dominating system alone. The procedure uses the bisection method. The parameter is assumed as having separable dominating and non-dominating intervals in its range.
        The tolerance in [0 1] is a percentage of the parameter interval width and is used to provide a termination condition for the
		bisection search.
        \return The intervals of safety and unsafety. */
	std::pair<Interval,Interval> parametric_dominance_1d_bisection(
			SystemVerificationInfo& dominating,
			SystemVerificationInfo& dominated,
			const RealConstant& parameter) const;
	/**
	 * \brief Performs a parametric dominance checking on two parameters \a xParam, \a yParam,
	 * discretizing with \a numPointsPerAxis points for each axis.
	 * \details The procedure uses the bisection method. Saves the results in a file called "<dominatingName>&<dominatedName>-<xName>-<yName>" and
	 * generates a "<dominatingName>&<dominatedName>-<xName>-<yName>.png" plot, where <systemName> is the name of the system,
	 * <xName> is the name of xParam and <yName> is the name of yParam.
	 */
	Parametric2DBisectionResults parametric_dominance_2d_bisection(
			SystemVerificationInfo& dominating,
			SystemVerificationInfo& dominated,
			const RealConstant& xParam,
			const RealConstant& yParam) const;

	/**
	 * \brief Performs a parametric dominance checking on a set of two parameters \a params, by using bisection.
	 */
	Parametric2DBisectionResults parametric_dominance_2d_bisection(
			SystemVerificationInfo& dominating,
			SystemVerificationInfo& dominated,
			const RealConstantSet& params) const;

	//@}

  private:

	//@{
	//! \name Safety methods

	/*! \brief Prove (once, i.e. for a given grid depth) that the the reachable set of \a system starting in \a initial_set
	 * definitely remains in the \a safe region.
	 * \details The \a constants are substituted into the system. */
	bool _safety_proving_once(HybridAutomaton& system,
							   const HybridImageSet& initial_set,
							   const HybridBoxes& safe_region,
							   const RealConstantSet& constants) const;
	/*! \brief Prove (once, i.e. for a given grid depth) that the reachable set of \a system starting in \a initial_set
	 * does definitely NOT remain in the \a safe region.
	 * \details The \a constants are substituted into the system. */
	bool _safety_disproving_once(HybridAutomaton& system,
							   const HybridImageSet& initial_set,
							   const HybridBoxes& safe_region,
							   const RealConstantSet& constants) const;

    /*! \brief Attempt (once, i.e. for a given grid depth) to verify that the reachable set of \a system starting in \a initial_set
     * remains in a safe_box.
     * \details The \a constants are substituted into the system. */
    tribool _safety_once(HybridAutomaton& system,
						 const HybridImageSet& initial_set,
						 const HybridBoxes& safe_region,
						 const RealConstantSet& constants) const;

	/*! \brief Performs iterative safety verification where \a constant is substituted into the system.
	 */
	tribool _safety(SystemVerificationInfo& verInfo,
					const RealConstant& constant) const;

	/*! \brief Performs iterative safety verification where the singleton \a value is substituted into the system for the given \a constant.
	 */
	tribool _safety(SystemVerificationInfo& verInfo,
					const RealConstant& constant,
					const Float& value) const;

	/*! \brief Performs iterative safety verification, with \a params_to_substitute substituted into the system.
	 * \details The \a constants are substituted in the system and are not allowed to be split */
	tribool _safety_nosplitting(SystemVerificationInfo& verInfo,
							 const RealConstantSet& constants) const;

	/*! \brief Performs one verification sweep along the X axis if \a sweepOnX is true, the Y axis otherwise. */
	void _parametric_safety_2d_bisection_sweep(Parametric2DBisectionResults& results,
						  	  	  	    SystemVerificationInfo& verInfo,
						  	  	  	    RealConstant xParam,
						  	  	  	    RealConstant yParam,
						  	  	  	    bool sweepOnX) const;

	//@}

	//@{
	//! \name Dominance methods

	/*! \brief Choose the settings for the next dominance iteration, given a bundle of information around a system and a set of constants
	 * that must be ignore when choosing the splitting factors of the system. */
	void _chooseDominanceSettings(SystemVerificationInfo& systemBundle,
								 const RealConstantSet& locked_constants,
								 const HybridGridTreeSet& domain_reach,
								 const HybridGridTreeSet& constraint_reach,
								 Semantics semantics) const;

	/**
	 * \brief Performs dominance checking with \a constants substituted into the \a dominating system.
	 * \details Verifies if the \a dominating system dominates the \a dominated system.
	 */
	tribool _dominance(SystemVerificationInfo& dominating,
					  SystemVerificationInfo& dominated,
					  const RealConstant& constants) const;

	/**
	 * \brief Performs dominance checking with \a constant substituted into the \a dominating system with a value of \a value.
	 * \details Verifies if the \a dominating system dominates the \a dominated system.
	 */
	tribool _dominance(SystemVerificationInfo& dominating,
					  SystemVerificationInfo& dominated,
					  const RealConstant& constant,
					  const Float& value) const;
public:
	/*! \brief Helper function to perform dominance in the more general case when some \a constants are substituted into the dominating system. */
	tribool _dominance(SystemVerificationInfo& dominating,
					   SystemVerificationInfo& dominated,
					   const RealConstantSet& constants) const;
private:
	/*! \brief Performs the positive part of dominance checking. */
	bool _dominance_proving(SystemVerificationInfo& dominating,
							 SystemVerificationInfo& dominated,
							 const RealConstantSet& constants) const;

	/*! \brief Performs the negative part of dominance checking. */
	bool _dominance_disproving(SystemVerificationInfo& dominating,
							 SystemVerificationInfo& dominated,
							 const RealConstantSet& constants) const;

	/*! \brief Performs one dominance sweep along the X axis if \a sweepOnX is true, the Y axis otherwise. */
	void _parametric_dominance_2d_bisection_sweep(Parametric2DBisectionResults& results,
						  	  	  	    		  SystemVerificationInfo& dominating,
						  	  	  	    		  SystemVerificationInfo& dominated,
						  	  	  	    		  RealConstant xParam,
						  	  	  	    		  RealConstant yParam,
						  	  	  	    		  bool sweepOnX) const;

	//@}

	//@{
	//! \name Other helper methods

	/* \brief Processes the \a result in order to update the \a positive_int interval, possibly updating \a negative_int too. */
	void _process_positive_bisection_result(const tribool& result,
											Interval& positive_int,
										    Interval& negative_int,
											const Float& current_value,
											const bool& safeOnBottom) const;

	/*! \brief Processes the \a result in order to update the \a negative_int interval, possibly updating \a positive_int too. */
	void _process_negative_bisection_result(const tribool& result,
											Interval& positive_int,
										    Interval& negative_int,
											const Float& current_value,
											const bool& safeOnBottom) const;

	/*! \brief Choose the initial evolution settings for safety verification of the proper analyser. */
	void _chooseInitialSafetySettings(
			HybridAutomaton& system,
			const HybridBoxes& domain,
			const HybridBoxes& safe,
			const RealConstantSet& locked_constants) const;

	/*! \brief Choose the initial evolution settings for safety verification of the proper analyser, given the \a semantics.*/
	void _chooseInitialSafetySettings(
			HybridAutomaton& system,
			const HybridBoxes& domain,
			const HybridBoxes& safe,
			const RealConstantSet& locked_constants,
			Semantics semantics) const;
	/*! \brief Choose the initial settings for dominance verification.
	 * \details Cannot set the analysers since they are used on different systems on each iteration.
	 */
	void _chooseInitialDominanceSettings() const;

	/*! \brief Tune the settings for the next iterative verification step. */
	void _tuneIterativeStepSettings(
			HybridAutomaton& system,
			const HybridGridTreeSet& hgts_domain,
			const HybridGridTreeSet& constraint_reach,
			Semantics semantics) const;

	/*! \brief Checks whether a grid depth value is allowed for use in iterative verification, based on the \a semantics. */
	bool _is_grid_depth_within_bounds(Semantics semantics) const;

	/*! \brief Updates the constraining information.
	 * \details It reads \a outer_approximation_cache for emptiness;
	 * then possibly sets \a outer_approximation_cache with \a new_outer_approximation. */
	void _update_safety_constraining(
			OuterApproximationCache& outer_approximation_cache,
			const HybridGridTreeSet& new_outer_approximation) const;

	// Reached region plotting methods
	void _plot_dirpath_init(std::string basename) const;
	void _plot_reach(
			const HybridGridTreeSet& reach,
			Semantics semantics) const;
	void _plot_dominance(
			const HybridGridTreeSet& reach,
			Semantics semantics,
			DominanceChecking dominance_checking) const;

	//@}

};

/* \brief Provides a better printing of a tribool verification result */
std::string pretty_print(tribool value);

/*! \brief Processes the \a positive_int and \a negative_int initial intervals based on the lower and upper results.
 *  \return A variable determining if we must proceed further with bisection refining, and another variable determining
 *  if positive values are found on the lower bound of \a positive_int (and consequently, negative values are found on the upper
 *  bound of \a negative_int ).
 */
std::pair<bool,bool> process_initial_bisection_results(
		Interval& positive_int,
		Interval& negative_int,
		const Interval& parameter_range,
		const tribool& lower_result,
		const tribool& upper_result);

/*! \brief Converts the positive/negative search intervals into positive/negative bounds.
 * \details The result is obtained by knowing the range of the parameter \a parameter_range and the side where
 * positive values hold, deduced from \a positiveOnBottom. */
std::pair<Interval,Interval> pos_neg_bounds_from_search_intervals(
		const Interval& positive_int,
		const Interval& negative_int,
		const Interval& parameter_range,
		bool positiveOnBottom);

/*! \brief Splits the parameters to the maximum based on the \a tolerance
 *  \details The \a numIntervalsPerParam is the number of intervals to split for each parameter.
 *  \return The resulting split parameters sets.
 */
std::list<RealConstantSet> maximally_split_parameters(const RealConstantSet& params,
									   	   	   	   	  const uint& maximum_parameter_depth);

}

#endif /* ARIADNE_VERIFIER_H_ */
