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
    mutable boost::shared_ptr< HybridReachabilityAnalyser > _outer_analyser, _lower_analyser;
  public:
	bool plot_results;
	HybridBoxes safe_region;
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
    //! \name Safety methods

	/*! \brief Attempt to verify that the reachable set of a system starting in an initial_set remains in a safe region inside a domain.
	 * \details This is done in an iterative way by tuning the evolution/analysis parameters. The \a verInfo contains all the information
	 * necessary for verification.
	 */
	tribool safety(SystemVerificationInfo& verInfo) const;

	/*! \brief Compute an underapproximation of the safety/unsafety intervals of \a parameter (defined as an interval) for the automaton
		\a system starting in \a initial_set, where the safe region is \a safe inside \a domain.
        \details The procedure uses the bisection method. The parameter is assumed as having separable safe and unsafe intervals in its range.
        The tolerance in [0 1] is a percentage of the parameter interval width and is used to provide a termination condition for the
		bisection search.
        \return The intervals of safety and unsafety. */
	std::pair<Interval,Interval> parametric_safety_1d_bisection(SystemVerificationInfo& verInfo,
										 							  const RealConstant& parameter,
										 							  const Float& tolerance) const;

	/**
	 * \brief Performs a parametric verification on two parameters \a xParam, \a yParam,
	 * discretizing with \a numPointsPerAxis points for each axis.
	 * \details The procedure uses the bisection method. Saves the results in a file called "<systemName>-<xName>-<yName>" and
	 * generates a "<systemName>-<xName>-<yName>.png" plot, where <systemName> is the name of the system,
	 * <xName> is the name of xParam and <yName> is the name of yParam.
	 */
	Parametric2DBisectionResults parametric_safety_2d_bisection(SystemVerificationInfo& verInfo,
								 const RealConstant& xParam,
								 const RealConstant& yParam,
								 const Float& tolerance,
								 const unsigned& numPointsPerAxis) const;

	/**
	 * \brief Performs a parametric verification on a set of two parameters \a params, by using bisection.
	 */
	Parametric2DBisectionResults parametric_safety_2d_bisection(SystemVerificationInfo& verInfo,
												   const RealConstantSet& params,
												   const Float& tolerance,
												   const unsigned& numPointsPerAxis) const;

	/**
	 * \brief Performs a parametric verification on a set of parameters \a params, by partitioning the parameters space.
	 * \details The \a logNumIntervalsPerParam variable determines how many times any parameter is split in two.
	 * The values in \a params are substituted into the system, the latter
	 * being restored to its initial conditions by the end of the method.
	 */
	ParametricPartitioningOutcomeList parametric_safety_partitioning(SystemVerificationInfo& verInfo,
													   const RealConstantSet& params,
													   const uint& logNumIntervalsPerParam) const;

	//@}

	//@{
	//! \name Dominance methods

	/**
	 * \brief Performs dominance checking.
	 * \details Verifies if the \a dominating system dominates the \a dominated system. Dominance is intended in terms of
	 * the outer reachability of a dominating system being inside the BOUNDING BOX of the lower reachability of a dominated
	 * system minus its epsilon.
	 */
	tribool dominance(SystemVerificationInfo& dominating,
					  SystemVerificationInfo& dominated) const;

	/*! \brief Compute an underapproximation of the dominating/non-dominating intervals of \a parameter for the dominance problem.
        \details The parameter is varied on the \a dominating system alone. The procedure uses the bisection method. The parameter is assumed as having separable dominating and non-dominating intervals in its range.
        The tolerance in [0 1] is a percentage of the parameter interval width and is used to provide a termination condition for the
		bisection search.
        \return The intervals of safety and unsafety. */
	std::pair<Interval,Interval> parametric_dominance_1d_bisection(SystemVerificationInfo& dominating,
			  	  	  	  	  	  	  	  	  	  	  	  	  	   SystemVerificationInfo& dominated,
										 						   const RealConstant& parameter,
										 						   const Float& tolerance) const;
	/**
	 * \brief Performs a parametric dominance checking on two parameters \a xParam, \a yParam,
	 * discretizing with \a numPointsPerAxis points for each axis.
	 * \details The procedure uses the bisection method. Saves the results in a file called "<dominatingName>&<dominatedName>-<xName>-<yName>" and
	 * generates a "<dominatingName>&<dominatedName>-<xName>-<yName>.png" plot, where <systemName> is the name of the system,
	 * <xName> is the name of xParam and <yName> is the name of yParam.
	 */
	Parametric2DBisectionResults parametric_dominance_2d_bisection(SystemVerificationInfo& dominating,
																   SystemVerificationInfo& dominated,
																   const RealConstant& xParam,
																   const RealConstant& yParam,
																   const Float& tolerance,
																   const unsigned& numPointsPerAxis) const;

	/**
	 * \brief Performs a parametric dominance checking on a set of two parameters \a params, by using bisection.
	 */
	Parametric2DBisectionResults parametric_dominance_2d_bisection(SystemVerificationInfo& dominating,
																   SystemVerificationInfo& dominated,
																   const RealConstantSet& params,
																   const Float& tolerance,
																   const unsigned& numPointsPerAxis) const;

	/**
	 * \brief Performs a parametric dominance checking on a set of parameters \a dominating_params of the \a dominating system, by partitioning
	 * the parameters space.
	 * \details The \a logNumIntervalsPerParam variable determines how many times any parameter is split in two.
     * The values in \a dominating_params are substituted into the \a dominating
	 * system alone, the latter being restored to its initial conditions by the end of the method.
	 */
	ParametricPartitioningOutcomeList parametric_dominance_partitioning(SystemVerificationInfo& dominating,
																	    SystemVerificationInfo& dominated,
																	    const RealConstantSet& dominating_params,
																	    const uint& logNumIntervalsPerParam) const;

	//@}

  private:

	//@{
	//! \name Basic verification methods

	/*! \brief Prove (once, i.e. for a given grid depth) that the the reachable set of \a system starting in \a initial_set
	 * definitely remains in the \a safe region. */
	bool _safety_positive_once(HybridAutomaton& system,
							   const HybridImageSet& initial_set,
							   const HybridBoxes& safe_region) const;
	/*! \brief Prove (once, i.e. for a given grid depth) that the reachable set of \a system starting in \a initial_set
	 * does definitely NOT remain in the \a safe region. */
	bool _safety_negative_once(HybridAutomaton& system,
							   const HybridImageSet& initial_set,
							   const HybridBoxes& safe_region) const;

    /*! \brief Attempt (once, i.e. for a given grid depth) to verify that the reachable set of \a system starting in \a initial_set
     * remains in a safe_box. */
    tribool _safety_once(HybridAutomaton& system,
						 const HybridImageSet& initial_set,
						 const HybridBoxes& safe_region) const;

	/*! \brief Performs iterative safety verification where \a parameter is substituted into the system.
	 */
	tribool _safety(SystemVerificationInfo& verInfo,
							 const RealConstant& parameter) const;

	/*! \brief Performs iterative safety verification where the singleton \a value is substituted into the system for the given \a parameter.
	 */
	tribool _safety(SystemVerificationInfo& verInfo,
						     const RealConstant& parameter, const Float& value) const;

	/*! \brief Performs iterative safety verification.
	 * \details The constants in the \a locked_constants set are not allowed to be split */
	tribool _safety_nosplitting(SystemVerificationInfo& verInfo,
							 const RealConstantSet& locked_constants) const;

	//@}

	/*! \brief Shows the verification \a result into the standard output. */
	void _log_verification_result(const tribool& result) const;

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

	/*! \brief Converts the positive/negative search intervals into positive/negative bounds.
	 * \details The result is obtained by knowing the range of the parameter \a parameter_range and the side where
	 * positive values hold, deduced from \a positiveOnBottom. */
	std::pair<Interval,Interval> _pos_neg_bounds_from_search_intervals(const Interval& positive_int,
												 	 	 	 	 	 	 	  const Interval& negative_int,
												 	 	 	 	 	 	 	  const Interval& parameter_range,
												 	 	 	 	 	 	 	  const bool& positiveOnBottom) const;

	/*! \brief Performs one verification sweep along the X axis if \a sweepOnX is true, the Y axis otherwise. */
	void _parametric_safety_2d_bisection_sweep(Parametric2DBisectionResults& results,
						  	  	  	    SystemVerificationInfo& verInfo,
						  	  	  	    RealConstant xParam,
						  	  	  	    RealConstant yParam,
						  	  	  	    const Float& tolerance,
						  	  	  	    const uint& numPointsPerAxis,
						  	  	  	    bool sweepOnX) const;

	/*! \brief Performs one dominance sweep along the X axis if \a sweepOnX is true, the Y axis otherwise. */
	void _parametric_dominance_2d_bisection_sweep(Parametric2DBisectionResults& results,
						  	  	  	    		  SystemVerificationInfo& dominating,
						  	  	  	    		  SystemVerificationInfo& dominated,
						  	  	  	    		  RealConstant xParam,
						  	  	  	    		  RealConstant yParam,
						  	  	  	    		  const Float& tolerance,
						  	  	  	    		  const uint& numPointsPerAxis,
						  	  	  	    		  bool sweepOnX) const;

	/*! \brief Set the parameters for the next dominance iteration, given a bundle of information around a system and a set of constants
	 * that must be ignore when choosing the splitting factors of the system. */
	void _setDominanceParameters(HybridReachabilityAnalyser& analyser,
								 SystemVerificationInfo& systemBundle,
								 const RealConstantSet& lockedConstants) const;

	/**
	 * \brief Performs dominance checking with \a parameter substituted into the \a dominating system.
	 * \details Verifies if the \a dominating system dominates the \a dominated system.
	 */
	tribool _dominance(SystemVerificationInfo& dominating,
					  SystemVerificationInfo& dominated,
					  const RealConstant& parameter) const;

	/**
	 * \brief Performs dominance checking with \a parameter substituted into the \a dominating system with a value of \a value.
	 * \details Verifies if the \a dominating system dominates the \a dominated system.
	 */
	tribool _dominance(SystemVerificationInfo& dominating,
					  SystemVerificationInfo& dominated,
					  const RealConstant& parameter,
					  const Float& value) const;

	/*! \brief Helper function to perform dominance in the more general case when some \a dominatingLockedConstants are enforced. */
	tribool _dominance(SystemVerificationInfo& dominating,
					   SystemVerificationInfo& dominated,
					   const RealConstantSet& dominatingLockedConstants) const;

	/*! \brief Performs the positive part of dominance checking. */
	bool _dominance_positive(SystemVerificationInfo& dominating,
							 SystemVerificationInfo& dominated,
							 const RealConstantSet& dominatingLockedConstants) const;

	/*! \brief Performs the negative part of dominance checking. */
	bool _dominance_negative(SystemVerificationInfo& dominating,
							 SystemVerificationInfo& dominated,
							 const RealConstantSet& dominatingLockedConstants) const;

};

/*! \brief Set the initial evolution parameters of the \a analyser and the grid of the automaton \a system, given the bounding domain \a domain and the safe region \a safe.*/
void setInitialParameters(HybridReachabilityAnalyser& analyser,
						  HybridAutomaton& system,
						  const HybridBoxes& domain,
						  const HybridBoxes& safe,
						  const RealConstantSet& locked_constants);

/* \brief Processes the \a positive_int and \a negative_int initial intervals based on the lower and upper results.
 * \return A variable determining if we must proceed further with bisection refining, and another variable determining
 * the bound where positive values are present.
 */
std::pair<bool,bool> process_initial_bisection_results(Interval& positive_int,
														Interval& negative_int,
														const Interval& parameter_range,
														const tribool& lower_result,
														const tribool& upper_result);

/*! \brief Splits the parameters to the maximum based on the \a tolerance
 * \details The \a numIntervalsPerParam is the number of intervals to split for each parameter.
 * \return The resulting split parameters sets.
 */
std::list<RealConstantSet> maximally_split_parameters(const RealConstantSet& params,
									   	   	   	   	  const uint& numIntervalsPerParam);

}

#endif /* ARIADNE_VERIFIER_H_ */
