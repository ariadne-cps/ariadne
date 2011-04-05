/***************************************************************************
 *            reachability_analyser.h
 *
 *  Copyright  2006-10  Alberto Casagrande, Pieter Collins, Luca Geretti
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

#ifndef ARIADNE_REACHABILITY_ANALYSER_H
#define ARIADNE_REACHABILITY_ANALYSER_H

#include <boost/smart_ptr.hpp>
#include <cstdarg>
#include <config.h>

#include "hybrid_set_interface.h"
#include "evolver_interface.h"
#include "discretiser_interface.h"
#include "reachability_analyser_interface.h"

#include "orbit.h"
#include "grid_set.h"
#include "hybrid_set.h"
#include "graphics.h"

#include "discretiser.h"
#include "hybrid_evolver.h"

#include "logging.h"
#include "parametric.h"

namespace Ariadne {
 
template<class ES> class Orbit;

class DiscreteState;

template<class BS> class HybridBasicSet;
typedef HybridBasicSet<Box> HybridBox;
typedef std::map<DiscreteState,Vector<Float> > HybridFloatVector;
typedef std::map<RealConstant,int,ConstantComparator<Real> > RealConstantIntMap;
typedef HybridAutomaton SystemType;
typedef HybridEvolver::EnclosureType EnclosureType;
typedef HybridEvolver::ContinuousEnclosureType ContinuousEnclosureType;

class HybridGrid;
class HybridGridCell;
class HybridGridTreeSet;

template<class ES> class HybridListSet;
template<class ES> class HybridDiscretiser;

/*! \brief A class for performing reachability analysis on a hybrid system-
 */
class HybridReachabilityAnalyser
    : public Loggable
{
  private:
    boost::shared_ptr< DiscreteEvolutionSettings > _settings;
    boost::shared_ptr< HybridDiscretiser<HybridEvolver::ContinuousEnclosureType> > _discretiser;
  public:
    typedef DiscreteEvolutionSettings EvolutionSettingsType;
    typedef HybridAutomaton SystemType;
    typedef SystemType::StateSpaceType StateSpaceType;
    typedef SystemType::TimeType TimeType;
    typedef HybridGridTreeSet SetApproximationType;
    typedef HybridEvolver::EnclosureType EnclosureType;
    typedef HybridEvolver::ContinuousEnclosureType ContinuousEnclosureType;
    typedef SetApproximationType (HybridReachabilityAnalyser::*OuterChainReachFuncPtr)(const SystemType&, const HybridImageSet&) const;
  public:
    //@{
    //! \name Constructors and destructors
    /*! \brief Virtual destructor */
    virtual ~HybridReachabilityAnalyser();

    /*! \brief Construct from a method for evolving basic sets. */
    HybridReachabilityAnalyser(const HybridDiscretiser<HybridEvolver::ContinuousEnclosureType>& discretiser);

    /*! \brief Construct from evolution parameters and a method for evolving basic sets. */
    template<class HybridEnclosureType>
    HybridReachabilityAnalyser(const EvolutionSettingsType& parameters,
                               const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver);

    template<class HybridEnclosureType>
    HybridReachabilityAnalyser(const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver);

    /*! \brief Make a dynamically-allocated copy. */
    virtual HybridReachabilityAnalyser* clone() const { return new HybridReachabilityAnalyser(*this); }
    //@}
  
    //@{ 
    //! \name Methods to set and get the settings controlling the accuracy
    /*! \brief The settings controlling the accuracy. */
    const EvolutionSettingsType& settings() const { return *this->_settings; }
    /*! \brief A reference to the settings controlling the accuracy. */
    EvolutionSettingsType& settings() { return *this->_settings; }
    //@}
  
    //@{

    //! \name Evaluation of systems on abstract sets
    /*! \brief Compute a lower-approximation to the set obtained by evolving \a system for \a time starting in \a initial_set. */
    virtual SetApproximationType lower_evolve(const SystemType& system,
                                              const HybridImageSet& initial_set, 
                                              const TimeType& time) const;
  
    /*! \brief Compute a lower-approximation to the reachable set of \a system starting in \a initial_set up to \a time (discrete part only). */
    virtual SetApproximationType
    lower_reach(const SystemType& system, 
                const HybridImageSet& initial_set, 
                const TimeType& time) const;
  
    /*! \brief Compute a lower-approximation to the reachable and evolved sets of \a system starting in \a initial_set up to \a time. */
    virtual std::pair<SetApproximationType,SetApproximationType>
    lower_reach_evolve(const SystemType& system, 
                       const HybridImageSet& initial_set, 
                       const TimeType& time) const;
  
    /*! \brief Compute an approximation to the set obtained by iterating \a time times \a system starting in \a initial_set. */
    virtual SetApproximationType upper_evolve(const SystemType& system,
                                              const HybridImageSet& initial_set, 
                                              const TimeType& time) const;
  
    /*! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a time times. */
    virtual SetApproximationType upper_reach(const SystemType& system,
                                             const HybridImageSet& initial_set, 
                                             const TimeType& timeType) const;
  
    /*! \brief Compute an approximation to the reachable and evolved sets of \a system starting in \a initial_set iterating at most \a time times. */
    virtual std::pair<SetApproximationType,SetApproximationType>
    upper_reach_evolve(const SystemType& system, 
                       const HybridImageSet& initial_set, 
                       const TimeType& time) const;
  
    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set. */
    virtual SetApproximationType chain_reach(const SystemType& system,
                                             const HybridImageSet& initial_set) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set and remaining in \a bounding_domain. \deprecated */
    virtual SetApproximationType chain_reach(const SystemType& system,
                                             const HybridImageSet& initial_set,
                                             const HybridBoxes& bounding_domain) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set, with
     * upper semantics; the method performs discretisation before transitions, then checks activations on the discretised cells.
     * \return The reach set */
    virtual SetApproximationType outer_chain_reach(SystemType& system,
												   const HybridImageSet& initial_set) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set, with
     * upper semantics.
     * \details The method performs discretisation before transitions, then checks activations on the discretised cells.
     * If \a skipIfOutOfTargetRegion is set, it checks online whether the reachable area is inside the \a target_region, skipping
     * further calculations in that case (useful for safety proving).
     * \return The reach set. */
    virtual SetApproximationType outer_chain_reach(SystemType& system,
												   const HybridImageSet& initial_set,
												   const HybridBoxes& target_region,
												   bool skipIfOutOfTargetRegion,
												   bool skipIfEnclosingTargetRegion) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set, with
     * lower semantics.
     * \details The resulting set is a subset of the outer-approximation of the whole evolution set.
     * \return The reach set and the falsification information (namely, only the reach bounds and the epsilon). */
    virtual std::pair<SetApproximationType,DisproveData> lower_chain_reach(SystemType& system,
																		   const HybridImageSet& initial_set) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set, with
     * lower semantics, checking for inclusion into a \a safe_region.
     * \details the method performs periodical discretisations and checks the new reached region for inclusion
     * The resulting set is a subset of the outer-approximation of the whole evolution set.
     * If \a enable_quick_safety_disproving is set, it checks online whether the reachable area is outside the \a safe_region, skipping
     * further calculations in that case (useful for safety disproving).
     * \return The reach set and the falsification information. */
    virtual std::pair<SetApproximationType,DisproveData> lower_chain_reach(SystemType& system,
																	   const HybridImageSet& initial_set,
																	   const HybridBoxes& safe_region,
																	   bool enable_quick_safety_disproving) const;
  
    /*! \brief Tunes the parameters of the internal evolver. */
    void tuneEvolverParameters(SystemType& system,
								 const HybridFloatVector& hmad,
								 uint maximum_grid_depth,
								 Semantics semantics);

	/*! \brief Tune the parameters for the next verification iteration, given a \a bounding_reach. */
	void tuneIterativeStepParameters(SystemType& system, const HybridGridTreeSet& bounding_reach, Semantics semantics);

	/*! \brief Whether to use a reachability restricting. */
	bool use_reachability_restricting() const;

	/*! \brief Whether to check the domain. */
	bool use_domain_checking() const;

    //@}

  public:

	// The reduction in the number of logical cores used in multithreading (down from the maximum concurrency of the machine) (zero by default)
	uint free_cores;
  
  public:

    typedef HybridTime T;
    typedef HybridAutomaton Sys;
    typedef HybridListSet<Box> BxLS;
    typedef HybridGrid Gr;
    typedef HybridGridCell GC;
    typedef HybridGridTreeSet GCLS;
    typedef HybridGridTreeSet GTS;
    typedef HybridOpenSetInterface OpSI;
    typedef HybridOvertSetInterface OvSI;
    typedef HybridCompactSetInterface CoSI;
    typedef std::map<RealConstant,int,ConstantComparator<Real> > RealConstantMap;

  private:

    /*! \brief Gets the calculus interface from the hybrid evolver. */
    const CalculusInterface<TaylorModel>& _getCalculusInterface() const;

    // Helper functions for operators on lists of sets.
    GTS _upper_reach(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
    GTS _upper_evolve(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
    std::pair<GTS,GTS> _upper_reach_evolve(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
    std::pair<GTS,GTS> _upper_reach_evolve_continuous(const Sys& sys, const list<EnclosureType>& initial_enclosures, const T& time, const int accuracy) const;
    SetApproximationType _outer_chain_reach_forward(const SystemType& system, const HybridImageSet& initial_set) const;

    /*! \brief Verifies whether the \a enclosure_bounds box fails the domain check in respect to \a domain.
     * \details Hides the fact that the check could not be performed at all in certain conditions.
     */
    bool _fails_domain_check(const Box& enclosure_bounds, const Box& domain) const;

    /*! \brief Pushes the enclosures from \a reachCells into \a result_enclosures.
     * \details Ignores enclosures that lie outside the domain.
     */
    void _outer_chain_reach_forward_pushTargetCells(const HybridGridTreeSet& reachCells,
    									   const SystemType& system,
    									   std::list<EnclosureType>& result_enclosures) const;

    /*! \brief Checks whether a box \a bx is outside any invariant from \a invariants. */
    bool _outer_chain_reach_isOutsideInvariants(const Box& bx,
    									   	    const std::map<DiscreteEvent,VectorFunction>& invariants) const;

    /*! \brief Pushes the enclosures from the \a source enclosure into the \a destination enclosure list, for all \a transitions.
     */
    void _outer_chain_reach_pushTargetEnclosures(const std::list<DiscreteTransition>& transitions,
												 const ContinuousEnclosureType& source,
												 const VectorFunction& dynamic,
												 const HybridGrid& grid,
												 std::list<EnclosureType>& result_enclosures) const;

    /*! \brief Pushes the enclosures from the \a source enclosure into the \a destination enclosure list for a specific transition \a trans.
     * \details Splits the \a source until the enclosure is definitely active for \a trans or the minimum allowed target cell widths \a minTargetCellWidths has been reached.
     */
    void _outer_chain_reach_pushTargetEnclosuresOfTransition(const DiscreteTransition& trans,
    														 const VectorFunction& dynamic,
    														 const ContinuousEnclosureType& source,
    														 const Vector<Float>& minTargetCellWidths,
    														 std::list<EnclosureType>& result_enclosures) const;

    /*! \brief Pushes the enclosures from the \a finalCells tree set into the \a result_enclosures list.
     */
    void _outer_chain_reach_pushLocalFinalCells(const HybridGridTreeSet& finalCells,
    									   std::list<EnclosureType>& result_enclosures) const;

    std::pair<SetApproximationType,DisproveData> _lower_chain_reach(const SystemType& system,
    																const HybridImageSet& initial_set,
    																const HybridBoxes& safe_region,
    																bool enable_quick_safety_disproving) const;

    /*! \brief Gets the set of all the split intervals from the \a system with a given \a tolerance.
     *  \details The calculation is performed over the domain with a splitting limit controlled by the \a tolerance, excluding
     *  those constants present in the locked_constants.
     *  Orders the list elements by first picking the leftmost subintervals, followed by the rightmost and then
     *  all the remaining from right to left. If no constants to split are available, returns the original constants.
     */
    std::list<RealConstantSet> _getSplitConstantsIntervalsSet(HybridAutomaton system,
    							  	  	  	  	  	  	  	  float tolerance) const;
};

template<class HybridEnclosureType>
HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver)
	: _settings(new EvolutionSettingsType())
	, _discretiser(new HybridDiscretiser<typename HybridEnclosureType::ContinuousStateSetType>(evolver))
	, free_cores(0)
{
}


template<class HybridEnclosureType>
HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const EvolutionSettingsType& parameters,
						   const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver)
	: _settings(new EvolutionSettingsType(parameters))
	, _discretiser(new HybridDiscretiser<typename HybridEnclosureType::ContinuousStateSetType>(evolver))
	, free_cores(0)
{
}


/*! \brief Generates a list of hybrid enclosures from the \a initial_set, depending on the minimum cell size
 * given by the \a grid. */
list<HybridBasicSet<TaylorSet> > split_initial_set(const HybridImageSet initial_set,
									   	   	   	    const HybridGrid grid,
									   	   	   	    int maximum_grid_depth,
									   	   	   	    Semantics semantics);

/*! \brief Gets for each non-singleton constant the factor determining the number of chunks its interval should be split into.
 *
 * \details Splits until the deviation of the derivatives is reasonably low in respect to the deviation calculated at the midpoint. This
 * limit value is expressed as a percentage using \a targetRatioPerc.
 *
 * @param system The system to get the accessible constants from.
 * @param targetRatioPerc The derivative widths ratio percentage to reach before termination.
 *
 * @return A split factor for each non-singleton accessible constant of the \a system.
 */
RealConstantIntMap getSplitFactorsOfConstants(HybridAutomaton& system, const RealConstantSet& locked_constants,
											   const Float& targetRatioPerc, const HybridBoxes& bounding_domain);

/*! \brief Gets the best constant among the \a working_constants of the \a system to split, in terms of
 * relative reduction of derivative widths compared to some \a referenceWidths.
 */
RealConstant getBestConstantToSplit(SystemType& system,
								    const RealConstantSet& working_constants,
							        const HybridFloatVector& referenceWidths,
							        const HybridBoxes& domain);

/*! \brief Helper function to get the maximum value of the derivative width ratio \f$ (w-w^r)/w_r \f$, where the \f$ w^r \f$ values
 * are stored in \a referenceWidths and the \f$ w \f$ values are obtained from the \a system.
 */
Float getMaxDerivativeWidthRatio(const HybridAutomaton& system,
								 const HybridFloatVector& referenceWidths,
								 const HybridBoxes& domain);

/*! \brief Helper function to get the widths of the derivatives from the \a system */
HybridFloatVector getDerivativeWidths(const HybridAutomaton& system,
									  const HybridBoxes& domain);

/*! \brief Gets the set of all the midpoints of the split intervals in \a intervals_set. */
std::list<RealConstantSet> getSplitConstantsMidpointsSet(const std::list<RealConstantSet>& intervals_set);

/*! \brief Splits the constant \a con into \a numParts parts. */
std::vector<RealConstant> split(const RealConstant& con, uint numParts);

/*! \brief Creates a set \a dest of all the possible combinations of split values from \a src. */
void fillSplitSet(const std::vector<std::vector<RealConstant> >& src,
				   std::vector<std::vector<RealConstant> >::iterator col_it,
				   std::vector<RealConstant>::iterator row_it,
				   RealConstantSet s,
				   std::list<RealConstantSet>& dest);

/*! \brief Splits \a target_encl for location \a target_loc, storing the result in \a initial_enclosures.
 * \details The function is recursive.
 */
void pushSplitTargetEnclosures(std::list<EnclosureType>& initial_enclosures,
						    const DiscreteState& target_loc,
						    const ContinuousEnclosureType& target_encl,
						    const Vector<Float>& minTargetCellWidths,
						    const Box& target_domain_constraint,
						    bool use_constraining);

/*! \brief Get the hybrid grid given the maximum derivative \a hmad and the \a bounding_domain parameter, where the grid is chosen differently for each location.
 * \details The grid is chosen so that each cell is included into the domains for all locations. */
HybridGrid getHybridGrid(const HybridFloatVector& hmad,
						 const HybridBoxes& domain);

/*! \brief Set the maximum enclosure cell from the hybrid grid \a hgrid and the \a maximum_grid_depth. */
Vector<Float> getMaximumEnclosureCell(const HybridGrid& hgrid, int maximum_grid_depth);

/*! \brief Get the hybrid maximum integration step size, under the assumption that given the maximum derivatives \a hmad,
	all variables in a step must cover a length greater than a length determined by the \a hgrid with a given \a maximum_grid_depth.
	\details The actual result is scaled based on the \a semantics. */
std::map<DiscreteState,Float> getHybridMaximumStepSize(
		const HybridFloatVector& hmad,
		const HybridGrid& hgrid,
		int maximum_grid_depth,
		Semantics semantics);

/*! \brief Set the lock to grid time of \system.
	\details The value is taken as the maximum over the times required by any variable on any location to cover a distance equal to
	the domain width of the location, moving at the maximum absolute derivative.
	ASSUMPTION: the continuous variables are preserved in order and quantity between discrete states. */
Float getLockToGridTime(const SystemType& system, const HybridBoxes& domain);

/*! \brief Get the hybrid maximum absolute derivatives of \system given a previously computed outer approximation
 *  \a outer_approx_constraint and a domain \a domain_constraint.
 * \details ASSUMPTION: the continuous variables are preserved in order and quantity between discrete states. */
HybridFloatVector getHybridMaximumAbsoluteDerivatives(const SystemType& system,
													  const HybridGridTreeSet& outer_approx_constraint,
													  const HybridBoxes& domain_constraint);

} // namespace Ariadne

#endif // ARIADNE_REACHABILITY_ANALYSER_H
