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
    boost::shared_ptr< DiscreteEvolutionParameters > _parameters;
	boost::shared_ptr< DiscreteEvolutionStatistics > _statistics;
    boost::shared_ptr< HybridDiscretiser<HybridEvolver::ContinuousEnclosureType> > _discretiser;
  public:
    typedef DiscreteEvolutionParameters EvolutionParametersType;
	typedef DiscreteEvolutionStatistics	EvolutionStatisticsType;
    typedef HybridAutomaton SystemType;
    typedef SystemType::StateSpaceType StateSpaceType;
    typedef SystemType::TimeType TimeType;
    typedef HybridGridTreeSet SetApproximationType;
    typedef HybridEvolver::EnclosureType EnclosureType;
    typedef HybridEvolver::ContinuousEnclosureType ContinuousEnclosureType;
    typedef SetApproximationType (HybridReachabilityAnalyser::*UpperChainReachFuncPtr)(const SystemType&, const HybridImageSet&) const;
  public:
    //@{
    //! \name Constructors and destructors
    /*! \brief Virtual destructor */
    virtual ~HybridReachabilityAnalyser();

    /*! \brief Construct from a method for evolving basic sets. */
    HybridReachabilityAnalyser(const HybridDiscretiser<HybridEvolver::ContinuousEnclosureType>& discretiser);

    /*! \brief Construct from evolution parameters and a method for evolving basic sets. */
    template<class HybridEnclosureType>
    HybridReachabilityAnalyser(const EvolutionParametersType& parameters,
                               const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver);

    template<class HybridEnclosureType>
    HybridReachabilityAnalyser(const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver);

    /*! \brief Make a dynamically-allocated copy. */
    virtual HybridReachabilityAnalyser* clone() const { return new HybridReachabilityAnalyser(*this); }
    //@}
  
    //@{ 
    //! \name Methods to set and get the parameters controlling the accuracy
    /*! \brief The parameters controlling the accuracy. */
    const EvolutionParametersType& parameters() const { return *this->_parameters; }
    /*! \brief A reference to the parameters controlling the accuracy. */
    EvolutionParametersType& parameters() { return *this->_parameters; }
    //@}

    //@{ 
    //! \name Methods to set and get the statistics related to the analyses
    /*! \brief The statistics stemming from the analyses. */
    const EvolutionStatisticsType& statistics() const { return *this->_statistics; }
    /*! \brief A reference to the statistics stemming from the analyses. */
    EvolutionStatisticsType& statistics() { return *this->_statistics; }
    //@}
  
    //@{

    //! \brief Resets the statistics
    void resetStatistics();

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
    virtual SetApproximationType upper_chain_reach_forward(SystemType& system,
														   const HybridImageSet& initial_set) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set, with
         * lower semantics; the method performs periodical discretisations and checks the new reached region for inclusion
         * The resulting set is a subset of the outer-approximation of the whole evolution set.
         * \return The reach set and the falsification information. */
    virtual std::pair<SetApproximationType,DisproveData> lower_chain_reach(SystemType& system,
																	   const HybridImageSet& initial_set) const;
  
	/*! \brief Tune the parameters for the next verification iteration, given a \a bounding_reach. */
	void tuneIterativeStepParameters(SystemType& system, const HybridGridTreeSet& bounding_reach);

    //@}

  public:

	// Determine if the upper/lower reached regions of verify_iterative (and of all the functions using it) must be saved into figures (false by default)
	bool plot_verify_results;
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

    /*! \brief Generates a list of hybrid enclosures from the \a initial_set, depending on the minimum cell size
     * given by the \a grid. */
    list<HybridBasicSet<TaylorSet> > _split_initial_set(const HybridImageSet initial_set,
										   	   	   	    const HybridGrid grid,
										   	   	   	    Semantics semantics) const;

    // Helper functions for operators on lists of sets.
    GTS _upper_reach(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
    GTS _upper_evolve(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
    std::pair<GTS,GTS> _upper_reach_evolve(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
    std::pair<GTS,GTS> _upper_reach_evolve_continuous(const Sys& sys, const list<EnclosureType>& initial_enclosures, const T& time, const int accuracy) const;
    SetApproximationType _upper_chain_reach_quick_proving(SystemType& system, const HybridImageSet& initial_set,UpperChainReachFuncPtr func) const;
    SetApproximationType _upper_chain_reach_forward(const SystemType& system, const HybridImageSet& initial_set) const;


    /*! \brief Pushes the enclosures from \a reachCells into \a destination.
     * \details Ignores enclosures that lie outside the domain.
     */
    void _upper_chain_reach_forward_pushTargetCells(const HybridGridTreeSet& reachCells,
    									   const SystemType& system,
    									   std::list<EnclosureType>& destination) const;

    /*! \brief Pushes the enclosures from the \a source enclosure into the \a destination enclosure list, for all \a transitions.
     */
    void _upper_chain_reach_pushTargetEnclosures(const std::list<DiscreteTransition>& transitions,
												 const ContinuousEnclosureType& source,
												 const HybridGrid& grid,
												 std::list<EnclosureType>& destination) const;

    /*! \brief Pushes the enclosures from the \a source enclosure into the \a destination enclosure list for a specific transition \a trans.
     * \details Splits the \a source until the enclosure is definitely active for \a trans or the minimum allowed target cell widths \a minTargetCellWidths has been reached.
     */
    void _upper_chain_reach_pushTargetEnclosuresOfTransition(const DiscreteTransition& trans,
    														 const ContinuousEnclosureType& source,
    														 const Vector<Float>& minTargetCellWidths,
    														 std::list<EnclosureType>& destination) const;

    /*! \brief Pushes the target enclosure from the \a source enclosure into the \a destination enclosure list for a specific transition \a trans.
     */
    void _upper_chain_reach_pushTransitioningEnclosure(const DiscreteTransition& trans,
    												   const ContinuousEnclosureType& source,
    												   std::list<EnclosureType>& destination) const;

    /*! \brief Pushes the enclosures from the \a finalCells tree set into the \a destination enclosure list.
     */
    void _upper_chain_reach_pushLocalFinalCells(const HybridGridTreeSet& finalCells,
    									   std::list<EnclosureType>& destination) const;

    /*! \brief Splits \a target_encl for location \a target_loc, storing the result in \a initial_enclosures.
     * \detail The function is recursive.
     */
    void _splitTargetEnclosures(std::list<EnclosureType>& initial_enclosures,
    						    const DiscreteState& target_loc,
    						    const ContinuousEnclosureType& target_encl,
    						    const Vector<Float>& minTargetCellWidths,
    						    const Box& target_bounding) const;

    /*! \brief Pushes the enclosures from the \a finalCells tree set into the \a destination enclosure list.
     *  \detail Does not check for inclusion into the domain.
     */
    void _upper_chain_reach_pushFinalCells_noDomainCheck(const HybridGridTreeSet& finalCells,
    									   	   	   	     std::list<EnclosureType>& destination) const;

    std::pair<SetApproximationType,DisproveData> _lower_chain_reach(const SystemType& system,
    																const HybridImageSet& initial_set) const;
};

template<class HybridEnclosureType>
HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver)
	: _parameters(new EvolutionParametersType())
	, _statistics(new EvolutionStatisticsType())
	, _discretiser(new HybridDiscretiser<typename HybridEnclosureType::ContinuousStateSetType>(evolver))
	, plot_verify_results(false)
	, free_cores(0)
{
}


template<class HybridEnclosureType>
HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const EvolutionParametersType& parameters,
						   const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver)
	: _parameters(new EvolutionParametersType(parameters))
	, _statistics(new EvolutionStatisticsType())
	, _discretiser(new HybridDiscretiser<typename HybridEnclosureType::ContinuousStateSetType>(evolver))
	, plot_verify_results(false)
	, free_cores(0)
{
}

/*! \brief Gets for each non-singleton constant the factor determining the number of chunks its interval should be split into.
 *
 * \details Splits until the deviation of the derivatives is reasonably low in respect to the deviation calculated at the midpoint. This
 * limit value is expressed as a percentage using \a tolerance.
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
							        const HybridBoxes& bounding_domain);

/*! \brief Helper function to get the maximum value of the derivative width ratio \f$ (w-w^r)/w_r \f$, where the \f$ w^r \f$ values
 * are stored in \a referenceWidths and the \f$ w \f$ values are obtained from the \a system.
 */
Float getMaxDerivativeWidthRatio(const HybridAutomaton& system,
								 const HybridFloatVector& referenceWidths,
								 const HybridBoxes& bounding_domain);

/*! \brief Helper function to get the widths of the derivatives from the \a system */
HybridFloatVector getDerivativeWidths(const HybridAutomaton& system,
									  const HybridBoxes& bounding_domain);

/*! \brief Splits the parameters to the maximum based on the \a tolerance
 * \details The \a numIntervalsPerParam is the number of intervals to split for each parameter.
 * \return The resulting split parameters sets.
 */
std::list<RealConstantSet> maximally_split_parameters(const RealConstantSet& params,
									   	   	   	   	  const uint& numIntervalsPerParam);

/*! \brief Gets the set of all the split intervals from the stored split factors.
 *  \details Orders the list elements by first picking the leftmost subintervals, followed by the rightmost and then
 *  all the remaining from right to left.
 */
std::list<RealConstantSet> getSplitConstantsIntervalsSet(const RealConstantIntMap& split_factors);

/*! \brief Gets the set of all the midpoints of the split intervals in \a intervals_set. */
std::list<RealConstantSet> getSplitConstantsMidpointsSet(const std::list<RealConstantSet>& intervals_set);

/*! \brief Get the hybrid grid given the maximum derivative \a hmad and the \a bounding_domain parameter, where the grid is chosen differently for each location.
 * \details The grid is chosen so that each cell is included into the domains for all locations. */
HybridGrid getHybridGrid(const HybridFloatVector& hmad,
						 const HybridBoxes& bounding_domain);

/*! \brief Set the maximum enclosure cell from the hybrid grid \a hgrid and the \a maximum_grid_depth. */
Vector<Float> getMaximumEnclosureCell(const HybridGrid& hgrid, int maximum_grid_depth);

/*! \brief Get the hybrid maximum integration step size, under the assumption that given the maximum derivatives \a hmad,
	all variables in a step must cover a length greater than a length determined by the \a hgrid with a given \a maximum_grid_depth. */
std::map<DiscreteState,Float> getHybridMaximumStepSize(const HybridFloatVector& hmad, const HybridGrid& hgrid, int maximum_grid_depth);

/*! \brief Set the lock to grid time of \system.
	\details The value is taken as the maximum over the times required by any variable on any location to cover a distance equal to
	the domain width of the location, moving at the maximum absolute derivative.
	ASSUMPTION: the continuous variables are preserved in order and quantity between discrete states. */
Float getLockToGridTime(const SystemType& system, const HybridBoxes& bounding_domain);

/*! \brief Get the hybrid maximum absolute derivatives of \system given a previously computed chain reach \a bounding_reach and
 * a domain \a bounding_domain.
 * \details ASSUMPTION: the continuous variables are preserved in order and quantity between discrete states. */
HybridFloatVector getHybridMaximumAbsoluteDerivatives(const SystemType& system,
													  const HybridGridTreeSet& bounding_reach,
													  const HybridBoxes& bounding_domain);

} // namespace Ariadne

#endif // ARIADNE_REACHABILITY_ANALYSER_H
