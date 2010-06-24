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

#include "discretiser.h"
#include "hybrid_evolver.h"

#include "logging.h"

namespace Ariadne {
 
template<class ES> class Orbit;

class DiscreteState;

template<class BS> class HybridBasicSet;
typedef HybridBasicSet<Box> HybridBox;
typedef std::map<DiscreteState,Vector<Float> > HybridFloatVector;

class HybridGrid;
class HybridGridCell;
class HybridGridTreeSet;

template<class ES> class HybridListSet;
template<class ES> class HybridDiscretiser;

/*! \brief A class for performing reachability analysis on a hybrid system.
	\details Log levels and the methods that are allowed to log output at such a level.
	<br>	1: Methods using the verify_iterative() method; 
	<br>	2: The verify_iterative() method;								
	<br>	3: The verify() and viable() methods;
	<br>	4: Internal methods for verify() and viable();
	<br>	5: Remaining analysis methods implemented from reachability_analyser_interface.h;
	<br>	6: Internal methods for the remaining analysis methods.
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
    //! \name Evaluation of systems on abstract sets
    /*! \brief Compute a lower-approximation to the set obtained by evolving \a system for \a time starting in \a initial_set. */
    virtual SetApproximationType lower_evolve(const SystemType& system,
                                              const HybridImageSet& initial_set, 
                                              const TimeType& time) const;
  
    /*! \brief Compute a lower-approximation to the reachable set of \a system starting in \a initial_set up to \a time. */
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
     * upper semantics; the method performs discretisation before transitions, then checks activations on the discretised cells. */
    virtual SetApproximationType upper_chain_reach(const SystemType& system,
												   const HybridImageSet& initial_set) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set, with
         * lower semantics; the method performs periodical discretisations and checks the new reached region for inclusion
         * The resulting set is a subset of the outer-approximation of the whole evolution set. */
        virtual SetApproximationType lower_chain_reach(const SystemType& system,
    												   const HybridImageSet& initial_set) const;
  
    /*! \brief Compute an outer-approximation to the viability kernel of \a system within \a bounding_set. */
    virtual SetApproximationType viable(const SystemType& system,
                                        const HybridImageSet& bounding_set) const;
  
    /*! \brief Attempt to verify that the reachable set of \a system starting in \a initial_set remains in \a safe_box. */
    virtual tribool verify(const SystemType& system, 
                           const HybridImageSet& initial_set, 
                           const HybridBoxes& safe_box);

	/*! \brief Attempt to verify that the reachable set of \a system starting in \a initial_set remains in \a safe_box inside \a domain, 
		in an iterative way by tuning the evolution/analysis parameters. */
	tribool verify_iterative(SystemType& system, 
							 const HybridImageSet& initial_set, 
							 const HybridBoxes& safe_box, 
							 const HybridBoxes& domain);

	/*! \brief Compute an underapproximation of the safety interval of \a parameter inside the \a parameter_interval for the automaton 
		\a system starting in \a initial_set, where the safe region is \a safe_box inside \a domain. 
        \details The procedure returns the interval of safety. */
	Interval safety_parametric(SystemType& system, 
							   const HybridImageSet& initial_set, 
							   const HybridBoxes& safe_box, 
			 				   const HybridBoxes& domain,
							   const RealConstant& parameter,
							   const Interval& parameter_interval,
							   const Float& tolerance);

	/*! \brief Compute an underapproximation of the safety/unsafety intervals of \a parameter inside the \a parameter_interval for the automaton 
		\a system starting in \a initial_set, where the safe region is \a safe_box inside \a domain. 
        \details The procedure returns the intervals of safety and unsafety. */
	std::pair<Interval,Interval> safety_unsafety_parametric(SystemType& system, 
										 					const HybridImageSet& initial_set, 
										 					const HybridBoxes& safe_box, 
			 							 					const HybridBoxes& domain,
										 					const RealConstant& parameter,
										 					const Interval& parameter_interval,
										 					const Float& tolerance);


    //@}

  public:

	// Determine if the upper/lower reached regions of verify_iterative (and of all the functions using it) must be saved into figures (false by default)
	bool plot_verify_results;
	// Determine if the reached region (resulting from the chain_reach_grid() procedure) is to be dumped into disk for those locations not involved in
	// the evolution, and reloaded when required again.
	bool chain_reach_dumping;
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

  private:

    /*! \brief Generates a list of hybrid enclosures from the \a initial_set, depending on the minimum cell size
     * given by the \a grid.
     */
    list<HybridBasicSet<TaylorSet> > _split_initial_set(const HybridImageSet initial_set,
										   const HybridGrid grid) const;

    // Helper functions for operators on lists of sets.
    GTS _upper_reach(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
    GTS _lower_reach(const Sys& system, std::list<EnclosureType>& initial_set,
    		         const HybridTime& time, const HybridTime& lock_time) const;
    GTS _upper_evolve(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
    std::pair<GTS,GTS> _upper_reach_evolve(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
    std::pair<GTS,GTS> _upper_reach_evolve_continuous(const Sys& sys, const list<EnclosureType>& initial_enclosures, const T& time, const int accuracy) const;

	/*! \brief Verify whether the automaton \a system starting in \a initial_set definitely remains in \a safe_box. */
	bool _safe(const SystemType& system, 
			   const HybridImageSet& initial_set, 
			   const HybridBoxes& safe_box);

	/*! \brief Verify whether the automaton \a system starting in \a initial_set possibly exits the \a safe_box. */	
	bool _unsafe(const SystemType& system, 
			 	 const HybridImageSet& initial_set, 
				 const HybridBoxes& safe_box);

	/*! \brief Get the hybrid maximum absolute derivatives of \system given a previously computed chain reach statistics. 
		\details ASSUMPTION: the continuous variables are preserved in order and quantity between discrete states. */
	HybridFloatVector _getHybridMaximumAbsoluteDerivatives(const SystemType& system) const;

	/*! \brief Set the hybrid maximum integration step size, under the assumption that given the maximum derivatives \a hmad,
		all variables in a step must cover a length greater than a length determined by the \a hgrid. */
	void _setHybridMaximumStepSize(const HybridFloatVector& hmad, const HybridGrid& hgrid);

	/*! \brief Set the hybrid maximum integration step size, under the assumption that given the maximum derivatives \a hmad,
		all variables in a step must cover a length greater than a length determined by the \a hgrid. The value is equal for all
		locations and corresponds to the largest integration step among the locations. */
	void _setEqualizedHybridMaximumStepSize(const HybridFloatVector& hmad, const HybridGrid& hgrid);

	/*! \brief Set the maximum enclosure cell from the hybrid grid \a hgrid. */
	void _setMaximumEnclosureCell(const HybridGrid& hgrid);

	/*! \brief Get the hybrid grid given the maximum derivative \a hmad and the bounding domain parameter, where the grid is chosen differently for each location.*/
	HybridGrid _getHybridGrid(const HybridFloatVector& hmad) const;

	/*! \brief Get the hybrid grid given the maximum derivative \a hmad and the bounding domain parameter, where the grid is chosen equally for all locations.*/
	HybridGrid _getEqualizedHybridGrid(const HybridFloatVector& hmad) const;

	/*! \brief Set the initial evolution parameters and the grid given the automaton \a system and the bounding domain \a domain.*/
	void _setInitialParameters(SystemType& system, 
							   const HybridBoxes& domain);

	/*! \brief Tune the parameters for the next verification iteration, given previous statistics. */
	void _tuneParameters(SystemType& system);
};

template<class HybridEnclosureType>
HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver)
    : _parameters(new EvolutionParametersType())
	, _statistics(new EvolutionStatisticsType())
    , _discretiser(new HybridDiscretiser<typename HybridEnclosureType::ContinuousStateSetType>(evolver))
{
	this->plot_verify_results = false;
	this->chain_reach_dumping = false;
	this->free_cores = 0;
}


template<class HybridEnclosureType>
HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const EvolutionParametersType& parameters,
                           const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver)
    : _parameters(new EvolutionParametersType(parameters))
	, _statistics(new EvolutionStatisticsType())
    , _discretiser(new HybridDiscretiser<typename HybridEnclosureType::ContinuousStateSetType>(evolver))
{
	this->plot_verify_results = false;
	this->chain_reach_dumping = false;
	this->free_cores = 0;
}


} // namespace Ariadne


#endif // ARIADNE_REACHABILITY_ANALYSER_H
