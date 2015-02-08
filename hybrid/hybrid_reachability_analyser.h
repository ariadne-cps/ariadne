/***************************************************************************
 *            hybrid_reachability_analyser.h
 *
 *  Copyright  2006-11  Alberto Casagrande, Pieter Collins
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

/*! \file hybrid_reachability_analyser.h
 *  \brief Methods for computing reachable sets of hybrid systems.
 */

#ifndef ARIADNE_HYBRID_REACHABILITY_ANALYSER_H
#define ARIADNE_HYBRID_REACHABILITY_ANALYSER_H

#include <boost/smart_ptr.hpp>

#include "solvers/configuration_interface.h"
#include "hybrid/hybrid_set_interface.h"
#include "hybrid/hybrid_evolver_interface.h"
#include "hybrid/hybrid_reachability_analyser_interface.h"

#include "hybrid/hybrid_orbit.h"
#include "hybrid/hybrid_grid.h"
#include "hybrid/hybrid_set.h"

#include "utility/logging.h"


namespace Ariadne {

template<class ES> class Orbit;

class DiscreteLocation;
class HybridAutomatonInterface;

template<class BS> class HybridBasicSet;
class HybridBox;

class HybridGrid;
class HybridGridCell;
class HybridGridTreeSet;
class HybridGridTreeSet;

template<class ES> class HybridListSet;

class HybridEvolverInterface;
class HybridReachabilityAnalyserConfiguration;

//! \ingroup AnalysisModule
//! \brief A class for performing reachability analysis on a hybrid system.
class HybridReachabilityAnalyser
    : public HybridReachabilityAnalyserInterface
{
  public:
    typedef HybridReachabilityAnalyserConfiguration ConfigurationType;
    typedef HybridAutomatonInterface SystemType;
    typedef SystemType::StateSpaceType StateSpaceType;
    typedef SystemType::TimeType TimeType;
    typedef HybridOvertSetInterface OvertSetInterfaceType;
    typedef HybridCompactSetInterface CompactSetInterfaceType;
    typedef HybridLocatedSetInterface LocatedSetInterfaceType;
    typedef HybridRegularSetInterface RegularSetInterfaceType;
    typedef HybridGridTreeSet SetApproximationType;
  private:
    std::shared_ptr< SystemType > _system;
    std::shared_ptr< HybridEvolverInterface > _evolver;
    std::shared_ptr< ConfigurationType > _configuration;
  public:
    //@{
    //! \name Constructors and destructors

    //! \brief Virtual destructor
    virtual ~HybridReachabilityAnalyser();

    //! \brief Construct from an evolver.
    HybridReachabilityAnalyser(
            const SystemType& system,
            const HybridEvolverInterface& evolver);

    //! \brief Make a dynamically-allocated copy.
    virtual HybridReachabilityAnalyser* clone() const;
    //@}

    //@{
    //! \name Methods to set and get the configuration properties of the class.

    //! \brief A reference to the analyser's configuration parameters.
    ConfigurationType& configuration() { return *this->_configuration; }
    const ConfigurationType& configuration() const { return *this->_configuration; }

    //! \brief Get the system associated with the analyser.
    virtual const SystemType& system() const { return *this->_system; }
    //@}



    //@{
    //! \name Evaluation of systems on abstract sets
    //! \brief Compute a lower-approximation to the set obtained by evolving the system for \a time starting in \a initial_set.
    virtual SetApproximationType
    lower_evolve(const OvertSetInterfaceType& initial_set,
                 const TimeType& time) const;

    //! \brief Compute a lower-approximation to the reachable and evolved sets of the system starting in \a initial_set up to \a time.
    virtual Pair<SetApproximationType,SetApproximationType>
    lower_reach_evolve(const OvertSetInterfaceType& initial_set,
                       const TimeType& time) const;

    //! \brief Compute a lower-approximation to the reachable set of the system starting in \a initial_set up to \a time.
    virtual SetApproximationType
    lower_reach(const OvertSetInterfaceType& initial_set,
                const TimeType& time) const;

    //! \brief Compute an infinite-time lower-approximation to the reachable set of the system starting in \a initial_set.
    virtual SetApproximationType
    lower_reach(const OvertSetInterfaceType& initial_set) const;

    //! \brief Compute an upper-approximation to the set obtained by iterating \a time times the system starting in \a initial_set.
    virtual SetApproximationType
    upper_evolve(const CompactSetInterfaceType& initial_set,
                 const TimeType& time) const;

    //! \brief Compute an upper-approximation to the reachable and evolved sets of the system starting in \a initial_set iterating at most \a time times.
    virtual Pair<SetApproximationType,SetApproximationType>
    upper_reach_evolve(const CompactSetInterfaceType& initial_set,
                       const TimeType& time) const;

    //! \brief Compute an upper-approximation to the reachable set of the system starting in \a initial_set iterating at most \a time times.
    virtual SetApproximationType
    upper_reach(const CompactSetInterfaceType& initial_set,
                const TimeType& timeType) const;

    //! \brief Compute a (possibly-restricted) approximation to the outer chain-reachable set of the system starting in \a initial_set.
    virtual SetApproximationType
    outer_chain_reach(const CompactSetInterfaceType& initial_set) const;

    //@}

  public:
    typedef HybridTime T;
    typedef HybridAutomatonInterface Sys;
    typedef HybridListSet<ExactBox> BxLS;
    typedef HybridGrid Gr;
    typedef HybridGridCell GC;
    typedef HybridGridTreeSet GCLS;
    typedef HybridGridTreeSet GTS;
    typedef HybridOpenSetInterface OpSI;
    typedef HybridOvertSetInterface OvSI;
    typedef HybridCompactSetInterface CoSI;
  private:
    // Helper functions for operators on lists of sets.
    Pair<HybridGridTreeSet,HybridGridTreeSet> _reach_evolve_resume(const ListSet<HybridEnclosure>& initial_enclosures,
            const HybridTime& time, const Int accuracy, ListSet<HybridEnclosure>& evolve_enclosures,
            Semantics semantics, const HybridEvolverInterface& evolver) const;
    HybridGridTreeSet _upper_reach(const HybridGridTreeSet& set, const HybridTime& time,
            const Int accuracy, const HybridEvolverInterface& evolver) const;
    HybridGridTreeSet _upper_evolve(const HybridGridTreeSet& set, const HybridTime& time,
            const Int accuracy, const HybridEvolverInterface& evolver) const;
    Void _adjoin_upper_reach_evolve(HybridGridTreeSet& reach_set, HybridGridTreeSet& final_set,
                                    const HybridGridTreeSet& set, const HybridTerminationCriterion& termination,
                                    const Int accuracy, const HybridEvolverInterface& evolver) const;
    //! \brief Perform restriction on \a set, using the overspill policy
    Void _checked_restriction(HybridGridTreeSet& set, const HybridGridTreeSet& bounding) const;
};


//! \brief Configuration for a HybridReachabilityAnalyser, essentially for controlling the
//!    accuracy of discretised evolution methods and reachability analysis.
class HybridReachabilityAnalyserConfiguration : public ConfigurationInterface {

  public:
    //! \brief The integer type.
    typedef Int IntType;
    //! \brief The unsigned integer type.
    typedef Nat UnsignedIntType;
    //! \brief The real type.
    typedef ExactNumber RealType;

    //! \brief Default constructor gives reasonable values.
    HybridReachabilityAnalyserConfiguration(HybridReachabilityAnalyser& analyser);

  private:

    //! \brief Reference back to the main class, for consistency checking.
    HybridReachabilityAnalyser& _analyser;

    //! \brief The time after which infinite-time evolution routines
    //! may approximate computed sets on a grid.
    //! \details
    //! This value should be set to the time after which the transient
    //! behaviour has mostly died out. If there are no transients (i.e. the system evolves
    //! for a certain time and then stops), then this parameter should be set to a value
    //! slightly higher than the maximum running time.
    //! <br>
    //! Setting this value to the time at which transients die out can improve the
    //! speed and accuracy of the computations.
    //!  <br>
    //! This property is only used by chain_reach() routines.
    RealType _transient_time;

    //! \brief The number of discrete steps after which infinite-time evolution
    //! routines may approximate computed sets on the grid.
    //! \details
    //! Note that the transients are assumed to have died out after <em>either</em>
    //! transient_time or transient_steps has been reached.
    //! <br>
    //! For example, if a hybrid system makes at most three steps before settling down
    //! into its limiting behaviour, then this parameter should be set to three.
    //! <br>
    //! Setting this value to the number of steps at which transients die out can improve the
    //! speed and accuracy of the computations.
    //! <br>
    //! This property is only used by chain_reach() routines.
    //! \sa #transient_time
    UnsignedIntType _transient_steps;
    // (See the #transient_time parameter.)

    //! \brief The time after which an upper evolution or reachability analysis routine
    //! may approximate computed sets on a grid, in order to use previously cached
    //! integration results for the grid.
    //! \details
    //! Increasing this parameter improves the accuracy of the computations.
    //! Setting this parameter too low usually results in meaningless computations.
    //! As a rule of thumb, a typical system trajectory should move at least four
    //! times the grid size between locking to the grid. <br>
    //! For forced oscillators, this parameter should be set to the forcing time,
    //! or a multiple or fraction thereof.
    //! <br>
    //! This property is only used for continuous-time computation.
    RealType _lock_to_grid_time;

    //! \brief A set of events after which an upper evolution or reachability analysis
    //! routine should approximate computed sets on a grid, in order to use previously
    //! cached integration results for the grid.
    //! \details
    //! Choosing regularly-occurring events may greatly improve the speed and
    //! accuracy of the computation.
    //! <br>
    //! This property is only used for continuous-time computation.
    Set<DiscreteEvent> _lock_to_grid_events;

    //! \brief The time after which an evolver may approximate computed sets on a grid,
    //! in order to use previously cached results for the grid.
    //! \details
    //! Increasing this parameter may improve the accuracy of the computations.
    //! If there is recurrence in the system, then this parameter should be set to
    //! the average recurrence time, if known.
    //!  <br>
    //! This property is only used for discrete-time computation.
    UnsignedIntType _lock_to_grid_steps;

    //! \brief The depth used for approximation on a grid for computations using upper semantics.
    //! \details
    //! Increasing this value increases the accuracy of the computation.
    //!  <br>
    //! This property is only used in upper_evolve(), upper_reach() and chain_reach() routines.
    IntType _maximum_grid_depth;

    //! \brief The maximum height used for approximation on a grid for chain reachability computations.
    //! \details
    //! Increasing this value increases the bounding domain over which computation is performed.
    //!  <br>
    //! This property is only used in the chain_reach() routines.
    IntType _maximum_grid_height;

    //! \brief The explicit bounding domain for approximation on a grid for chain reachability computations.
    //! \details
    //! If defined, this property combines with _maximum_grid_height to define the actual bounding domain.
    //!  <br>
    //! This property is only used in the chain_reach() routines.
    std::shared_ptr<HybridBoxes> _bounding_domain_ptr;

    //! \brief The grid used for approximation.
    //! \details Provided as a pointer in order to avoid double construction. It must always be defined.
    std::shared_ptr<HybridGrid> _grid_ptr;

    //! \brief The policy for overspill in outer chain reach computation.
    //! \details A OVERSPILL_IGNORE simply bypasses checks, OVERSPILL_WARNING issues a warning and OVERSPILL_ERROR
    //! issues a OuterChainOverspill exception.
    ChainOverspillPolicy _outer_overspill_policy;

  public:

    const RealType& transient_time() const { return _transient_time; }
    Void set_transient_time(const RawFloat64 value) { _transient_time = RealType(value); }

    const UnsignedIntType& transient_steps() const { return _transient_steps; }
    Void set_transient_steps(const UnsignedIntType value) { _transient_steps = value; }

    const RealType& lock_to_grid_time() const { return _lock_to_grid_time; }
    Void set_lock_to_grid_time(const RawFloat64 value) { _lock_to_grid_time = RealType(value); }

    const UnsignedIntType& lock_to_grid_steps() const { return _lock_to_grid_steps; }
    Void set_lock_to_grid_steps(const UnsignedIntType value) { _lock_to_grid_steps = value; }

    const Set<DiscreteEvent>& lock_to_grid_events() const { return _lock_to_grid_events; }
    Void set_lock_to_grid_events(const Set<DiscreteEvent> value) { _lock_to_grid_events = value; }

    const IntType& maximum_grid_depth() const { return _maximum_grid_depth; }
    Void set_maximum_grid_depth(const IntType value) { _maximum_grid_depth = value; }

    const IntType& maximum_grid_height() const { return _maximum_grid_height; }
    Void set_maximum_grid_height(const IntType value) { _maximum_grid_height = value; }

    const std::shared_ptr<HybridBoxes>& bounding_domain_ptr() const { return _bounding_domain_ptr; }
    //! \brief Check the consistency in respect to the system space, then set the bounding domain.
    Void set_bounding_domain_ptr(const std::shared_ptr<HybridBoxes> value);

    const HybridGrid& grid() const { return *_grid_ptr; }
    HybridGrid& grid() { return *_grid_ptr; }
    //! \brief Check the consistency in respect to the system space, then set the grid.
    Void set_grid(const std::shared_ptr<HybridGrid> value_ptr);

    Void set_scaling(const RealVariable& v, RawFloat64 s) {
        dynamic_cast<SimpleHybridScaling&>(this->grid().scalings()).set_scaling(v,ExactFloat64(s)); }

    const ChainOverspillPolicy& outer_overspill_policy() const { return _outer_overspill_policy; }
    Void set_outer_overspill_policy(const ChainOverspillPolicy value) { _outer_overspill_policy = value; }

  public:

    virtual OutputStream& write(OutputStream& os) const;
};


} // namespace Ariadne


#endif // ARIADNE_HYBRID_REACHABILITY_ANALYSER_H
