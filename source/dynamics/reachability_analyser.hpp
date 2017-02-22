/***************************************************************************
 *            reachability_analyser.hpp
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

/*! \file reachability_analyser.hpp
 *  \brief Methods for computing reachable sets of hybrid systems.
 */

#ifndef ARIADNE_REACHABILITY_ANALYSER_HPP
#define ARIADNE_REACHABILITY_ANALYSER_HPP

#include "solvers/configuration_interface.hpp"
#include "geometry/set.hpp"
#include "geometry/grid_set.hpp"
#include "geometry/function_set.hpp"

#include "dynamics/orbit.hpp"
#include "dynamics/vector_field_evolver.hpp"
#include "dynamics/reachability_analyser_interface.hpp"

#include "utility/logging.hpp"


namespace Ariadne {

template<class ES> class Orbit;

class DiscreteLocation;
class AutomatonInterface;

class Grid;
class GridCell;
class PavingType;
class PavingType;

template<class ES> class ListSet;

class VectorFieldEvolverType;
class ReachabilityAnalyserConfiguration;

enum class ChainOverspillPolicy : char;

struct SafetyCertificate {
    ValidatedSierpinskian is_safe;
    GridTreeSet chain_reach_set;
    GridTreeSet safe_set;
};

//! \ingroup AnalysisModule
//! \brief A class for performing reachability analysis on a hybrid system.
class ReachabilityAnalyser
    : public ReachabilityAnalyserInterface<VectorField>
    , public Loggable
{
  public:
    typedef ReachabilityAnalyserConfiguration ConfigurationType;
    typedef VectorField SystemType;
    typedef SystemType::TimeType TimeType;
    typedef VectorFieldEvolver EvolverType;
    typedef Enclosure EnclosureType;
    typedef Grid GridType;
    typedef GridTreeSet PavingType;
    typedef OvertSetInterface OvertSetInterfaceType;
    typedef OpenSetInterface OpenSetInterfaceType;
    typedef CompactSetInterface CompactSetInterfaceType;
    typedef LocatedSetInterface LocatedSetInterfaceType;
    typedef RegularSetInterface RegularSetInterfaceType;
    typedef SafetyCertificate SafetyCertificateType;
  private:
    std::shared_ptr< VectorField > _system;
    std::shared_ptr< VectorFieldEvolver > _evolver;
    std::shared_ptr< ConfigurationType > _configuration;
  public:
    //@{
    //! \name Constructors and destructors

    //! \brief Virtual destructor
    virtual ~ReachabilityAnalyser();

    //! \brief Construct from an evolver.
    ReachabilityAnalyser(
            const SystemType& system,
            const VectorFieldEvolver& evolver);

    //! \brief Make a dynamically-allocated copy.
    virtual ReachabilityAnalyser* clone() const;
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
    //! \brief Compute a lower-approximation to the set obtained by evolving the system starting in \a initial_set until \a time.
    virtual PavingType
    lower_evolve(const OvertSetInterfaceType& initial_set,
                 const TimeType& time) const;

    //! \brief Compute a lower-approximation to the reachable and evolved sets of the system starting in \a initial_set up to \a time.
    virtual Pair<PavingType,PavingType>
    lower_reach_evolve(const OvertSetInterfaceType& initial_set,
                       const TimeType& time) const;

    //! \brief Compute a lower-approximation to the reachable set of the system starting in \a initial_set up to \a time.
    virtual PavingType
    lower_reach(const OvertSetInterfaceType& initial_set,
                const TimeType& time) const;

    //! \brief Compute an infinite-time lower-approximation to the reachable set of the system starting in \a initial_set.
    virtual PavingType
    lower_reach(const OvertSetInterfaceType& initial_set) const;

    //! \brief Compute an upper-approximation to the set obtained by evolving the system starting in \a initial_set up to \a time.
    virtual PavingType
    upper_evolve(const CompactSetInterfaceType& initial_set,
                 const TimeType& time) const;

    //! \brief Compute upper-approximations to the reachable and evolved sets of the system starting in \a initial_set up to \a time.
    virtual Pair<PavingType,PavingType>
    upper_reach_evolve(const CompactSetInterfaceType& initial_set,
                       const TimeType& time) const;

    //! \brief Compute an upper-approximation to the reachable set of the system starting in \a initial_set up to \a time.
    virtual PavingType
    upper_reach(const CompactSetInterfaceType& initial_set,
                const TimeType& time) const;

    //! \brief Compute a (possibly-restricted) approximation to the outer chain-reachable set of the system starting in \a initial_set.
    virtual PavingType
    outer_chain_reach(const CompactSetInterfaceType& initial_set) const;

    //! \brief Test if the system is safe.
    virtual SafetyCertificate
    verify_safety(const CompactSetInterfaceType& initial_set, const OpenSetInterfaceType& safe_set) const;

    //@}

  private:
    // Helper functions for operators on lists of sets.
    Pair<PavingType,PavingType> _reach_evolve_resume(const ListSet<EnclosureType>& initial_enclosures,
            const TimeType& time, const Int accuracy, ListSet<EnclosureType>& evolve_enclosures,
            Semantics semantics, const EvolverType& evolver) const;
    PavingType _upper_reach(const PavingType& set, const TimeType& time,
            const Int accuracy, const EvolverType& evolver) const;
    PavingType _upper_evolve(const PavingType& set, const TimeType& time,
            const Int accuracy, const EvolverType& evolver) const;
    Void _adjoin_upper_reach_evolve(PavingType& reach_set, PavingType& final_set,
                                    const PavingType& set, const TimeType& time,
                                    const Int accuracy, const EvolverType& evolver) const;
    //! \brief Perform restriction on \a set, using the overspill policy
    Void _checked_restriction(PavingType& set, const PavingType& bounding) const;
};


//! \brief Configuration for a ReachabilityAnalyser, essentially for controlling the
//!    accuracy of discretised evolution methods and reachability analysis.
class ReachabilityAnalyserConfiguration : public ConfigurationInterface {

  public:
    //! \brief The integer type.
    typedef Int IntType;
    //! \brief The unsigned integer type.
    typedef Nat UnsignedIntType;
    //! \brief The type representing time.
    typedef typename VectorField::TimeType TimeType;

    //! \brief Default constructor gives reasonable values.
    ReachabilityAnalyserConfiguration(ReachabilityAnalyser& analyser);

  private:

    //! \brief Reference back to the main class, for consistency checking.
    ReachabilityAnalyser& _analyser;

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
    TimeType _transient_time;

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
    TimeType _lock_to_grid_time;

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
    std::shared_ptr<ExactBoxType> _bounding_domain_ptr;

    //! \brief The grid used for approximation.
    //! \details Provided as a pointer in order to avoid double construction. It must always be defined.
    std::shared_ptr<Grid> _grid_ptr;

    //! \brief The policy for overspill in outer chain reach computation.
    //! \details A OVERSPILL_IGNORE simply bypasses checks, OVERSPILL_WARNING issues a warning and OVERSPILL_ERROR
    //! issues a OuterChainOverspill exception.
    ChainOverspillPolicy _outer_overspill_policy;

  public:

    const TimeType& transient_time() const { return _transient_time; }
    Void set_transient_time(const TimeType value) { _transient_time = TimeType(value); }

    const TimeType& lock_to_grid_time() const { return _lock_to_grid_time; }
    Void set_lock_to_grid_time(const TimeType value) { _lock_to_grid_time = TimeType(value); }

    const IntType& maximum_grid_depth() const { return _maximum_grid_depth; }
    Void set_maximum_grid_depth(const IntType value) { _maximum_grid_depth = value; }

    const IntType& maximum_grid_height() const { return _maximum_grid_height; }
    Void set_maximum_grid_height(const IntType value) { _maximum_grid_height = value; }

    const std::shared_ptr<ExactBoxType>& bounding_domain_ptr() const { return _bounding_domain_ptr; }
    Void set_bounding_domain_ptr(const std::shared_ptr<ExactBoxSet> value);

    const Grid& grid() const { return *_grid_ptr; }
    Void set_grid(const std::shared_ptr<Grid> grid_ptr);

    const ChainOverspillPolicy& outer_overspill_policy() const { return _outer_overspill_policy; }
    Void set_outer_overspill_policy(const ChainOverspillPolicy value) { _outer_overspill_policy = value; }

  public:

    virtual OutputStream& write(OutputStream& os) const;
};


} // namespace Ariadne


#endif // ARIADNE_REACHABILITY_ANALYSER_HPP
