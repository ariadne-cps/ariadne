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


#include "hybrid_set_interface.h"
#include "evolver_interface.h"
#include "discretiser_interface.h"
#include "reachability_analyser_interface.h"

#include "orbit.h"
#include "hybrid_orbit.h"
#include "grid_set.h"
#include "hybrid_grid.h"
#include "hybrid_set.h"

#include "logging.h"


namespace Ariadne {

template<class ES> class Orbit;

class DiscreteLocation;
class HybridAutomatonInterface;

template<class BS> class HybridBasicSet;
typedef HybridBasicSet<Box> HybridBox;

class HybridGrid;
class HybridGridCell;
class HybridGridTreeSet;
class HybridGridTreeSet;

template<class ES> class HybridListSet;
template<class ES> class HybridDiscretiser;

class HybridEvolverInterface;


/*! \brief A class for performing reachability analysis on a hybrid system.
 */
class HybridReachabilityAnalyser
    : public ReachabilityAnalyserInterface<HybridAutomatonInterface>
    , public Loggable
{
  private:
    boost::shared_ptr< DiscreteEvolutionParameters > _parameters;
    boost::shared_ptr< HybridEvolverInterface > _evolver;
    boost::shared_ptr< HybridScalingInterface > _scaling;
  public:
    typedef DiscreteEvolutionParameters EvolutionParametersType;
    typedef HybridAutomatonInterface SystemType;
    typedef SystemType::StateSpaceType StateSpaceType;
    typedef SystemType::TimeType TimeType;
    typedef HybridOvertSetInterface OvertSetInterfaceType;
    typedef HybridCompactSetInterface CompactSetInterfaceType;
    typedef HybridLocatedSetInterface LocatedSetInterfaceType;
    typedef HybridRegularSetInterface RegularSetInterfaceType;
    typedef HybridGridTreeSet SetApproximationType;
    typedef HybridBoxes BoundingSetType;
  public:
    //@{
    //! \name Constructors and destructors
    /*! \brief Virtual destructor */
    virtual ~HybridReachabilityAnalyser();

    /*! \brief Construct from a method for evolving basic sets. */
    HybridReachabilityAnalyser(const HybridEvolverInterface& evolver);

    /*! \brief Construct from evolution parameters and a method for evolving basic sets. */
    HybridReachabilityAnalyser(const EvolutionParametersType& parameters,
                               const HybridEvolverInterface& evolver);

    /*! \brief Make a dynamically-allocated copy. */
    virtual HybridReachabilityAnalyser* clone() const;
    //@}

    //@{
    //! \name Methods to set and get the parameters controlling the accuracy.
    /*! \brief The parameters controlling the accuracy. */
    const EvolutionParametersType& parameters() const { return *this->_parameters; }
    /*! \brief A reference to the parameters controlling the accuracy. */
    EvolutionParametersType& parameters() { return *this->_parameters; }
    //! \brief Set the length scales for the variables in each locations.
    void set_scaling(const HybridScalingInterface& hsc);
    //@}


    //@{
    //! \name Evaluation of systems on abstract sets
    /*! \brief Compute a lower-approximation to the set obtained by evolving \a system for \a time starting in \a initial_set. */
    virtual HybridGridTreeSet
    lower_evolve(const SystemType& system,
                 const OvertSetInterfaceType& initial_set,
                 const TimeType& time) const;

    /*! \brief Compute a lower-approximation to the reachable set of \a system starting in \a initial_set up to \a time. */
    virtual HybridGridTreeSet
    lower_reach(const SystemType& system,
                const OvertSetInterfaceType& initial_set,
                const TimeType& time) const;

    /*! \brief Compute a lower-approximation to the reachable and evolved sets of \a system starting in \a initial_set up to \a time. */
    virtual Pair<SetApproximationType,SetApproximationType>
    lower_reach_evolve(const SystemType& system,
                       const OvertSetInterfaceType& initial_set,
                       const TimeType& time) const;

    /*! \brief Compute an approximation to the set obtained by iterating \a time times \a system starting in \a initial_set. */
    virtual SetApproximationType
    upper_evolve(const SystemType& system,
                 const CompactSetInterfaceType& initial_set,
                 const TimeType& time) const;

    /*! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a time times. */
    virtual SetApproximationType
    upper_reach(const SystemType& system,
                const CompactSetInterfaceType& initial_set,
                const TimeType& timeType) const;

    /*! \brief Compute an approximation to the reachable and evolved sets of \a system starting in \a initial_set iterating at most \a time times. */
    virtual Pair<SetApproximationType,SetApproximationType>
    upper_reach_evolve(const SystemType& system,
                       const CompactSetInterfaceType& initial_set,
                       const TimeType& time) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set. */
    virtual SetApproximationType
    chain_reach(const SystemType& system,
                const CompactSetInterfaceType& initial_set) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set and remaining in \a bounding_domain. \deprecated */
    virtual SetApproximationType
    chain_reach(const SystemType& system,
                const CompactSetInterfaceType& initial_set,
                const BoundingSetType& bounding_domain) const;

    /*! \brief Compute an outer-approximation to the viability kernel of \a system within \a bounding_set. */
    virtual SetApproximationType
    viable(const HybridAutomatonInterface& system,
           const HybridCompactSetInterface& bounding_set) const;

    /*! \brief Attempt to verify that the reachable set of \a system starting in \a initial_set remains in \a safe_set. */
    virtual tribool
    verify(const HybridAutomatonInterface& system,
           const HybridLocatedSetInterface& initial_set,
           const HybridRegularSetInterface& safe_set) const;
    //@}

  public:
    typedef HybridTime T;
    typedef HybridAutomatonInterface Sys;
    typedef HybridListSet<Box> BxLS;
    typedef HybridGrid Gr;
    typedef HybridGridCell GC;
    typedef HybridGridTreeSet GCLS;
    typedef HybridGridTreeSet GTS;
    typedef HybridOpenSetInterface OpSI;
    typedef HybridOvertSetInterface OvSI;
    typedef HybridCompactSetInterface CoSI;
  public:
    // Helper functions for operators on lists of sets.
    HybridGridTreeSet _upper_reach(const HybridAutomatonInterface& sys, const HybridGridTreeSet& set, const HybridTime& time, const int accuracy) const;
    HybridGridTreeSet _upper_evolve(const HybridAutomatonInterface& sys, const HybridGridTreeSet& set, const HybridTime& time, const int accuracy) const;
    Pair<HybridGridTreeSet,HybridGridTreeSet> _upper_reach_evolve(const HybridAutomatonInterface& sys, const HybridGridTreeSet& set, const HybridTime& time, const int accuracy) const;
  private:
    // Helper functions for approximating sets
    HybridGrid _hybrid_grid(const Sys& sys) const;
};



} // namespace Ariadne


#endif // ARIADNE_HYBRID_REACHABILITY_ANALYSER_H
