/***************************************************************************
 *      reachability_analyser.h
 *
 *  Copyright  2006-8  Alberto Casagrande, Pieter Collins
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


#include "hybrid_set_interface.h"
#include "evolver_interface.h"
#include "discretiser_interface.h"
#include "reachability_analyser_interface.h"

#include "orbit.h"
#include "grid_set.h"
#include "hybrid_set.h"

#include "discretiser.h"

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



/*! \brief A class for performing reachability analysis on a hybrid system.
 */
class HybridReachabilityAnalyser
    : public ReachabilityAnalyserInterface<HybridAutomatonInterface>
    , public Loggable
{
  private:
    boost::shared_ptr< DiscreteEvolutionParameters > _parameters;
    boost::shared_ptr< DiscretiserInterface<HybridAutomatonInterface,HybridGridCell> > _discretiser;
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
    HybridReachabilityAnalyser(const DiscretiserInterface<HybridAutomatonInterface,HybridGridCell>& discretiser);

    /*! \brief Construct from evolution parameters and a method for evolving basic sets. */
    template<class HybridEnclosureType>
    HybridReachabilityAnalyser(const EvolutionParametersType& parameters,
    const EvolverInterface<HybridAutomatonInterface,HybridEnclosureType>& evolver);

    template<class HybridEnclosureType>
    HybridReachabilityAnalyser(const EvolverInterface<HybridAutomatonInterface,HybridEnclosureType>& evolver);

    /*! \brief Make a dynamically-allocated copy. */
    virtual HybridReachabilityAnalyser* clone() const { return new HybridReachabilityAnalyser(*this); }
    //@}

    //@{
    //! \name Methods to set and get the parameters controlling the accuracy.
    /*! \brief The parameters controlling the accuracy. */
    const EvolutionParametersType& parameters() const { return *this->_parameters; }
    /*! \brief A reference to the parameters controlling the accuracy. */
    EvolutionParametersType& parameters() { return *this->_parameters; }
    //@}


    //@{
    //! \name Evaluation of systems on abstract sets
    /*! \brief Compute a lower-approximation to the set obtained by evolving \a system for \a time starting in \a initial_set. */
    virtual SetApproximationType lower_evolve(const SystemType& system,
     const OvertSetInterfaceType& initial_set,
     const TimeType& time) const;

    /*! \brief Compute a lower-approximation to the reachable set of \a system starting in \a initial_set up to \a time. */
    virtual SetApproximationType
    lower_reach(const SystemType& system,
    const OvertSetInterfaceType& initial_set,
    const TimeType& time) const;

    /*! \brief Compute a lower-approximation to the reachable and evolved sets of \a system starting in \a initial_set up to \a time. */
    virtual std::pair<SetApproximationType,SetApproximationType>
    lower_reach_evolve(const SystemType& system,
     const OvertSetInterfaceType& initial_set,
     const TimeType& time) const;

    /*! \brief Compute an approximation to the set obtained by iterating \a time times \a system starting in \a initial_set. */
    virtual SetApproximationType upper_evolve(const SystemType& system,
     const CompactSetInterfaceType& initial_set,
     const TimeType& time) const;

    /*! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a time times. */
    virtual SetApproximationType upper_reach(const SystemType& system,
    const CompactSetInterfaceType& initial_set,
    const TimeType& timeType) const;

    /*! \brief Compute an approximation to the reachable and evolved sets of \a system starting in \a initial_set iterating at most \a time times. */
    virtual std::pair<SetApproximationType,SetApproximationType>
    upper_reach_evolve(const SystemType& system,
     const CompactSetInterfaceType& initial_set,
     const TimeType& time) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set. */
    virtual SetApproximationType chain_reach(const SystemType& system,
      const CompactSetInterfaceType& initial_set) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set and remaining in \a bounding_domain. \deprecated */
    virtual SetApproximationType chain_reach(const SystemType& system,
      const CompactSetInterfaceType& initial_set,
      const BoundingSetType& bounding_domain) const;

    /*! \brief Compute an outer-approximation to the viability kernel of \a system within \a bounding_set. */
    virtual SetApproximationType viable(const HybridAutomatonInterface& system,
    const HybridCompactSetInterface& bounding_set) const;

    /*! \brief Attempt to verify that the reachable set of \a system starting in \a initial_set remains in \a safe_set. */
    virtual tribool verify(const HybridAutomatonInterface& system,
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
    GTS _upper_reach(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
    GTS _upper_evolve(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
    std::pair<GTS,GTS> _upper_reach_evolve(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
  private:
    // Helper functions for approximating sets
};


template<class HybridEnclosureType>
HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const EvolverInterface<HybridAutomatonInterface,HybridEnclosureType>& evolver)
    : _parameters(new EvolutionParametersType())
    , _discretiser(new HybridDiscretiser<HybridEnclosureType>(evolver))
{
}


template<class HybridEnclosureType>
HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const EvolutionParametersType& parameters,
      const EvolverInterface<HybridAutomatonInterface,HybridEnclosureType>& evolver)
    : _parameters(new EvolutionParametersType(parameters))
    , _discretiser(new HybridDiscretiser<HybridEnclosureType>(evolver))
{
}


} // namespace Ariadne


#endif // ARIADNE_REACHABILITY_ANALYSER_H
