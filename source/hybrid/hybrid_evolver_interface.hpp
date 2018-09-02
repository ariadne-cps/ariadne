/***************************************************************************
 *            hybrid_evolver_interface.hpp
 *
 *  Copyright  2011  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file hybrid_evolver_interface.hpp
 *  \brief Interface for evolver classes in the hybrid space.
 */

#ifndef ARIADNE_HYBRID_EVOLVER_INTERFACE_HPP
#define ARIADNE_HYBRID_EVOLVER_INTERFACE_HPP

#include "../dynamics/evolver_interface.hpp"

#include "../hybrid/hybrid_time.hpp"
#include "../hybrid/hybrid_set.decl.hpp"
#include "../hybrid/hybrid_orbit.hpp"

#include "../hybrid/discrete_event.hpp"

#include "../output/logging.hpp"

namespace Ariadne {

class HybridAutomatonInterface;
class HybridEnclosure;

class HybridTerminationCriterion {
  public:
    typedef HybridTime::ContinuousTimeType ContinuousTimeType;
    typedef HybridTime::DiscreteTimeType DiscreteTimeType;
  private:
    ContinuousTimeType _maximum_time;
    DiscreteTimeType _maximum_steps;
    Set<DiscreteEvent> _terminating_events;
  public:
    HybridTerminationCriterion(ContinuousTimeType tmax, DiscreteTimeType nmax, Set<DiscreteEvent> evnts)
        : _maximum_time(tmax), _maximum_steps(nmax), _terminating_events(evnts) { }
    HybridTerminationCriterion(ContinuousTimeType tmax, DiscreteTimeType nmax)
        : HybridTerminationCriterion(tmax,nmax,Set<DiscreteEvent>()) { }
    HybridTerminationCriterion(const HybridTime& maximum_time)
        : HybridTerminationCriterion(maximum_time.continuous_time(),maximum_time.discrete_time()) { }
    //! \brief The maximum continuous (real, physical) time.
    const ContinuousTimeType& maximum_time() const { return this->_maximum_time; }
    //! \brief The maximum number of discrete steps taken.
    const DiscreteTimeType& maximum_steps() const { return this->_maximum_steps; }
    //! \brief The maximum number of discrete steps taken.
    const Set<DiscreteEvent>& terminating_events() const { return this->_terminating_events; }
};
OutputStream& operator<<(OutputStream& os, const HybridTerminationCriterion& termination);

//! \brief Interface for hybrid evolvers using HybridEnclosure as the enclosure type.
//! \details The class is loggable in order to allow verbosity tuning at the analyser layer.
class HybridEvolverInterface
    : public EvolverInterface<HybridAutomatonInterface,HybridEnclosure,HybridTerminationCriterion>
    , public Loggable
{
  public:
    //! \brief Make a dynamically-allocated copy.
    virtual HybridEvolverInterface* clone() const = 0;

    //@{
    //! \name Main evolution functions.

    //! \brief Compute an approximation to the orbit set using the given semantics, starting from an initial enclosure.
    //!   Useful for continuing a partially-computed orbit.
    virtual Orbit<EnclosureType> orbit(const EnclosureType& initial_enclosure,const TerminationType& termination,Semantics semantics) const = 0;

    //! \brief Compute an approximation to the orbit set using the given semantics, starting from a box.
    //!   Useful for computing the evolution starting from a cell of a grid.
    virtual Orbit<EnclosureType> orbit(const HybridExactBoxType& initial_box,const TerminationType& termination,Semantics semantics) const = 0;
    //! \brief Compute an approximation to the orbit set using the given semantics, starting from a set described by bounds and constraints.
    //!   Useful for computing the evolution starting from user-provided set.
    virtual Orbit<EnclosureType> orbit(const HybridBoundedConstraintSet& initial_set,const TerminationType& termination,Semantics semantics) const = 0;

    //! \brief Compute an approximation to the evolution set under the given semantics.
    virtual Pair<EnclosureListType,EnclosureListType> reach_evolve(const EnclosureType& initial_set, const TerminationType& termination, Semantics semantics=Semantics::UPPER) const = 0;
    //@}

    //@{
    //! \name Auxiliary set conversion functionality

    //! \brief Set construct an enclosure from a box, such as one obtained from a grid.
    virtual EnclosureType enclosure(const HybridExactBoxType& initial_box) const = 0;
    //! \brief Set construct an enclosure from a user-provided set.
    virtual EnclosureType enclosure(const HybridBoundedConstraintSet& initial_set) const = 0;

    //@}

};


//! \brief Factory for hybrid evolver interface classes.
class HybridEvolverFactoryInterface
    : public EvolverFactoryInterface<HybridAutomatonInterface,HybridEnclosure,HybridTerminationCriterion>
{
  public:

    //! \brief Create a copy of the factory.
    virtual HybridEvolverFactoryInterface* clone() const = 0;
    //! \brief Create a hybrid evolver interface object around a \a system.
    virtual HybridEvolverInterface* create(const HybridAutomatonInterface& system) const = 0;
};


} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLVER_INTERFACE_HPP
