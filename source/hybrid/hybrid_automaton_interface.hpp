/***************************************************************************
 *            hybrid_automaton_interface.hpp
 *
 *  Copyright  2010  Pieter Collins
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

/*! \file hybrid_automaton_interface.hpp
 *  \brief Interface for hybrid automata to be used by evolver classes.
 */

#ifndef ARIADNE_HYBRID_AUTOMATON_INTERFACE_HPP
#define ARIADNE_HYBRID_AUTOMATON_INTERFACE_HPP

#include <cassert>
#include "../function/function.hpp"
#include "../hybrid/discrete_event.hpp"
#include "../hybrid/discrete_location.hpp"

namespace Ariadne {


class HybridTime;
class HybridSpace;
class HybridGrid;
class DiscreteEvent;

class DiscreteLocation;

template<class T> class Space;
typedef Space<Real> RealSpace;

enum class EventKind : std::uint8_t { INVARIANT, PROGRESS, PERMISSIVE, URGENT, IMPACT };
inline OutputStream& operator<<(OutputStream&, const EventKind& evk);

class HybridEvolverInterface;
class HybridEnclosure;

class SystemSpecificationError : public std::runtime_error {
  public:
    SystemSpecificationError(const StringType& what) : std::runtime_error(what) { }
};

class NonExistentModeError : public SystemSpecificationError {
  public:
    NonExistentModeError(const StringType& what) : SystemSpecificationError(what) { }
};

class IndistinguishableModeError : public SystemSpecificationError {
  public:
    IndistinguishableModeError(const StringType& what) : SystemSpecificationError(what) { }
};

class DuplicateEventError : public SystemSpecificationError {
  public:
    DuplicateEventError(const StringType& what) : SystemSpecificationError(what) { }
};

class MultipleGuardError : public SystemSpecificationError {
  public:
    MultipleGuardError(const StringType& what) : SystemSpecificationError(what) { }
};

class MultipleTransitionError : public SystemSpecificationError {
  public:
    MultipleTransitionError(const StringType& what) : SystemSpecificationError(what) { }
};

class AlgebraicLoopError : public SystemSpecificationError {
  public:
    AlgebraicLoopError(const StringType& what) : SystemSpecificationError(what) { }
};

class OverspecifiedSystemError : public SystemSpecificationError {
  public:
    OverspecifiedSystemError(const StringType& what) : SystemSpecificationError(what) { }
};

class OverspecifiedDynamicError : public OverspecifiedSystemError {
  public:
    OverspecifiedDynamicError(const StringType& what) : OverspecifiedSystemError(what) { }
};

class OverspecifiedResetError : public OverspecifiedSystemError {
  public:
    OverspecifiedResetError(const StringType& what) : OverspecifiedSystemError(what) { }
};

class UnderspecifiedSystemError : public SystemSpecificationError {
  public:
    UnderspecifiedSystemError(const StringType& what) : SystemSpecificationError(what) { }
};

class UnderspecifiedDynamicError : public UnderspecifiedSystemError {
  public:
    UnderspecifiedDynamicError(const StringType& what) : UnderspecifiedSystemError(what) { }
};

class UnderspecifiedConstraintError : public UnderspecifiedSystemError {
  public:
    UnderspecifiedConstraintError(const StringType& what) : UnderspecifiedSystemError(what) { }
};

class UnderspecifiedResetError : public UnderspecifiedSystemError {
  public:
    UnderspecifiedResetError(const StringType& what) : UnderspecifiedSystemError(what) { }
};


class StateSpecificationError : public std::runtime_error {
  public:
    StateSpecificationError(const StringType& what) : std::runtime_error(what) { }
};

class IndistinguishableLocationsError : public StateSpecificationError {
  public:
    IndistinguishableLocationsError(const StringType& what) : StateSpecificationError(what) { }
};

class OverspecifiedLocationException : public IndistinguishableLocationsError {
  public:
    OverspecifiedLocationException(const StringType& what) : IndistinguishableLocationsError(what) { }
};


//! \ingroup SystemModule
//! \brief Base interface for hybrid systems, to allow different types to be used in evolution routines.
class HybridAutomatonInterface {
  public:
    //! \brief The type used to represent time.
    typedef HybridTime TimeType;
    //! \brief The type used to describe the state space.
    typedef HybridSpace StateSpaceType;
    //! \brief The type used to describe the state space.
    typedef HybridEvolverInterface EvolverType;
    typedef HybridEnclosure EnclosureType;

  public:
    //! \brief Virtual destructor.
    virtual ~HybridAutomatonInterface() = default;

    //! \brief Cloning operator.
    virtual HybridAutomatonInterface* clone() const = 0;

    //@{
    //! \name Data access and queries.

    //! \brief The name of the automaton.
    virtual const Identifier& name() const = 0;

    //@{
    //! \name Methods for testing and extracting the discrete dynamics.

    //! \brief Test if the hybrid automaton has a valid discrete mode \a location.
    virtual Bool has_mode(DiscreteLocation location) const = 0;

    //! \brief The set of all events possible in the given \a location.
    virtual Set<DiscreteEvent> events(DiscreteLocation location) const = 0;

    //! \brief The kind (permissive, urgent etc) of the event.
    virtual EventKind event_kind(DiscreteLocation location, DiscreteEvent event) const = 0;

    //! \brief Test if the hybrid automaton has an invariant or urgent guard constraint in the \a location labelled by \a event.
    virtual Bool has_invariant(DiscreteLocation location, DiscreteEvent event) const = 0;

    //! \brief Test if the hybrid automaton has a guard constraint in the \a location labelled by \a event.
    virtual Bool has_guard(DiscreteLocation location, DiscreteEvent event) const = 0;

    //! \brief Test if the hybrid automaton has a discrete transition in \a source due to \a event.
    virtual Bool has_transition(DiscreteLocation source, DiscreteEvent event) const = 0;

    //! \brief The target location of \a event starting in the \a source location.
    virtual DiscreteLocation target(DiscreteLocation source, DiscreteEvent event) const = 0;

    //@}

    //@{
    //! \name Methods for extracting the continuous dynamics.

    //! \brief The dimension of the state space in the given \a location.
    virtual DimensionType dimension(DiscreteLocation location) const = 0;

    //! \brief The function giving the auxiliary variables in the mode \a location.
    virtual EffectiveVectorMultivariateFunction auxiliary_function(DiscreteLocation location) const = 0;

    //! \brief The dynamic valid in the mode \a location.
    virtual EffectiveVectorMultivariateFunction dynamic_function(DiscreteLocation location) const = 0;

    //! \brief The constraint function defining the invariant or time-can-progress predicate \f$p(x)\leq0\f$.
    virtual EffectiveScalarMultivariateFunction invariant_function(DiscreteLocation location, DiscreteEvent event) const = 0;

    //! \brief The constraint function defining the condition \f$c(x)\geq0\f$ under which a transition occurs.
    virtual EffectiveScalarMultivariateFunction guard_function(DiscreteLocation location, DiscreteEvent event) const = 0;

    //! \brief The dynamic valid in the mode \a location.
    virtual EffectiveVectorMultivariateFunction reset_function(DiscreteLocation location, DiscreteEvent event) const = 0;

    //! \brief The hybrid state space.
    virtual HybridSpace state_space() const = 0;

    //! \brief The hybrid space of state and auxiliary variables.
    virtual HybridSpace state_auxiliary_space() const = 0;

    //! \brief The continuous state space in the \a location.
    virtual RealSpace continuous_state_space(DiscreteLocation location) const = 0;

    //! \brief The continuous state space in the \a location.
    virtual RealSpace continuous_auxiliary_space(DiscreteLocation location) const = 0;

    //! \brief The continuous state space in the \a location.
    virtual List<RealExpression> auxiliary_expressions(DiscreteLocation location) const { assert(false); return List<RealExpression>(); }


    //@}

    //@{
    //! \name Input/output methods

    //! \brief Write to an output stream.
    virtual OutputStream& write(OutputStream& os) const = 0;
    //@}

};

inline OutputStream& operator<<(OutputStream& os, const HybridAutomatonInterface& ha) {
    ha.write(os);
    return os;
}

inline OutputStream& operator<<(OutputStream& os, const EventKind& evk) {
    switch(evk) {
        case EventKind::INVARIANT: os<<"invariant"; break;
        case EventKind::PROGRESS: os<<"progress"; break;
        case EventKind::PERMISSIVE: os<<"permissive"; break;
        case EventKind::URGENT: os<<"urgent"; break;
        case EventKind::IMPACT: os<<"impact"; break;
        default: ARIADNE_FAIL_MSG("Unhandled EventKind for output streaming.\n");
    } return os;
}

} //namespace Ariadne

#endif // ARIADNE_HYBRID_AUTOMATON_INTERFACE_HPP
