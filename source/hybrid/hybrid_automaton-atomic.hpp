/***************************************************************************
 *            hybrid_automaton-atomic.hpp
 *
 *  Copyright  2004-11  Alberto Casagrande, Pieter Collins
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

/*! \file hybrid_automaton-atomic.hpp
 *  \brief Atomic hybrid system class, where each location is defined by a single string variable.
 */

#ifndef ARIADNE_ATOMIC_HYBRID_AUTOMATON_HPP
#define ARIADNE_ATOMIC_HYBRID_AUTOMATON_HPP

#include "../hybrid/hybrid_automaton_interface.hpp"
#include "../hybrid/hybrid_automaton-composite.hpp"

namespace Ariadne {

class AtomicHybridAutomaton;

typedef StringVariable AtomicDiscreteVariable;

class AtomicDiscreteLocation : public StringConstant {
    using StringConstant::StringConstant;
};

//! \ingroup SystemModule
//! \brief A hybrid automaton, comprising continuous-time behaviour
//! at each discrete mode, coupled by instantaneous discrete transitions.
//!  The state space is given by a hybrid set.
//!
//! A hybrid automaton is a dynamic system with evolution in both
//! continuous time and discrete time.
//! The state space is a product \f$X=\bigcup\{q\}\times X_q\f$
//! where \f$q\f$ is the <em>discrete state</em> and \f$X_q\f$
//! is the <em>continuous state space</em> of corresponding to
//! each discrete state.
//!
//! For each %AtomicDiscreteMode, the dynamics is given by a
//! %VectorField describing the continuous dynamics,
//! and a %Set giving an invariants which must be satisified at
//! all times.
//!
//! The discrete time behaviour is specified by %AtomicDiscreteTransition
//! objects.
//! Each discrete transition represents an jump from a \a source
//! mode to a \a target mode.
//! There can be at most one discrete transition in an automaton
//! with the same event and source.
//!
//! A discrete transision can either be \em forced or \em unforced.
//! A forced transition much occur as soon as it is activated.
//! An unforced transition may occur at any time it is activated,
//! but is only forced to occur if the continuous evolution is
//! blocked by an invariant.
//!
//! \sa \ref AtomicDiscreteMode, \ref AtomicDiscreteTransition, \ref CompositeHybridAutomaton
class AtomicHybridAutomaton
    : public HybridAutomaton
{
  public:
    //! \brief The type used to represent time.
    typedef HybridTime TimeType;
    //! \brief The type used to represent real numbers.
    typedef double RealType ;
    //! \brief The type used to describe the state space.
    typedef HybridSpace StateSpaceType;

  private:
    StringVariable _variable;
  public:
    //@{
    //! \name Constructors and destructors

    //! \brief Construct an empty automaton with the given name
    AtomicHybridAutomaton(const Identifier& name) : HybridAutomaton(name), _variable(name) { }

    //! \brief Construct dynamically-allocated copy. (Not currently implemented)
    AtomicHybridAutomaton* clone() const { return new AtomicHybridAutomaton(*this); }

    //! \brief  Destructor.
    ~AtomicHybridAutomaton() { };
    //@}

    //@{
    //! \name Methods for building the automaton.

    //! \brief Adds a discrete mode to the automaton.
    Void new_mode(AtomicDiscreteLocation location,
                  List<RealAssignment> const& auxiliary,
                  List<DottedRealAssignment> const& dynamic) {
        this->HybridAutomaton::new_mode(this->_variable|location,auxiliary,dynamic);
    }

    //! \brief Adds a discrete mode to the automaton.
    Void new_mode(AtomicDiscreteLocation location,
                  List<DottedRealAssignment> const& dynamic) {
        this->HybridAutomaton::new_mode(this->_variable|location,dynamic);
    }

    //! \brief Adds a discrete mode to the automaton.
    Void new_mode(AtomicDiscreteLocation location,
                  List<RealAssignment> const& auxiliary) {
        this->HybridAutomaton::new_mode(this->_variable|location,auxiliary);
    }


    //! \brief Adds a new internal/output event with a given enabling \a activation condition and triggering \a invariant.
    Void new_action(AtomicDiscreteLocation location,
                    ContinuousPredicate invariant,
                    DiscreteEvent event,
                    ContinuousPredicate activation,
                    EventKind kind=EventKind::PERMISSIVE) {
        ARIADNE_ASSERT(kind==EventKind::PERMISSIVE);
        this->HybridAutomaton::new_action(this->_variable|location,invariant,event,activation,kind);
    }

    //! \brief Adds a new urgent internal/output event with a given enabling \a guard.
    Void new_action(AtomicDiscreteLocation location,
                    DiscreteEvent event,
                    ContinuousPredicate guard,
                    EventKind kind=EventKind::URGENT) {
        ARIADNE_ASSERT(kind==EventKind::URGENT || kind==EventKind::IMPACT);
        this->HybridAutomaton::new_guard(this->_variable|location,event,guard,kind);
    }

    //! \brief Adds a discrete mode to the automaton.
    Void new_invariant(AtomicDiscreteLocation location,
                       ContinuousPredicate const& invariant,
                       DiscreteEvent event) {
        this->HybridAutomaton::new_invariant(this->_variable|location,invariant,event);
    }

    //! \brief Adds a guard to the automaton.
    Void new_guard(AtomicDiscreteLocation location,
                   DiscreteEvent event,
                   ContinuousPredicate const& guard,
                   EventKind kind) {
        this->HybridAutomaton::new_guard(this->_variable|location,event,guard,kind);
    }

    //! \brief Adds an urgent guard to the automaton.
    Void new_urgent_guard(AtomicDiscreteLocation location,
                          DiscreteEvent event,
                          ContinuousPredicate const& guard) {
        this->HybridAutomaton::new_guard(this->_variable|location,event,guard,EventKind::URGENT);
    }

    //! \brief Adds a permissive guard to the automaton.
    Void new_permissive_guard(AtomicDiscreteLocation location,
                              DiscreteEvent event,
                              ContinuousPredicate const& guard) {
        this->HybridAutomaton::new_guard(this->_variable|location,event,guard,EventKind::PERMISSIVE);
    }

    //! \brief Adds an update/reset to the automaton.
    Void new_update(AtomicDiscreteLocation source,
                    DiscreteEvent event,
                    AtomicDiscreteLocation target,
                    List<PrimedRealAssignment> const& reset) {
        this->HybridAutomaton::new_update(this->_variable|source,event,this->_variable|target,reset);
    }

    //! \brief Adds a reset to the automaton. (Same as new_transition.)
    Void new_reset(AtomicDiscreteLocation source,
                   DiscreteEvent event,
                   AtomicDiscreteLocation target,
                   List<PrimedRealAssignment> const& reset) {
        this->HybridAutomaton::new_reset(this->_variable|source,event,this->_variable|target,reset);
    }


    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    Void new_transition(AtomicDiscreteLocation source,
                        DiscreteEvent event,
                        ContinuousPredicate const& guard,
                        AtomicDiscreteLocation target,
                        List<PrimedRealAssignment> const& reset,
                        EventKind urgency=EventKind::URGENT) {
        this->HybridAutomaton::new_transition(this->_variable|source,event,this->_variable|target,reset,guard,urgency);
    }

    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    Void new_transition(AtomicDiscreteLocation source,
                        DiscreteEvent event,
                        AtomicDiscreteLocation target,
                        List<PrimedRealAssignment> const& reset,
                        ContinuousPredicate const& guard,
                        EventKind urgency=EventKind::URGENT) {
        this->HybridAutomaton::new_transition(this->_variable|source,event,this->_variable|target,reset,guard,urgency);
    }

    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //! The reset is trivial. This form is for the case that there are no continuous state variables in the new location.
    Void new_transition(AtomicDiscreteLocation source,
                        DiscreteEvent event,
                        AtomicDiscreteLocation target,
                        ContinuousPredicate const& guard,
                        EventKind urgency=EventKind::URGENT) {
        this->HybridAutomaton::new_transition(this->_variable|source,event,this->_variable|target,guard,urgency);
    };

    //! \brief Adds an unguarded transition to the automaton.
    //! The guard is the constant "True" i.e. the event is an input event.
    Void new_transition(AtomicDiscreteLocation source,
                        DiscreteEvent event,
                        AtomicDiscreteLocation target,
                        List<PrimedRealAssignment> const& reset) {
        this->HybridAutomaton::new_transition(this->_variable|source,event,this->_variable|target,reset);
    }

    //! \brief Adds an unguarded transition to the automaton.
    //! The guard is the constant "True" i.e. the event is an input event.
    Void new_transition(AtomicDiscreteLocation source,
                        DiscreteEvent event,
                        AtomicDiscreteLocation target) {
        this->HybridAutomaton::new_transition(this->_variable|source,event,this->_variable|target);
    }


    //@}

    //! \brief The variable used to define the discrete component.
    const StringVariable& variable() const {
        return this->_variable; }

    //! \brief The discrete mode with given discrete state.
    const DiscreteMode& mode(AtomicDiscreteLocation location) const {
        return this->HybridAutomaton::mode(this->_variable|location); }

    //! \brief The target for \a event from location \a source. Returns \a source if \a event is not present.
    AtomicDiscreteLocation get_target(AtomicDiscreteLocation source, DiscreteEvent event) const {
        DiscreteLocation source_mode(this->_variable|source);
        DiscreteLocation target_mode=this->HybridAutomaton::target(source_mode,event);
        return AtomicDiscreteLocation(target_mode[this->_variable]); }

    //@}

    //! \brief Write to an output stream.
    OutputStream& write(OutputStream& os) const {
        return this->HybridAutomaton::write(os); }
};

inline Pair<StringVariable,String> operator|(const AtomicHybridAutomaton& ha, const AtomicDiscreteLocation& q) {
    return ha.variable() | q;
}

inline OutputStream& operator<<(OutputStream& os, const AtomicHybridAutomaton& ha) {
    return ha.write(os);
}


} // namespace Ariadne

#endif // ARIADNE_ATOMIC_HYBRID_AUTOMATON_HPP
