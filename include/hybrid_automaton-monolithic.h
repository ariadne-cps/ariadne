/***************************************************************************
 *            hybrid_automaton-monolithic.h
 *
 *  Copyright  2004-8  Alberto Casagrande, Pieter Collins
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

/*! \file hybrid_automaton.h
 *  \brief Main hybrid system class.
 */

#ifndef ARIADNE_MONOLITHIC_HYBRID_AUTOMATON_H
#define ARIADNE_MONOLITHIC_HYBRID_AUTOMATON_H

#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "function.h"
#include "discrete_location.h"
#include "discrete_event.h"

#include "hybrid_automaton_interface.h"

namespace Ariadne {



class HybridTime;
class HybridSpace;

class MonolithicHybridAutomaton;



/*
 * //! \brief A discrete mode of a hybrid automaton, comprising continuous evolution given by a vector field
//! within and invariant constraint set.
//!
//! A %DiscreteMode can only be created using the new_mode() method in
//! the %MonolithicHybridAutomaton class.
//!
//! \sa MonolithicHybridAutomaton, DiscreteTransition
class DiscreteMode {
    friend class MonolithicHybridAutomaton;
  private:

    // The discrete mode's discrete state.
    DiscreteLocation _location;

    // The discrete mode's vector field.
    RealVectorFunction _dynamic;

    // The discrete mode's invariants.
    Map< DiscreteEvent, Pair<EventKind,RealScalarFunction> _invariants;

  public:
    //! \brief The mode's discrete state.
    DiscreteLocation location() const {
        return this->_location; }

    //! \brief The discrete mode's dynamic (a vector field).
    const RealVectorFunction& dynamic() const {
        return this->_dynamic; }

    //! \brief The discrete mode's invariants.
    const Map< DiscreteEvent, RealScalarFunction >& invariants() const {
        return this->_invariants; }

    const RealScalarFunction& invariant(const DiscreteEvent& event) const {
        return this->_invariants.find(event)->second; }

    //! \brief The dimension of the discrete mode.
    uint dimension() const;

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream& os) const;

  private:
    // Construct discrete mode.
    //
    // \param id is the identifier of the mode.
    // \param dynamic is the mode's vector field.
    // \param invariants is the mode's invariants.
    DiscreteMode(DiscreteLocation location,
                 const RealVectorFunction& dynamic);

};


std::ostream& operator<<(std::ostream& os, const DiscreteMode& dm);

*/



//! \ingroup SystemModule
//! \brief A hybrid automaton, comprising continuous-time behaviour
//! at each discrete mode, coupled by instantaneous discrete transitions.
//! The state space is given by a hybrid set.
//!
//! A hybrid automaton is a dynamic system with evolution in both
//! continuous time and discrete time.
//! The state space is a product \f$X=\bigcup\{q\}\times X_q\f$
//! where \f$q\f$ is the <em>discrete state</em> and \f$X_q\f$
//! is the <em>continuous state space</em> of corresponding to
//! each discrete state.
//!
//! For each %DiscreteMode, the dynamics is given by a
//! %VectorField describing the continuous dynamics,
//! and a %Set giving an invariants which must be satisified at
//! all times.
//!
//! The discrete time behaviour is specified by %DiscreteTransition
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
//! \sa \link Ariadne::DiscreteMode \c DiscreteMode \endlink, \link Ariadne::DiscreteTransition \c DiscreteTransition \endlink
class MonolithicHybridAutomaton
    : public HybridAutomatonInterface
{
    struct Invariant {
        Invariant(DiscreteLocation q, DiscreteEvent e, RealScalarFunction g, EventKind k)
            : _location(q), _event(e), _guard(g), _kind(k) { }
        DiscreteLocation _location;
        DiscreteEvent _event;
        RealScalarFunction _guard;
        EventKind _kind;
    };
    struct Transition {
        Transition(DiscreteLocation s, DiscreteEvent e, DiscreteLocation t, RealVectorFunction r, RealScalarFunction g, EventKind k)
            : _source(s), _event(e), _target(t), _reset(r), _guard(g), _kind(k) { }
        DiscreteLocation _source;
        DiscreteEvent _event;
        DiscreteLocation _target;
        RealVectorFunction _reset;
        RealScalarFunction _guard;
        EventKind _kind;
    };
    struct Mode {
        Mode(DiscreteLocation q, RealSpace s, RealVectorFunction f);
        DiscreteLocation _location;
        List<Identifier> _variable_names;
        RealVectorFunction _dynamic;
        Map< DiscreteEvent, Invariant >  _invariants;
        Map< DiscreteEvent, Transition >  _transitions;
        uint dimension() const { return _dynamic.result_size(); }
    };
    friend std::ostream& operator<<(std::ostream&, const MonolithicHybridAutomaton&);
  public:
    //! \brief The type used to represent time.
    typedef HybridTime TimeType;
    //! \brief The type used to describe the state space.
    typedef HybridSpace StateSpaceType;
  private:
    //! \brief The list of the hybrid automaton's discrete modes.
    Map< DiscreteLocation, Mode > _modes;

  public:
    //@{
    //! \name Constructors and destructors

    //! \brief Construct an empty automaton with no name
    MonolithicHybridAutomaton();

    //! \brief Construct an empty automaton with the given name (Deprecated)
    MonolithicHybridAutomaton(const std::string& name);

    //! \brief  Destructor.
    ~MonolithicHybridAutomaton();
    //@}

    //@{
    //! \name Methods for building the automaton.

    //! \brief Adds a discrete mode to the automaton.
    //!
    //!   \param state is the mode's discrete state.
    //!   \param dynamic is the mode's vector field.
    //!
    //! The variables are given default names x0, x1 etc.
    void new_mode(DiscreteLocation state,
                  RealVectorFunction dynamic);

    //! \brief Adds a discrete mode to the automaton.
    //!
    //!   \param state is the mode's discrete state.
    //!   \param space is a list of state variables.
    //!   \param dynamic is the mode's vector field.
    void new_mode(DiscreteLocation state,
                  RealSpace space,
                  RealVectorFunction dynamic);

    //! \brief Adds an invariant to a mode of the automaton.
    //!
    //!   \param state is the mode's discrete state.
    //!   \param label is a discrete event labelling the invariant;
    //!            must be different from all transition events
    //!   \param invariant is the new invariant condition, in the form \f$g(x)<0\f$.
    //!   \param kind determines whether the constraint is a true invariant, or a progress predicate.

    void new_invariant(DiscreteLocation state,
                       DiscreteEvent label,
                       RealScalarFunction invariant,
                       EventKind kind=PROGRESS
                      );


    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //!
    //!    \param source is the transition's source location.
    //!    \param event is the transition's event.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param guard is the transition's guard function. The transition may occur when \f$g(x)\geq0\f$.
    //!    \param kind determines whether the transision is an impact, urgent or permissive.
    void new_transition(DiscreteLocation source,
                        DiscreteEvent event,
                        DiscreteLocation target,
                        RealVectorFunction reset,
                        RealScalarFunction guard,
                        EventKind kind);

    //@}

    //@{
    //! \name Data access and queries.

    //! \brief Test if the hybrid automaton has a discrete mode \a location.
    bool has_mode(DiscreteLocation location) const;

    //! \brief Test if the hybrid automaton has an action (corresponding to an invariant or guard) in \a location labelled by \a event.
    bool has_guard(DiscreteLocation location, DiscreteEvent event) const;

    //! \brief Test if the hybrid automaton has an invariant in \a location labelled by \a event.
    bool has_invariant(DiscreteLocation location, DiscreteEvent event) const;

    //! \brief Test if the hybrid automaton has a discrete transition with \a event_id and \a source_id.
    bool has_transition(DiscreteLocation source, DiscreteEvent event) const;

    //! \brief The dimension of the state space of in the mode \a location.
    uint dimension(DiscreteLocation location) const;

    //! \brief The dynamic valid in the mode \a location.
    virtual RealVectorFunction dynamic_function(DiscreteLocation location) const;

    //! \brief The set of all events possible in the given \a location.
    virtual Set<DiscreteEvent> events(DiscreteLocation location) const;

    //! \brief The kind (permissive, urgent etc) of the event.
    virtual EventKind event_kind(DiscreteLocation location, DiscreteEvent event) const;

    //! \brief The constraint function defining the invariant or time-can-progress predicate \f$p(x)\leq0\f$.
    virtual RealScalarFunction invariant_function(DiscreteLocation location, DiscreteEvent event) const;

    //! \brief The constraint function defining the condition \f$c(x)\geq0\f$ under which a transition occurs.
    virtual RealScalarFunction guard_function(DiscreteLocation location, DiscreteEvent event) const;

    //! \brief The target location of \a event starting in the \a source location.
    virtual DiscreteLocation target(DiscreteLocation source, DiscreteEvent event) const;

    //! \brief The dynamic valid in the mode \a location.
    virtual RealVectorFunction reset_function(DiscreteLocation location, DiscreteEvent event) const;

    //! \brief The discrete mode for the given \a location.
    const Mode& mode(DiscreteLocation source) const;

    //! \brief The invariant in the given \a location with the given \a event label.
    const Invariant& invariant(DiscreteLocation location, DiscreteEvent event) const;

    //! \brief The discrete transition with given \a source location and \a event.
    const Transition& transition(DiscreteLocation source, DiscreteEvent event) const;

    //! \brief The set of all locations.
    Set<DiscreteLocation> locations() const;

    //! \brief A hybrid space, comprising a continuous space for each (reachable) location of the automaton.
    virtual HybridSpace state_space() const;

    //! \brief The continuous state space in the given \a location.
    virtual RealSpace continuous_state_space(DiscreteLocation location) const;


    //@}

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream&) const;

};

inline std::ostream& operator<<(std::ostream& os, const MonolithicHybridAutomaton& ha) {
    return ha.write(os);
}


} // namespace Ariadne

#endif // ARIADNE_HYBRID_AUTOMATON_H
