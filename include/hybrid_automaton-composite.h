/***************************************************************************
 *            hybrid_automaton-composite.h
 *
 *  Copyright  2004-11  Alberto Casagrande, Pieter Collins
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

/*! \file hybrid_automaton-composite.h
 *  \brief Main composite hybrid system class.
 */

#ifndef ARIADNE_COMPOSITE_HYBRID_AUTOMATON_H
#define ARIADNE_COMPOSITE_HYBRID_AUTOMATON_H

#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "function.h"
#include "discrete_location.h"
#include "discrete_event.h"
#include "assignment.h"
#include "expression.h"
#include "logging.h"

#include "hybrid_automaton_interface.h"

namespace Ariadne {

class HybridTime;
class HybridSpace;
class HybridSet;

class DiscreteLocation;
class DiscreteMode;
class DiscreteTransition;
class HybridAutomaton;
class CompositeHybridAutomaton;


//! \related HybridAutomaton
//! \brief A discrete transition of a hybrid automaton, representing an instantaneous
//! jump from one discrete mode to another, governed by an activation set and a reset map.
//!
//! A %DiscreteTransition can only be created using the new_transition() method in
//! the %HybridAutomaton class.
//!
//! An invariant is modelled by a discrete transition with negative event id and null reset pointer.
//!
//! \sa \ref HybridAutomaton, \ref DiscreteMode
class DiscreteTransition
{
    friend class DiscreteMode;
    friend class HybridAutomaton;
    friend class CompositeHybridAutomaton;
  private:
    //  The source of the discrete transition.
    DiscreteLocation _source;

    //  The discrete transition's identificator.
    DiscreteEvent _event;

    //  The target of the discrete transition.
    DiscreteLocation _target;

    //  The guard predicate of the discrete transition.
    ContinuousPredicate _guard;

    //  The reset assignments of the discrete transition.
    List<PrimedRealAssignment> _reset;

    //  The kind of the event.
    EventKind _kind;

  public:

    //! \brief The source mode of the discrete transition.
    DiscreteLocation source() const { return this->_source; }

    //! \brief The discrete event associated with the discrete transition.
    DiscreteEvent event() const { return this->_event; }

    //! \brief The target mode of the discrete transition.
    DiscreteLocation target() const { return this->_target; };

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream& os) const;
  private:
    DiscreteTransition(DiscreteLocation source,
                       DiscreteEvent event,
                       DiscreteLocation target,
                       const List<PrimedRealAssignment>& reset,
                       const ContinuousPredicate& guard,
                       EventKind kind);
};

inline std::ostream& operator<<(std::ostream& os, const DiscreteTransition& dt) {
    return dt.write(os); }

inline bool operator<(const DiscreteTransition& transition1, const DiscreteTransition& transition2) {
    return transition1.event() < transition2.event()
        || (transition1.event() == transition2.event()
            && transition1.source() < transition2.source());
}

//! \related HybridAutomaton
//! \brief A discrete mode of a hybrid automaton, comprising continuous evolution given by a vector field
//! within and invariant constraint set.
//!
//! A %DiscreteMode can only be created using the new_mode() method in
//! the %HybridAutomaton class.
//!
//! \sa \ref HybridAutomaton, \ref DiscreteTransition
class DiscreteMode {
    friend class DiscreteTransition;
    friend class HybridAutomaton;
    friend class CompositeHybridAutomaton;
  private:
    // The discrete mode's discrete state.
    DiscreteLocation _location;

    // The algebraic equations
    List<RealAssignment> _auxiliary;

    // The algebraic equations
    List<DottedRealAssignment> _dynamic;

    // The discrete mode's invariants.
    Map<DiscreteEvent,ContinuousPredicate> _invariants;

    // The discrete mode's urgent guards.
    Map<DiscreteEvent, ContinuousPredicate> _guards;

    // The kind of each event
    Map<DiscreteEvent, EventKind> _kinds;

    // The target locations for the discrete transitions.
    Map<DiscreteEvent,DiscreteLocation> _targets;

    // The update assignments for the doscrete transitions.
    Map<DiscreteEvent, List<PrimedRealAssignment> > _resets;
  private:
    DiscreteMode(DiscreteLocation, List<RealAssignment>const&, List<DottedRealAssignment>const&);
    DiscreteMode(DiscreteLocation);
    DiscreteMode();
  public:
    //! \brief The mode's discrete state.
    DiscreteLocation location() const { return this->_location; }

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream& os) const;
};


inline std::ostream& operator<<(std::ostream& os, const DiscreteMode& dm) {
    return dm.write(os); }

inline bool operator<(const DiscreteMode& mode1, const DiscreteMode& mode2) {
    return mode1.location() < mode2.location(); }







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
//! \sa \ref DiscreteMode, \ref DiscreteTransition, \ref CompositeHybridAutomaton
class HybridAutomaton
    : public HybridAutomatonInterface
    , public Loggable
{
    friend class CompositeHybridAutomaton;
  public:
    //! \brief The type used to represent time.
    typedef HybridTime TimeType;
    //! \brief The type used to represent real numbers.
    typedef double RealType ;
    //! \brief The type used to describe the state space.
    typedef HybridSpace StateSpaceType;

    static List<RealAssignment> sort(const List<RealAssignment>& auxiliary);

    //! \brief The list of the hybrid automaton's discrete modes.
    Map< DiscreteLocation, DiscreteMode > _modes;

  protected:
    void _new_mode(DiscreteLocation location,
                   List<RealAssignment> const& auxilary,
                   List<DottedRealAssignment> const& dynamic);

    void _new_invariant_(DiscreteLocation location,
                         ContinuousPredicate invariant,
                         DiscreteEvent event);

    void _new_guard_(DiscreteLocation location,
                     DiscreteEvent event,
                     ContinuousPredicate guard,
                     EventKind kind);

    void _new_action(DiscreteLocation location,
                     ContinuousPredicate invariant,
                     DiscreteEvent event,
                     ContinuousPredicate guard,
                     EventKind kind);

    void _new_update(DiscreteLocation location,
                     DiscreteEvent event,
                     DiscreteLocation target,
                     List<PrimedRealAssignment> const& reset);

  public:
    //@{
    //! \name Constructors and destructors

    //! \brief Construct an empty automaton with the given name
    HybridAutomaton();

    //! \brief Construct an empty automaton with the given name
    HybridAutomaton(const List<StringVariable>& discrete_variables);

    //! \brief Construct dynamically-allocated copy. (Not currently implemented)
    HybridAutomaton* clone() const;

    //! \brief  Destructor.
    ~HybridAutomaton();
    //@}

    //@{
    //! \name Methods for building the automaton.

    //! \brief Adds a discrete mode to the automaton.
    void new_mode(DiscreteLocation location,
                  List<RealAssignment> const& auxiliary,
                  List<DottedRealAssignment> const& dynamic) {
        this->_new_mode(location,auxiliary,dynamic);
    }

    //! \brief Adds a discrete mode to the automaton.
    void new_mode(DiscreteLocation location,
                  List<DottedRealAssignment> const& dynamic) {
        this->_new_mode(location,List<RealAssignment>(),dynamic);
    }

    //! \brief Adds a discrete mode to the automaton.
    void new_mode(DiscreteLocation location,
                  List<RealAssignment> const& auxiliary) {
        this->_new_mode(location,auxiliary,List<DottedRealAssignment>());
    }

    //! \brief Adds a discrete mode to the automaton without any dynamics.
    void new_mode(DiscreteLocation location) {
        this->_new_mode(location,List<RealAssignment>(),List<DottedRealAssignment>());
    }

    //! \brief Adds a discrete mode to the automaton.
    void new_mode(List<RealAssignment> const& auxiliary, List<DottedRealAssignment> const& dynamic) {
        this->_new_mode(DiscreteLocation(),auxiliary,dynamic);
    }

    //! \brief Adds a discrete mode to the automaton.
    void new_mode(List<DottedRealAssignment> const& dynamic) {
        this->_new_mode(DiscreteLocation(),List<RealAssignment>(),dynamic);
    }

    //! \brief Adds a discrete mode to the automaton.
    void new_mode(List<RealAssignment> const& auxiliary) {
        this->_new_mode(DiscreteLocation(),auxiliary,List<DottedRealAssignment>());
    }

    //! \brief Adds a new internal/output event with a given enabling \a guard condition and triggering \a invariant.
    void new_action(DiscreteLocation location,
                    ContinuousPredicate invariant,
                    DiscreteEvent event,
                    ContinuousPredicate guard,
                    EventKind kind=PERMISSIVE) {
        this->_new_action(location,invariant,event,guard,kind);
    }

    //! \brief Adds a new internal/output event with a given enabling \a guard condition and triggering \a invariant.
    void new_action(DiscreteLocation location,
                    DiscreteEvent event,
                    ContinuousPredicate guard,
                    EventKind kind=URGENT) {
        ARIADNE_ASSERT(kind==URGENT);
        this->_new_action(location,!guard,event,guard,kind);
    }

    //! \brief Adds an invariant to the automaton.
    void new_invariant(DiscreteLocation location,
                       ContinuousPredicate const& invariant,
                       DiscreteEvent event) {
        //this->_new_action(location,event,invariant,ContinuousPredicate(tribool(false)),INVARIANT);
        this->_new_invariant_(location,invariant,event);
    }

    //! \brief Adds a guard to the automaton.
    void new_guard(DiscreteLocation location,
                   DiscreteEvent event,
                   ContinuousPredicate const& guard,
                   EventKind kind) {
        //this->_new_action(location,event,ContinuousPredicate(true),guard,kind);
        this->_new_guard_(location,event,guard,kind);
    }

    //! \brief Adds an update/reset to the automaton.
    void new_update(DiscreteLocation source,
                    DiscreteEvent event,
                    DiscreteLocation target,
                    List<PrimedRealAssignment> const& reset) {
        this->_new_update(source,event,target,reset);
    }

    //! \brief Adds a discrete update to the automaton to a mode with no continuous state variables.
    void new_update(DiscreteLocation source,
                    DiscreteEvent event,
                    DiscreteLocation target) {
        this->_new_update(source,event,target,List<PrimedRealAssignment>());
    }

    //! \brief Adds a reset to an automaton with a single mode.
    void new_update(DiscreteEvent event,
                    List<PrimedRealAssignment> const& reset) {
        this->_new_update(DiscreteLocation(),event,DiscreteLocation(),List<PrimedRealAssignment>());
    }

    //! \brief Adds a reset to the automaton. (Same as new_update.)
    void new_reset(DiscreteLocation source,
                   DiscreteEvent event,
                   DiscreteLocation target,
                   List<PrimedRealAssignment> const& reset) {
        this->_new_update(source,event,target,reset);
    }

    //! \brief Adds a reset to an automaton with a single mode.
    void new_reset(DiscreteEvent event,
                    List<PrimedRealAssignment> const& reset) {
        this->_new_update(DiscreteLocation(),event,DiscreteLocation(),List<PrimedRealAssignment>());
    }

    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //!
    //!    \param source is the transition's source location.
    //!    \param event is the transition's event.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param guard is the transition's activation region.
    //!    \param urgency is a flag indicating whether the transition is urgent i.e. occurs as soon as it is activated.
    void new_transition(DiscreteLocation source,
                        DiscreteEvent event,
                        DiscreteLocation target,
                        const List<PrimedRealAssignment>& reset,
                        const ContinuousPredicate& guard,
                        EventKind kind=urgent) {
        if(kind==urgent || kind==impact) { this->_new_action(source,!guard,event,guard,kind); }
        else if(kind==permissive) { this->_new_guard_(source,event,guard,kind); }
        else { ARIADNE_FAIL_MSG("Unhandled event kind "<<kind); }
        this->_new_update(source,event,target,reset);
    }

    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //! The reset is trivial. This form is for the case that there are no continuous state variables in the new location.
    void new_transition(DiscreteLocation source,
                        DiscreteEvent event,
                        DiscreteLocation target,
                        ContinuousPredicate const& guard,
                        EventKind kind=urgent) {
        this->new_transition(source,event,target,List<PrimedRealAssignment>(),guard,kind);
    }

    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //! The reset is trivial. This form is for the case that there are no continuous state variables in the new location.
    void new_transition(DiscreteLocation source,
                        DiscreteEvent event,
                        ContinuousPredicate const& guard,
                        DiscreteLocation target,
                        EventKind kind=urgent) {
        this->new_transition(source,event,target,List<PrimedRealAssignment>(),guard,kind);
    }

    //! \brief Adds an unguarded transition to the automaton.
    //! The guard is the constant "True" i.e. the event is an input event.
    void new_transition(DiscreteLocation source,
                        DiscreteEvent event,
                        DiscreteLocation target,
                        List<PrimedRealAssignment> const & reset) {
        this->_new_update(source,event,target,reset);
    }

    //! \brief Adds an unguarded transition to the automaton.
    //! The guard is the constant "True" i.e. the event is an input event.
    void new_transition(DiscreteLocation source,
                        DiscreteEvent event,
                        DiscreteLocation target) {
        this->_new_update(source,event,target,List<PrimedRealAssignment>());
    }

    //! \brief Adds a discrete transition for an automaton with a single mode.
    void new_transition(DiscreteEvent event,
                        List<PrimedRealAssignment> const & reset,
                        ContinuousPredicate const& guard,
                        EventKind kind=urgent) {
        this->new_transition(DiscreteLocation(),event,DiscreteLocation(),reset,guard,kind);
    }

    //@}

    //@{
    //! \name Data access and queries.

    //! \brief The modes of the automaton.
    Map<DiscreteLocation,DiscreteMode> const& modes() const;

    //! \brief The set of discrete locations.
    Set<DiscreteLocation> locations() const;

    //! \brief Test if the hybrid automaton has a discrete mode with the given exact \a location.
    bool has_location(DiscreteLocation location) const;


    //! \brief The discrete mode with given discrete state.
    const DiscreteMode& mode(DiscreteLocation location) const;
    DiscreteMode& mode(DiscreteLocation location);

    //! \brief Checks validity of the mode for the given \a location.
    //! Only checks for underspecified dynamics; overspecification is determined at build-time.
    void check_mode(DiscreteLocation location) const;
    //@}

    //@{
    //! \name Access for compositional hybrid automata.

    //! \brief The state (dotted) variables in the given location.
    List<RealVariable> state_variables(DiscreteLocation location) const;
    //! \brief The auxiliary (algebraic/output) variables in the given location.
    List<RealVariable> auxiliary_variables(DiscreteLocation location) const;
    //! \brief The input variables in the given location.
    Set<RealVariable> input_variables(DiscreteLocation location) const;
    //! \brief The target location when the \a event occurs in the \a source location.
    DiscreteLocation target_location(DiscreteLocation source, DiscreteEvent event) const;
    //! \brief The algebraic equations valid in the given location.
    List<RealAssignment> auxiliary_assignments(DiscreteLocation location) const;
    //! \brief The differential equations valid in the given location.
    List<DottedRealAssignment> dynamic_assignments(DiscreteLocation location) const;
    //! \brief The invariant predicates valid in the given location.
    ContinuousPredicate invariant_predicate(DiscreteLocation location, DiscreteEvent action) const;
    //! \brief The guard predicate for the given event in the given location.
    ContinuousPredicate guard_predicate(DiscreteLocation location, DiscreteEvent event) const;
    //! \brief The differential equations valid in the given location.
    List<PrimedRealAssignment> reset_assignments(DiscreteLocation source, DiscreteEvent event) const;
    //@}

    //@{
    //! \name Functions for conformance to HybridAutomatonInterface

    virtual HybridSpace state_space() const;
    //! \brief The continuous state space in the given location.
    virtual RealSpace continuous_state_space(DiscreteLocation) const;
    //! \brief The dimension of the continuous state space in the given location.
    virtual uint dimension(DiscreteLocation) const;

    //! \brief Test if the hybrid automaton has a discrete mode corresponding to the given location.
    virtual bool has_mode(DiscreteLocation location) const;
    //! \brief Test if the hybrid automaton has an invariant (either explicit or from an urgent transition) with the given \a event label in \a location.
    virtual bool has_invariant(DiscreteLocation location, DiscreteEvent event) const;
    //! \brief Tests if the automaton has an invariant or transition corresponding to the given location and event.
    virtual bool has_guard(DiscreteLocation, DiscreteEvent) const;
    //! \brief Test if the hybrid automaton has a discrete transition starting from the given location with the given event.
    virtual bool has_transition(DiscreteLocation source, DiscreteEvent event) const;
    //! \brief The events which are active in the given location.
    virtual Set<DiscreteEvent> events(DiscreteLocation) const;


    //! \brief The type of the event (urgent, permissive, impact etc).
    virtual EventKind event_kind(DiscreteLocation location, DiscreteEvent event) const;
    //! \brief The target for \a event from location \a source. Returns \a source if \a event is not present.
    virtual DiscreteLocation target(DiscreteLocation location, DiscreteEvent event) const;

    //! \brief The function outputting the auxiliary variables \f$y=h(x)\f$ in the location.
    virtual RealVectorFunction auxiliary_function(DiscreteLocation location) const;
    //! \brief The function outputting the differential equations \f$\dot{x}=f(x)\f$ in the location.
    virtual RealVectorFunction dynamic_function(DiscreteLocation location) const;
    //! \brief The invariant function \f$i(x)\leq 0\f$ corresponding to the given event.
    virtual RealScalarFunction invariant_function(DiscreteLocation location, DiscreteEvent event) const;
    //! \brief The guard function \f$g(x)\geq 0\f$ corresponding to the given event.
    virtual RealScalarFunction guard_function(DiscreteLocation location, DiscreteEvent event) const;
    //! \brief The reset function \f$x'=r(x)\f$ for the given event.
    virtual RealVectorFunction reset_function(DiscreteLocation location, DiscreteEvent event) const;

    //@}

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream&) const;
};

inline std::ostream& operator<<(std::ostream& os, const HybridAutomaton& ha) {
    return ha.write(os);
}



//! \ingroup SystemModule
//! \brief A hybrid automaton, formed by running a finite number of HybridAutomaton classes in parallel.
//!
//! \sa \ref HybridAutomaton
class CompositeHybridAutomaton
    : public HybridAutomatonInterface
    , public Loggable
{
    mutable DiscreteMode _cached_mode;
    void _cache_mode(DiscreteLocation) const;
  public:
    typedef HybridTime TimeType;
  public:
    //@{
    //! \name Constructors.

    //! \brief Default constructor creates empty automaton.
    CompositeHybridAutomaton();
    //! \brief Convert a single atomic hybrid automaton to a composite hybrid automaton for analysis.
    CompositeHybridAutomaton(const HybridAutomaton&);
    //! \brief Create the parallel composition of a list of atomic hybrid automata.
    CompositeHybridAutomaton(const List<HybridAutomaton>&);
    //@}

    //@{
    //! \name Methods for querying the component automata.

    //! \brief The number of component automata.
    uint number_of_components() const;
    //! \brief The \a i<sup>th</sup> component automaton.
    const HybridAutomaton& component(uint i) const;
    //@}

    //@{
    //! \name Methods for finding the modes, transitions and events of the composite automaton.

    //! The mode corresponding to the given location.
    DiscreteMode const& mode(DiscreteLocation) const;
    //! \brief Tests if the automaton has a mode corresponding to the given location.
    bool has_mode(DiscreteLocation) const;
    //! \brief Tests if the automaton has an invariant corresponding to the given location and event.
    bool has_invariant(DiscreteLocation, DiscreteEvent) const;
    //! \brief Tests if the automaton has an invariant or transition corresponding to the given location and event.
    bool has_guard(DiscreteLocation, DiscreteEvent) const;
    //! \brief Tests if the automaton has a transition corresponding to the given location and event.
    bool has_transition(DiscreteLocation, DiscreteEvent) const;

    //@}

    //@{
    //! \name Methods fo rextracting the continuous dynamics of the automaton.


    //! \brief The target location when \a event occurs in the \a source location.
    DiscreteLocation target_location(DiscreteLocation source, DiscreteEvent event) const;

    //! \brief The continuous variables which are defined in the location.
    List<RealVariable> variables(DiscreteLocation) const;
    //! \brief The continuous variables which are state variables in the location.
    //! The state variables are those defined by a differential equation i.e. the dotted variables.
    List<RealVariable> state_variables(DiscreteLocation) const;
    //! \brief The dependent continuous variables which are not state variables in the location.
    //! These variables are defined by algebraic equations.
    List<RealVariable> auxiliary_variables(DiscreteLocation) const;
    //! \brief The continuous variables which are not state or auxiliary variables in the location,
    //! but are required in one of the dynamic equations, constraints or resets.
    Set<RealVariable> input_variables(DiscreteLocation) const;

    //! \brief The algebraic equations valid in the location, ordered so that the defining equation for a variable
    //! occurs before any equation using that variable.
    List<RealAssignment> auxiliary_assignments(DiscreteLocation location) const;
    //! \brief The differential equations valid in the location.
    List<DottedRealAssignment> dynamic_assignments(DiscreteLocation location) const;
    //! \brief The reset equations used when the \a event occurs in the \a source location.
    List<PrimedRealAssignment> reset_assignments(DiscreteLocation source, DiscreteEvent event) const;
    //! \brief The invariant (time-can-progress predicates) corresponding to the given \a event.
    ContinuousPredicate invariant_predicate(DiscreteLocation location, DiscreteEvent event) const;
    //! \brief The guard (activation predicate) corresponding to the given \a event.
    ContinuousPredicate guard_predicate( DiscreteLocation location, DiscreteEvent event) const;
    //@}


    //@{
    //! \name Functions for conformance to HybridAutomatonInterface

    HybridSpace state_space() const;
    //! \brief The continuous state space in the given location.
    RealSpace continuous_state_space(DiscreteLocation) const;
    //! \brief The dimension of the continuous state space in the given location.
    uint dimension(DiscreteLocation) const;

    //! \brief The events which are active in the given location.
    Set<DiscreteEvent> events(DiscreteLocation) const;
    //! \brief The target for \a event from location \a source. Returns \a source if \a event is not present.
    DiscreteLocation target(DiscreteLocation location, DiscreteEvent event) const;

    //! \brief The function outputting the auxiliary variables \f$y=h(x)\f$ in the location.
    RealVectorFunction auxiliary_function(DiscreteLocation location) const;
    //! \brief The function outputting the differential equations \f$\dot{x}=f(x)\f$ in the location.
    RealVectorFunction dynamic_function(DiscreteLocation location) const;
    //! \brief The reset function \f$x'=r(x)\f$ for the given event.
    RealVectorFunction reset_function(DiscreteLocation location, DiscreteEvent event) const;
    //! \brief The invariant function \f$i(x)\leq 0\f$ corresponding to the given event.
    RealScalarFunction invariant_function(DiscreteLocation location, DiscreteEvent event) const;
    //! \brief The guard function \f$g(x)\geq 0\f$ corresponding to the given event.
    RealScalarFunction guard_function(DiscreteLocation location, DiscreteEvent event) const;
    //! \brief The type of the event (urgent, permissive, impact etc).
    EventKind event_kind(DiscreteLocation location, DiscreteEvent event) const;

    //@}


    //@{
    //! \name Methods for checking the validity of the automaton.

    //! \brief Checks whether the equations valid in the location are valid.
    //!
    //! Includes a check for algebraic dependencies, over-defined variables, under-defined variables, and
    //! variables which should be defined in a reset but are not.
    void check_mode(DiscreteLocation) const;
    //! \brief Runs check_mode() in any mode reachable under the discrete dynamics from the given initial location.
    void check_reachable_modes(DiscreteLocation) const;
    //! \brief Runs check_mode() in any mode reachable under the discrete dynamics from the given initial locations.
    void check_reachable_modes(const Set<DiscreteLocation>&) const;
    //@}

    //@{
    //! \name Discrete reachability analysis.

    //! \brief Performs a discrete reachability analysis from the given initial location.
    Set<DiscreteLocation> discrete_reachability(DiscreteLocation) const;
    //! \brief Performs a discrete reachability analysis from the given initial locations.
    Set<DiscreteLocation> discrete_reachability(const Set<DiscreteLocation>&) const;
    //@}

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream&) const;
  private:
    List<HybridAutomaton> _components;
};

CompositeHybridAutomaton parallel_composition(const List<HybridAutomaton>& components);

inline std::ostream& operator<<(std::ostream& os, const CompositeHybridAutomaton& ha) {
    return ha.write(os); }

} // namespace Ariadne

#endif // ARIADNE_HYBRID_AUTOMATON_H
