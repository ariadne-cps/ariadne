/***************************************************************************
 *            hybrid_automaton.h
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
#include "formula.h"
#include "logging.h"

#include "hybrid_automaton_interface.h"

namespace Ariadne {

class HybridTime;
class HybridSpace;
class HybridSet;
class HybridGrid;

class AtomicDiscreteMode;
class AtomicDiscreteTransition;
class AtomicHybridAutomaton;

class ScalarFunction;
class VectorFunction;
class Grid;


/*! \brief A discrete transition of a hybrid automaton, representing an instantaneous
 * jump from one discrete mode to another, governed by an activation set and a reset map.
 *
 * A %AtomicDiscreteTransition can only be created using the new_transition() method in
 * the %AtomicHybridAutomaton class.
 *
 * An invariant is modelled by a discrete transition with negative event id and null reset pointer.
 *
 * \sa \link Ariadne::AtomicHybridAutomaton \c AtomicHybridAutomaton \endlink, \link Ariadne::AtomicDiscreteMode \c AtomicDiscreteMode \endlink
 */
class AtomicDiscreteTransition
{
    friend class AtomicDiscreteMode;
    friend class AtomicHybridAutomaton;
  private:
    //  The discrete transition's identificator.
    DiscreteEvent _event;

    // IMPORTANT: Can't use pointers here since they break default copy
    // used to create composite automaton

    //  The source of the discrete transition.
    AtomicDiscreteLocation _source;

    //  The target of the discrete transition.
    AtomicDiscreteLocation _target;

    //  The guard predicate of the discrete transition.
    ContinuousPredicate _guard_predicate;

    //  The reset assignments of the discrete transition.
    List<RealUpdateAssignment> _update_assignments;

    //  Whether or not the transition is forced (Deprecated)
    bool _forced;

  public:

    //! \brief The discrete event associated with the discrete transition.
    DiscreteEvent event() const { return this->_event; }

    //! \brief The source mode of the discrete transition.
    AtomicDiscreteLocation source() const { return this->_source; }

    //! \brief The target mode of the discrete transition.
    AtomicDiscreteLocation target() const { return this->_target; };

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream& os) const;
  private:
    AtomicDiscreteTransition(DiscreteEvent event,
                             const AtomicDiscreteMode& source_mode,
                             const AtomicDiscreteMode& target_mode,
                             const List<RealUpdateAssignment>& reset,
                             const ContinuousPredicate& guard);

    AtomicDiscreteTransition(DiscreteEvent event,
                             AtomicDiscreteLocation source,
                             AtomicDiscreteLocation target,
                             const List<RealUpdateAssignment>& reset,
                             const ContinuousPredicate& guard);
};

inline std::ostream& operator<<(std::ostream& os, const AtomicDiscreteTransition& dt) {
    return dt.write(os); }




/*! \brief A discrete mode of a hybrid automaton, comprising continuous evolution given by a vector field
 * within and invariant constraint set.
 *
 * A %AtomicDiscreteMode can only be created using the new_mode() method in
 * the %AtomicHybridAutomaton class.
 *
 * \sa \link Ariadne::AtomicHybridAutomaton \c AtomicHybridAutomaton \endlink, \link Ariadne::AtomicDiscreteTransition \c AtomicDiscreteTransition \endlink
 */
class AtomicDiscreteMode {
    friend class AtomicDiscreteTransition;
    friend class AtomicHybridAutomaton;
  private:
    // The discrete mode's discrete state.
    AtomicDiscreteLocation _location;

    // The algebraic equations
    List<RealAlgebraicAssignment> _algebraic_assignments;

    // The algebraic equations
    List<RealDifferentialAssignment> _differential_assignments;

    // The discrete mode's invariants.
    Map<DiscreteEvent,ContinuousPredicate> _invariant_predicates;

    // The discrete mode's urgent guards.
    Map<DiscreteEvent,ContinuousPredicate> _urgent_guard_predicates;

    // The discrete mode's urgent guards.
    Map<DiscreteEvent,ContinuousPredicate> _permissive_guard_predicates;

    // The target locations for the discrete transitions.
    Map<DiscreteEvent,AtomicDiscreteLocation> _targets;

    // The update assignments for the doscrete transitions.
    Map<DiscreteEvent, List<RealUpdateAssignment> > _update_assignments;
  private:
    AtomicDiscreteMode(AtomicDiscreteLocation, const List<RealAlgebraicAssignment>&,const List<RealDifferentialAssignment>&);
  public:
    //! \brief The mode's discrete state.
    AtomicDiscreteLocation location() const { return this->_location; }

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream& os) const;
};


inline std::ostream& operator<<(std::ostream& os, const AtomicDiscreteMode& dm) {
    return dm.write(os); }

inline bool operator<(const AtomicDiscreteMode& mode1, const AtomicDiscreteMode& mode2) {
    return mode1.location() < mode2.location(); }

inline bool operator<(const AtomicDiscreteTransition& transition1, const AtomicDiscreteTransition& transition2) {
    return transition1.event() < transition2.event()
        || (transition1.event() == transition2.event()
            && transition1.source() < transition2.source());
}





/*! \brief A hybrid automaton, comprising continuous-time behaviour
 *  at each discrete mode, coupled by instantaneous discrete transitions.
 *  The state space is given by a hybrid set.
 *
 * A hybrid automaton is a dynamic system with evolution in both
 * continuous time and discrete time.
 * The state space is a product \f$X=\bigcup\{q\}\times X_q\f$
 * where \f$q\f$ is the <em>discrete state</em> and \f$X_q\f$
 * is the <em>continuous state space</em> of corresponding to
 * each discrete state.
 *
 * For each %AtomicDiscreteMode, the dynamics is given by a
 * %VectorField describing the continuous dynamics,
 * and a %Set giving an invariants which must be satisified at
 * all times.
 *
 * The discrete time behaviour is specified by %AtomicDiscreteTransition
 * objects.
 * Each discrete transition represents an jump from a \a source
 * mode to a \a target mode.
 * There can be at most one discrete transition in an automaton
 * with the same event and source.
 *
 * A discrete transision can either be \em forced or \em unforced.
 * A forced transition much occur as soon as it is activated.
 * An unforced transition may occur at any time it is activated,
 * but is only forced to occur if the continuous evolution is
 * blocked by an invariant.
 *
 * \sa \link Ariadne::AtomicDiscreteMode \c AtomicDiscreteMode \endlink, \link Ariadne::AtomicDiscreteTransition \c AtomicDiscreteTransition \endlink, \link Ariadne::CompositeHybridAutomaton \c CompositeHybridAutomaton \endlink

 */
class AtomicHybridAutomaton
    : public Loggable
{
  public:
    //! \brief The type used to represent time.
    typedef HybridTime TimeType;
    //! \brief The type used to represent real numbers.
    typedef double RealType ;
    //! \brief The type used to describe the state space.
    typedef HybridSpace StateSpaceType;

  private:
    //! \brief The hybrid automaton's name.
    String _name;

    //! \brief The list of the hybrid automaton's discrete modes.
    Map< AtomicDiscreteLocation, AtomicDiscreteMode > _modes;

  public:
    //@{
    //! \name Constructors and destructors

    //! \brief Construct an empty automaton with no name
    AtomicHybridAutomaton();

    //! \brief Construct an empty automaton with the given name
    AtomicHybridAutomaton(const String& name);

    //! \brief Construct dynamically-allocated copy. (Not currently implemented)
    AtomicHybridAutomaton* clone() const;

    //! \brief  Destructor.
    ~AtomicHybridAutomaton();
    //@}

    //@{
    //! \name Methods for building the automaton.

    //! \brief Adds a discrete mode to the automaton.
    void new_mode(AtomicDiscreteLocation state,
                  const List<RealAlgebraicAssignment>& equations,
                  const List<RealDifferentialAssignment>& dynamic);

    //! \brief Adds a discrete mode to the automaton.
    void new_mode(AtomicDiscreteLocation state,
                  const List<RealDifferentialAssignment>& dynamic);

    //! \brief Adds a discrete mode to the automaton.
    void new_mode(AtomicDiscreteLocation state,
                  const List<RealAlgebraicAssignment>& equations);


    //! \brief Adds a discrete mode to the automaton.
    void new_invariant(AtomicDiscreteLocation state,
                       DiscreteEvent event,
                       const ContinuousPredicate& constraint);

    //! \brief Adds an urgent guard to the automaton.
    void new_urgent_guard(AtomicDiscreteLocation state,
                          DiscreteEvent event,
                          const ContinuousPredicate& constraint);

    //! \brief Adds a permissive guard to the automaton.
    void new_permissive_guard(AtomicDiscreteLocation state,
                              DiscreteEvent event,
                              const ContinuousPredicate& constraint);

    //! \brief Adds a reset to the automaton. (Same as new_transition.)
    void new_reset(AtomicDiscreteLocation state,
                   DiscreteEvent event,
                   AtomicDiscreteLocation target,
                   const List<RealUpdateAssignment>& reset);


    //! \brief Adds an unguarded transition to the automaton.
    //! The guard is the constant "True" i.e. the event is an input event.
    void new_transition(AtomicDiscreteLocation state,
                        DiscreteEvent event,
                        AtomicDiscreteLocation target,
                        const List<RealUpdateAssignment>& reset);

    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //!
    //!    \param source is the transition's source location.
    //!    \param event is the transition's event.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param guard is the transition's activation region.
    void new_transition(AtomicDiscreteLocation source,
                        DiscreteEvent event,
                        const ContinuousPredicate& guard,
                        AtomicDiscreteLocation target,
                        const List<RealUpdateAssignment>& reset,
                        Urgency urgency=urgent);

    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //!
    //!    \param source is the transition's source location.
    //!    \param event is the transition's event.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param guard is the transition's activation region.
    void new_transition(AtomicDiscreteLocation source,
                        DiscreteEvent event,
                        AtomicDiscreteLocation target,
                        const List<RealUpdateAssignment>& reset,
                        const ContinuousPredicate& guard,
                        Urgency urgency=urgent);

    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //! The reset is trivial. This form is for the case that there are no continuous state variables in the new location.
    //!
    //!    \param source is the transition's source location.
    //!    \param event is the transition's event.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param guard is the transition's activation region.
    void new_transition(AtomicDiscreteLocation source,
                        DiscreteEvent event,
                        AtomicDiscreteLocation target,
                        const ContinuousPredicate& guard);


    //! \brief Set the grid controlling relative scaling in the mode. \deprecated
    void set_grid(AtomicDiscreteLocation location, const Grid& grid);

    //! \brief Set the grid controlling relative scaling. This method sets the same grid for every mode. \deprecated
    void set_grid(const Grid& grid);

    //! \brief Set the hybrid grid controlling relative scaling. \deprecated
    void set_grid(const HybridGrid& hgrid);

    //@}

    //@{
    //! \name Data access and queries.

    //! \brief Returns the hybrid automaton's name.
    const String& name() const;

    //! \brief The set of discrete locations.
    Set<AtomicDiscreteLocation> locations() const;

    //! \brief The discrete events corresponding to an invariant in \a source.
    Set<DiscreteEvent> invariant_events(AtomicDiscreteLocation source) const;

    //! \brief The discrete events corresponding to a transition \a source.
    Set<DiscreteEvent> transition_events(AtomicDiscreteLocation source) const;

    //! \brief The discrete events corresponding to an invariant \a source.
    Set<DiscreteEvent> blocking_events(AtomicDiscreteLocation source) const;

    //! \brief The discrete events corresponding to an urgent transition in \a source.
    Set<DiscreteEvent> urgent_events(AtomicDiscreteLocation source) const;

    //! \brief The discrete events corresponding to an permissive transitions in \a source.
    Set<DiscreteEvent> permissive_events(AtomicDiscreteLocation source) const;


    //! \brief Test if the hybrid automaton has a discrete mode with the given \a location.
    bool has_mode(AtomicDiscreteLocation location) const;

    //! \brief Test if the hybrid automaton has a discrete transition starting from the given location with the given event.
    bool has_transition(AtomicDiscreteLocation source, DiscreteEvent event) const;

    //! \brief Test if the hybrid automaton has an invariant (either explicit or from an urgent transition) with the given \a event label in \a location.
    bool has_invariant(AtomicDiscreteLocation location, DiscreteEvent event) const;

    //! \brief Test if the hybrid automaton has a guard (urgent or permissive) with \a event_id and \a source_id.
    bool has_guard(AtomicDiscreteLocation location, DiscreteEvent event) const;


    //! \brief The discrete modes of the automaton.
    Set<AtomicDiscreteMode> modes() const;

    //! \brief The discrete mode with given discrete state.
    const AtomicDiscreteMode& mode(AtomicDiscreteLocation location) const;

    //! \brief The discrete mode with given discrete state.
    AtomicDiscreteMode& mode(AtomicDiscreteLocation location);


    //! \brief The target for \a event from location \a source. Returns \a source if \a event is not present.
    AtomicDiscreteLocation target(AtomicDiscreteLocation location, DiscreteEvent event) const;

    //! \brief The target for \a event from location \a source. Returns the identity if \a event is not present.
    List<RealUpdateAssignment> reset(AtomicDiscreteLocation location, DiscreteEvent event) const;

    //! \brief The invariant for action label \a event in \a location.
    const ContinuousPredicate& invariant(AtomicDiscreteLocation location, DiscreteEvent event) const;

    //! \brief The guard for the discrete transition with given \a event and \a source location.
    const ContinuousPredicate& guard(AtomicDiscreteLocation location, DiscreteEvent event) const;


    //@}

    //@{
    //! \name New-style access for compositional hybrid automata.


    //! \brief The state (dotted) variables in the given location.
    List<RealVariable> state_variables(AtomicDiscreteLocation location) const;
    //! \brief The auxiliary (algebraic/output) variables in the given location.
    List<RealVariable> auxiliary_variables(AtomicDiscreteLocation location) const;
    //! \brief The algebraic equations valid in the given location.
    List<RealAssignment> algebraic_assignments(AtomicDiscreteLocation location) const;
    //! \brief The differential equations valid in the given location.
    List<DottedRealAssignment> differential_assignments(AtomicDiscreteLocation location) const;
    //! \brief The differential equations valid in the given location.
    List<PrimedRealAssignment> update_assignments(AtomicDiscreteLocation source, DiscreteEvent event) const;
    //! \brief The invariant predicates valid in the given location.
    ContinuousPredicate invariant_predicate(AtomicDiscreteLocation location, DiscreteEvent action) const;
    //! \brief The guard predicate for the given event in the given location.
    ContinuousPredicate guard_predicate(AtomicDiscreteLocation location, DiscreteEvent event) const;
    //@}

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream&) const;
};

inline std::ostream& operator<<(std::ostream& os, const AtomicHybridAutomaton& ha) {
    return ha.write(os); }



/*! \brief A hybrid automaton, formed by running a finite number of AtomicHybridAutomaton classes in parallel.
 *
 * \sa \link Ariadne::AtomicHybridAutomaton \c AtomicHybridAutomaton \endlink
 */
class CompositeHybridAutomaton
    : public HybridAutomatonInterface
    , public Loggable
{
  public:
    typedef HybridTime TimeType;
  public:
    //@{
    //! \name Constructors.

    //! \brief Default constructor creates empty automaton.
    CompositeHybridAutomaton();
    //! \brief Convert a single atomic hybrid automaton to a composite hybrid automaton for analysis.
    CompositeHybridAutomaton(const AtomicHybridAutomaton&);
    //! \brief Create the parallel composition of a list of atomic hybrid automata.
    CompositeHybridAutomaton(const List<AtomicHybridAutomaton>&);
    //@}

    //@{
    //! \name Methods for querying the component automata.

    //! \brief The number of component automata.
    uint number_of_components() const;
    //! \brief The \a i<sup>th</sup> component automaton.
    const AtomicHybridAutomaton& component(uint i) const;
    //@}

    // The signatures of these methods are wrong.
    //const Set<AtomicDiscreteMode> modes() const;
    //const Set<AtomicDiscreteTransition> transitions(DiscreteLocation location) const;

    //@{
    //! \name Methods for finding the modes, transitions and events of the composite automaton.

    //! \brief Tests if the automaton has a mode corresponding to the given location.
    bool has_mode(DiscreteLocation) const;
    //! \brief Tests if the automaton has a transition corresponding to the given location and event.
    bool has_transition(DiscreteLocation, DiscreteEvent) const;

    //! \brief The dimension of the state spacec in the given \a location.
    uint dimension(DiscreteLocation location) const;

    //! \brief The set of all events possible in the given \a location.
    Set<DiscreteEvent> events(DiscreteLocation location) const;
    //! \brief The set of events corresponding to a discrete transition.
    Set<DiscreteEvent> transition_events(DiscreteLocation) const;
    //! \brief The set of events corresponding to an invariant or time-can-progress predicate.
    Set<DiscreteEvent> invariant_events(DiscreteLocation) const;
    //! \brief The set of events corresponding to a non-urgent action.
    Set<DiscreteEvent> permissive_events(DiscreteLocation) const;
    //! \brief The set of events corresponding to an urgent action.
    Set<DiscreteEvent> urgent_events(DiscreteLocation) const;
    //! \brief The set of events which prohibit further continuous evolution.
    //! Equal to the union of invariant_events and urgent_events.
    Set<DiscreteEvent> blocking_events(DiscreteLocation) const;

    //! \brief The target location when \a event occurs in the \a source location.
    DiscreteLocation target(DiscreteLocation source, DiscreteEvent event) const;
    //@}

    //@{
    //! \name Methods for extracting the continuous dynamics of the automaton.


    //! \brief The continuous variables which are defined in the location.
    List<RealVariable> variables(DiscreteLocation) const;
    //! \brief The continuous variables which are state variables in the location.
    //! The state variables are those defined by a differential equation i.e. the dotted variables.
    List<RealVariable> state_variables(DiscreteLocation) const;
    //! \brief The dependent continuous variables which are not state variables in the location.
    //! These variables are defined by algebraic equations.
    List<RealVariable> auxiliary_variables(DiscreteLocation) const;

    //! \brief The algebraic equations valid in the location, ordered so that the defining equation for a variable
    //! occurs before any equation using that variable.
    List<RealAssignment> algebraic_assignments(DiscreteLocation location) const;
    //! \brief The differential equations valid in the location.
    List<DottedRealAssignment> differential_assignments(DiscreteLocation location) const;
    //! \brief The reset equations used when the \a event occurs in the \a source location.
    List<PrimedRealAssignment> update_assignments(DiscreteLocation source, DiscreteEvent event) const;
    //! \brief The invariant (time-can-progress predicates) corresponding to the given \a event.
    ContinuousPredicate invariant_predicate(DiscreteLocation location, DiscreteEvent event) const;
    //! \brief The guard (activation predicate) corresponding to the given \a event.
    ContinuousPredicate guard_predicate( DiscreteLocation location, DiscreteEvent event) const;
    //@}


    EventKind event_kind(DiscreteLocation location, DiscreteEvent event) const;
    VectorFunction output_function(DiscreteLocation) const;
    VectorFunction dynamic_function(DiscreteLocation) const;
    VectorFunction reset_function(DiscreteLocation, DiscreteEvent) const;
    ScalarFunction invariant_function(DiscreteLocation, DiscreteEvent) const;
    ScalarFunction guard_function(DiscreteLocation, DiscreteEvent) const;

    Grid grid(DiscreteLocation) const;

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
    List<AtomicHybridAutomaton> _components;
};

CompositeHybridAutomaton parallel_composition(const List<AtomicHybridAutomaton>& components);

inline std::ostream& operator<<(std::ostream& os, const CompositeHybridAutomaton& ha) {
    return ha.write(os); }


} // namespace Ariadne

#endif // ARIADNE_HYBRID_AUTOMATON_H
