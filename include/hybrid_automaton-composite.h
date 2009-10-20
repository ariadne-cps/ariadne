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
    // \brief The discrete transition's identificator.
    DiscreteEvent _event;

    // \brief The source of the discrete transition.
    const AtomicDiscreteMode* _source_mode;

    // \brief The target of the discrete transition.
    const AtomicDiscreteMode* _target_mode;

    // \brief The guard predicate of the discrete transition.
    ContinuousPredicate _guard_predicate;

    // \brief The reset assignments of the discrete transition.
    List<RealUpdateAssignment> _update_assignments;

    // \brief Whether or not the transition is forced (Deprecated)
    bool _forced;

  public:

    //! \brief The discrete event associated with the discrete transition.
    DiscreteEvent event() const { return this->_event; }

    //! \brief The source mode of the discrete transition.
    const AtomicDiscreteMode& source_mode() const { return *this->_source_mode; }

    //! \brief The target mode of the discrete transition.
    const AtomicDiscreteMode& target_mode() const { return *this->_target_mode; }

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream& os) const;
  private:
    // Construct from shared pointers (for internal use).
    AtomicDiscreteTransition(DiscreteEvent event,
                             const AtomicDiscreteMode& source,
                             const AtomicDiscreteMode& target,
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

    // The discrete mode's grid for reachability analysis.
    Map<DiscreteEvent,AtomicDiscreteTransition> _transitions;
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
            && transition1.source_mode().location() < transition2.source_mode().location());
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
 * \sa \link Ariadne::AtomicDiscreteMode \c AtomicDiscreteMode \endlink, \link Ariadne::AtomicDiscreteTransition \c AtomicDiscreteTransition \endlink

 */
class AtomicHybridAutomaton
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
    const AtomicDiscreteMode&
    new_mode(AtomicDiscreteLocation state,
             const List<RealAlgebraicAssignment>& equations,
             const List<RealDifferentialAssignment>& dynamic);

    //! \brief Adds a discrete mode to the automaton.
    const AtomicDiscreteMode&
    new_mode(AtomicDiscreteLocation state,
             const List<RealDifferentialAssignment>& dynamic);

    //! \brief Adds a discrete mode to the automaton.
    const AtomicDiscreteMode&
    new_mode(AtomicDiscreteLocation state,
             const List<RealAlgebraicAssignment>& equations);


    //! \brief Adds a discrete mode to the automaton.
    const AtomicDiscreteMode&
    new_invariant(AtomicDiscreteLocation state,
                  DiscreteEvent event,
                  const ContinuousPredicate& constraint);

    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //!
    //!    \param source is the transition's source location.
    //!    \param event is the transition's event.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param guard is the transition's activation region.
    const AtomicDiscreteTransition& new_transition(AtomicDiscreteLocation source,
                                             DiscreteEvent event,
                                             AtomicDiscreteLocation target,
                                             const List<RealUpdateAssignment>& reset,
                                             const ContinuousPredicate& guard);

    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //! The guard is the constant "True" i.e. the event is an input event.
    //!
    //!    \param source is the transition's source location.
    //!    \param event is the transition's event.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    const AtomicDiscreteTransition& new_transition(AtomicDiscreteLocation source,
                                             DiscreteEvent event,
                                             AtomicDiscreteLocation target,
                                             const List<RealUpdateAssignment>& reset);

     //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //! The reset is trivial. This form is for the case that there are no continuous state variables in the new location.
    //!
    //!    \param source is the transition's source location.
    //!    \param event is the transition's event.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param guard is the transition's activation region.
    const AtomicDiscreteTransition& new_transition(AtomicDiscreteLocation source,
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

    //! \brief The discrete events possible in location \a source.
    Set<DiscreteEvent> events(AtomicDiscreteLocation source) const;

    //! \brief Test if the hybrid automaton has a discrete mode with the given \a location.
    bool has_mode(AtomicDiscreteLocation location) const;

    //! \brief Test if the hybrid automaton has a discrete transition with the given \a event label in \a location.
    bool has_invariant(AtomicDiscreteLocation source, DiscreteEvent event) const;

    //! \brief Test if the hybrid automaton has a discrete transition with \a event_id and \a source_id.
    bool has_transition(AtomicDiscreteLocation source, DiscreteEvent event) const;

    //! \brief The discrete mode with given discrete state.
    const AtomicDiscreteMode& mode(AtomicDiscreteLocation location) const;

    //! \brief The discrete mode with given discrete state.
    AtomicDiscreteMode& mode(AtomicDiscreteLocation location);

    //! \brief The discrete transition with given \a event and \a source location.
    const ContinuousPredicate& invariant(AtomicDiscreteLocation location, DiscreteEvent event) const;

    //! \brief The discrete transition with given \a event and \a source location.
    const AtomicDiscreteTransition& transition(AtomicDiscreteLocation source, DiscreteEvent event) const;

    //! \brief The discrete transitions from location \a source.
    Set<AtomicDiscreteTransition> transitions(AtomicDiscreteLocation source) const;

    //@}

    //@{
    //! \name New-style access for compositional hybrid automata.

    //! \brief The target location of the discrete event from the given discrete location.
    AtomicDiscreteLocation target(const AtomicDiscreteLocation& source, const DiscreteEvent& event);
    //! \brief The state (dotted) variables in the given location.
    List<RealVariable> state_variables(AtomicDiscreteLocation location) const;
    //! \brief The auxiliary (algebraic/output) variables in the given location.
    List<RealVariable> auxiliary_variables(AtomicDiscreteLocation location) const;
    //! \brief The algebraic equations valid in the given location.
    List<RealAssignment> algebraic_assignments(const AtomicDiscreteLocation& location) const;
    //! \brief The differential equations valid in the given location.
    List<DottedRealAssignment> differential_assignments(const AtomicDiscreteLocation& location) const;
    //! \brief The differential equations valid in the given location.
    List<PrimedRealAssignment> update_assignments(const AtomicDiscreteLocation& source, const DiscreteEvent& event) const;
    //! \brief The invariant predicates valid in the given location.
    Map<DiscreteEvent,ContinuousPredicate> invariant_predicates(const AtomicDiscreteLocation& location) const;
    ContinuousPredicate invariant_predicate(const AtomicDiscreteLocation& location, const DiscreteEvent& action) const;
    //! \brief The guard predicate for the given event in the given location.
    ContinuousPredicate guard_predicate(const AtomicDiscreteLocation& location, const DiscreteEvent& event) const;
    //@}

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream&) const;
};

inline std::ostream& operator<<(std::ostream& os, const AtomicHybridAutomaton& ha) {
    return ha.write(os); }



class CompositeHybridAutomaton {
  public:
    CompositeHybridAutomaton();
    CompositeHybridAutomaton(const AtomicHybridAutomaton&);
    CompositeHybridAutomaton(const List<AtomicHybridAutomaton>&);
    uint number_of_components() const;
    const AtomicHybridAutomaton& component(uint) const;

    const Set<AtomicDiscreteMode> modes() const;
    const Set<AtomicDiscreteTransition> transitions(const DiscreteLocation& location) const;

    bool has_mode(const DiscreteLocation&) const;
    bool has_transition(const DiscreteLocation&, const DiscreteEvent&) const;
    Set<DiscreteEvent> events(const DiscreteLocation&) const;
    DiscreteLocation target(const DiscreteLocation&, const DiscreteEvent&) const;

    List<RealVariable> variables(const DiscreteLocation&) const;
    List<RealVariable> state_variables(const DiscreteLocation&) const;
    List<RealVariable> auxiliary_variables(const DiscreteLocation&) const;
    List<RealAssignment> algebraic_assignments(const DiscreteLocation&) const;
    List<DottedRealAssignment> differential_assignments(const DiscreteLocation&) const;
    List<PrimedRealAssignment> update_assignments(const DiscreteLocation&, const DiscreteEvent&) const;

    Map<DiscreteEvent,ContinuousPredicate> invariant_predicates(const DiscreteLocation&) const;
    ContinuousPredicate invariant_predicate(const DiscreteLocation&, const DiscreteEvent&) const;
    Map<DiscreteEvent,ContinuousPredicate> guard_predicates(const DiscreteLocation&) const;
    ContinuousPredicate guard_predicate( const DiscreteLocation&, const DiscreteEvent&) const;

    VectorFunction output_function(const DiscreteLocation&) const;
    VectorFunction dynamic_function(const DiscreteLocation&) const;
    VectorFunction reset_function(const DiscreteLocation&, const DiscreteEvent&) const;

    Map<DiscreteEvent,ScalarFunction> invariant_functions(const DiscreteLocation&) const;
    Map<DiscreteEvent,ScalarFunction> guard_functions(const DiscreteLocation&) const;
    ScalarFunction invariant_function(const DiscreteLocation&, const DiscreteEvent&) const;
    ScalarFunction guard_function(const DiscreteLocation&, const DiscreteEvent&) const;

    void check_mode(const DiscreteLocation&) const;

    std::ostream& write(std::ostream&) const;
  private:
    List<AtomicHybridAutomaton> _components;
};

inline CompositeHybridAutomaton parallel_composition(const List<AtomicHybridAutomaton>& components) {
    return CompositeHybridAutomaton(components); }

inline std::ostream& operator<<(std::ostream& os, const CompositeHybridAutomaton& ha) {
    return ha.write(os); }

Set<DiscreteLocation> discrete_reachability(const CompositeHybridAutomaton&, const DiscreteLocation&);

} // namespace Ariadne

#endif // ARIADNE_HYBRID_AUTOMATON_H
