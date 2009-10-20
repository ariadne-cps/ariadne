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
#include "discrete_state.h"
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


/*! \brief A discrete mode of a hybrid automaton, comprising continuous evolution given by a vector field
 * within and invariant constraint set.
 *
 * A %AtomicDiscreteMode can only be created using the new_mode() method in
 * the %AtomicHybridAutomaton class.
 *
 * \sa \link Ariadne::AtomicHybridAutomaton \c AtomicHybridAutomaton \endlink, \link Ariadne::AtomicDiscreteTransition \c AtomicDiscreteTransition \endlink
 */
class AtomicDiscreteMode {
    friend class AtomicHybridAutomaton;
  private:

    // The discrete mode's discrete state.
    DiscreteState _location;

    // The discrete mode's state variables
    RealSpace _state_space;

    // The discrete mode's auxiliary/output variables
    RealSpace _auxiliary_space;

    // The discrete mode's input variables
    RealSpace _input_space;

    // The algebraic equations
    List<RealAlgebraicAssignment> _algebraic_assignments;

    // The algebraic equations
    List<RealDifferentialAssignment> _differential_assignments;

    // The discrete mode's invariants.
    Map< DiscreteEvent, ContinuousPredicate > _invariant_predicates;

    // The discrete mode's vector field.
    VectorFunction _dynamic;

    // The discrete mode's invariants.
    std::map< DiscreteEvent, ScalarFunction > _invariants;

    // The discrete mode's grid for reachability analysis.
    boost::shared_ptr< const Grid > _grid;
  public:
    //! \brief The mode's discrete state.
    DiscreteState location() const {
        return this->_location; }

    //! \brief The discrete mode's dynamic (a vector field).
    const VectorFunction& dynamic() const {
        return this->_dynamic; }

    //! \brief The discrete mode's invariants, converted to vector functions.
    const std::map< DiscreteEvent, VectorFunction > vector_invariants() const {
        std::map<DiscreteEvent,VectorFunction> result;
        for(std::map<DiscreteEvent,ScalarFunction>::const_iterator iter=this->_invariants.begin();
            iter!=this->_invariants.end(); ++iter)
        {
            result[iter->first]=VectorFunction(1u,iter->second);
        }
        return result;
    }

    //! \brief The discrete mode's invariants.
    const std::map< DiscreteEvent, ScalarFunction >& invariants() const {
        return this->_invariants; }

    //! \brief The discrete mode's default spacial grid.
    const Grid& grid() const {
        return *this->_grid; }

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
    AtomicDiscreteMode(DiscreteState location,
                 const VectorFunction& dynamic);

    // Construct from objects managed by shared pointers (for internal use)
    AtomicDiscreteMode(DiscreteState location,
                 const VectorFunction dynamic,
                 const std::vector< ScalarFunction >& invariants);

};


inline std::ostream& operator<<(std::ostream& os, const AtomicDiscreteMode& dm) {
    return dm.write(os); }

inline bool operator<(const AtomicDiscreteMode& mode1, const AtomicDiscreteMode& mode2) {
    return mode1.location() < mode2.location(); }




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
    friend class AtomicHybridAutomaton;
  private:
    // \brief The discrete transition's identificator.
    DiscreteEvent _event;

    // \brief The source of the discrete transition.
    const AtomicDiscreteMode* _source;

    // \brief The target of the discrete transition.
    const AtomicDiscreteMode* _target;

    // \brief The activation region of the discrete transition.
    ScalarFunction _activation;

    // \brief The reset of the discrete transition.
    VectorFunction _reset;

    // \brief Whether or not the transition is forced.
    bool _forced;

    List<RealUpdateAssignment> _update_assignments;

    ContinuousPredicate _guard_predicate;
  public:

    //! \brief The discrete event associated with the discrete transition.
    DiscreteEvent event() const {
        return this->_event; }

    //! \brief The source mode of the discrete transition.
    const AtomicDiscreteMode& source() const {
        return *this->_source; }

    //! \brief The target of the discrete transition.
    const AtomicDiscreteMode& target() const {
        return *this->_target; }


    //! \brief The activation region of the discrete transition.
    const ScalarFunction& activation() const {
        return this->_activation;
    }

    //! \brief The activation region of the discrete transition.
    const VectorFunction vector_activation() const {
        return VectorFunction(1u,this->_activation);
    }


    //! \brief The reset map of the discrete transition.
    const VectorFunction& reset() const {
        return this->_reset;
    }

    //! \brief True if the transition is forced (occurs as soon as it is activated).
    bool forced() const {
        return this->_forced;
    }

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream& os) const;
  private:


    // Construct from shared pointers (for internal use).
    AtomicDiscreteTransition(DiscreteEvent event,
                       const AtomicDiscreteMode& source,
                       const AtomicDiscreteMode& target,
                       const VectorFunction& reset,
                       const ScalarFunction& activation,
                       bool forced=false);

};

inline std::ostream& operator<<(std::ostream& os, const AtomicDiscreteTransition& dt) {
    return dt.write(os); }

inline bool operator<(const AtomicDiscreteTransition& transition1, const AtomicDiscreteTransition& transition2) {
    return transition1.event() < transition2.event()
        || (transition1.event() == transition2.event()
            && transition1.source().location() < transition2.source().location());
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


    typedef std::map<DiscreteEvent,ScalarFunction>::const_iterator invariant_const_iterator;
    typedef std::set<AtomicDiscreteTransition>::const_iterator discrete_transition_const_iterator;
    typedef std::set<AtomicDiscreteMode>::const_iterator discrete_mode_const_iterator;
    typedef std::set<AtomicDiscreteTransition>::const_iterator transition_const_iterator;
    typedef std::set<AtomicDiscreteMode>::const_iterator mode_const_iterator;
  private:
    //! \brief The hybrid automaton's name.
    std::string _name;

    //! \brief The list of the hybrid automaton's discrete modes.
    std::set< AtomicDiscreteMode > _modes;

    //! \brief The hybrid automaton's transitions.
    std::set< AtomicDiscreteTransition > _transitions;

  public:
    //@{
    //! \name Constructors and destructors

    //! \brief Construct an empty automaton with no name
    AtomicHybridAutomaton();

    //! \brief Construct an empty automaton with the given name
    AtomicHybridAutomaton(const std::string& name);

    //! \brief Construct dynamically-allocated copy. (Not currently implemented)
    AtomicHybridAutomaton* clone() const;

    //! \brief  Destructor.
    ~AtomicHybridAutomaton();
    //@}

    //@{
    //! \name Methods for building the automaton.

    //! \brief Adds a discrete mode to the automaton.
    const AtomicDiscreteMode& new_mode(DiscreteState state,
                                 const List<RealAlgebraicAssignment>& equations,
                                 const List<RealDifferentialAssignment>& dynamic);

    //! \brief Adds a discrete mode to the automaton.
    const AtomicDiscreteMode& new_mode(DiscreteState state,
                                 const List<RealDifferentialAssignment>& dynamic);

    //! \brief Adds a discrete mode to the automaton.
    const AtomicDiscreteMode& new_mode(DiscreteState state,
                                 const List<RealAlgebraicAssignment>& equations);

    //! \brief Adds a discrete mode to the automaton.
    const AtomicDiscreteMode& new_invariant(DiscreteState state,
                                      DiscreteEvent event,
                                      const ContinuousPredicate& constraint);

     //! \brief Adds a discrete mode to the automaton.
    //!
    //!   \param state is the mode's discrete state.
    //!   \param dynamic is the mode's vector field.
    const AtomicDiscreteMode& new_mode(DiscreteState state,
                                 const VectorFunction& dynamic);

    //! \brief Adds an invariant to a mode of the automaton.
    //!
    //!   \param state is the mode's discrete state.
    //!   \param invariant is the new invariant condition, in the form \f$g(x)<0\f$.

    const AtomicDiscreteMode& new_invariant(DiscreteState state,
                                      const ScalarFunction& invariant);

    //! \brief Adds an invariants to a mode of the automaton.
    //!
    //!   \param state is the mode's discrete state.
    //!   \param invariants is the new invariants condition.

    const AtomicDiscreteMode& new_invariant(DiscreteState state,
                                      const VectorFunction& invariants);

    //! \brief Adds an invariants to a mode of the automaton.
    //!
    //!    \param mode is the discrete mode.
    //!    \param invariants is the new invariants condition.

    const AtomicDiscreteMode& new_invariant(const AtomicDiscreteMode& mode,
                                      const VectorFunction& invariants);


    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    //!    \param forced determines whether the transision is forced (urgent) or unforced (permissive).
    const AtomicDiscreteTransition& new_transition(DiscreteEvent event,
                                             DiscreteState source,
                                             DiscreteState target,
                                             const VectorFunction& reset,
                                             const ScalarFunction& activation,
                                             bool forced);

    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region; must be a one-dimensional valued function.
    //!    \param forced determines whether the transision is forced (urgent) or unforced (permissive).
    const AtomicDiscreteTransition& new_transition(DiscreteEvent event,
                                             DiscreteState source,
                                             DiscreteState target,
                                             const VectorFunction& reset,
                                             const VectorFunction& activation,
                                             bool forced);


    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //!
    //!    \param source is the transition's source location.
    //!    \param event is the transition's event.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param guard is the transition's activation region.
    const AtomicDiscreteTransition& new_transition(DiscreteState source,
                                             DiscreteEvent event,
                                             DiscreteState target,
                                             const List<RealUpdateAssignment>& reset,
                                             const ContinuousPredicate& guard);

    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //! The guard is the constant "True" i.e. the event is an input event.
    //!
    //!    \param source is the transition's source location.
    //!    \param event is the transition's event.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    const AtomicDiscreteTransition& new_transition(DiscreteState source,
                                             DiscreteEvent event,
                                             DiscreteState target,
                                             const List<RealUpdateAssignment>& reset);

     //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //! The reset is trivial. This form is for the case that there are no continuous state variables in the new location.
    //!
    //!    \param source is the transition's source location.
    //!    \param event is the transition's event.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param guard is the transition's activation region.
    const AtomicDiscreteTransition& new_transition(DiscreteState source,
                                             DiscreteEvent event,
                                             DiscreteState target,
                                             const ContinuousPredicate& guard);

   //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param activation is the transition's activation region.
    //!    \param forced determines whether the transision is forced (urgent) or unforced (permissive).
    //! As there is no reset parameter, the target location cannot have any state variables.
    const AtomicDiscreteTransition& new_transition(DiscreteEvent event,
                                             DiscreteState source,
                                             DiscreteState target,
                                             const ContinuousPredicate& guard,
                                             bool urgent);

    //! \brief Adds a forced (urgent) discrete transition to the automaton
    //! using the discrete states to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    const AtomicDiscreteTransition& new_forced_transition(DiscreteEvent event,
                                                    DiscreteState source,
                                                    DiscreteState target,
                                                    const VectorFunction& reset,
                                                    const VectorFunction& activation);

    //! \brief Adds an unforced (non-urgent) discrete transition to the automaton
    //! using the discrete states to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    const AtomicDiscreteTransition& new_unforced_transition(DiscreteEvent event,
                                                      DiscreteState source,
                                                      DiscreteState target,
                                                      const VectorFunction& reset,
                                                      const VectorFunction& activation);

    //! \brief Adds a discrete transition to the automaton using the discrete modes to specify the source and target.
    //!
    //!    \param event is the discrete transition's discrete event.
    //!    \param source is the discrete transition's source mode.
    //!    \param target is the discrete transition's target mode.
    //!    \param reset is the discrete transition's reset.
    //!    \param activation is the discrete transition's activation region.
    //!    \param forced determines whether the transition is forced or unforced.
    const AtomicDiscreteTransition& new_transition(DiscreteEvent event,
                                             const AtomicDiscreteMode& source,
                                             const AtomicDiscreteMode& target,
                                             const VectorFunction& reset,
                                             const VectorFunction& activation,
                                             bool forced);

    //! \brief Set the grid controlling relative scaling in the mode.
    void set_grid(DiscreteState location, const Grid& grid);

    //! \brief Set the grid controlling relative scaling. This method sets the same grid for every mode.
    void set_grid(const Grid& grid);

    //! \brief Set the hybrid grid controlling relative scaling.
    void set_grid(const HybridGrid& hgrid);

    //@}

    //@{
    //! \name Data access and queries.

    //! \brief Returns the hybrid automaton's name.
    const std::string& name() const;

    //! \brief Test if the hybrid automaton has a discrete mode with discrete state \a state.
    bool has_mode(DiscreteState state) const;

    //! \brief Test if the hybrid automaton has a discrete transition with \a event_id and \a source_id.
    bool has_transition(DiscreteEvent event, DiscreteState source) const;
    bool has_transition(DiscreteState source, DiscreteEvent event) const;

    //! \brief The discrete mode with given discrete state.
    const AtomicDiscreteMode& mode(DiscreteState state) const;

    //! \brief The discrete transition with given \a event and \a source location.
    const AtomicDiscreteTransition& transition(DiscreteEvent event, DiscreteState source) const;
    const AtomicDiscreteTransition& transition(DiscreteState source, DiscreteEvent event) const;

    //! \brief The set of discrete modes. (Not available in Python interface)
    const std::set< AtomicDiscreteMode >& modes() const;

    //! \brief The set of discrete transitions. (Not available in Python interface)
    const std::set< AtomicDiscreteTransition >& transitions() const;

    //! \brief The discrete transitions from location \a source.
    std::set< AtomicDiscreteTransition > transitions(DiscreteState source) const;

    //! \brief The blocking events (invariants and urgent transitions) in \a location.
    std::map<DiscreteEvent,ScalarFunction> blocking_guards(DiscreteState location) const;

    //! \brief The permissive events (invariants and urgent transitions) in \a location.
    std::map<DiscreteEvent,ScalarFunction> permissive_guards(DiscreteState location) const;

    //! \brief The state space of the system.
    HybridSpace state_space() const;

    //! \brief The hybrid set giving the invariants for each discrete location.
    HybridSet invariant() const;

    //! \brief The natural grid to use in the specified location.
    Grid grid(DiscreteState location) const;

    //! \brief The natural grid to use in the over all locations.
    HybridGrid grid() const;
    //@}

    //@{
    //! \name New-style access for compositional hybrid automata.

    //! \brief The target location of the discrete event from the given discrete location.
    DiscreteState target(const DiscreteState& source, const DiscreteEvent& event);
    // Alternative form for backwards compatibility
    DiscreteState target(const DiscreteEvent& event, const DiscreteState& source);
    //! \brief The state (dotted) variables in the given location.
    List<RealVariable> state_variables(DiscreteState location) const;
    //! \brief The auxiliary (algebraic/output) variables in the given location.
    List<RealVariable> auxiliary_variables(DiscreteState location) const;
    //! \brief The algebraic equations valid in the given location.
    List<RealAssignment> algebraic_assignments(const DiscreteState& location) const;
    //! \brief The differential equations valid in the given location.
    List<DottedRealAssignment> differential_assignments(const DiscreteState& location) const;
    //! \brief The differential equations valid in the given location.
    List<PrimedRealAssignment> update_assignments(const DiscreteState& source, const DiscreteEvent& event) const;
    //! \brief The invariant predicates valid in the given location.
    Map<DiscreteEvent,ContinuousPredicate> invariant_predicates(const DiscreteState& location) const;
    ContinuousPredicate invariant_predicate(const DiscreteState& location, const DiscreteEvent& action) const;
    //! \brief The guard predicate for the given event in the given location.
    ContinuousPredicate guard_predicate(const DiscreteState& location, const DiscreteEvent& event) const;
    //@}

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream&) const;
};

inline std::ostream& operator<<(std::ostream& os, const AtomicHybridAutomaton& ha) {
    return ha.write(os); }



class CompositeHybridAutomaton {
    typedef List<DiscreteState> DiscreteLocation;
  public:
    CompositeHybridAutomaton(const AtomicHybridAutomaton&);
    CompositeHybridAutomaton(const List<AtomicHybridAutomaton>&);
    uint number_of_components() const;
    const AtomicHybridAutomaton& component(uint) const;

    const Set<AtomicDiscreteMode> modes() const;
    const Set<AtomicDiscreteTransition> transitions(const DiscreteLocation& location) const;

    bool has_mode(const DiscreteLocation&) const;
    bool has_transition(const DiscreteLocation&, const DiscreteEvent&) const;
    DiscreteLocation target(const DiscreteLocation&, const DiscreteEvent&) const;

    List<RealVariable> variables(const DiscreteLocation&) const;
    List<RealVariable> state_variables(const DiscreteLocation&) const;
    List<RealVariable> auxiliary_variables(const DiscreteLocation&) const;
    List<RealAssignment> algebraic_assignments(const DiscreteLocation&) const;
    List<DottedRealAssignment> differential_assignments(const DiscreteLocation&) const;
    List<PrimedRealAssignment> update_assignments(const DiscreteLocation&, const DiscreteEvent&) const;

    Map<DiscreteEvent,ContinuousPredicate> invariant_predicates(const DiscreteLocation&) const;
    ContinuousPredicate invariant_predicate(const DiscreteLocation&, const DiscreteEvent&) const;
    ContinuousPredicate guard_predicate( const DiscreteLocation&, const DiscreteEvent&) const;

    VectorFunction output(const DiscreteLocation&) const;
    VectorFunction dynamic(const DiscreteLocation&) const;
    VectorFunction reset(const DiscreteLocation&, const DiscreteEvent&) const;

    Map<DiscreteEvent,ScalarFunction> invariants(const DiscreteLocation&) const;
    ScalarFunction invariant(const DiscreteLocation&, const DiscreteEvent&) const;
    ScalarFunction guard(const DiscreteLocation&, const DiscreteEvent&) const;

    std::ostream& write(std::ostream&) const;
  private:
    List<AtomicHybridAutomaton> _components;
};

inline CompositeHybridAutomaton parallel_composition(const List<AtomicHybridAutomaton>& components) {
    return CompositeHybridAutomaton(components); }

inline std::ostream& operator<<(std::ostream& os, const CompositeHybridAutomaton& ha) {
    return ha.write(os); }



} // namespace Ariadne

#endif // ARIADNE_HYBRID_AUTOMATON_H
