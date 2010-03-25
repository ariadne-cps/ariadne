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
class HybridSet;
class HybridGrid;

class DiscreteMode;
class DiscreteTransition;
class MonolithicHybridAutomaton;

class ScalarFunction;
class VectorFunction;
class Grid;


/*! \brief A discrete mode of a hybrid automaton, comprising continuous evolution given by a vector field
 * within and invariant constraint set.
 *
 * A %DiscreteMode can only be created using the new_mode() method in
 * the %MonolithicHybridAutomaton class.
 *
 * \sa \link Ariadne::MonolithicHybridAutomaton \c MonolithicHybridAutomaton \endlink, \link Ariadne::DiscreteTransition \c DiscreteTransition \endlink
 */
class DiscreteMode {
    friend class MonolithicHybridAutomaton;
  private:

    // The discrete mode's discrete state.
    DiscreteLocation _location;

    // The discrete mode's vector field.
    VectorFunction _dynamic;

    // The discrete mode's invariants.
    Map< DiscreteEvent, ScalarFunction > _invariants;

    // The discrete mode's grid for reachability analysis.
    boost::shared_ptr< const Grid > _grid;
  public:
    //! \brief The mode's discrete state.
    DiscreteLocation location() const {
        return this->_location; }

    //! \brief The discrete mode's dynamic (a vector field).
    const VectorFunction& dynamic() const {
        return this->_dynamic; }

    //! \brief The discrete mode's invariants.
    const Map< DiscreteEvent, ScalarFunction >& invariants() const {
        return this->_invariants; }

    const ScalarFunction& invariant(const DiscreteEvent& event) const {
        return this->_invariants.find(event)->second; }

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
    DiscreteMode(DiscreteLocation location,
                 const VectorFunction& dynamic);

    // Construct from objects managed by shared pointers (for internal use)
    DiscreteMode(DiscreteLocation location,
                 const VectorFunction dynamic,
                 const std::vector< ScalarFunction >& invariants);

};


std::ostream& operator<<(std::ostream& os, const DiscreteMode& dm);

inline bool operator<(const DiscreteMode& mode1, const DiscreteMode& mode2) {
    return mode1.location() < mode2.location(); }




/*! \brief A discrete transition of a hybrid automaton, representing an instantaneous
 * jump from one discrete mode to another, governed by an activation set and a reset map.
 *
 * A %DiscreteTransition can only be created using the new_transition() method in
 * the %MonolithicHybridAutomaton class.
 *
 * An invariant is modelled by a discrete transition with negative event id and null reset pointer.
 *
 * \sa \link Ariadne::MonolithicHybridAutomaton \c MonolithicHybridAutomaton \endlink, \link Ariadne::DiscreteMode \c DiscreteMode \endlink
 */
class DiscreteTransition
{
    friend class MonolithicHybridAutomaton;
  private:
    // \brief The discrete transition's identificator.
    DiscreteEvent _event;

    // \brief The source of the discrete transition.
    const DiscreteMode* _source;

    // \brief The target of the discrete transition.
    const DiscreteMode* _target;

    // \brief The activation region of the discrete transition.
    ScalarFunction _activation;

    // \brief The reset of the discrete transition.
    VectorFunction _reset;

    // \brief Whether or not the transition is urgent or permissive.
    Urgency _urgency;

  public:

    //! \brief The discrete event associated with the discrete transition.
    DiscreteEvent event() const {
        return this->_event; }

    //! \brief The source mode of the discrete transition.
    const DiscreteMode& source_mode() const {
        return *this->_source; }

    //! \brief The target of the discrete transition.
    const DiscreteMode& target_mode() const {
        return *this->_target; }


    //! \brief The source mode of the discrete transition.
    DiscreteLocation source() const {
        return this->_source->location(); }

    //! \brief The target of the discrete transition.
    DiscreteLocation target() const {
        return this->_target->location(); }


    //! \brief The activation constraint function of the discrete transition.
    const ScalarFunction& activation() const {
        return this->_activation;
    }

    //! \brief The guard/activation constraint function of the discrete transition.
    const ScalarFunction& guard() const {
        return this->_activation;
    }

    //! \brief The reset map of the discrete transition.
    const VectorFunction& reset() const {
        return this->_reset;
    }

    //! \brief True if the transition is forced (occurs as soon as it is activated).
    bool forced() const {
        return this->_urgency;
    }

    //! \brief True if the transition is forced (occurs as soon as it is activated).
    Urgency urgency() const {
        return this->_urgency;
    }

  private:


    // Construct from shared pointers (for internal use).
    DiscreteTransition(DiscreteEvent event,
                       const DiscreteMode& source,
                       const DiscreteMode& target,
                       const VectorFunction& reset,
                       const ScalarFunction& activation,
                       Urgency urgency);

};

std::ostream& operator<<(std::ostream& os, const DiscreteTransition& dt);

inline bool operator<(const DiscreteTransition& transition1, const DiscreteTransition& transition2) {
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
 * For each %DiscreteMode, the dynamics is given by a
 * %VectorField describing the continuous dynamics,
 * and a %Set giving an invariants which must be satisified at
 * all times.
 *
 * The discrete time behaviour is specified by %DiscreteTransition
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
 * \sa \link Ariadne::DiscreteMode \c DiscreteMode \endlink, \link Ariadne::DiscreteTransition \c DiscreteTransition \endlink

 */
class MonolithicHybridAutomaton
    : public HybridAutomatonInterface
{
  public:
    //! \brief The type used to represent time.
    typedef HybridTime TimeType;
    //! \brief The type used to represent real numbers.
    typedef double RealType ;
    //! \brief The type used to describe the state space.
    typedef HybridSpace StateSpaceType;


    typedef Map<DiscreteEvent,VectorFunction>::const_iterator invariant_const_iterator;
    typedef Set<DiscreteTransition>::const_iterator discrete_transition_const_iterator;
    typedef Set<DiscreteMode>::const_iterator discrete_mode_const_iterator;
  private:
    //! \brief The hybrid automaton's name.
    std::string _name;

    //! \brief The list of the hybrid automaton's discrete modes.
    Set< DiscreteMode > _modes;

    //! \brief The hybrid automaton's transitions.
    Set< DiscreteTransition > _transitions;

  public:
    //@{
    //! \name Constructors and destructors

    //! \brief Construct an empty automaton with no name
    MonolithicHybridAutomaton();

    //! \brief Construct an empty automaton with the given name
    MonolithicHybridAutomaton(const std::string& name);

    //! \brief Construct dynamically-allocated copy. (Not currently implemented)
    MonolithicHybridAutomaton* clone() const;

    //! \brief  Destructor.
    ~MonolithicHybridAutomaton();
    //@}

    //@{
    //! \name Methods for building the automaton.

    //! \brief Adds a discrete mode to the automaton.
    //!
    //!   \param state is the mode's discrete state.
    //!   \param dynamic is the mode's vector field.
    const DiscreteMode& new_mode(DiscreteLocation state,
                                 const VectorFunction& dynamic);

    //! \brief Adds an invariant to a mode of the automaton.
    //!
    //!   \param state is the mode's discrete state.
    //!   \param invariant is the new invariant condition, in the form \f$g(x)<0\f$.

    const DiscreteMode& new_invariant(DiscreteLocation state,
                                      const ScalarFunction& invariant);

    //! \brief Adds an invariants to a mode of the automaton.
    //!
    //!   \param state is the mode's discrete state.
    //!   \param invariants is the new invariants condition.

    const DiscreteMode& new_invariant(DiscreteLocation state,
                                      const VectorFunction& invariants);

    //! \brief Adds an invariants to a mode of the automaton.
    //!
    //!    \param mode is the discrete mode.
    //!    \param invariants is the new invariants condition.

    const DiscreteMode& new_invariant(const DiscreteMode& mode,
                                      const VectorFunction& invariants);


    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    //!    \param forced determines whether the transision is forced (urgent) or unforced (permissive).
    const DiscreteTransition& new_transition(DiscreteEvent event,
                                             DiscreteLocation source,
                                             DiscreteLocation target,
                                             const VectorFunction& reset,
                                             const ScalarFunction& activation,
                                             Urgency urgency);

    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    //!    \param forced determines whether the transision is forced (urgent) or unforced (permissive).
    const DiscreteTransition& new_transition(DiscreteEvent event,
                                             DiscreteLocation source,
                                             DiscreteLocation target,
                                             const VectorFunction& reset,
                                             const VectorFunction& activation,
                                             Urgency urgency);

    //! \brief Adds a forced (urgent) discrete transition to the automaton
    //! using the discrete states to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    const DiscreteTransition& new_forced_transition(DiscreteEvent event,
                                                    DiscreteLocation source,
                                                    DiscreteLocation target,
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
    const DiscreteTransition& new_unforced_transition(DiscreteEvent event,
                                                      DiscreteLocation source,
                                                      DiscreteLocation target,
                                                      const VectorFunction& reset,
                                                      const VectorFunction& activation);
/*
    //! \brief Adds a discrete transition to the automaton using the discrete modes to specify the source and target.
    //!
    //!    \param event is the discrete transition's discrete event.
    //!    \param source is the discrete transition's source mode.
    //!    \param target is the discrete transition's target mode.
    //!    \param reset is the discrete transition's reset.
    //!    \param activation is the discrete transition's activation region.
    //!    \param forced determines whether the transition is forced or unforced.
    const DiscreteTransition& new_transition(DiscreteEvent event,
                                             const DiscreteMode& source,
                                             const DiscreteMode& target,
                                             const VectorFunction& reset,
                                             const VectorFunction& activation,
                                             bool forced);
*/

    //! \brief Set the grid controlling relative scaling in the mode.
    void set_grid(DiscreteLocation location, const Grid& grid);

    //! \brief Set the grid controlling relative scaling. This method sets the same grid for every mode.
    void set_grid(const Grid& grid);

    //! \brief Set the hybrid grid controlling relative scaling.
    void set_grid(const HybridGrid& hgrid);
    //@}

    //@{
    //! \name Data access and queries.

    //! \brief Returns the hybrid automaton's name.
    const std::string& name() const;

    //! \brief Test if the hybrid automaton has a discrete mode \a location.
    bool has_mode(DiscreteLocation location) const;

    //! \brief Test if the hybrid automaton has a discrete transition with \a event_id and \a source_id.
    bool has_transition(DiscreteLocation source, DiscreteEvent event) const;
    bool has_transition(DiscreteEvent event, DiscreteLocation source) const;

    //! \brief The discrete mode with location \a location. \deprecated
    const DiscreteMode& mode(DiscreteLocation location) const;

    //! \brief The discrete transition with given \a event and \a source location. \deprecated
    const DiscreteTransition& transition(DiscreteLocation source, DiscreteEvent event) const;
    const DiscreteTransition& transition(DiscreteEvent event, DiscreteLocation source) const;

    //! \brief The dimension of the state space of in the mode \a location.
    uint dimension(DiscreteLocation location) const;

    //! \brief The dynamic valid in the mode \a location.
    VectorFunction dynamic(DiscreteLocation location) const;

    //! \brief The invariants valid in the mode \a location.
    Map<DiscreteEvent,ScalarFunction> invariants(DiscreteLocation location) const;

    //! \brief The guards active in the mode \a location.
    Map<DiscreteEvent,ScalarFunction> guards(DiscreteLocation location) const;
    Map<DiscreteEvent,ScalarFunction> activations(DiscreteLocation location) const;

    //! \brief The target location of \a event starting in the \a source location.
    DiscreteLocation target(DiscreteLocation source, DiscreteEvent event) const;
    Map<DiscreteEvent,DiscreteLocation> targets(DiscreteLocation source) const;

    //! \brief The dynamic valid in the mode \a location.
    VectorFunction reset(DiscreteLocation source, DiscreteEvent event) const;
    Map<DiscreteEvent,VectorFunction> resets(DiscreteLocation source) const;

    //! \brief The set of discrete modes. (Not available in Python interface)
    const Set< DiscreteMode >& modes() const;

    //! \brief The set of discrete transitions. (Not available in Python interface)
    const Set< DiscreteTransition >& transitions() const;

    //! \brief The discrete transitions from location \a source.
    Set< DiscreteTransition > transitions(DiscreteLocation source) const;

    //! \brief The blocking events (invariants and urgent transitions) in \a location.
    Map<DiscreteEvent,ScalarFunction> blocking_guards(DiscreteLocation location) const;

    //! \brief The permissive events (invariants and urgent transitions) in \a location.
    Map<DiscreteEvent,ScalarFunction> permissive_guards(DiscreteLocation location) const;

    //! \brief The state space of the system.
    HybridSpace state_space() const;

    //! \brief The hybrid set giving the invariants for each discrete location.
    HybridSet invariant() const;

    //! \brief The natural grid to use in the specified location.
    Grid grid(DiscreteLocation location) const;

    //! \brief The natural grid to use in the over all locations.
    HybridGrid grid() const;

    virtual Set<DiscreteEvent> urgent_events(DiscreteLocation) const;
    virtual Set<DiscreteEvent> permissive_events(DiscreteLocation) const;
    virtual Set<DiscreteEvent> blocking_events(DiscreteLocation) const;
    virtual Set<DiscreteEvent> invariant_events(DiscreteLocation) const;
    virtual Set<DiscreteEvent> transition_events(DiscreteLocation) const;
    virtual VectorFunction dynamic_function(DiscreteLocation) const;
    virtual VectorFunction reset_function(DiscreteLocation, DiscreteEvent) const;
    virtual ScalarFunction guard_function(DiscreteLocation, DiscreteEvent) const;
    virtual ScalarFunction invariant_function(DiscreteLocation, DiscreteEvent) const;
    //@}

};

std::ostream& operator<<(std::ostream& os, const MonolithicHybridAutomaton& ha);


} // namespace Ariadne

#endif // ARIADNE_HYBRID_AUTOMATON_H
