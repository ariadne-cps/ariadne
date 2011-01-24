/***************************************************************************
 *            hybrid_io_automaton.h
 *
 *  Copyright  2010  Davide Bresolin
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

/*! \file hybrid_io_automaton.h
 *  \brief Main hybrid system class.
 */

#ifndef ARIADNE_HYBRID_IO_AUTOMATON_H
#define ARIADNE_HYBRID_IO_AUTOMATON_H

#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "expression.h"
#include "function.h"
#include "discrete_state.h"
#include "discrete_event.h"
#include "variables.h"

namespace Ariadne {

class HybridAutomaton;

/*! \brief A discrete mode of a hybrid I/O automaton, comprising continuous evolution given by a vector field
 * within and invariant constraint set.
 *
 * A %DiscreteIOMode can only be created using the new_mode() method in
 * the %HybridIOAutomaton class.
 *
 * \sa \link Ariadne::HybridIOAutomaton \c HybridIOAutomaton \endlink, \link Ariadne::DiscreteIOTransition \c DiscreteIOTransition \endlink
 */
class DiscreteIOMode {
    friend class HybridIOAutomaton;
  private:

    // The discrete mode's discrete state.
    DiscreteState _location;

    // The discrete mode's vector field, described by a RealExpression for each variable.
    std::map< RealVariable, RealExpression > _dynamics;
    
    // The discrete mode's invariants.
    std::list< RealExpression > _invariants;

  public:
    //! \brief The mode's discrete state.
    DiscreteState location() const {
        return this->_location; }

    //! \brief Returns true if the dynamic for a given variable is defined in the mode.
    bool has_dynamics(const RealVariable& var) {
        return (this->_dynamics.find(var) != this->_dynamics.end());
    }

    //! \brief The discrete mode's dynamic for a given variable.
    const RealExpression& dynamics(const RealVariable& var);

    //! \brief The discrete mode's dynamic for all variables.
    const std::map< RealVariable, RealExpression >& dynamics() const {
        return this->_dynamics;
    }

    //! \brief The discrete mode's invariants.
    const std::list< RealExpression >& invariants() const {
        return this->_invariants; }
   
	//! \brief Substitute the constant \a c into the corresponding Constant \a con, if present, on the invariants and dynamic functions.
	void substitute(const Constant<Real>& con, const Real& c);

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream& os) const;

    // Construct discrete mode.
    //
    // \param id is the identifier of the mode.
    // \param dynamic is the mode's dynamic.
    // \param invariants is the mode's invariants.
    DiscreteIOMode(DiscreteState location);

    DiscreteIOMode(DiscreteState location,
                 const std::map< RealVariable, RealExpression >& dynamic);

    DiscreteIOMode(DiscreteState location,
                 const std::map< RealVariable, RealExpression >& dynamic,
                 const std::list< RealExpression >& invariants);

    // Set the dynamic of variable var.
    void set_dynamics(const RealVariable& var,
                     const RealExpression& dyn);
   
    // Add a new invariant to the discrete mode.
    void add_invariant(const RealExpression& inv);

};


std::ostream& operator<<(std::ostream& os, const DiscreteIOMode& dm);



/*! \brief A discrete transition of a hybrid I/O automaton, representing an instantaneous
 * jump from one discrete mode to another, governed by an activation set and a reset map.
 *
 * A %DiscreteIOTransition can only be created using the new_transition() method in
 * the %HybridIOAutomaton class.
 *
 * An invariant is modelled by a discrete transition with negative event id and null reset pointer.
 *
 * \sa \link Ariadne::HybridIOAutomaton \c HybridIOAutomaton \endlink, \link Ariadne::DiscreteIOMode \c DiscreteIOMode \endlink
 */
class DiscreteIOTransition
{
    friend class HybridIOAutomaton;
  private:
    // \brief The discrete transition's identificator.
    DiscreteEvent _event;

    // \brief The source of the discrete transition.
    DiscreteState _source;

    // \brief The target of the discrete transition.
    DiscreteState _target;

    // \brief The activation region of the discrete transition.
    RealExpression _activation;

    // \brief The reset of the discrete transition, defined by a RealExpression for each variable.
    std::map< RealVariable, RealExpression > _reset;

    // \brief Whether or not the transition is forced.
    bool _forced;

  public:
    //! \brief The default constructor.
    DiscreteIOTransition();

    //! \brief The discrete event associated with the discrete transition.
    DiscreteEvent event() const {
        return this->_event; }

    //! \brief The source mode of the discrete transition.
    DiscreteState source() const {
        return this->_source; }

    //! \brief The target of the discrete transition.
    DiscreteState target() const {
        return this->_target; }

	//! \brief Substitute the constant \a c into the corresponding Constant \a con, if present, on the reset and activation functions.
	void substitute(const Constant<Real>& con, const Real& c);

    //! \brief The activation region of the discrete transition.
    const RealExpression& activation() const {
        return this->_activation;
    }

    //! \brief Returns true if the reset map is defined for a given variable
    bool has_reset(const RealVariable& var) {
        return ( this->_reset.find(var) != this->_reset.end() );
    }

    //! \brief The reset map of the discrete transition for a given variable
    const RealExpression& reset(const RealVariable& var); 
    
    //! \brief The reset map of the discrete transition for all variables.
    const std::map< RealVariable, RealExpression >& reset() const {
        return this->_reset;
    }

    //! \brief True if the transition is forced (occurs as soon as it is activated).
    bool forced() const {
        return this->_forced;
    }

    // Construct a new transition (for internal use).    
    DiscreteIOTransition(DiscreteEvent event,
                         DiscreteState source,
                         DiscreteState target,
                         bool forced=false);

    DiscreteIOTransition(DiscreteEvent event,
                         DiscreteState source,
                         DiscreteState target,
                         const RealExpression& activation,
                         bool forced=false);

    DiscreteIOTransition(DiscreteEvent event,
                         DiscreteState source,
                         DiscreteState target,
                         const std::map< RealVariable, RealExpression >& reset,
                         bool forced=false);

    DiscreteIOTransition(DiscreteEvent event,
                         DiscreteState source,
                         DiscreteState target,
                         const std::map< RealVariable, RealExpression >& reset,
                         const RealExpression& activation,
                         bool forced=false);
                         
    void set_event(DiscreteEvent event);
    
    void set_source(DiscreteState source);
    
    void set_target(DiscreteState target);

    // Set the reset map for a given variable
    void set_reset(const RealVariable& var,
                   const RealExpression& res);

    void set_reset(const std::map< RealVariable, RealExpression >& reset);

    // Set the activation region of the transition
    void set_activation(const RealExpression& act);    
};

std::ostream& operator<<(std::ostream& os, const DiscreteIOTransition& dt);


/*! \brief A hybrid I/O automaton, comprising continuous-time behaviour
 *  at each discrete mode, coupled by instantaneous discrete transitions.
 *
 * A hybrid automaton is a dynamic system with evolution in both
 * continuous time and discrete time, and is defined by:
 *
 * A set of input, output, and internal %RealVariable,
 *
 * A set of input, output, and internal %DiscreteEvent,
 *
 * A set of %DiscreteIOMode, 
 * 
 * A set of %DiscreteIOTransition.
 *
 * For each %DiscreteIOMode, the dynamics is given by a %RealExpression for each controlled 
 * (i.e. output or internal) variable describing the continuous dynamics,
 * and a set of %RealExpression giving an invariants which must be satisified at
 * all times.
 *
 * The discrete time behaviour is specified by %DiscreteIOTransition
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
 * \sa \link Ariadne::DiscreteIOMode \c DiscreteIOMode \endlink, \link Ariadne::DiscreteIOTransition \c DiscreteIOTransition \endlink

 */
class HybridIOAutomaton
{
  public:
    typedef std::list<DiscreteIOTransition>::const_iterator discrete_transition_const_iterator;
    typedef std::list<DiscreteIOMode>::const_iterator discrete_mode_const_iterator;
  private:
    //! \brief The hybrid automaton's name.
    std::string _name;

    //! \brief The sets of input, output, and internal variables.
    std::set< RealVariable > _input_vars;
    std::set< RealVariable > _output_vars;
    std::set< RealVariable > _internal_vars;

    //! \brief The sets of input, output, and internal events.
    std::set< DiscreteEvent > _input_events;
    std::set< DiscreteEvent > _output_events;
    std::set< DiscreteEvent > _internal_events;

    //! \brief The list of the hybrid automaton's discrete modes.
    std::list< DiscreteIOMode > _modes;

    //! \brief The hybrid automaton's transitions.
    std::list< DiscreteIOTransition > _transitions;
    
    //! \brief The accessible constants.
    //! \details This set does not necessarily reflect the whole set of constants in all functions.
    std::list< RealConstant > _accessible_constants;

    //! \brief The grid for the controlled variables of the automaton.
    std::map< RealVariable, Float > _grid;
    
    //! \brief Access to a discrete mode (for internal use only)
    DiscreteIOMode& _mode(DiscreteState state);

    //! \brief Access to a discrete transition (for internal use only)   
    DiscreteIOTransition& _transition(DiscreteEvent event, DiscreteState state);

  public:
    //@{
    //! \name Constructors and destructors

    //! \brief Construct an empty automaton with no name
    HybridIOAutomaton();

    //! \brief Construct an empty automaton with the given name
    HybridIOAutomaton(const std::string& name);

    //! \brief Construct an empty automaton with the given name and with the given sets of variables
    HybridIOAutomaton(const std::string& name,
                      const std::set< RealVariable >& input_vars,
                      const std::set< RealVariable >& output_vars,
                      const std::set< RealVariable >& internal_vars);


    //! \brief Construct an empty automaton with the given name and with the given sets of variables and events
    HybridIOAutomaton(const std::string& name,
                      const std::set< RealVariable >& input_vars,
                      const std::set< RealVariable >& output_vars,
                      const std::set< RealVariable >& internal_vars,
                      const std::set< DiscreteEvent >& input_events,
                      const std::set< DiscreteEvent >& output_events,
                      const std::set< DiscreteEvent >& internal_events);

    //@}

    //@{
    //! \name Methods for building the automaton.

    //! \brief Change the automaton's name.
    const std::string& set_name(const std::string& name);

    //! \brief Adds an input variable \param u to the automaton.
    const std::set< RealVariable >& add_input_var(const RealVariable& u);

    //! \brief Adds an output variable \param y to the automaton.
    const std::set< RealVariable >& add_output_var(const RealVariable& y);

    //! \brief Adds an internal variable \param x to the automaton.
    const std::set< RealVariable >& add_internal_var(const RealVariable& x);

    //! \brief Adds an input event \param e to the automaton.
    const std::set< DiscreteEvent >& add_input_event(const DiscreteEvent& e);

    //! \brief Adds an output event \param e to the automaton.
    const std::set< DiscreteEvent >& add_output_event(const DiscreteEvent& e);

    //! \brief Adds an internal event \param e to the automaton.
    const std::set< DiscreteEvent >& add_internal_event(const DiscreteEvent& e);
    
    //! \brief Registers a constant into the list of accessible constants.
    //! \details If the constant is already registered, no action is performed.
    void register_accessible_constant(RealConstant c);

    //! \brief Adds a discrete mode to the automaton.
    //!
    //!   \param state is the mode's discrete state.
    //!   \param dynamic is the mode's vector field.
    const DiscreteIOMode& new_mode(DiscreteState state);
    
    const DiscreteIOMode& new_mode(DiscreteState state,
                                   const std::map< RealVariable, RealExpression >& dynamic);

    const DiscreteIOMode& new_mode(const DiscreteIOMode& mode);
    
    //! \brief Sets the dynamic for a variable in a mode of the automaton.
    //!
    //!   \param state is the mode's discrete state.
    //|   \param var is the automaton's controlled variable.
    //!   \param dynamic is the dynamics for var.
    const DiscreteIOMode& set_dynamics(DiscreteState state,
                                      const RealVariable& var,
                                      const RealExpression& dynamics);

    //! \brief Adds an invariant to a mode of the automaton.
    //!
    //!   \param state is the mode's discrete state.
    //!   \param invariant is the new invariant condition, in the form \f$g(x)<0\f$.
    const DiscreteIOMode& new_invariant(DiscreteState state,
                                        const RealExpression& invariant);

    //! \brief Adds a discrete transition to the automaton.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    //!    \param forced determines whether the transision is forced (urgent) or unforced (permissive).
    const DiscreteIOTransition& new_transition(DiscreteEvent event,
                                               DiscreteState source,
                                               DiscreteState target,
                                               const std::map< RealVariable, RealExpression >& reset,
                                               const RealExpression& activation,
                                               bool forced);

    const DiscreteIOTransition& new_transition(DiscreteEvent event,
                                               DiscreteState source,
                                               DiscreteState target,
                                               const std::map< RealVariable, RealExpression >& reset,
                                               bool forced);

    const DiscreteIOTransition& new_transition(DiscreteEvent event,
                                               DiscreteState source,
                                               DiscreteState target,
                                               const RealExpression& activation,
                                               bool forced);

    const DiscreteIOTransition& new_transition(DiscreteEvent event,
                                               DiscreteState source,
                                               DiscreteState target,
                                               bool forced);

    const DiscreteIOTransition& new_transition(const DiscreteIOTransition& trans);
    
    //! \brief Adds a forced (urgent) discrete transition to the automaton.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    const DiscreteIOTransition& new_forced_transition(DiscreteEvent event,
                                                      DiscreteState source,
                                                      DiscreteState target,
                                                      const std::map< RealVariable, RealExpression >& reset,
                                                      const RealExpression& activation) 
    {
        return this->new_transition(event, source, target, reset, activation, true);
    }
 
    const DiscreteIOTransition& new_forced_transition(DiscreteEvent event,
                                                      DiscreteState source,
                                                      DiscreteState target,
                                                      const std::map< RealVariable, RealExpression >& reset) 
    {
        return this->new_transition(event, source, target, reset, true);
    }
 
    const DiscreteIOTransition& new_forced_transition(DiscreteEvent event,
                                                      DiscreteState source,
                                                      DiscreteState target,
                                                      const RealExpression& activation)
    {
        return this->new_transition(event, source, target, activation, true);
    }

    const DiscreteIOTransition& new_forced_transition(DiscreteEvent event,
                                                      DiscreteState source,
                                                      DiscreteState target)
    {
        return this->new_transition(event, source, target, true);
    }

    //! \brief Adds an unforced (non-urgent) discrete transition to the automaton.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    const DiscreteIOTransition& new_unforced_transition(DiscreteEvent event,
                                                        DiscreteState source,
                                                        DiscreteState target,
                                                        const std::map< RealVariable, RealExpression >& reset,
                                                        const RealExpression& activation)
    {
        return this->new_transition(event, source, target, reset, activation, false);
    }

    const DiscreteIOTransition& new_unforced_transition(DiscreteEvent event,
                                                        DiscreteState source,
                                                        DiscreteState target,
                                                        const std::map< RealVariable, RealExpression >& reset)
    {
        return this->new_transition(event, source, target, reset, false);
    }

    const DiscreteIOTransition& new_unforced_transition(DiscreteEvent event,
                                                        DiscreteState source,
                                                        DiscreteState target,
                                                        const RealExpression& activation)
    {
        return this->new_transition(event, source, target, activation, false);
    }

    const DiscreteIOTransition& new_unforced_transition(DiscreteEvent event,
                                                        DiscreteState source,
                                                        DiscreteState target)
    {
        return this->new_transition(event, source, target, false);
    }

    //! \brief Set the reset function for a discrete transition in the automaton.
    //!
    //!    \param event is the transition's discrete event.
    //!    \param source is the transition's source mode.
    //!    \param reset is a map that gives the reset function for each variable.
    //!    
    const DiscreteIOTransition& set_reset(DiscreteEvent event,
                                          DiscreteState source,
                                          const std::map< RealVariable, RealExpression >& reset);

    //! \brief Set the reset function for a variable in a discrete transition of the automaton.
    //!
    //!    \param event is the transition's discrete event.
    //!    \param source is the transition's source mode.
    //!    \param var is the variable
    //!    \param reset is the reset function for the given variable.
    //!    
    const DiscreteIOTransition& set_reset(DiscreteEvent event,
                                          DiscreteState source,
                                          const RealVariable& var,
                                          const RealExpression& reset);
                                          
    //! \brief Set the activation region for a discrete transition in the automaton.
    //!
    //!    \param event is the transition's discrete event.
    //!    \param source is the transition's source mode.
    //!    \param activation is the activation region of the transition.
    //!    
    const DiscreteIOTransition& set_activation(DiscreteEvent event,
                                               DiscreteState source,
                                               const RealExpression& activation);                                              
                                              
    //! \brief Set the grid controlling relative scaling. This method sets the same grid for all controlled variables.
    void set_grid(const std::map< RealVariable, Float>& grid);

    //! \brief Set the grid controlling relative scaling for a given variable.
    void set_grid(const RealVariable& var, Float scaling);

	//! \brief Substitute the constant \a c into the corresponding Constant \a con, if present, on all the functions of modes and transitions.
	void substitute(const Constant<Real>& con, const Real& c);

	//! \brief Substitute the value of the Constant \a con into the corresponding Constant on all the functions of modes and transitions.
	void substitute(const Constant<Real>& con) { this->substitute(con,con.value()); }

	//@}

    //@{
    //! \name Data access and queries.

    //! \brief Returns the hybrid automaton's name.
    const std::string& name() const;

    //! \brief The sets of input, output, and internal variables.
    const std::set< RealVariable >& input_vars() const;
    const std::set< RealVariable >& output_vars() const;
    const std::set< RealVariable >& internal_vars() const;

    //! \brief The sets of controlled variables of the automaton (i.e., internal and output variables).
    std::set< RealVariable > controlled_vars() const;        

    //! \brief The sets of input, output, and internal events.
    const std::set< DiscreteEvent >& input_events() const;
    const std::set< DiscreteEvent >& output_events() const;
    const std::set< DiscreteEvent >& internal_events() const;

    //! \brief The sets of controlled events of the automaton (i.e., internal and output events).
    std::set< DiscreteEvent > controlled_events() const;

    //! \brief Test if the hybrid automaton has a input variable named \a u.
    bool has_input_var(const RealVariable& u) const;
    
    //! \brief Test if the hybrid automaton has a output variable named \a y.
    bool has_output_var(const RealVariable& y) const;
    
    //! \brief Test if the hybrid automaton has a internal variable named \a x.
    bool has_internal_var(const RealVariable& x) const;

    //! \brief Test if the hybrid automaton has a input event named \a e.
    bool has_input_event(const DiscreteEvent& e) const;
    
    //! \brief Test if the hybrid automaton has a output event named \a e.
    bool has_output_event(const DiscreteEvent& e) const;
    
    //! \brief Test if the hybrid automaton has a internal event named \a e.
    bool has_internal_event(const DiscreteEvent& e) const;
    
    //! \brief Test if the hybrid automaton has a discrete mode with discrete state \a state.
    bool has_mode(DiscreteState state) const;

    //! \brief Test if the hybrid automaton has a discrete transition with \a event_id and \a source_id.
    bool has_transition(DiscreteEvent event, DiscreteState source) const;

    //! \brief Test if the hybrid automaton has a constant with the same label.
    bool has_accessible_constant(const RealConstant& con) const;

    //! \brief Get the value of an accessible constant.
    const Real& accessible_constant_value(const String& name) const;

    //! \brief The discrete mode with given discrete state.
    const DiscreteIOMode& mode(DiscreteState state) const;

    //! \brief The discrete transition with given \a event and \a source location.
    const DiscreteIOTransition& transition(DiscreteEvent event, DiscreteState source) const;

    //! \brief The set of discrete modes. 
    const std::list< DiscreteIOMode >& modes() const;

    //! \brief The set of discrete transitions. 
    const std::list< DiscreteIOTransition >& transitions() const;

    //! \brief The set of accessible constants.
    const std::list< RealConstant >& accessible_constants() const;

    //! \brief The discrete transitions from location \a source.
    std::list< DiscreteIOTransition > transitions(DiscreteState source) const;

    //! \brief The natural grid to use over all variables.
    const std::map< RealVariable, Float >& grid() const;
    
    //! \brief The grid scaling for a given variable.
    Float grid(const RealVariable& var) const;

    //@}

};

//@{
//! \name Auxiliary functions for HybridIOAutomaton.
//!
//! \brief Convert an HybridIOAutomaton to a monolithic HybridAutomaton
//!     The input HybridIOAutomaton should be closed (no input variables or events) and
//!     with all dynamics, activations, and resets specified. Otherwise an exception is thrown.
//!
std::pair< HybridAutomaton, RealSpace > make_monolithic_automaton(const HybridIOAutomaton& hioa);

//! \brief Checks whether an HybridIOAutomaton is an elastic controller under the AASAP methodology.
bool is_elastic_controller(const HybridIOAutomaton& hioa);

//! \brief Convert an HybridIOAutomaton to a relaxed HybridIOAutomaton under the AASAP semantics. 
//! The input automaton must be an elastic controller.
HybridIOAutomaton aasap_relaxation(const HybridIOAutomaton& hioa);

//! \brief Compose two HybridIOAutomaton.
//!     The two I/O automata should be compatible. No conflict between variables and events.
HybridIOAutomaton compose(const std::string& name, 
                           const HybridIOAutomaton& ha1, 
                           const HybridIOAutomaton& ha2,
                           const DiscreteState& init1,
                           const DiscreteState& init2);

//! \brief Output an HybridIOAutomaton to a stream.
std::ostream& operator<<(std::ostream& os, const HybridIOAutomaton& ha);

//@}

} // namespace Ariadne

#endif // ARIADNE_HYBRID_IO_AUTOMATON_H
