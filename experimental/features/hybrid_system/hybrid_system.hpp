/***************************************************************************
 *            hybrid_system.hpp
 *
 *  Copyright  2009  Pieter Collins
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

/*! \file hybrid_system.hpp
 *  \brief Main compositional hybrid system class.
 */

#ifndef ARIADNE_HYBRID_SYSTEM_HPP
#define ARIADNE_HYBRID_SYSTEM_HPP

#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <map>

#include <memory>

#include "../utility/declarations.hpp"

#include "../output/logging.hpp"

#include "../utility/container.hpp"
#include "../symbolic/valuation.hpp"
#include "../symbolic/assignment.hpp"
#include "../symbolic/expression.hpp"

#include "../hybrid/discrete_location.hpp"

namespace Ariadne {


class Event;
class EventSet;

class HybridTime;
class HybridSpace;
class HybridGrid;

class HybridRealExpressionBoundedConstraintSet;
typedef HybridRealExpressionBoundedConstraintSet HybridExpressionSet;

class DiscreteLocation;
class DiscreteSpace;

template<class R> class ExpressionInterface;
class Grid;

template<class T> class List;
template<class T> class Set;
template<class k, class V> class Map;



//! \brief A hybrid system, comprising continuous-time behaviour
//! at each discrete mode, coupled by instantaneous discrete transitions.
//! The state space is given by a hybrid set.
//! UNDER DEVELOPMENT
//! \sa \ref HybridAutomaton.
class HybridSystem
    : public Loggable
{
  public:
    //! \brief The type used to represent time.
    typedef HybridTime TimeType;
    //! \brief The type used to represent real numbers.
    typedef double RealType ;
    //! \brief The type used to describe the state space.
    typedef HybridSpace StateSpaceType;

    typedef std::shared_ptr<const RealScalarFunction> ScalarFunctionPtr;
    typedef std::shared_ptr<const RealVectorFunction> VectorFunctionPtr;


/*
    typedef std::map< Event, std::shared_ptr<const VectorFunction> >::ConstIterator invariant_const_iterator;
    typedef std::set< DiscreteTransition >::ConstIterator discrete_transition_const_iterator;
    typedef std::set< DiscreteMode >::ConstIterator discrete_mode_const_iterator;
*/
  private:
  public:


    //! \brief The list of the hybrid automaton's discrete modes.
    //std::set< DiscreteMode > _modes;

    //struct DiscreteEquation { DiscretePredicate loc; DiscreteVariable lhs; DiscreteExpression rhs; };
    struct DifferentialEquation { DiscretePredicate loc; RealVariable lhs; RealExpression rhs; };
    struct AlgebraicEquation { DiscretePredicate loc; RealVariable lhs; RealExpression rhs; };
    struct InvariantPredicate { DiscretePredicate loc; ContinuousPredicate pred; };
    struct DiscreteUpdate { EventSet evnts; DiscretePredicate loc; StringVariable lhs; StringExpression rhs;  };
    struct ContinuousUpdate { EventSet evnts; DiscretePredicate loc; RealVariable lhs; RealExpression rhs;  };
    struct GuardPredicate { EventSet evnts; DiscretePredicate loc; ContinuousPredicate pred; };
    struct DisabledEvents { EventSet evnts; DiscretePredicate loc; };

    List<DifferentialEquation> _differential_equations;
    List<AlgebraicEquation> _algebraic_equations;
    List<DiscreteUpdate> _discrete_updates;
    List<ContinuousUpdate> _continuous_updates;
    List<GuardPredicate> _guard_predicates;
    List<InvariantPredicate> _invariant_predicates;
    List<DisabledEvents> _disabled_events;

    typedef List<DifferentialEquation>::ConstIterator dynamic_const_iterator;
    typedef List<AlgebraicEquation>::ConstIterator relation_const_iterator;
    typedef List<DiscreteUpdate>::ConstIterator switch_const_iterator;
    typedef List<ContinuousUpdate>::ConstIterator jump_const_iterator;
    typedef List<InvariantPredicate>::ConstIterator invariant_const_iterator;
    typedef List<GuardPredicate>::ConstIterator guard_const_iterator;
    typedef List<DisabledEvents>::ConstIterator disabled_const_iterator;
  public:
    //@{
    //! \name Constructors and destructors

    //! \brief Construct an empty automaton with no name
    HybridSystem();

    //! \brief Construct dynamically-allocated copy. (Not currently implemented)
    HybridSystem* clone() const;

    //! \brief  Destructor.
    ~HybridSystem();
    //@}

    //@{
    //! \name Methods for building the automaton.

    // Methods for rules valid in certain modes
    //! \brief Adds a algebraic equation to the system.
    Void new_equation(DiscretePredicate q, RealAssignment a) {
        AlgebraicEquation eqn={q,a.lhs,a.rhs}; _algebraic_equations.push_back(eqn); };
    //! \brief Adds a differential equation to the system.
    Void new_dynamic(DiscretePredicate q, RealDynamic d) {
        DifferentialEquation eqn={q,d.lhs.base(),d.rhs}; _differential_equations.push_back(eqn); };
    //! \brief Adds a discrete reset to the system.
    Void new_transition(EventSet e, DiscretePredicate q, StringUpdate a) {
        DiscreteUpdate eqn={e,q,a.lhs.base(),a.rhs}; _discrete_updates.push_back(eqn); }
    Void new_transition(DiscretePredicate q, EventSet e, StringUpdate a) {
        DiscreteUpdate eqn={e,q,a.lhs.base(),a.rhs}; _discrete_updates.push_back(eqn); }
    //! \brief Adds a reset equation to the system.
    Void new_reset(EventSet e, DiscretePredicate q, RealUpdate a) {
        ContinuousUpdate eqn={e,q,a.lhs.base(),a.rhs}; _continuous_updates.push_back(eqn); }
    Void new_reset(DiscretePredicate q, EventSet e, RealUpdate a) {
        ContinuousUpdate eqn={e,q,a.lhs.base(),a.rhs}; _continuous_updates.push_back(eqn); }
    //! \brief Adds a guard predicate to the system.
    Void new_guard(EventSet e, DiscretePredicate q, ContinuousPredicate p) {
        GuardPredicate eqn={e,q,p}; _guard_predicates.push_back(eqn); }
    Void new_guard(DiscretePredicate q, EventSet e, ContinuousPredicate p) {
        GuardPredicate eqn={e,q,p}; _guard_predicates.push_back(eqn); }
    //! \brief Adds a guard predicate and invariant to the system.
    Void new_guard(EventSet e, DiscretePredicate q,  ContinuousPredicate a, ContinuousPredicate i) {
        GuardPredicate aeqn={e,q,a}; _guard_predicates.push_back(aeqn);
        InvariantPredicate ieqn={q,i}; _invariant_predicates.push_back(ieqn); }
    Void new_guard(DiscretePredicate q, EventSet e, ContinuousPredicate a, ContinuousPredicate i) {
        GuardPredicate aeqn={e,q,a}; _guard_predicates.push_back(aeqn);
        InvariantPredicate ieqn={q,i}; _invariant_predicates.push_back(ieqn); }
    //! \brief Adds a guard predicate to the system.
    Void new_guard(EventSet e, DiscretePredicate q, Bool p) {
        GuardPredicate eqn={e,q,ContinuousPredicate(ValidatedKleenean(p))}; _guard_predicates.push_back(eqn); }
    Void new_guard(DiscretePredicate q, EventSet e, Bool p) {
        GuardPredicate eqn={e,q,ContinuousPredicate(ValidatedKleenean(p))}; _guard_predicates.push_back(eqn); }
    //! \brief Adds an invariant to the system.
    Void new_invariant(DiscretePredicate q, ContinuousPredicate p) {
        InvariantPredicate eqn={q,p}; _invariant_predicates.push_back(eqn); }
    //! \brief Disables events in a given set of locations.
    Void new_disabled_events(EventSet e, DiscretePredicate q) {
        DisabledEvents dis={e,q}; _disabled_events.push_back(dis); }
    Void new_disabled_events(DiscretePredicate q, EventSet e) {
        DisabledEvents dis={e,q}; _disabled_events.push_back(dis); }

    // Methods for rules valid in all modes.
    //! \brief Adds a algebraic equation to the system, valid in all modes.
    Void new_equation(RealAssignment a) { this->new_equation(DiscretePredicate(true),a); }
    //! \brief Adds a differential equation to the system.
    Void new_dynamic(RealDynamic d) { this->new_dynamic(DiscretePredicate(true),d); }
    //! \brief Adds a discrete reset to the system, valid in all modes.
    Void new_transition(EventSet e, StringUpdate du) { this->new_transition(e,DiscretePredicate(true),du); }
    //! \brief Adds a reset equation to the system, valid in all modes.
    Void new_reset(EventSet e, RealUpdate u) { this->new_reset(e,DiscretePredicate(true),u); }
    //! \brief Adds a guard predicate to the system, valid in all modes.
    Void new_guard(EventSet e, ContinuousPredicate p) { this->new_guard(e,DiscretePredicate(true),p); }
    //! \brief Adds a guard predicate to the system, valid in all modes.
    Void new_guard(EventSet e, Bool p) { this->new_guard(e,DiscretePredicate(true),ContinuousPredicate(ValidatedKleenean(p))); }
    //! \brief Adds an invariant to the system, valid in all modes.
    Void new_invariant(ContinuousPredicate p) { this->new_invariant(DiscretePredicate(true),p); }
    //! \brief Disables events in all locations.
    Void new_disabled_events(EventSet e) {
    DisabledEvents dis={e,DiscretePredicate(true)}; _disabled_events.push_back(dis); }

    // Methods for rules valid for all events.
    //! \brief Adds a discrete reset to the system, valid in all modes and for all events.
    Void new_transition(StringUpdate du) { this->new_transition(EventSet::all(),DiscretePredicate(true),du); }
    //! \brief Adds a reset equation to the system, valid in all modes and for all events.
    Void new_reset(RealUpdate u) { this->new_reset(EventSet::all(),DiscretePredicate(true),u); }

    //@}

    //@{
    //! \name Data access and queries.

    //! \brief The set of all discrete events.
    EventSet events() const;
    //! \brief The set of all discrete (string and integer) variables.
    VariableSet discrete_variables() const;
    //! \brief The set of all variables which are specified in the given discrete \a state.
    VariableSet result_variables(const DiscreteValuation& state) const;
    //! \brief The set of all variables which are used in the given discrete \a state, but have no defining equation.
    VariableSet argument_variables(const DiscreteValuation& state) const;
    //! \brief The set of all continuous variables which occur as result or argument variables in the discrete \a state.
    VariableSet continuous_variables(const DiscreteValuation& state) const;
    //! \brief The set of all <em>state</em> variables in the discrete \a state,
    //!   which are those continuous variables defined by.a differential equation.
    VariableSet state_variables(const DiscreteValuation& location) const;
    //! \brief The set of all continuous variables defined by algebraic equations.
    VariableSet algebraic_variables(const DiscreteValuation& state) const;
    //! \brief The set of all continuous variables defined by algebraic equations which are not output variables.
    VariableSet auxiliary_variables(const DiscreteValuation& state) const;
    //! \brief The set of all external continuous variables which are
    //!   not defined by an algebraic or differential equation.
    VariableSet input_variables(const DiscreteValuation& state) const;
    //! \brief The set of external continuous variables with a defining equatioin.
    VariableSet output_variables(const DiscreteValuation& state) const;

    //! \brief Check that the continuous dynamics in \a location is valid, which
    //!   means that every continuous variable has exactly one defining equation.
    Bool check_dynamic(const DiscreteValuation& location) const;
    //! \brief Check that the guards in \a location are well-defined, which means
    //!   that they only use continuous variables defined in the given \a location.
    Bool check_guards(const DiscreteValuation& location) const;
    //! \brief Check that the reset map for the transition
    //!   from the \a source location to the \a target location when the \a event occurs
    //!   only uses variables defined in the source location,
    //!   and specifies exactly the state variables in the target location.
    Bool check_reset(const Event& event, const DiscreteValuation& source, const DiscreteValuation& target) const;


    //! \brief The set of all events which can occur in the discrete \a location.
    EventSet events(const DiscreteValuation& location) const;
    //! \brief The vector field describing the continuous dynamic in the discrete \a location, defined as
    //!   a function on the space of state variables.
    VectorFunction dynamic(const DiscreteValuation& location) const;
    //! \brief The target location for the transition when \a event occurs in the \a source location.
    DiscreteValuation target(const Event& event, const DiscreteValuation& source) const;
    //! \brief The reset map for the transition when \a event occurs in the \a source location,
    //!  defined as a function between the spaces of state variables.
    VectorFunction reset(const Event& event, const DiscreteValuation& source) const;
    //! \brief The guard \f$g(x)\geq0\f$ for the for the transition when \a event occurs in the \a source location,
    //!  defined as a function on the space of state variables.
    ScalarFunction guard(const Event& event, const DiscreteValuation& source) const;

    //! \brief The algebraic equations valid in the given discrete \a state in arbitrary order.
    Set<RealAssignment> unordered_equations(const DiscreteValuation& state) const;

    //! \brief The algebraic equations valid in the given discrete \a state, ordered so that
    //!   variables defined in an equation are only used in subsequent equations
    List<RealAssignment> equations(const DiscreteValuation& state) const;
    //! \brief The differential equations valid in the given discrete \a state.
    List<RealDynamic> dynamics(const DiscreteValuation& state) const;
    //! \brief The discrete part of the transition when \a event occurs in the given \a source location.
    List<StringUpdate> switching(const Event& event, const DiscreteValuation& source) const;
    //! \brief The continuous part of the transition when \a event occurs in the given \a source location.
    List<RealUpdate> resets(const Event& event, const DiscreteValuation& state) const;
    //! \brief The guard conditions valid in the given discrete \a state.
    Map<Event,ContinuousPredicate> guards(const DiscreteValuation& state) const;
    //! \brief The guard condition for \a event to occur in the discrete \a state.
    ContinuousPredicate guard_predicate(const Event& event, const DiscreteValuation& state) const;

    //@}


    //@{
    //! \name Operations on systems.

    //! \brief The parallel composition of two systems.
    friend HybridSystem parallel_composition(const HybridSystem&, const HybridSystem&);

    //@}
};

OutputStream& operator<<(OutputStream& os, const HybridSystem& hs);
OutputStream& operator<<(OutputStream& os, const HybridSystem::AlgebraicEquation& ae);
OutputStream& operator<<(OutputStream& os, const HybridSystem::DifferentialEquation& de);
OutputStream& operator<<(OutputStream& os, const HybridSystem::DiscreteUpdate& da);
OutputStream& operator<<(OutputStream& os, const HybridSystem::ContinuousUpdate& re);
OutputStream& operator<<(OutputStream& os, const HybridSystem::GuardPredicate& g);
OutputStream& operator<<(OutputStream& os, const HybridSystem::InvariantPredicate& inv);




} // namespace Ariadne

#endif // ARIADNE_HYBRID_SYSTEM_HPP
