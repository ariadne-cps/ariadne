/***************************************************************************
 *            hybrid_system.h
 *
 *  Copyright  2009  Pieter Collins
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

/*! \file hybrid_system.h
 *  \brief Main compositional hybrid system class.
 */

#ifndef ARIADNE_HYBRID_SYSTEM_H
#define ARIADNE_HYBRID_SYSTEM_H

#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

#include "logging.h"

#include "formula.h"
#include "container.h"
#include "valuation.h"

#include "discrete_state.h"

namespace Ariadne {


class Event;
class EventSet;

class HybridTime;
class HybridSpace;
class HybridSet;
class HybridGrid;

class DiscreteState;
class DiscreteSpace;

class ScalarFunction;
class VectorFunction;

template<class R> class ExpressionInterface;
class Grid;

template<class T> class List;
template<class T> class Set;
template<class k, class V> class Map;



/*! \brief A hybrid system, comprising continuous-time behaviour
 *  at each discrete mode, coupled by instantaneous discrete transitions.
 *  The state space is given by a hybrid set.
 * \sa \link Ariadne::HybridAutomaton \c HybridAutomaton \endlink.

 */
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

    typedef boost::shared_ptr<const ScalarFunction> ScalarFunctionPtr;
    typedef boost::shared_ptr<const VectorFunction> VectorFunctionPtr;


/*
    typedef std::map< Event, boost::shared_ptr<const VectorFunction> >::const_iterator invariant_const_iterator;
    typedef std::set< DiscreteTransition >::const_iterator discrete_transition_const_iterator;
    typedef std::set< DiscreteMode >::const_iterator discrete_mode_const_iterator;
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

    typedef List<DifferentialEquation>::const_iterator dynamic_const_iterator;
    typedef List<AlgebraicEquation>::const_iterator relation_const_iterator;
    typedef List<DiscreteUpdate>::const_iterator switch_const_iterator;
    typedef List<ContinuousUpdate>::const_iterator jump_const_iterator;
    typedef List<InvariantPredicate>::const_iterator invariant_const_iterator;
    typedef List<GuardPredicate>::const_iterator guard_const_iterator;
    typedef List<DisabledEvents>::const_iterator disabled_const_iterator;
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
    void new_equation(DiscretePredicate q, RealAssignment a) {
        AlgebraicEquation eqn={q,a.lhs,a.rhs}; _algebraic_equations.push_back(eqn); };
    //! \brief Adds a differential equation to the system.
    void new_dynamic(DiscretePredicate q, RealDynamic d) {
        DifferentialEquation eqn={q,d.lhs.base,d.rhs}; _differential_equations.push_back(eqn); };
    //! \brief Adds a discrete reset to the system.
    void new_transition(EventSet e, DiscretePredicate q, StringUpdate a) {
        DiscreteUpdate eqn={e,q,a.lhs.base,a.rhs}; _discrete_updates.push_back(eqn); }
    //! \brief Adds a reset equation to the system.
    void new_reset(EventSet e, DiscretePredicate q, RealUpdate a) {
        ContinuousUpdate eqn={e,q,a.lhs.base,a.rhs}; _continuous_updates.push_back(eqn); }
    //! \brief Adds a guard predicate to the system.
    void new_guard(EventSet e, DiscretePredicate q, ContinuousPredicate p) {
        GuardPredicate eqn={e,q,p}; _guard_predicates.push_back(eqn); }
    //! \brief Adds a guard predicate and invariant to the system.
    void new_guard(EventSet e, DiscretePredicate q, ContinuousPredicate a, ContinuousPredicate i) {
        GuardPredicate aeqn={e,q,a}; _guard_predicates.push_back(aeqn);
        InvariantPredicate ieqn={q,i}; _invariant_predicates.push_back(ieqn); }
    //! \brief Adds a guard predicate to the system.
    void new_guard(EventSet e, DiscretePredicate q, bool p) {
        GuardPredicate eqn={e,q,ContinuousPredicate(tribool(p))}; _guard_predicates.push_back(eqn); }
    //! \brief Adds an invariant to the system.
    void new_invariant(DiscretePredicate q, ContinuousPredicate p) {
        InvariantPredicate eqn={q,p}; _invariant_predicates.push_back(eqn); }
    //! \brief Disables events in a given set of locations.
    void new_disabled_events(EventSet e, DiscretePredicate q) {
    DisabledEvents dis={e,q}; _disabled_events.push_back(dis); }

    // Methods for rules valid in all modes.
    //! \brief Adds a algebraic equation to the system, valid in all modes.
    void new_equation(RealAssignment a) { this->new_equation(DiscretePredicate(true),a); }
    //! \brief Adds a differential equation to the system.
    void new_dynamic(RealDynamic d) { this->new_dynamic(DiscretePredicate(true),d); }
    //! \brief Adds a discrete reset to the system, valid in all modes.
    void new_transition(EventSet e, StringUpdate du) { this->new_transition(e,DiscretePredicate(true),du); }
    //! \brief Adds a reset equation to the system, valid in all modes.
    void new_reset(EventSet e, RealUpdate u) { this->new_reset(e,DiscretePredicate(true),u); }
    //! \brief Adds a guard predicate to the system, valid in all modes.
    void new_guard(EventSet e, ContinuousPredicate p) { this->new_guard(e,DiscretePredicate(true),p); }
    //! \brief Adds a guard predicate to the system, valid in all modes.
    void new_guard(EventSet e, bool p) { this->new_guard(e,DiscretePredicate(true),ContinuousPredicate(tribool(p))); }
    //! \brief Adds an invariant to the system, valid in all modes.
    void new_invariant(ContinuousPredicate p) { this->new_invariant(DiscretePredicate(true),p); }
    //! \brief Disables events in all locations.
    void new_disabled_events(EventSet e) {
    DisabledEvents dis={e,DiscretePredicate(true)}; _disabled_events.push_back(dis); }

    // Methods for rules valid for all events.
    //! \brief Adds a discrete reset to the system, valid in all modes and for all events.
    void new_transition(StringUpdate du) { this->new_transition(EventSet::all(),DiscretePredicate(true),du); }
    //! \brief Adds a reset equation to the system, valid in all modes and for all events.
    void new_reset(RealUpdate u) { this->new_reset(EventSet::all(),DiscretePredicate(true),u); }

    //@}

    //@{
    //! \name Data access and queries.

    //! \brief .
    EventSet events() const;
    //! \brief .
    VariableSet discrete_variables() const;
    //! \brief .
    VariableSet result_variables(const DiscreteValuation& state) const;
    //! \brief .
    VariableSet argument_variables(const DiscreteValuation& state) const;
    //! \brief .
    VariableSet continuous_variables(const DiscreteValuation& state) const;
    //! \brief .
    VariableSet state_variables(const DiscreteValuation& state) const;
    //! \brief .
    VariableSet algebraic_variables(const DiscreteValuation& state) const;
    //! \brief .
    VariableSet auxiliary_variables(const DiscreteValuation& state) const;
    //! \brief .
    VariableSet input_variables(const DiscreteValuation& state) const;
    //! \brief .
    VariableSet output_variables(const DiscreteValuation& state) const;

    //! \brief .
    bool check_dynamic(const DiscreteValuation& location) const;
    //! \brief .
    bool check_guards(const DiscreteValuation& location) const;
    //! \brief .
    bool check_reset(const Event& event, const DiscreteValuation& source, const DiscreteValuation& target) const;


    //! \brief .
    EventSet events(const DiscreteValuation& location) const;
    //! \brief .
    VectorFunction dynamic(const DiscreteValuation& location) const;
    //! \brief .
    DiscreteValuation target(const Event& event, const DiscreteValuation& source) const;
    //! \brief .
    VectorFunction reset(const Event& event, const DiscreteValuation& source) const;
    //! \brief .
    ScalarFunction guard(const Event& event, const DiscreteValuation& source) const;

    //! \brief .
    Set<RealAssignment> unordered_equations(const DiscreteValuation& state) const;

    //! \brief .
    List<RealAssignment> equations(const DiscreteValuation& state) const;
    //! \brief .
    List<RealDynamic> dynamics(const DiscreteValuation& state) const;
    //! \brief .
    List<StringUpdate> switching(const Event& event, const DiscreteValuation& state) const;
    //! \brief .
    List<RealUpdate> resets(const Event& event, const DiscreteValuation& state) const;
    //! \brief .
    Map<Event,ContinuousPredicate> guards(const DiscreteValuation& state) const;
    //! \brief .
    ContinuousPredicate guard_predicate(const Event& event, const DiscreteValuation& state) const;

    //@}


    //@{
    //! \name Operations on systems.

    //! \brief The parallel composition of two systems.
    friend HybridSystem parallel_composition(const HybridSystem&, const HybridSystem&);

    //@}
};

std::ostream& operator<<(std::ostream& os, const HybridSystem& hs);
std::ostream& operator<<(std::ostream& os, const HybridSystem::AlgebraicEquation& ae);
std::ostream& operator<<(std::ostream& os, const HybridSystem::DifferentialEquation& de);
std::ostream& operator<<(std::ostream& os, const HybridSystem::DiscreteUpdate& da);
std::ostream& operator<<(std::ostream& os, const HybridSystem::ContinuousUpdate& re);
std::ostream& operator<<(std::ostream& os, const HybridSystem::GuardPredicate& g);
std::ostream& operator<<(std::ostream& os, const HybridSystem::InvariantPredicate& inv);




} // namespace Ariadne

#endif // ARIADNE_HYBRID_SYSTEM_H
