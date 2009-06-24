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

#include "discrete_automaton.h"
#include "formula.h"

namespace Ariadne {


class DiscreteState;
class DiscreteEvent;

class HybridTime;
class HybridSpace;
class HybridSet;
class HybridGrid;

class DiscreteMode;
class DiscreteTransition;
class HybridAutomaton;

class ExpressionInterface;
class FunctionInterface;
class FormulaInterface;
class Grid;

/*! \brief A discrete mode of a hybrid automaton, comprising continuous evolution given by a vector field
 * within and invariant constraint set.
 *
 * A %DiscreteMode can only be created using the new_mode() method in
 * the %HybridAutomaton class.
 *
 * \sa \link Ariadne::HybridAutomaton \c HybridAutomaton \endlink, \link Ariadne::DiscreteTransition \c DiscreteTransition \endlink
 */
class DiscreteMode {
    friend class HybridAutomaton;
    typedef boost::shared_ptr<const FunctionInterface> FunctionPtr;
    typedef boost::shared_ptr<const FunctionInterface> ExpressionPtr;
  private:

    // The discrete mode's discrete state.
    DiscreteState _location;

    // The discrete mode's vector field.
    FunctionPtr _dynamic;
    // The discrete mode's invariants.
    std::map< DiscreteEvent, FunctionPtr > _invariants;

    // The discrete mode's grid for reachability analysis.
    boost::shared_ptr< const Grid > _grid;
  public:
    //! \brief The mode's discrete state.
    DiscreteState location() const {
        return this->_location; }

    //! \brief The discrete mode's dynamic (a vector field).
    const FunctionInterface& dynamic() const {
        return *this->_dynamic; }

    //! \brief The discrete mode's dynamic (a vector field).
    FunctionPtr dynamic_ptr() const {
        return this->_dynamic; }

    //! \brief The discrete mode's invariants.
    const std::map< DiscreteEvent, FunctionPtr >& invariants() const {
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
    DiscreteMode(DiscreteState location,
                 const FunctionInterface& dynamic);

    // Construct from objects managed by shared pointers (for internal use)
    DiscreteMode(DiscreteState location,
                 const FunctionPtr dynamic,
                 const std::vector< FunctionPtr >& invariants);

};


std::ostream& operator<<(std::ostream& os, const DiscreteMode& dm);

inline bool operator<(const DiscreteMode& mode1, const DiscreteMode& mode2) {
    return mode1.location() < mode2.location(); }




/*! \brief A discrete transition of a hybrid automaton, representing an instantaneous
 * jump from one discrete mode to another, governed by an activation set and a reset map.
 *
 * A %DiscreteTransition can only be created using the new_transition() method in
 * the %HybridAutomaton class.
 *
 * An invariant is modelled by a discrete transition with negative event id and null reset pointer.
 *
 * \sa \link Ariadne::HybridAutomaton \c HybridAutomaton \endlink, \link Ariadne::DiscreteMode \c DiscreteMode \endlink
 */
class DiscreteTransition
{
    friend class HybridAutomaton;
  private:
    // \brief The discrete transition's identificator.
    DiscreteEvent _event;

    // \brief The source of the discrete transition.
    const DiscreteMode* _source;

    // \brief The target of the discrete transition.
    const DiscreteMode* _target;

    // \brief The activation region of the discrete transition.
    boost::shared_ptr< const FunctionInterface > _activation;

    // \brief The reset of the discrete transition.
    boost::shared_ptr< const FunctionInterface > _reset;

    // \brief Whether or not the transition is forced.
    bool _forced;

  public:

    //! \brief The discrete event associated with the discrete transition.
    DiscreteEvent event() const {
        return this->_event; }

    //! \brief The source mode of the discrete transition.
    const DiscreteMode& source() const {
        return *this->_source; }

    //! \brief The target of the discrete transition.
    const DiscreteMode& target() const {
        return *this->_target; }


    //! \brief The activation region of the discrete transition.
    boost::shared_ptr<const FunctionInterface> activation_ptr() const {
        return this->_activation;
    }

    //! \brief The activation region of the discrete transition.
    const FunctionInterface& activation() const {
        return *this->_activation;
    }

    //! \brief The reset map of the discrete transition.
    const FunctionInterface& reset() const {
        return *this->_reset;
    }

    //! \brief The reset map of the discrete transition.
    boost::shared_ptr<const FunctionInterface> reset_ptr() const {
        return this->_reset;
    }

    //! \brief True if the transition is forced (occurs as soon as it is activated).
    bool forced() const {
        return this->_forced;
    }

  private:


    // Construct from shared pointers (for internal use).
    DiscreteTransition(DiscreteEvent event,
                       const DiscreteMode& source,
                       const DiscreteMode& target,
                       const FunctionInterface& reset,
                       const FunctionInterface& activation,
                       bool forced=false);

    // Construct from shared pointers (for internal use). */
    DiscreteTransition(DiscreteEvent event,
                       const DiscreteMode& source,
                       const DiscreteMode& target,
                       const boost::shared_ptr< FunctionInterface > reset,
                       const boost::shared_ptr< FunctionInterface > activation,
                       bool forced=false);
};

std::ostream& operator<<(std::ostream& os, const DiscreteTransition& dt);

inline bool operator<(const DiscreteTransition& transition1, const DiscreteTransition& transition2) {
    return transition1.event() < transition2.event()
        || (transition1.event() == transition2.event()
            && transition1.source().location() < transition2.source().location());
}




/*! \brief A hybrid system, comprising continuous-time behaviour
 *  at each discrete mode, coupled by instantaneous discrete transitions.
 *  The state space is given by a hybrid set.
 * \sa \link Ariadne::HybridAutomaton \c HybridAutomaton \endlink.

 */
class HybridSystem
{
  public:
    //! \brief The type used to represent time.
    typedef HybridTime TimeType;
    //! \brief The type used to represent real numbers.
    typedef double RealType ;
    //! \brief The type used to describe the state space.
    typedef HybridSpace StateSpaceType;

    typedef boost::shared_ptr<const ExpressionInterface> ExpressionPtr;
    typedef boost::shared_ptr<const FunctionInterface> FunctionPtr;
    typedef boost::shared_ptr<const FormulaInterface> FormulaPtr;


/*
    typedef std::map< DiscreteEvent, boost::shared_ptr<const FunctionInterface> >::const_iterator invariant_const_iterator;
    typedef std::set< DiscreteTransition >::const_iterator discrete_transition_const_iterator;
    typedef std::set< DiscreteMode >::const_iterator discrete_mode_const_iterator;
*/
  private:
  public:


    //! \brief The list of the hybrid automaton's discrete modes.
    //std::set< DiscreteMode > _modes;

    //struct DiscreteEquation { DiscretePredicate loc; DiscreteVariable lhs; DiscreteFormula rhs; };
    struct DifferentialEquation { DiscretePredicate loc; RealVariable lhs; RealFormula rhs; };
    struct AlgebraicEquation { DiscretePredicate loc; RealVariable lhs; RealFormula rhs; };
    struct DiscreteAssignment { DiscreteEventSet e; DiscretePredicate loc; DiscreteVariable lhs; DiscreteFormula rhs;  };
    struct UpdateEquation { DiscreteEventSet e; DiscretePredicate loc; RealVariable lhs; RealFormula rhs;  };
    struct GuardPredicate { DiscreteEventSet e; DiscretePredicate loc; RealPredicate pred; };
    struct InvariantPredicate { DiscretePredicate loc; RealPredicate pred; };

    std::vector<DifferentialEquation> _differential_equations;
    std::vector<AlgebraicEquation> _algebraic_equations;
    std::vector<DiscreteAssignment> _discrete_assignments;
    std::vector<UpdateEquation> _update_equations;
    std::vector<GuardPredicate> _guard_predicates;
    std::vector<InvariantPredicate> _invariant_predicates;

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
    void new_equation(DiscretePredicate q, Assignment a) {
        AlgebraicEquation eqn={q,a.lhs,a.rhs}; _algebraic_equations.push_back(eqn); };
    void new_dynamic(DiscretePredicate q, Dynamic d) {
        DifferentialEquation eqn={q,d.lhs,d.rhs}; _differential_equations.push_back(eqn); };
    void new_reset(DiscreteEventSet e, DiscretePredicate q, Ariadne::DiscreteUpdate a) {
        DiscreteAssignment eqn={e,q,a.lhs,a.rhs}; _discrete_assignments.push_back(eqn); }
    void new_reset(DiscreteEventSet e, DiscretePredicate q, Ariadne::Update a) {
        UpdateEquation eqn={e,q,a.lhs,a.rhs}; _update_equations.push_back(eqn); }
    void new_guard(DiscreteEventSet e, DiscretePredicate q, RealPredicate p) {
        GuardPredicate eqn={e,q,p}; _guard_predicates.push_back(eqn); }
    void new_invariant(DiscretePredicate q, RealPredicate p) {
        InvariantPredicate eqn={q,p}; _invariant_predicates.push_back(eqn); }

    // Methods for rules valid in all modes.
    void new_invariant(RealPredicate p) { this->new_invariant(DiscretePredicate(true),p); }
    void new_equation(Assignment a) { this->new_equation(DiscretePredicate(true),a); }
    void new_dynamic(Dynamic d) { this->new_dynamic(DiscretePredicate(true),d); }
    void new_reset(DiscreteEventSet e, Ariadne::DiscreteUpdate du) { this->new_reset(e,DiscretePredicate(true),du); }
    void new_reset(DiscreteEventSet e, Ariadne::Update u) { this->new_reset(e,DiscretePredicate(true),u); }
    void new_guard(DiscreteEventSet e, RealPredicate p) { this->new_guard(e,DiscretePredicate(true),p); }

    // Methods for rules valid for all events.
    void new_reset(Ariadne::DiscreteUpdate du) { this->new_reset(DiscreteEventSet::all(),DiscretePredicate(true),du); }
    void new_reset(Ariadne::Update u) { this->new_reset(DiscreteEventSet::all(),DiscretePredicate(true),u); }

    Space state_variables(const DiscreteState& state);
    Space algebraic_variables(const DiscreteState& state);
    Space input_variables(const DiscreteState& state);

    //! \brief Adds a discrete mode to the automaton.
    //!
    //!   \param state is the mode's discrete state.
    //!   \param dynamic is the mode's vector field.
    const DiscreteMode& new_mode(DiscreteState state,
                                 const FunctionInterface& dynamic);


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

    //! \brief The discrete mode with given discrete state.
    const DiscreteMode& mode(DiscreteState state) const;

    //! \brief The discrete transition with given \a event and \a source location.
    const DiscreteTransition& transition(DiscreteEvent event, DiscreteState source) const;

    //! \brief The set of discrete modes. (Not available in Python interface)
    const std::set< DiscreteMode >& modes() const;

    //! \brief The set of discrete transitions. (Not available in Python interface)
    const std::set< DiscreteTransition >& transitions() const;

    //! \brief The discrete transitions from location \a source.
    std::set< DiscreteTransition > transitions(DiscreteState source) const;

    //! \brief The blocking events (invariants and urgent transitions) in \a location.
    std::map<DiscreteEvent,FunctionPtr> blocking_guards(DiscreteState location) const;

    //! \brief The permissive events (invariants and urgent transitions) in \a location.
    std::map<DiscreteEvent,FunctionPtr> permissive_guards(DiscreteState location) const;

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
    //! \name Operations on systems.

    //! \brief The parallel composition of two systems.
    HybridSystem parallel_composition(const HybridSystem&, const HybridSystem&);

    //@}
};

std::ostream& operator<<(std::ostream& os, const HybridSystem& hs);
std::ostream& operator<<(std::ostream& os, const HybridSystem::AlgebraicEquation& ae);
std::ostream& operator<<(std::ostream& os, const HybridSystem::DifferentialEquation& de);
std::ostream& operator<<(std::ostream& os, const HybridSystem::DiscreteAssignment& da);
std::ostream& operator<<(std::ostream& os, const HybridSystem::UpdateEquation& re);
std::ostream& operator<<(std::ostream& os, const HybridSystem::GuardPredicate& g);
std::ostream& operator<<(std::ostream& os, const HybridSystem::InvariantPredicate& inv);




} // namespace Ariadne

#endif // ARIADNE_HYBRID_SYSTEM_H
