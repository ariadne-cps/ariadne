/***************************************************************************
 *            hybrid_automaton-restrictive.hpp
 *
 *  Copyright  2004-20  Alberto Casagrande, Pieter Collins
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

/*! \file hybrid_automaton-restrictive.hpp
 *  \brief Class for hybrid systems defined by a list of constraints on the dynamics.
 */

#ifndef ARIADNE_RESTRICTIVE_HYBRID_AUTOMATON_HPP
#define ARIADNE_RESTRICTIVE_HYBRID_AUTOMATON_HPP

#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "../utility/tuple.hpp"
#include "../function/function.hpp"
#include "../hybrid/discrete_location.hpp"
#include "../hybrid/discrete_event.hpp"
#include "../symbolic/assignment.hpp"
#include "../symbolic/expression.hpp"
#include "../symbolic/valuation.hpp"
#include "../output/logging.hpp"

#include "../hybrid/hybrid_automaton_interface.hpp"

namespace Ariadne {

class Grid;

class HybridTime;
class HybridSpace;

class HybridRealExpressionBoundedConstraintSet;
typedef HybridRealExpressionBoundedConstraintSet HybridExpressionSet;

class DiscreteMode;
class DiscreteTransition;

typedef ContinuousPredicate RealPredicate;

EffectiveVectorMultivariateFunction dynamic_function(Space<Real>& space, const List<RealAssignment>& algebraic, const List<DottedRealAssignment>& differential);
EffectiveScalarMultivariateFunction constraint_function(Space<Real>& space, const List<RealAssignment>& algebraic, const RealPredicate& constraint);


template<class T> class FinitarySet
{
    Bool _is_infinite;
    Set<T> _set_or_complement;
  public:
    FinitarySet() : _is_infinite(false), _set_or_complement() { }
    FinitarySet(const T& t) : _is_infinite(false), _set_or_complement() { _set_or_complement.insert(t); }
    FinitarySet(const InitializerList<T>& s) : _is_infinite(false), _set_or_complement(s) { }
    FinitarySet(const Set<T>& s) : _is_infinite(false), _set_or_complement(s) { }
    FinitarySet(Bool c, const Set<T>& s) : _is_infinite(c), _set_or_complement(s) { }
    Bool is_finite() const { return !this->_is_infinite; }
    Bool is_infinite() const { return this->_is_infinite; }
    Bool contains(const T& t) const { return this->_set_or_complement.contains(t) xor _is_infinite; }
    const Set<T>& _underlying_set() const { return this->_set_or_complement; }
    Void adjoin(const FinitarySet<T>& set) { *this=join(*this,set); }
    Void restrict(const FinitarySet<T>& set) { *this=intersection(*this,set); }
};
template<class T> FinitarySet<T> complement(const FinitarySet<T>& s) {
    return FinitarySet<T>(!s.is_infinite(),s._underlying_set()); }
template<class T> FinitarySet<T> difference(const FinitarySet<T>& s1,const FinitarySet<T>& s2) {
    return intersection(s1,complement(s2)); }
template<class T> FinitarySet<T> join(const FinitarySet<T>& s1,const FinitarySet<T>& s2) {
    return complement(intersection(complement(s1),complement(s2))); }
template<class T> FinitarySet<T> intersection(const FinitarySet<T>& s1,const FinitarySet<T>& s2) {
    if(s1.is_finite()) {
        if(s2.is_finite()) { return FinitarySet<T>(false,intersection(s1._underlying_set(),s2._underlying_set())); }
        else { return FinitarySet<T>(false,difference(s1._underlying_set(),s2._underlying_set())); }
    } else {
        if(s2.is_finite()) { return FinitarySet<T>(false,difference(s2._underlying_set(),s1._underlying_set())); }
        else { return FinitarySet<T>(true,join(s1._underlying_set(),s2._underlying_set())); }
    }
}
template<class T> OutputStream& operator<<(OutputStream& os, const FinitarySet<T>& s) {
    if(s.is_infinite()) { os << "~"; }
    return os << s._underlying_set();
}

typedef FinitarySet<DiscreteEvent> EventSet;

class EventPredicate {
   public:
     Bool operator() (const DiscreteEvent& e);
};

typedef ContinuousPredicate RealPredicate;

typedef List<PrimedStringAssignment> DiscreteUpdate;
typedef List<PrimedRealAssignment> ContinuousUpdate;
DiscretePredicate predicate(const DiscreteLocation& location);
DiscreteUpdate update(const DiscreteLocation& location);


class RestrictiveDiscreteMode {
  public:
    DiscreteLocation location;
    List<RealAssignment> auxiliary;
    List<DottedRealAssignment> dynamic;

    Map<DiscreteEvent, RealPredicate> invariants;
    Map<DiscreteEvent, Tuple< RealPredicate,DiscreteLocation,List<PrimedRealAssignment> > > transitions;
};
OutputStream& operator<<(OutputStream& os, const RestrictiveDiscreteMode& mode);

Set<Identifier> names(const Set<RealVariable>& v);

//! \ingroup SystemModule
//! \ingroup HybridAutomataSubModule
//! \brief A hybrid system defined using purely restrictive operations.
//! \details Trivially compositional by definition.
class RestrictiveHybridAutomaton
//    : public virtual HybridAutomatonInterface
//    , public Loggable
: public Loggable
{
  public:
    //! \brief The type used to represent time.
    typedef HybridTime TimeType;
    //! \brief The type used to represent real numbers.
    typedef Real RealType ;
    //! \brief The type used to describe the state space.
    typedef HybridSpace StateSpaceType;

  private:
    List< Tuple< DiscretePredicate, List<RealAssignment> > > _assignments;
    List< Tuple< DiscretePredicate, List<DottedRealAssignment> > > _dotted_assignments;

    List< Tuple< DiscretePredicate, EventSet > > _disabled_events;
    List< Tuple< DiscretePredicate, EventSet, Set<Identifier> > > _nonjumping_discrete_variables;
    List< Tuple< DiscretePredicate, EventSet, Set<Identifier> > > _nonjumping_continuous_state_variables;

    List< Tuple< DiscretePredicate, DiscreteEvent, List<PrimedStringAssignment> > > _discrete_updates;
    List< Tuple< DiscretePredicate, DiscreteEvent, List<PrimedRealAssignment> > > _primed_assignments;
    List< Tuple< DiscretePredicate, DiscreteEvent, RealPredicate> > _invariant_predicates;
    List< Tuple< DiscretePredicate, DiscreteEvent, RealPredicate> > _guard_predicates;

    mutable Map< DiscreteLocation, RestrictiveDiscreteMode > _cached_modes;
  protected:

    // Construct a system only using the given discrete variables and events
    RestrictiveDiscreteMode compute_mode(const DiscreteLocation& q) const;

  public:
    //! \brief  Constructor.
    RestrictiveHybridAutomaton();

    //! \brief  Destructor.
    virtual ~RestrictiveHybridAutomaton() = default;
    //! \brief Construct dynamically-allocated copy.
    virtual RestrictiveHybridAutomaton* clone() const;

    friend RestrictiveHybridAutomaton compose(const List<RestrictiveHybridAutomaton>& components);

    //@{
    //! \name Methods for building the automaton.

    Void disable_events(DiscretePredicate q, EventSet e);
    Void nonjumping_variables(DiscretePredicate q, EventSet e, Set<StringVariable> v);
    Void nonjumping_variables(DiscretePredicate q, EventSet e, Set<RealVariable> v);

    Void new_update(DiscretePredicate q, DiscreteEvent e, List<PrimedStringAssignment> u);
    Void new_reset(DiscretePredicate q, DiscreteEvent e, List<PrimedRealAssignment> r);
    Void new_dynamic(DiscretePredicate q, List<DottedRealAssignment> d);
    Void new_dynamic(List<DottedRealAssignment> d);
    Void new_auxiliary(DiscretePredicate q, List<RealAssignment> a);
    Void new_auxiliary(List<RealAssignment> a);
    Void new_invariant(DiscretePredicate q, DiscreteEvent e, RealPredicate i);
    Void new_invariant(DiscretePredicate q, RealPredicate i);
    Void new_guard(DiscretePredicate q, DiscreteEvent e, RealPredicate g);

    Void new_transition(DiscretePredicate q, DiscreteEvent e, PrimedStringAssignment u, RealPredicate g);
    Void new_transition(DiscretePredicate q, DiscreteEvent e, PrimedStringAssignment u, List<PrimedRealAssignment> r, RealPredicate g);
    Void new_transition(DiscretePredicate q, DiscreteEvent e, List<PrimedStringAssignment> u, List<PrimedRealAssignment> r, RealPredicate g);

    Void new_mode(DiscreteLocation, List<RealAssignment>, List<DottedRealAssignment>);
    Void new_mode(DiscreteLocation, List<DottedRealAssignment>, List<RealAssignment>);
    Void new_transition(DiscreteLocation q, DiscreteEvent e, DiscreteLocation t, List<PrimedRealAssignment> r, RealPredicate g);
    Void new_transition(DiscreteLocation q, DiscreteEvent e, DiscreteLocation t, PrimedRealAssignment r, RealPredicate g);
    Void new_transition(DiscreteLocation q, DiscreteEvent e, DiscreteLocation t, RealPredicate g);

    //@}

    //@{
    //! \name Data access and queries.

    //! \brief The discrete mode with given discrete state.
    const RestrictiveDiscreteMode& mode(DiscreteLocation location) const;

    //! \brief The discrete events which are active in \a source.
    Set<DiscreteEvent> events(DiscreteLocation source) const;


    //! \brief Test if the hybrid automaton has an invariant (either explicit or from an urgent transition) with the given \a event label in \a location.
    Bool has_invariant(DiscreteLocation location, DiscreteEvent event) const;

    //! \brief Tests if the automaton has a guard predicate corresponding to the given location and event.
    Bool has_guard(DiscreteLocation, DiscreteEvent) const;


    //! \brief The discrete events corresponding to an invariant in \a source.
    Set<DiscreteEvent> invariant_events(DiscreteLocation source) const;

    //! \brief The discrete events corresponding to guards in \a source.
    Set<DiscreteEvent> guard_events(DiscreteLocation source) const;

    //! \brief The kind of response that needs to be taken
    //! when the guard corresponding to the \a event becomes true in the \a location.
    EventKind event_kind(DiscreteLocation location, DiscreteEvent event) const;



    //@}

    //@{
    //! \name Functions for conformance to HybridAutomatonInterface


    HybridSpace state_space() const;
    DimensionType dimension(DiscreteLocation) const;
    RealSpace continuous_state_space(DiscreteLocation) const;
    EffectiveVectorMultivariateFunction output_function(DiscreteLocation) const;
    EffectiveVectorMultivariateFunction auxiliary_function(DiscreteLocation) const;
    EffectiveVectorMultivariateFunction dynamic_function(DiscreteLocation) const;
    EffectiveVectorMultivariateFunction reset_function(DiscreteLocation, DiscreteEvent) const;
    EffectiveScalarMultivariateFunction constraint_function(DiscreteLocation, DiscreteEvent) const;
    EffectiveScalarMultivariateFunction invariant_function(DiscreteLocation, DiscreteEvent) const;
    EffectiveScalarMultivariateFunction guard_function(DiscreteLocation, DiscreteEvent) const;

    //@}


    //@{
    //! \name Methods for extracting the discrete dynamics of the automaton.

    //! \brief Test if the hybrid automaton has a valid discrete mode with the given \a location.
    Bool has_mode(DiscreteLocation location) const;

    //! \brief Test if the hybrid automaton has a discrete transition starting from the given location with the given event.
    Bool has_transition(DiscreteLocation source, DiscreteEvent event) const;

    //! \brief The discrete events corresponding to a discrete transition in \a source.
    Set<DiscreteEvent> transition_events(DiscreteLocation source) const;

    //! \brief The target for \a event from location \a source. Returns \a source if \a event is not present.
    DiscreteLocation target(DiscreteLocation location, DiscreteEvent event) const;

    //@}


    //@{
    //! \name Methods for extracting the continuous dynamics of the automaton.

    //! \brief The continuous variables which are required in the location.
    Set<UntypedVariable> argument_variables(DiscreteLocation) const;
    //! \brief The continuous variables which are state variables in the \a location.
    //! The state variables are those defined by a differential equation i.e. the dotted variables.
    List<RealVariable> state_variables(DiscreteLocation location) const;
    //! \brief The dependent continuous variables which are not state variables in the \a location.
    //! These variables are defined by algebraic equations.
    List<RealVariable> auxiliary_variables(DiscreteLocation location) const;
    //! \brief The continuous variables which are needed to compute the dynamics in the \a location
    //! but which are not state or auxiliary variables.
    Set<UntypedVariable> input_variables(DiscreteLocation location) const;

    //! \brief The algebraic equations valid in the location, ordered so that the defining equation for a variable
    //! occurs before any equation using that variable.
    List<RealAssignment> algebraic_assignments(DiscreteLocation location) const;
    //! \brief The differential equations valid in the location.
    List<DottedRealAssignment> differential_assignments(DiscreteLocation location) const;
    //! \brief The reset equations used when the \a event occurs in the \a source location.
    ContinuousUpdate update_assignments(DiscreteLocation source, DiscreteEvent event) const;
    //! \brief The invariant (time-can-progress predicates) corresponding to the given \a event.
    ContinuousPredicate invariant_predicate(DiscreteLocation location, DiscreteEvent event) const;
    //! \brief The guard (activation predicate) corresponding to the given \a event.
    ContinuousPredicate guard_predicate( DiscreteLocation location, DiscreteEvent event) const;
    //@}


/*
     //@{
    //! \name Methods for extracting the continuous dynamics of the automaton for conformance with HybridAutomatonInterface.

    //! \brief The dimension of the state spacec in the given \a location.
    virtual DimensionType dimension(DiscreteLocation location) const;
    //! \brief The output function on Euclidean state space. Used for outputting auxiliary variables.
    virtual EffectiveVectorMultivariateFunction output_function(DiscreteLocation location) const;
    //! \brief The function defining the differential equation \f$\der{x}=f(x)\f$ valid in the \a location.
    virtual EffectiveVectorMultivariateFunction dynamic_function(DiscreteLocation location) const;
    //! \brief The function defining the reset \f$x'=r(x)\f$ when the \a event occurs in the \a source location.
    virtual EffectiveVectorMultivariateFunction reset_function(DiscreteLocation source, DiscreteEvent event) const;
    //! \brief The function defining the guard \f$g(x) \geq 0\f$ for the given \a event to occur in \a location.
    virtual EffectiveScalarMultivariateFunction guard_function(DiscreteLocation location, DiscreteEvent event) const;
    //! \brief The function defining the invariant \f$p(x)\leq 0\f$ for continuous evolution to be blocked in \a location.
    virtual EffectiveScalarMultivariateFunction invariant_function(DiscreteLocation location, DiscreteEvent event) const;

    //! \brief The grid for the state variables in \a location.
    virtual Grid grid(DiscreteLocation location) const;

    virtual HybridGrid grid() const; // Deprecated

    virtual HybridSpace state_space() const;
    virtual RealSpace continuous_state_space(DiscreteLocation loc) const;

   //@}
*/
    //@{
    //! \name Methods for checking the validity of the automaton.

    //! \brief Checks whether the equations valid in the location are valid.
    //!
    //! Includes a check for algebraic dependencies, over-defined variables, under-defined variables, and
    //! variables which should be defined in a reset but are not.
    Void check_mode(DiscreteLocation) const;
    Void check_overspecification(DiscreteLocation) const;
    //! \brief Runs check_mode() in any mode reachable under the discrete dynamics from the given initial location(s).
    Void check_reachable_modes(const Set<DiscreteLocation>&) const;
    Void check_reachable_modes(DiscreteLocation) const;
    //@}

    //@{
    //! \name Discrete reachability analysis.

    //! \brief Performs a discrete reachability analysis from the given initial location.
    Set<DiscreteLocation> reachable_locations(const Set<DiscreteLocation>&) const;
    Set<DiscreteLocation> reachable_locations(const DiscreteLocation&) const;
    //@}

    //! \brief Write to an output stream.
    OutputStream& _write(OutputStream&) const;

  private:
    DiscreteLocation _compute_target(DiscreteLocation source, DiscreteEvent target) const;
    Map<DiscreteEvent,DiscreteLocation> _compute_targets(DiscreteLocation source) const;
    static List<RealAssignment> _order_algebraic_assignments(const List<RealAssignment>&);
};

inline OutputStream& operator<<(OutputStream& os, const RestrictiveHybridAutomaton& hs) {
    return hs._write(os);
}

RestrictiveHybridAutomaton compose(const List<RestrictiveHybridAutomaton>& components);

} // namespace Ariadne

#endif // ARIADNE_RESTRICTIVE_HYBRID_AUTOMATON_HPP
