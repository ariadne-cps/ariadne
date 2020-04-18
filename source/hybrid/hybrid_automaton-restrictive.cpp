/***************************************************************************
 *            hybrid_automaton-restrictive.cpp
 *
 *  Copyright  2010-20  Alberto Casagrande, Pieter Collins
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

#include "../function/functional.hpp"
#include "../config.hpp"

#include <map>

#include "../utility/macros.hpp"
#include "../utility/stlio.hpp"
#include "../utility/tuple.hpp"
#include "../symbolic/expression.hpp"
#include "../symbolic/assignment.hpp"
#include "../symbolic/space.hpp"
#include "../function/function.hpp"
#include "../hybrid/hybrid_time.hpp"
#include "../hybrid/hybrid_space.hpp"

#include "../hybrid/hybrid_automaton-restrictive.hpp"

namespace Ariadne {


Bool is_blocking(EventKind);
Bool is_activating(EventKind);

namespace {

// Shorthand for testing whether a location satisfies a predicate. For internal use only.
Bool operator==(const DiscretePredicate& predicate, const DiscreteLocation& location) {
    return evaluate(predicate,location);
}

Bool operator==(const EventSet& set, const DiscreteEvent& event) {
    return set.contains(event);
}

template<class T>
List<T> filter(const List<Tuple<DiscretePredicate,List<T> > >& rules, const DiscreteLocation& q) {
    List<T> result;
    for(SizeType i=0; i!=rules.size(); ++i) {
        if(get_first(rules[i])==q) { result.append(get_second(rules[i])); }
    }
    return result;
}

template<class T1, class T2>
List< Tuple<T1,T2> > filter(const List<Tuple<DiscretePredicate,T1,T2> >& rules, const DiscreteLocation& q) {
    List< Tuple<T1,T2> > result;
    for(SizeType i=0; i!=rules.size(); ++i) {
        if(get_first(rules[i])==q) { result.append(make_tuple(get_second(rules[i]),get_third(rules[i]))); }
    }
    return result;
}

template<class T>
List<T> filter(const List<Tuple<DiscretePredicate,EventSet,List<T> > >& rules, const DiscreteLocation& q, const DiscreteEvent& e) {
    List<T> result;
    for(SizeType i=0; i!=rules.size(); ++i) {
        if(get_first(rules[i])==q && get_second(rules[i]).contains(e)) { result.append(get_third(rules[i])); }
    }
    return result;
}

template<class T>
T select(const List<Tuple<DiscretePredicate,EventSet,T> >& rules, const DiscreteLocation& q, const DiscreteEvent& e) {
    ARIADNE_ASSERT(rules.size()>0);
    Bool found=false;
    T result=get_third(rules[0]);
    for(SizeType i=0; i!=rules.size(); ++i) {
        if(get_first(rules[i])==q && get_second(rules[i]).contains(e)) {
            if(found) { ARIADNE_FAIL_MSG("More than one rule matching ("<<q<<","<<e<<") in "<<rules); }
            found=true; result=get_third(rules[i]);
        }
    }
    if(found==false) { ARIADNE_FAIL_MSG("No rule matching ("<<q<<","<<e<<") in "<<rules); }
    return result;
}

template<class T>
T select(const List<Tuple<DiscretePredicate,DiscreteEvent,T> >& rules, const DiscreteLocation& q, const DiscreteEvent& e) {
    ARIADNE_ASSERT(rules.size()>0);
    Bool found=false;
    T result=get_third(rules[0]);
    for(SizeType i=0; i!=rules.size(); ++i) {
        if(get_first(rules[i])==q && get_second(rules[i])==e) {
            if(found) { ARIADNE_FAIL_MSG("More than one rule matching ("<<q<<","<<e<<") in "<<rules); }
            found=true; result=get_third(rules[i]);
        }
    }
    if(found==false) { ARIADNE_FAIL_MSG("No rule matching ("<<q<<","<<e<<") in "<<rules); }
    return result;
}

EventSet filter_join(const List< Tuple<DiscretePredicate,EventSet> >& events, const DiscreteLocation& q) {
    EventSet result;
    for(SizeType i=0; i!=events.size(); ++i) {
        if(get_first(events[i])==q) { result=join(result,get_second(events[i])); }
    }
    return result;
}

template<class Assignment> List<RealVariable> left_hand_sides(const List<Assignment>& assignments) {
    List<RealVariable> variables;
    for(SizeType i=0; i!=assignments.size(); ++i) {
        variables.append(assignments[i].variable().base());
    }
    return variables;
}

template<class Assignment> List<RealExpression> right_hand_sides(const List<Assignment>& assignments) {
    List<RealExpression> expressions;
    for(SizeType i=0; i!=assignments.size(); ++i) {
        expressions.append(assignments[i].expression());
    }
    return expressions;
}

template<class X>
List<Identifier> names(const List<Variable<X> >& variables) {
    List<Identifier> result;
    for(SizeType i=0; i!=variables.size(); ++i) {
        result.append(variables[i].name());
    }
    return result;
}

template<class Tp, template<class> class Dec>
List< Variable<Tp> > base(const List< Dec<Tp> >& variables) {
    List< Variable<Tp> > result;
    for(SizeType i=0; i!=variables.size(); ++i) {
        result.append(variables[i].base());
    }
    return result;
}

template<class T>
Set<T> make_set(const List<T>& list) {
    return Set<T>(list.begin(),list.end());
}

template<class Var>
Bool unique(const List<Var>& variables) {
    Set<Var> found;
    for(SizeType i=0; i!=variables.size(); ++i) {
        if(found.contains(variables[i])) {
            return false;
        } else {
            found.insert(variables[i]);
        }
    }
    return true;
}

template<class Var>
Set<Var> duplicates(const List<Var>& variables) {
    Set<Var> result;
    Set<Var> found;
    for(SizeType i=0; i!=variables.size(); ++i) {
        if(found.contains(variables[i])) {
            result.insert(variables[i]);
        } else {
            found.insert(variables[i]);
        }
    }
    return result;
}


} // namespace


// Construct the set of events such that for each variable v in vars , there is a rule saying that v is nonjumping
EventSet ignorable(const List< Tuple<EventSet,Set<Identifier> > >& rules, const Set<Identifier>& vars);
// Construct the set of events such that for each variable v in vars , there is a rule saying that v is nonjumping
EventSet ignorable(const List< Tuple<DiscretePredicate, EventSet, Set<Identifier> > >& nonjumping, const DiscreteLocation& q, const Set<Identifier>& vars);

DiscretePredicate make_predicate(const DiscreteLocation& q);
List<PrimedStringAssignment> make_update(const DiscreteLocation& q);

List<PrimedRealAssignment> primed_real_assignments(const Set<Identifier>& nonjumping);

Set<Identifier> names(const Set<RealVariable>& v) {
    Set<Identifier> r;
    for(Set<RealVariable>::ConstIterator viter=v.begin(); viter!=v.end(); ++viter) {
        r.insert(viter->name());
    }
    return r;
}

// Construct the set of events such that for each variable v in vars , there is a rule saying that v is nonjumping
EventSet ignorable(const List< Tuple<EventSet,Set<Identifier> > >& rules, const Set<Identifier>& vars) {
    EventSet result;
    for(Set<Identifier>::ConstIterator var_iter=vars.begin(); var_iter!=vars.end(); ++var_iter) {
        EventSet variable_does_not_jump;
        for(SizeType i=0; i!=rules.size(); ++i) {
            if(get_second(rules[i]).contains(*var_iter)) {
                variable_does_not_jump.adjoin(get_first(rules[i]));
            }
        }
        result.restrict(variable_does_not_jump);
    }

    return result;
}

// Construct the set of events such that for each variable v in vars , there is a rule saying that v is nonjumping
EventSet ignorable(const List< Tuple<DiscretePredicate, EventSet, Set<Identifier> > >& nonjumping, const DiscreteLocation& q, const Set<Identifier>& vars) {
    return ignorable(filter(nonjumping,q),vars);
}


List<PrimedRealAssignment> primed_real_assignments(const Set<Identifier>& nonjumping) {
    List<PrimedRealAssignment> result;
    for(Set<Identifier>::ConstIterator var_iter=nonjumping.begin(); var_iter!=nonjumping.end(); ++var_iter) {
        RealVariable var(*var_iter);
        result.append(next(var)=var);
    }
    return result;
}



RestrictiveHybridAutomaton::RestrictiveHybridAutomaton()
{
}

Void RestrictiveHybridAutomaton::disable_events(DiscretePredicate q, EventSet e) {
    this->_disabled_events.append(make_tuple(q,e)); }
Void RestrictiveHybridAutomaton::nonjumping_variables(DiscretePredicate q, EventSet e, Set<RealVariable> v) {
    this->_nonjumping_continuous_state_variables.append(make_tuple(q,e,names(v))); }
Void RestrictiveHybridAutomaton::new_update(DiscretePredicate q, DiscreteEvent e, List<PrimedStringAssignment> u) {
    this->_discrete_updates.append(make_tuple(q,e,u)); }
Void RestrictiveHybridAutomaton::new_reset(DiscretePredicate q, DiscreteEvent e, List<PrimedRealAssignment> r) {
    this->_primed_assignments.append(make_tuple(q,e,r)); }
Void RestrictiveHybridAutomaton::new_dynamic(DiscretePredicate q, List<DottedRealAssignment> d) {
    this->_dotted_assignments.append(make_tuple(q,d)); }
Void RestrictiveHybridAutomaton::new_dynamic(List<DottedRealAssignment> d) {
    this->_dotted_assignments.append(make_tuple(DiscretePredicate(true),d)); }
Void RestrictiveHybridAutomaton::new_auxiliary(DiscretePredicate q, List<RealAssignment> a) {
    this->_assignments.append(make_tuple(q,a)); }
Void RestrictiveHybridAutomaton::new_auxiliary(List<RealAssignment> a) {
    this->_assignments.append(make_tuple(DiscretePredicate(true),a)); }
Void RestrictiveHybridAutomaton::new_invariant(DiscretePredicate q, DiscreteEvent e, RealPredicate i) {
    this->_invariant_predicates.append(make_tuple(q,e,i)); }
Void RestrictiveHybridAutomaton::new_invariant(DiscretePredicate q, RealPredicate i) {
    this->_invariant_predicates.append(make_tuple(q,DiscreteEvent(),i)); }
Void RestrictiveHybridAutomaton::new_guard(DiscretePredicate q, DiscreteEvent e, RealPredicate g) {
    this->_guard_predicates.append(make_tuple(q,e,g)); }
/*
Void RestrictiveHybridAutomaton::new_transition(DiscretePredicate q, DiscreteEvent e, List<PrimedStringAssignment> u,
                                                List<PrimedRealAssignment> r, RealPredicate g) {
    this->new_update(q,e,u); this->new_reset(q,e,r); this->new_guard(q,e,g); }
Void RestrictiveHybridAutomaton::new_transition(DiscretePredicate q, DiscreteEvent e, PrimedStringAssignment u, PrimedRealAssignment r, RealPredicate g) {
    this->new_reset(q,e,make_list(r)); this->new_guard(q,e,g); }
Void RestrictiveHybridAutomaton::new_transition(DiscretePredicate q, DiscreteEvent e, PrimedStringAssignment u, RealPredicate g) {
    this->new_update(q,e,make_list(u)); this->new_guard(q,e,g); }
*/

RestrictiveDiscreteMode RestrictiveHybridAutomaton::compute_mode(const DiscreteLocation& location) const
{
    RestrictiveDiscreteMode mode;

    mode.location=location;

    // Compute the continuous dynamics
    mode.auxiliary = filter(this->_assignments,location);
    mode.dynamic = filter(this->_dotted_assignments,location);

    Set<Identifier> discrete_variables = location.defined();
    Set<Identifier> state_variables = make_set(names(left_hand_sides(mode.dynamic)));

    // Extract the disabled events
    FinitarySet<DiscreteEvent> disabled_events = filter_join(this->_disabled_events,location);

    // If an event does not affect any discrete or continuous state variables, then it can be safely ignored
    FinitarySet<DiscreteEvent> discrete_ignorable_events = ignorable(this->_nonjumping_discrete_variables,location,discrete_variables);
    FinitarySet<DiscreteEvent> continuous_ignorable_events = ignorable(this->_nonjumping_continuous_state_variables,location,state_variables);
    FinitarySet<DiscreteEvent> ignorable_events = intersection(discrete_ignorable_events,continuous_ignorable_events);

    FinitarySet<DiscreteEvent> nonjumping_events = join(disabled_events,ignorable_events);
    ARIADNE_ASSERT(nonjumping_events.is_infinite());

    // The set of events leading to discrete transitions
    Set<DiscreteEvent> transition_events = nonjumping_events._underlying_set();

    // Compute the discrete transitions
    for(Set<DiscreteEvent>::ConstIterator event_iter=transition_events.begin(); event_iter!=transition_events.end(); ++event_iter) {
        DiscreteEvent event=*event_iter;

        // Compute the targets
        DiscreteLocation target;
        for(SizeType i=0; i!=this->_nonjumping_discrete_variables.size(); ++i) {
            if(get_first(this->_nonjumping_discrete_variables[i])==location && get_second(this->_nonjumping_discrete_variables[i])==event) {
                Set<Identifier> const& nonjumping_variables=get_third(this->_nonjumping_discrete_variables[i]);
                for(Set<Identifier>::ConstIterator var_iter=nonjumping_variables.begin(); var_iter!=nonjumping_variables.end(); ++var_iter) {
                    target[StringVariable(*var_iter)]=location[StringVariable(*var_iter)];
                }
            }
        }
        for(SizeType i=0; i!=this->_discrete_updates.size(); ++i) {
            if(get_first(this->_discrete_updates[i])==location && get_second(this->_discrete_updates[i])==event) {
                List<PrimedStringAssignment> const& updates=get_third(this->_discrete_updates[i]);
                for(List<PrimedStringAssignment>::ConstIterator asn_iter=updates.begin(); asn_iter!=updates.end(); ++asn_iter) {
                    const StringVariable lhs_var=asn_iter->lhs.base();
                    target[asn_iter->lhs.base()]=evaluate(asn_iter->rhs,location);
                }
            }
        }

        // Compute the resets
        List<PrimedRealAssignment> reset;
        for(SizeType i=0; i!=this->_nonjumping_continuous_state_variables.size(); ++i) {
            if(get_first(this->_nonjumping_continuous_state_variables[i])==location && get_second(this->_nonjumping_continuous_state_variables[i])==event) {
                reset.append(primed_real_assignments(get_third(this->_nonjumping_continuous_state_variables[i])));
            }
        }
        for(SizeType i=0; i!=this->_primed_assignments.size(); ++i) {
            if(get_first(this->_primed_assignments[i])==location && get_second(this->_primed_assignments[i])==event) {
                reset.append(get_third(this->_primed_assignments[i]));
            }
        }

        // Compute the guards
        RealPredicate guard = select(this->_guard_predicates,location,event);

        mode.transitions.insert(event,make_tuple(guard,target,reset));
    }

    return mode;
}

RestrictiveDiscreteMode const& RestrictiveHybridAutomaton::mode(DiscreteLocation location) const {
    if(!this->_cached_modes.has_key(location)) {
        this->_cached_modes.insert(location,this->compute_mode(location));
    }
    return this->_cached_modes[location];
}

Set<DiscreteLocation> RestrictiveHybridAutomaton::reachable_locations(DiscreteLocation const& location) const {
    ARIADNE_NOT_IMPLEMENTED;
}

Void RestrictiveHybridAutomaton::check_mode(DiscreteLocation) const {
    ARIADNE_NOT_IMPLEMENTED;
}


RestrictiveHybridAutomaton* RestrictiveHybridAutomaton::clone() const {
    return new RestrictiveHybridAutomaton(*this);
}


OutputStream& RestrictiveHybridAutomaton::_write(OutputStream& os) const {
    return os << "RestrictiveHybridAutomaton"
              << "(\n  updates=" << _discrete_updates
              << ",\n  auxiliary=" << _assignments
              << ",\n  dynamic=" << _dotted_assignments
              << ",\n  invariants=" << _invariant_predicates
              << ",\n  guards=" << _guard_predicates
              << ",\n  reset=" << _primed_assignments
              << "\n)\n";
}


OutputStream& operator<<(OutputStream& os, const RestrictiveDiscreteMode& mode) {
    return os << "  DiscreteMode"
              << "(\n    location=" << mode.location
              << ",\n    auxiliary=" << mode.auxiliary
              << ",\n    dynamic=" << mode.dynamic
              << ",\n    invariants=" << mode.invariants
              << ",\n    transitions=" << mode.transitions
              << "\n  )\n";
}

RestrictiveHybridAutomaton compose(const List<RestrictiveHybridAutomaton>&) {
    ARIADNE_NOT_IMPLEMENTED;
}


DiscretePredicate make_predicate(const DiscreteLocation& q) {
    DiscretePredicate p(true);
    if(q.begin()==q.end()) {
        return p;
    }
    DiscreteLocation::ConstIterator iter = q.begin();
    p=(StringVariable(iter->first)==iter->second);
    ++iter;
    for( ; iter!=q.end(); ++iter) {
        p=p && (StringVariable(iter->first)==iter->second);
    }
    return p;
}

List<PrimedStringAssignment> make_update(const DiscreteLocation& q) {
    List<PrimedStringAssignment> u;
    for(DiscreteLocation::ConstIterator iter=q.begin(); iter!=q.end(); ++iter) {
        u.append(next(StringVariable(iter->first))=iter->second);
    }
    return u;
}



Void RestrictiveHybridAutomaton::new_transition(DiscretePredicate s, DiscreteEvent e, PrimedStringAssignment u, RealPredicate g) {
    this->new_transition(s,e,make_list(u),List<PrimedRealAssignment>(),g); }
Void RestrictiveHybridAutomaton::new_transition(DiscretePredicate s, DiscreteEvent e, PrimedStringAssignment u, List<PrimedRealAssignment> r, RealPredicate g) {
    this->new_transition(s,e,make_list(u),r,g); }
Void RestrictiveHybridAutomaton::new_transition(DiscretePredicate s, DiscreteEvent e, List<PrimedStringAssignment> u,
                                                List<PrimedRealAssignment> r, RealPredicate g) {
    this->new_guard(s,e,g);
    this->new_update(s,e,u);
    this->new_reset(s,e,r);
}


Void RestrictiveHybridAutomaton::new_mode(DiscreteLocation q, List<RealAssignment> a, List<DottedRealAssignment> d) {
    DiscretePredicate p=make_predicate(q);
    this->new_auxiliary(p,a);
    this->new_dynamic(p,d);
}
Void RestrictiveHybridAutomaton::new_mode(DiscreteLocation q, List<DottedRealAssignment> d, List<RealAssignment> a) {
    new_mode(q,a,d);
}


Void RestrictiveHybridAutomaton::new_transition(DiscreteLocation s, DiscreteEvent e, DiscreteLocation t, List<PrimedRealAssignment> r, RealPredicate g) {
    this->new_transition(make_predicate(s),e,make_update(t),r,g); }
Void RestrictiveHybridAutomaton::new_transition(DiscreteLocation s, DiscreteEvent e, DiscreteLocation t, PrimedRealAssignment r, RealPredicate g) {
    this->new_transition(s,e,t,make_list(r),g); }
Void RestrictiveHybridAutomaton::new_transition(DiscreteLocation s, DiscreteEvent e, DiscreteLocation t, RealPredicate g) {
    this->new_transition(s,e,t,List<PrimedRealAssignment>(),g); }




EffectiveVectorMultivariateFunction
dynamic_function(Space<Real>& space, const List<RealAssignment>& algebraic, const List<DottedRealAssignment>& differential)
{
    RealExpression default_expression;
    Vector<RealExpression> results(differential.size(),default_expression);
    for(SizeType i=0; i!=differential.size(); ++i) {
        results[space.index(differential[i].lhs.base())]=substitute(differential[i].rhs,algebraic);
    }

    return make_function(space,results);
}

















#ifdef ARIADNE_DISABLE


inline Map<RealVariable,RealInterval> make_map(const List<RealVariableInterval>& b) {
    Map<RealVariable,RealInterval> res;
    for(Nat i=0; i!=b.size(); ++i) {
        res.insert(b[i].variable(),RealInterval(b[i].lower(),b[i].upper()));
    }
    return res;
}

DiscreteLocation evaluate(const DiscreteUpdate& update, const DiscreteLocation& location) {
    DiscreteLocation result;
    for(SizeType i=0; i!=update.size(); ++i) {
        result.insert(update[i].variable().base(),evaluate(update[i].expression(),location));
    }
    return result;
}

Map<DiscreteEvent,DiscreteLocation> evaluate(const Map<DiscreteEvent,DiscreteUpdate>& updates, const DiscreteLocation& location) {
    Map<DiscreteEvent,DiscreteLocation> result;
    for(Map<DiscreteEvent,DiscreteUpdate>::ConstIterator iter=updates.begin(); iter!=updates.end(); ++iter) {
        result.insert(iter->first,evaluate(iter->second,location));
    }
    return result;
}



List<DottedRealVariable> dot(const List<RealVariable>& v) {
    List<DottedRealVariable> result;
    for(SizeType i=0; i!=v.size(); ++i) { result.append(dot(v[i])); }
    return result;
}

List<PrimedRealVariable> next(const List<RealVariable>& v) {
    List<PrimedRealVariable> result;
    for(SizeType i=0; i!=v.size(); ++i) { result.append(next(v[i])); }
    return result;
}

List<RealVariable> base(const List<PrimedRealVariable>& v) {
    List<RealVariable> result;
    for(SizeType i=0; i!=v.size(); ++i) { result.append(v[i].base()); }
    return result;
}

List<RealVariable> base(const List<DottedRealVariable>& v) {
    List<RealVariable> result;
    for(SizeType i=0; i!=v.size(); ++i) { result.append(v[i].base()); }
    return result;
}






DiscreteMode::
DiscreteMode(const DiscreteLocation& location)
    : _location(location)
{
}

DiscreteMode::
DiscreteMode(const DiscreteLocation& location,
             const List<RealAssignment>& equations,
             const List<DottedRealAssignment>& dynamic)
    : _location(location), _assignments(equations), _dotted_assignments(dynamic)
{
}

OutputStream&
operator<<(OutputStream& os, const DiscreteMode& mode)
{
    os << "DiscreteMode( "
       << "location=" << mode._location;
    if(mode._assignments.size()>0) {
        os << ", algebraic_equations="<<mode._assignments; }
    if(mode._dotted_assignments.size()>0) {
        os << ", differential_equations="<<mode._dotted_assignments; }
    if(mode._invariant_predicates.size()>0) {
        os << ", invariants="<<mode._invariant_predicates; }
    if(mode._guard_predicates.size()>0) {
        os << ", guards="<<mode._guard_predicates; }
    if(mode._targets.size()>0) {
        os << ", targets="<<mode._targets; }
    if(mode._primed_assignments.size()>0) {
        os << ", resets="<<mode._primed_assignments; }
    return os << " )";
}

Bool
operator<(const DiscreteMode& mode1, const DiscreteMode& mode2)
{
    return mode1.location() < mode2.location();
}


DiscreteTransition::
DiscreteTransition(const DiscreteEvent& event,
                   const DiscreteLocation& source,
                   const DiscreteLocation& target,
                   const List<PrimedRealAssignment>& reset,
                   const ContinuousPredicate& guard,
                   const EventKind& kind)
    : _event(event), _source(source), _target(target),
      _guard_predicate(guard), _primed_assignments(reset), _kind(kind)
{
}

OutputStream&
operator<<(OutputStream& os, const DiscreteTransition& transition)
{
    return os << "DiscreteTransition( "
              << "event=" << transition._event << ", "
              << "source=" << transition._source << ", "
              << "target=" << transition._target << ", "
              << "reset=" << transition._primed_assignments << ", "
              << "guard=" << transition._guard_predicate << ", "
              << "kind=" << transition._kind << " )";
}




CompositionalHybridAutomaton::~CompositionalHybridAutomaton()
{
}

CompositionalHybridAutomaton*
CompositionalHybridAutomaton::clone() const
{
    return new CompositionalHybridAutomaton(*this);
}

CompositionalHybridAutomaton::
CompositionalHybridAutomaton(const List<StringVariable>& discrete_variables, const List<DiscreteEvent>& discrete_event_list)
{
    Set<DiscreteEvent> discrete_events(discrete_event_list);
    for(List<StringVariable>::ConstIterator variable_iter=discrete_variables.begin();
        variable_iter!=discrete_variables.end(); ++variable_iter)
    {
        this->_jumping_events.insert(*variable_iter,discrete_events);
    }
}


const DiscreteMode&
CompositionalHybridAutomaton::mode(DiscreteLocation location) const
{
    if(!this->_cached_modes.has_key(location)) {
        DiscreteMode mode(location);
        mode._assignments = this->_order_algebraic_assignments(filter(this->_assignments,location));
        mode._dotted_assignments = filter(this->_dotted_assignments,location);
        mode._invariant_predicates = filter1(this->_invariant_predicates,location);
        mode._guard_predicates = filter1(this->_guard_predicates,location);
        mode._primed_assignments = filter(this->_primed_assignments,location);
        mode._targets = evaluate(filter(this->_discrete_updates,location),location);
        this->_cached_modes.insert(location,mode);
    }
    return this->_cached_modes.find(location)->second;
}


Bool
CompositionalHybridAutomaton::has_mode(DiscreteLocation location) const
{
    try {
        this->check_mode(location);
        return true;
    }
    catch(const SystemSpecificationError& e) {
        return false;
    }
}

Bool
CompositionalHybridAutomaton::has_invariant(DiscreteLocation q, DiscreteEvent e) const
{
    return this->mode(q)._invariant_predicates.has_key(e);
}

Bool
CompositionalHybridAutomaton::has_guard(DiscreteLocation q, DiscreteEvent e) const
{
    return this->mode(q)._guard_predicates.has_key(e);
}

Bool
CompositionalHybridAutomaton::has_transition(DiscreteLocation q, DiscreteEvent e) const
{
    return this->mode(q)._targets.has_key(e) || this->mode(q)._primed_assignments.has_key(e);
}

Set<DiscreteEvent>
CompositionalHybridAutomaton::guard_events(DiscreteLocation q) const
{
    Set<DiscreteEvent> result;
    for(SizeType i=0; i!=this->_guard_predicates.size(); ++i) {
        if(get_first(this->_guard_predicates[i]) == q) {
            result.insert(get_second(this->_guard_predicates[i]));
        }
    }
    return result;
}

Set<DiscreteEvent>
CompositionalHybridAutomaton::invariant_events(DiscreteLocation q) const
{
    Set<DiscreteEvent> result;
    for(SizeType i=0; i!=this->_invariant_predicates.size(); ++i) {
        if(get_first(this->_invariant_predicates[i]) == q) {
            result.insert(get_second(this->_invariant_predicates[i]));
        }
    }
    return result;
}







OutputStream&
CompositionalHybridAutomaton::_write(OutputStream& os) const
{
    os << "\nHybridAutomaton( \n";
    os << "  discrete_transitions="<<this->_discrete_updates<<"\n";
    os << "  algebraic_equations="<<this->_assignments<<"\n";
    os << "  invariant_predicates="<<this->_invariant_predicates<<"\n";
    os << "  guard_conditions="<<this->_guard_predicates<<"\n";
    os << "  differential_equations="<<this->_dotted_assignments<<"\n";
    os << "  reset_equations="<<this->_primed_assignments<<"\n";
    return os << ")\n";
}



Set<UntypedVariable>
CompositionalHybridAutomaton::argument_variables(DiscreteLocation location) const {
    return join(arguments(this->differential_assignments(location)),arguments(this->algebraic_assignments(location)));
}

List<RealVariable>
CompositionalHybridAutomaton::state_variables(DiscreteLocation location) const {
    return left_hand_sides(this->differential_assignments(location));
}

List<RealVariable>
CompositionalHybridAutomaton::auxiliary_variables(DiscreteLocation location) const {
    return left_hand_sides(this->algebraic_assignments(location));
}

Set<UntypedVariable>
CompositionalHybridAutomaton::input_variables(DiscreteLocation location) const {
    Set<UntypedVariable> result = this->argument_variables(location);
    result.remove(this->auxiliary_variables(location));
    result.remove(this->state_variables(location));
    return result;
}



List<RealAssignment>
CompositionalHybridAutomaton::algebraic_assignments(DiscreteLocation location) const {
    return this->mode(location)._assignments;
}

List<DottedRealAssignment>
CompositionalHybridAutomaton::differential_assignments(DiscreteLocation location) const {
    return this->mode(location)._dotted_assignments;
}

List<PrimedRealAssignment>
CompositionalHybridAutomaton::update_assignments(DiscreteLocation source, DiscreteEvent event) const {
    return this->mode(source)._primed_assignments[event];
}

ContinuousPredicate
CompositionalHybridAutomaton::invariant_predicate(DiscreteLocation location, DiscreteEvent action) const {
    return this->mode(location)._invariant_predicates[action];
}

ContinuousPredicate
CompositionalHybridAutomaton::guard_predicate(DiscreteLocation source, DiscreteEvent event) const {
    return this->mode(source)._guard_predicates[event];
}




class CompositeHybridSpace
    : public HybridSpaceInterface
{
  public:
    ~CompositeHybridSpace() { _system_ptr = 0; }
    CompositeHybridSpace(const CompositeHybridAutomaton& ha) : _system_ptr(&ha) { }
    virtual CompositeHybridSpace* clone() const { return new CompositeHybridSpace(*this); }
    virtual Bool has_location(const DiscreteLocation& q) const { return this->_system_ptr->has_mode(q); }
    virtual RealSpace operator[](const DiscreteLocation& q) const { return this->_system_ptr->continuous_state_space(q); }
    virtual OutputStream& _write(OutputStream& os) const { return os << "CompositeHybridSpace( " << *this->_system_ptr << " )"; }
  private:
    const CompositeHybridAutomaton* _system_ptr;
};


CompositeHybridAutomaton::CompositeHybridAutomaton()
    : _components() { }





Set<DiscreteEvent>
CompositionalHybridAutomaton::events(DiscreteLocation location) const
{
    // What should the semantics be here?
    //return join(join(this->invariant_events(location),this->guard_events(location)),this->transition_events(location));
    return join(this->invariant_events(location),this->guard_events(location));
}


DiscreteLocation
CompositionalHybridAutomaton::_compute_target(DiscreteLocation source, DiscreteEvent event) const {
    DiscreteLocation result;
    for(SizeType i=0; i!=this->_discrete_updates.size(); ++i) {
        if(get_second(this->_discrete_updates[i]) == event && get_first(this->_discrete_updates[i]) == source) {
            DiscreteUpdate const& update = get_third(_discrete_updates[i]);
            for(Nat j=0; j!=update.size(); ++j) {
                result.insert(update[j].variable().base(), evaluate(update[j].expression(),source));
            }
        }
    }
    return result;
}


Map<DiscreteEvent, DiscreteLocation>
CompositionalHybridAutomaton::_compute_targets(DiscreteLocation source) const {
    Map<DiscreteEvent,DiscreteLocation> result;
    for(SizeType i=0; i!=this->_discrete_updates.size(); ++i) {
        if(get_first(this->_discrete_updates[i]) == source) {
            DiscreteEvent const& event = get_second(_discrete_updates[i]);
            DiscreteUpdate const& update = get_third(_discrete_updates[i]);
            DiscreteLocation& target = result[event];
            for(Nat j=0; j!=update.size(); ++j) {
                target.insert(update[j].variable().base(), evaluate(update[j].expression(),source));
            }
        }
    }
    return result;
}


DimensionType
CompositeHybridAutomaton::dimension(DiscreteLocation location) const {
    return this->state_variables(location).size();
}

HybridSpace
CompositeHybridAutomaton::state_space() const {
    return new CompositeHybridSpace(*this);
}

RealSpace
CompositeHybridAutomaton::continuous_state_space(DiscreteLocation location) const {
    return RealSpace(this->state_variables(location));
}

List<RealVariable>
CompositeHybridAutomaton::variables(DiscreteLocation location) const {
    return catenate(this->state_variables(location),this->auxiliary_variables(location));
}

List<RealVariable>
CompositeHybridAutomaton::state_variables(DiscreteLocation location) const {
    List<RealVariable> result;
    for(SizeType i=0; i!=this->_components.size(); ++i) {
        ARIADNE_ASSERT(this->_components[i].has_mode(location[i]));
        // Append state space in mode i; note that this automatically checks whether a variable is already present
        result.append(this->_components[i].state_variables(location[i]));
    }
    return result;
}

DiscreteLocation
CompositionalHybridAutomaton::target(DiscreteLocation source, DiscreteEvent event) const {
    if(this->_cached_modes.has_key(source)) {
        return _cached_modes.at(source)._targets[event];
    } else {
        return this->_compute_target(source,event);
    }
}

EventKind
CompositionalHybridAutomaton::event_kind(DiscreteLocation location, DiscreteEvent event) const
{
    const DiscreteMode& mode = this->mode(location);
    if(mode._invariant_predicates.has_key(event)) {
        return PROGRESS;
    } else if(mode._guard_predicates.has_key(event)) {
        assert(false); // FIXME; return URGENT;
    } else {
        ARIADNE_THROW(std::runtime_error,"CompositionalHybridAutomaton::event_kind(...)",
                      "No event "<<event<<" defined in location "<<location);
    }
}





DimensionType
CompositionalHybridAutomaton::dimension(DiscreteLocation location) const {
    return this->state_variables(location).size();
}


EffectiveVectorMultivariateFunction
CompositionalHybridAutomaton::output_function(DiscreteLocation location) const {
    Space<Real> space=this->state_variables(location);
    List<RealAssignment> algebraic=this->algebraic_assignments(location);
    List<RealExpression> results(algebraic.size(),RealExpression(0.0));
    for(SizeType i=0; i!=algebraic.size(); ++i) { results[i]=algebraic[i].variable(); }
    return EffectiveVectorMultivariateFunction(Ariadne::dimension(space),Ariadne::formula(results,algebraic,space));
}

EffectiveVectorMultivariateFunction
CompositionalHybridAutomaton::dynamic_function(DiscreteLocation location) const {
    Space<Real> space=this->state_variables(location);
    List<RealAssignment> algebraic=this->algebraic_assignments(location);
    List<DottedRealAssignment> differential=this->differential_assignments(location);
    List<RealExpression> results(differential.size(),RealExpression(0.0));
    for(SizeType i=0; i!=differential.size(); ++i) { results[space.index(differential[i].variable().base())]=substitute(differential[i].expression(),algebraic); }
    return EffectiveVectorMultivariateFunction(Ariadne::dimension(space),Ariadne::formula(results,algebraic,space));
}

EffectiveVectorMultivariateFunction
CompositionalHybridAutomaton::reset_function(DiscreteLocation source, DiscreteEvent event) const {
    DiscreteLocation target=this->target(source,event);
    Space<Real> source_space=this->state_variables(source);
    Space<Real> target_space=this->state_variables(target);
    List<RealAssignment> algebraic=this->algebraic_assignments(source);
    List<PrimedRealAssignment> update=this->update_assignments(source,event);
    List<RealExpression> results(update.size(),RealExpression(0.0));
    for(SizeType i=0; i!=update.size(); ++i) { results[target_space.index(update[i].variable().base())]=update[i].expression(); }
    return EffectiveVectorMultivariateFunction(Ariadne::dimension(source_space),Ariadne::formula(results,algebraic,source_space));
}

EffectiveScalarMultivariateFunction
CompositionalHybridAutomaton::invariant_function(DiscreteLocation location, DiscreteEvent event) const {
    Space<Real> space=this->state_variables(location);
    List<RealAssignment> algebraic=this->algebraic_assignments(location);
    RealExpression invariant=indicator(invariant_predicate(location,event),Sign::NEGATIVE);
    return EffectiveScalarMultivariateFunction(Ariadne::dimension(space),Ariadne::formula(invariant,algebraic,space));
}

EffectiveScalarMultivariateFunction
CompositionalHybridAutomaton::guard_function(DiscreteLocation location, DiscreteEvent event) const {
    Space<Real> space=this->state_variables(location);
    List<RealAssignment> algebraic=this->algebraic_assignments(location);
    RealExpression guard=indicator(guard_predicate(location,event),Sign::POSITIVE);
    return EffectiveScalarMultivariateFunction(Ariadne::dimension(space),Ariadne::formula(guard,algebraic,space));
}



Grid
CompositionalHybridAutomaton::grid(DiscreteLocation location) const
{
    List<RealVariable> variables=this->state_variables(location);
    Vector<FloatDP> lengths(variables.size());
    for(SizeType i=0; i!=variables.size(); ++i) {
        lengths[i]=variables[i].resolution();
    }
    return Grid(lengths);
}

HybridGrid
CompositionalHybridAutomaton::grid() const
{
    return HybridGrid(*this);
}

Set<DiscreteEvent>
CompositeHybridAutomaton::urgent_events(DiscreteLocation location) const
{
    Set<DiscreteEvent> result;
    Set<DiscreteEvent> duplicates;
    for(SizeType i=0; i!=this->_components.size(); ++i) {
        Set<DiscreteEvent> found=this->_components[i].urgent_events(location[i]);
        duplicates.adjoin(intersection(result,found));
        result.adjoin(found);
    }
    if(!duplicates.empty()) {
        ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::urgent_events",
                      "Events "<<duplicates<<" are urgent in more than one component in location "<<location);
    }
    return result;
}


Set<DiscreteEvent>
CompositeHybridAutomaton::permissive_events(DiscreteLocation location) const
{
    Set<DiscreteEvent> result;
    for(SizeType i=0; i!=this->_components.size(); ++i) {
        result.adjoin(this->_components[i].permissive_events(location[i]));
    }
    return result;
}


>>>>>>> .r1509






// Order the assignments so that each variable depends only on previous ones.
List<RealAssignment>
CompositionalHybridAutomaton::_order_algebraic_assignments(const List<RealAssignment>& assignments)
{
    List<RealAssignment> result;

    Set< RealVariable > result_variables;
    Map< RealVariable, Set<RealVariable> > unresolved_dependencies;
    Map< RealVariable, Nat > assignment_table;

    for(SizeType i=0; i!=assignments.size(); ++i) {
        result_variables.insert(assignments[i].variable());
        assignment_table.insert(assignments[i].variable(),i);
    }

    for(SizeType i=0; i!=assignments.size(); ++i) {
        unresolved_dependencies.insert(assignments[i].variable(),intersection(real_arguments(assignments[i].expression()),result_variables));
    }

    while(!unresolved_dependencies.empty()) {
        // Check for an which has been resolved and add to list of ok equations
        Bool found=false;
        for(Map< RealVariable, Set<RealVariable> >::Iterator iter=unresolved_dependencies.begin();
            iter!=unresolved_dependencies.end(); ++iter)
        {
            if(iter->second.empty()) {
                RealVariable variable=iter->first;
                result.append(assignments[assignment_table[variable]]);
                for(Map< RealVariable, Set<RealVariable> >::Iterator rem_iter=unresolved_dependencies.begin();
                        rem_iter!=unresolved_dependencies.end(); ++rem_iter) {
                    rem_iter->second.erase(variable);
                }
                unresolved_dependencies.erase(iter);
                found=true;
                break;
            }
        }

        if(!found) {
            ARIADNE_THROW(SystemSpecificationError,"CompositionalHybridAutomaton::order_algebraic_assignments(List<Assignment>)",
                          "Cannot order assignments "<<assignments<<" due to algebraic loops.");
        }
    }

    assert(result.size()==assignments.size());
    return result;
}



Void
CompositionalHybridAutomaton::check_overspecification(DiscreteLocation location) const
{
    // Check for overspecification of state and auxiliary variables
    List<RealVariable> auxiliary_variables=filter_left_hand_sides(this->_assignments,location);
    List<RealVariable> state_variables=filter_left_hand_sides(this->_dotted_assignments,location);
    List<RealVariable> location_variables=catenate(state_variables,auxiliary_variables);

    if(!unique(location_variables)) {
        ARIADNE_THROW(OverspecifiedSystemError,"CompositionalHybridAutomaton::check_mode(DiscreteLocation)",
                      "Variables "<<duplicates(location_variables)<<" in location "<<location<<" have more than one defining equation.");
    }

    // Check overspecification of invariants and guards
    List<DiscreteEvent> invariant_predicates = filter_events(this->_invariant_predicates,location);
    List<DiscreteEvent> guard_predicates = filter_events(this->_guard_predicates,location);

    if(!unique(invariant_predicates)) {
        ARIADNE_THROW(OverspecifiedSystemError,"CompositionalHybridAutomaton::check_mode(DiscreteLocation)",
                      "Events "<<duplicates(invariant_predicates)<<" in location "<<location<<" are labels of more than one invariant.")
    }

    if(!unique(guard_predicates)) {
        ARIADNE_THROW(OverspecifiedSystemError,"CompositionalHybridAutomaton::check_mode(DiscreteLocation)",
                      "Events "<<duplicates(guard_predicates)<<" in location "<<location<<" are labels of more than one guard.")
    }

    // NOTE: Should we check duplication of invariant and guard labels?

    // Check overspecification in reset
    Map<DiscreteEvent,DiscreteLocation> targets = this->_compute_targets(location);
    for(Map<DiscreteEvent,DiscreteLocation>::ConstIterator target_iter=targets.begin(); target_iter!=targets.end(); ++target_iter) {

        // Compute reset variables and target auxiliary variables
        const DiscreteEvent& event = target_iter->first;
        const DiscreteLocation& target = target_iter->second;
        List<RealVariable> reset_variables = filter_left_hand_sides(this->_primed_assignments,location,event);
        Set<RealVariable> target_auxiliary_variables = Set<RealVariable>(filter_left_hand_sides(this->_assignments,target));

        // Check duplication of reset equations
        if(!unique(reset_variables)) {
        ARIADNE_THROW(OverspecifiedSystemError,"CompositionalHybridAutomaton::check_mode(DiscreteLocation)",
                      "Variables "<<duplicates(reset_variables)<<" have more than one defining equation in reset for event \""<<event<<"\" in location "<<location<<".");
        }

        // Check overspecification of target auxiliary variables in reset.
        if(!target_auxiliary_variables.disjoint(reset_variables)) {
            ARIADNE_THROW(OverspecifiedSystemError,"CompositionalHybridAutomaton::check_mode(DiscreteLocation)",
                          "Variables "<<intersection(Set<RealVariable>(reset_variables),target_auxiliary_variables)<<" in location "<<target<<
                          " are defined in reset "<<event<<" in location "<<location);
        }

    }

}


Void
CompositionalHybridAutomaton::check_mode(DiscreteLocation location) const {
    this->check_overspecification(location);

    List<RealAssignment> equations=filter(this->_assignments,location);
    List<DottedRealAssignment> dynamics=filter(this->_dotted_assignments,location);

    // Check for overspecification of state and auxiliary variables
    List<RealVariable> state_variables=left_hand_sides(dynamics);
    List<RealVariable> auxiliary_variables=left_hand_sides(equations);
    List<RealVariable> location_variables=catenate(state_variables,auxiliary_variables);
    Set<RealVariable> result_variables(location_variables);

    Set<RealVariable> undefined_variables=difference(real_arguments(equations),result_variables);
    if(!undefined_variables.empty()) {
        ARIADNE_THROW(UnderspecifiedSystemError,"CompositionalHybridAutomaton::check_mode(DiscreteLocation)",
                      "Variables "<<undefined_variables<<" are used in the definition of the algebraic equations "<<equations<<" in location "<<location<<" with state variables "<<state_variables<<", but are not defined.");
    }

    undefined_variables=difference(real_arguments(dynamics),result_variables);
    if(!undefined_variables.empty()) {
        ARIADNE_THROW(UnderspecifiedSystemError,"CompositionalHybridAutomaton::check_mode(DiscreteLocation)",
                      "Variables "<<undefined_variables<<" are used in the definition of the differential equations "<<dynamics<<" in location "<<location<<" with state variables "<<state_variables<<", but are not defined.");
    }

    List<RealAssignment> ordered_equations=this->_order_algebraic_assignments(equations);

    Map<DiscreteEvent,ContinuousPredicate> const& invariants=this->mode(location)._invariant_predicates;
    for(Map<DiscreteEvent,ContinuousPredicate>::ConstIterator invariant_iter=invariants.begin();
        invariant_iter!=invariants.end(); ++invariant_iter)
    {
        DiscreteEvent action=invariant_iter->first;
        const ContinuousPredicate& invariant=invariant_iter->second;
        if(!subset(real_arguments(invariant),result_variables)) {
            undefined_variables=difference(real_arguments(invariant),result_variables);
            ARIADNE_THROW(UnderspecifiedSystemError,"CompositionalHybridAutomaton::check_mode(DiscreteLocation)",
                          "Variables "<<undefined_variables<<" are used in the invariant "<<invariant<<" with label \""<<action<<"\" in location "<<location<<" with variables "<<location_variables<<", but are not defined.");
        }
    }

    Set<DiscreteEvent> events=this->events(location);
    for(Set<DiscreteEvent>::ConstIterator event_iter=events.begin(); event_iter!=events.end(); ++event_iter) {

        // Get transition information
        DiscreteEvent event=*event_iter;
        DiscreteLocation target=this->target(location,event);
        ContinuousPredicate guard=this->guard_predicate(location,event);
        List<PrimedRealAssignment> reset=this->update_assignments(location,event);

        Set<RealVariable> target_state_variables=make_set(this->state_variables(target));
        List<RealVariable> reset_variables=base(assigned(reset));

        // Check arguments of guard
        if(!subset(real_arguments(guard),result_variables)) {
            undefined_variables=difference(real_arguments(guard),result_variables);
            ARIADNE_THROW(UnderspecifiedSystemError,"CompositionalHybridAutomaton::check_mode(DiscreteLocation)",
                          "Variables "<<undefined_variables<<" are used in the guard "<<guard<<" for event \""<<event<<"\" in location "<<location<<" with variables "<<location_variables<<", but are not defined.");
        }

        // Check arguments of reset
        if(!subset(real_arguments(reset),result_variables)) {
            undefined_variables=difference(real_arguments(reset),result_variables);
            ARIADNE_THROW(UnderspecifiedSystemError,"CompositionalHybridAutomaton::check_mode(DiscreteLocation)",
                          "Variables "<<undefined_variables<<" are used in the reset "<<reset<<" for event \""<<event<<"\" in location "<<location<<" with variables "<<location_variables<<", but are not defined.");
        }

        Set<RealVariable> updated_variables=make_set(reset_variables);

        if(!subset(updated_variables,target_state_variables)) {
            Set<RealVariable> extra_variables=difference(updated_variables,target_state_variables);
            ARIADNE_THROW(OverspecifiedSystemError,"CompositionalHybridAutomaton::check_mode(DiscreteLocation)",
                          "Variables "<<extra_variables<<" are defined in the reset "<<reset<<" for event \""<<event<<"\" in location "<<location<<", but are not state variables "<<target_state_variables<<" in the target "<<target<<".");
        }

        if(!subset(target_state_variables,updated_variables)) {
            undefined_variables=difference(target_state_variables,updated_variables);
            ARIADNE_THROW(UnderspecifiedSystemError,"CompositionalHybridAutomaton::check_mode(DiscreteLocation)",
                          "Variables "<<undefined_variables<<" are state variables in location "<<target<<", but are not defined in the reset "<<reset<<" for the transition \""<<event<<"\" from location "<<location<<".");
        }

    } // Finished checking transitions

    return;
}


Void
CompositionalHybridAutomaton::check_reachable_modes(DiscreteLocation initial_location) const
{
    Set<DiscreteLocation> initial_locations;
    initial_locations.insert(initial_location);
    this->check_reachable_modes(initial_locations);
}


Void
CompositionalHybridAutomaton::check_reachable_modes(const Set<DiscreteLocation>& initial_locations) const
{
    Set<DiscreteLocation> reachable_locations=this->reachable_locations(initial_locations);
    for(Set<DiscreteLocation>::ConstIterator iter=reachable_locations.begin(); iter!=reachable_locations.end(); ++iter) {
        this->check_mode(*iter);
    }
}

Set<DiscreteLocation>
CompositionalHybridAutomaton::reachable_locations(const DiscreteLocation& initial_location) const
{
    return this->reachable_locations(Set<DiscreteLocation>(initial_location));
}

Set<DiscreteLocation>
CompositionalHybridAutomaton::reachable_locations(const Set<DiscreteLocation>& initial_locations) const
{
    const CompositionalHybridAutomaton& automaton=*this;
    typedef Nat Nat;

    Set<DiscreteLocation> reached=initial_locations;
    Set<DiscreteLocation> working=initial_locations;
    Set<DiscreteLocation> found;

    Nat step=0u;

    while(!working.empty()) {
        ++step;
        for(Set<DiscreteLocation>::ConstIterator source_iter=working.begin(); source_iter!=working.end(); ++source_iter) {
            DiscreteLocation location=*source_iter;
            ARIADNE_LOG(5,"  mode: "<<location<<":\n");

            Set<DiscreteEvent> events=automaton.guard_events(location);
            ARIADNE_LOG(5,"    events: "<<events<<"\n");
            for(Set<DiscreteEvent>::ConstIterator event_iter=events.begin(); event_iter!=events.end(); ++event_iter) {
                DiscreteEvent event=*event_iter;
                ARIADNE_LOG(5,"      transition: "<<event<<" -> ");
                DiscreteLocation target=automaton.target(location,event);
                if(verbosity>=5) { std::clog<<target<<"\n"; }
                if(!reached.contains(target)) {
                    found.insert(target);
                    reached.insert(target);
                }
           }

        }
        ARIADNE_LOG(5,"\nstep "<<step<<" found: "<<found<<"\n\n");
        working.clear();
        working.swap(found);
    }

    return reached;
}


CompositionalHybridAutomaton compose(const List<CompositionalHybridAutomaton>& components)
{
    CompositionalHybridAutomaton result;
    for(SizeType i=0; i!=components.size(); ++i) {
        const CompositionalHybridAutomaton& component = components[i];
        result._discrete_events.adjoin(component._discrete_events);
        result._discrete_variables.adjoin(component._discrete_variables);
        result._continuous_variables.adjoin(component._continuous_variables);

        result._discrete_updates.append(component._discrete_updates);
        result._discrete_updates.append(component._assignments);
        result._discrete_updates.append(component._dotted_assignments);
        result._discrete_updates.append(component._primed_assignments);
        result._discrete_updates.append(component._invariant_predicates);
        result._discrete_updates.append(component._guard_predicates);

    }

    // Introduce nonjumping constraints for discrete variables
    for(SizeType i=0; i!=components.size(); ++i) {
        for(Nat j=0; j!=component[i]._dotted_assignments.size(); ++j) {
            Tuple<DiscretePredicate,DottedRealAssignment> dynamic = _components[i]._dotted_assignments[j];
            this->_nonjumping[i]._append(dynamic.first,dynamic.second.left_hand_side());
        }
    }

    List<Pair<Int,double> > lst;
    make_map(lst);
    return result;

}


ComponentHybridAutomaton::
~ComponentHybridAutomaton()
{
}


ComponentHybridAutomaton::
ComponentHybridAutomaton(const List<StringVariable>& discrete_variables, const List<DiscreteEvent>& discrete_events)
    : CompositionalHybridAutomaton(discrete_variables,discrete_events), _variables(discrete_variables), _events(discrete_events)
{
}


Void
ComponentHybridAutomaton::
new_equation(const DiscretePredicate& locations,
             const List<RealAssignment>& equations)
{
    for(SizeType i=0; i!=equations.size(); ++i) {
        this->_assignments.append(make_tuple(locations,equations[i]));
    }
}

Void
ComponentHybridAutomaton::
new_dynamic(const DiscretePredicate& locations,
            const List<DottedRealAssignment>& dynamic)
{
    for(SizeType i=0; i!=dynamic.size(); ++i) {
        this->_dotted_assignments.append(make_tuple(locations,dynamic[i]));
    }
}

Void
ComponentHybridAutomaton::
new_invariant(const DiscretePredicate& locations,
              DiscreteEvent event,
              const ContinuousPredicate& invariant,
              EventKind kind)
{
    this->_invariant_predicates.append(make_tuple(locations,event,invariant));
}


Void
ComponentHybridAutomaton::
new_guard(const DiscretePredicate& locations,
          DiscreteEvent event,
          const ContinuousPredicate& guard,
          EventKind kind)
{
    this->_guard_predicates.append(make_tuple(locations,event,guard));
}


Void
ComponentHybridAutomaton::
new_reset(const DiscretePredicate& sources,
          DiscreteEvent event,
          const DiscreteUpdate& update,
          const List<PrimedRealAssignment>& reset)
{
    this->_discrete_updates.append(make_tuple(sources,event,update));
    for(SizeType i=0; i!=reset.size(); ++i) {
        this->_primed_assignments.append(make_tuple(sources,event,reset[i]));
    }
}



Void
ComponentHybridAutomaton::
new_transition(const DiscretePredicate& sources,
               DiscreteEvent event,
               const DiscreteUpdate& update,
               const List<PrimedRealAssignment>& reset,
               const ContinuousPredicate& guard,
               EventKind kind)
{
    this->_guard_predicates.append(make_tuple(sources,event,guard));
    this->_discrete_updates.append(make_tuple(sources,event,update));
    for(SizeType i=0; i!=reset.size(); ++i) {
        this->_primed_assignments.append(make_tuple(sources,event,reset[i]));
    }
}






OutputStream&
CompositeHybridAutomaton::_write(OutputStream& os) const
{
    return os << "CompositeHybridAutomaton(\n" << this->_components << "\n)\n";
}


#endif // ARIADNE_DISABLE


} // namespace Ariadne
