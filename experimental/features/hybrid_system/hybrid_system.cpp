/***************************************************************************
 *            hybrid_system.cpp
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

#include <map>

#include "../function/functional.hpp"

#include "../config.hpp"

#include "../utility/macros.hpp"
#include "../output/logging.hpp"

#include "../utility/stlio.hpp"
#include "../symbolic/expression.hpp"
#include "../function/function.hpp"
#include "../hybrid/hybrid_time.hpp"
#include "hybrid_system.hpp"
#include "../geometry/grid.hpp"

namespace Ariadne {

template<class R, class Op=OperatorCode, class A1=R, class A2=A1> class BinaryExpression
    : public ExpressionInterface<R>
{
  public:
    BinaryExpression(Op op, const ExpressionInterface<A1>& expr1, const ExpressionInterface<A2>& expr2)
        : _op(op), _arg1(expr1.clone()), _arg2(expr2.clone()) { }
    BinaryExpression(Op op, const ExpressionInterface<A1>* expr1, const ExpressionInterface<A2>* expr2)
        : _op(op), _arg1(expr1), _arg2(expr2) { }
    BinaryExpression(Op op, std::shared_ptr< const ExpressionInterface<A1> > expr1, std::shared_ptr< const ExpressionInterface<A2> > expr2)
        : _op(op), _arg1(expr1), _arg2(expr2)  { }
    virtual String operator_name() const { return name(_op); }
    virtual OperatorCode type() const { return static_cast<OperatorCode>(_op); }
    virtual BinaryExpression<R,Op,A1,A2>* clone() const { return new BinaryExpression<R,Op,A1,A2>(_op,_arg1._ptr,_arg2._ptr); }
    virtual Set<String> arguments() const { return join(this->_arg1.arguments(),this->_arg2.arguments()); }
    virtual OutputStream& write(OutputStream& os) const {
        return os << "(" << _arg1 << symbol(_op) << _arg2 << ")"; }
  protected:
    virtual ExpressionInterface<R>* simplify() const { return this->clone(); }
  public:
    Op _op;
    Expression<A1> _arg1;
    Expression<A2> _arg2;
};


using std::make_pair;

std::vector<StringType> Event::_names=std::vector<StringType>();


HybridSystem::~HybridSystem()
{
}

HybridSystem::HybridSystem()
{
}

EventSet
HybridSystem::events() const
{
    EventSet result;
    for(guard_const_iterator iter=this->_guard_predicates.begin(); iter!=this->_guard_predicates.end(); ++iter) {
        if(iter->evnts.finite()) { result.adjoin(iter->evnts); }
    }
    return result;
}

VariableSet
HybridSystem::discrete_variables() const
{
    VariableSet variables;
    for(guard_const_iterator iter=this->_guard_predicates.begin(); iter!=this->_guard_predicates.end(); ++iter) {
        variables.adjoin(iter->loc.arguments());
    }
    return variables;
}

VariableSet
HybridSystem::result_variables(const DiscreteValuation& discrete_state) const
{
    VariableSet result;
    for(dynamic_const_iterator iter=this->_differential_equations.begin(); iter!=this->_differential_equations.end(); ++iter) {
        if(evaluate(iter->loc,discrete_state)) { result.insert(iter->lhs); }
    }
    for(relation_const_iterator iter=this->_algebraic_equations.begin(); iter!=this->_algebraic_equations.end(); ++iter) {
        if(evaluate(iter->loc,discrete_state)) { result.insert(iter->lhs); }
    }
    return result;
}

VariableSet
HybridSystem::argument_variables(const DiscreteValuation& discrete_state) const
{
    VariableSet result;
    for(dynamic_const_iterator iter=this->_differential_equations.begin(); iter!=this->_differential_equations.end(); ++iter) {
        if(evaluate(iter->loc,discrete_state)) { adjoin(result,iter->rhs.arguments()); }
    }
    for(relation_const_iterator iter=this->_algebraic_equations.begin(); iter!=this->_algebraic_equations.end(); ++iter) {
        if(evaluate(iter->loc,discrete_state)) { adjoin(result,iter->rhs.arguments()); }
    }
    return result;
}

VariableSet
HybridSystem::state_variables(const DiscreteValuation& discrete_state) const
{
    VariableSet result;
    for(dynamic_const_iterator iter=this->_differential_equations.begin(); iter!=this->_differential_equations.end(); ++iter) {
        if(evaluate(iter->loc,discrete_state)) { result.insert(iter->lhs); }
    }
    return result;
}

VariableSet
HybridSystem::algebraic_variables(const DiscreteValuation& discrete_state) const
{
    VariableSet result;
    for(relation_const_iterator iter=this->_algebraic_equations.begin(); iter!=this->_algebraic_equations.end(); ++iter) {
        if(evaluate(iter->loc,discrete_state)) { result.insert(iter->lhs); }
    }
    return result;
}

VariableSet
HybridSystem::auxiliary_variables(const DiscreteValuation& discrete_state) const
{
    return intersection(this->algebraic_variables(discrete_state),this->argument_variables(discrete_state));
}

VariableSet
HybridSystem::input_variables(const DiscreteValuation& discrete_state) const
{
    return difference(this->argument_variables(discrete_state),this->result_variables(discrete_state));
}

VariableSet
HybridSystem::output_variables(const DiscreteValuation& discrete_state) const
{
    return difference(this->argument_variables(discrete_state),this->result_variables(discrete_state));
}

template<class T>
List<T>
flatten(Map< T, Set<T> >& dag)
{
    List<T> result;
    while(!dag.empty()) {
        Set<T> visited;
        T t=dag.begin()->first;
        while(!dag[t].empty()) {
            visited.insert(t);
            t=*dag[t].begin();
            if(contains(visited,t)) {
                ARIADNE_THROW(std::runtime_error,"flatten","Algebraic loop involving variable \""<<t<<"\" with dependencies "<<dag);
            }
        }
        result.push_back(t);
        dag.erase(t);
        for(typename Map<T,Set<T> >::Iterator iter=dag.begin(); iter!=dag.end(); ++iter) {
            iter->second.erase(t);
        }
    }
    return result;
}

DiscreteValuation HybridSystem::target(const Event& event, const DiscreteValuation& source) const {
    DiscreteValuation target;
    for(switch_const_iterator iter=_discrete_updates.begin(); iter!=this->_discrete_updates.end(); ++iter) {
        if(iter->evnts.contains(event)) {
            if(evaluate(iter->loc,source)) {
                String val=evaluate(iter->rhs,source);
                target.set(iter->lhs,val);
            }
        }
    }
    return target;
}


Bool HybridSystem::check_dynamic(const DiscreteValuation& location) const
{
    typedef Map<UntypedVariable,Set<UntypedVariable> >::Iterator dependencies_iterator;

    Set<UntypedVariable> independent_variables;
    Set<UntypedVariable> differential_variables;
    Set<UntypedVariable> algebraic_variables;
    Set<UntypedVariable> auxiliary_variables;
    Set<UntypedVariable> input_variables;
    Set<UntypedVariable> output_variables;

    Map<UntypedVariable,Set<UntypedVariable> > differential_dependencies;
    Map<UntypedVariable,Set<UntypedVariable> > algebraic_dependencies;

    for(dynamic_const_iterator iter=this->_differential_equations.begin(); iter!=this->_differential_equations.end(); ++iter) {
        Bool active;
        try {
            active=evaluate(iter->loc,location);
        } catch(...) {
            ARIADNE_ASSERT_MSG(false,"cannot determine activation of "<<*iter<<" in location "<<location);
        }
        if(active) {
            ARIADNE_ASSERT_MSG(!has_key(differential_dependencies,iter->lhs),"Variable "<<iter->lhs.name()<<" is assigned in two different rules");
            differential_variables.insert(iter->lhs);
            differential_dependencies.insert(make_pair(iter->lhs,iter->rhs.arguments()));
            adjoin(independent_variables,iter->rhs.arguments());
        }
    }

    for(relation_const_iterator iter=this->_algebraic_equations.begin(); iter!=this->_algebraic_equations.end(); ++iter) {
        if(evaluate(iter->loc,location)) {
            ARIADNE_ASSERT_MSG(!has_key(differential_dependencies,iter->lhs),"");
            ARIADNE_ASSERT_MSG(!has_key(algebraic_dependencies,iter->lhs),"");
            algebraic_variables.insert(iter->lhs);
            algebraic_dependencies.insert(make_pair(iter->lhs,iter->rhs.arguments()));
            adjoin(independent_variables,iter->rhs.arguments());
       }
    }

    ARIADNE_LOG(3,"differential_variables="<<differential_variables<<"\n");
    ARIADNE_LOG(3,"algebraic_variables="<<algebraic_variables<<"\n");
    ARIADNE_LOG(3,"independent_variables="<<independent_variables<<"\n");

    input_variables=difference(independent_variables,join(differential_variables,algebraic_variables));
    ARIADNE_LOG(3,"input_variables="<<independent_variables<<"\n");
    ARIADNE_ASSERT_MSG(input_variables.empty(),"Variables "<<input_variables<<" are used, but have no defining rules");

    ARIADNE_LOG(3,"algebraic_dependencies="<<algebraic_dependencies<<"\n");
    for(dependencies_iterator iter=algebraic_dependencies.begin(); iter!=algebraic_dependencies.end(); ++iter) {
        restrict(iter->second,algebraic_variables);
    }
    std::vector<Identifier> ordered_algebraic_variables=flatten(algebraic_dependencies);
    ARIADNE_LOG(3,"ordered_algebraic_variables="<<ordered_algebraic_variables<<"\n");

    return true;
}

Bool HybridSystem::check_reset(const Event& event, const DiscreteValuation& source, const DiscreteValuation& target) const
{
    Set<Identifier> source_variables=join(this->state_variables(source),this->auxiliary_variables(source));
    Set<Identifier> target_state_variables=this->state_variables(source);
    Set<Identifier> updated_variables;

    for(jump_const_iterator iter=this->_continuous_updates.begin(); iter!=this->_continuous_updates.end(); ++iter) {
        if(iter->evnts.contains(event) && evaluate(iter->loc,source)) {
            ARIADNE_ASSERT_MSG(subset(iter->rhs.arguments(),source_variables),"");
            updated_variables.insert(iter->lhs);
        }
    }

    ARIADNE_LOG(3,"source_variables="<<source_variables<<"\n");
    ARIADNE_LOG(3,"target_state_variables="<<target_state_variables<<"\n");
    ARIADNE_LOG(3,"updated_variables="<<updated_variables<<"\n");

    ARIADNE_ASSERT_MSG(subset(updated_variables,target_state_variables),"");
    ARIADNE_ASSERT_MSG(subset(target_state_variables,updated_variables),"");

    return true;
}


Bool HybridSystem::check_guards(const DiscreteValuation& location) const
{
    Set<Identifier> variables=join(this->state_variables(location),this->auxiliary_variables(location));

    EventSet unguarded_events=EventSet::all();

    for(guard_const_iterator iter=this->_guard_predicates.begin(); iter!=this->_guard_predicates.end(); ++iter) {
        if(evaluate(iter->loc,location)) {
            if(iter->evnts.finite()) {
                ARIADNE_ASSERT_MSG(subset(iter->pred.arguments(),variables),"");
            } else {
                const Expression<Kleenean>& pred=iter->pred;
                if(!(pred==false)) { ARIADNE_ASSERT_MSG(true,"Nontrivial guard applied to an infinite set of events"); }
            }
            unguarded_events.remove(iter->evnts);
        }
    }
    ARIADNE_ASSERT_MSG(unguarded_events.size()==0,"Events "<<unguarded_events<<" have no active guard in location "<<location);
    return true;
}



/* Functions used to create a hybrid automaton. */
EventSet
HybridSystem::events(const DiscreteValuation& location) const
{
    EventSet events;
    for(jump_const_iterator iter=this->_continuous_updates.begin(); iter!=this->_continuous_updates.end(); ++iter) {
        if(evaluate(iter->loc,location)) {
            const RealVariable& lhs=iter->lhs;
            const RealExpression& rhs=iter->rhs;
            const RealExpression lhse=RealExpression(lhs);
            if(!identical(lhse,rhs)) {
                events.adjoin(iter->evnts);
            }
            else {
                ARIADNE_WARN("Reset "<<*iter<<" is implicitly trivial.");
            }
        }
    }
    ARIADNE_ASSERT(events.finite());
    return events;
}


VectorFunction
HybridSystem::dynamic(const DiscreteValuation& location) const
{
    Set<RealAssignment> algebraic_equations=this->unordered_equations(location);
    ARIADNE_ASSERT_MSG(algebraic_equations.empty(),"Current implementation can only compute dynamic if there are no algebraic equations.");
    List<RealDynamic> differential_equations=this->dynamics(location);
    RealSpace space=this->state_variables(location);

    List<RealExpression> expressions;
    for(Nat i=0; i!=differential_equations.size(); ++i) {
        const RealExpression& rhs=differential_equations[i].rhs;
        expressions.push_back(rhs);
        //expressions.push_back(differential_equations[i].rhs());
    }

    return VectorFunction(expressions,space);
}


VectorFunction
HybridSystem::reset(const Event& event, const DiscreteValuation& source) const
{
    DiscreteValuation target = this->target(event,source);
    Set<RealAssignment> source_algebraic_equations=this->unordered_equations(source);
    ARIADNE_ASSERT_MSG(source_algebraic_equations.empty(),"Current implementation can only compute reset if source has no algebraic equations.");
    List<RealUpdate> update_equations=this->resets(event,source);
    RealSpace source_space=this->state_variables(source);
    RealSpace target_space=this->state_variables(source);

    Map<RealVariable,RealExpression> assignments;
    for(List<RealUpdate>::ConstIterator update_iter=update_equations.begin(); update_iter!=update_equations.end(); ++update_iter) {
        assignments.insert(update_iter->lhs.base(),update_iter->rhs);
    }

    ARIADNE_NOT_IMPLEMENTED
    //return VectorFunction(target_space,assignments,source_space);
}


ScalarFunction
HybridSystem::guard(const Event& event, const DiscreteValuation& location) const
{
    ContinuousPredicate guard_predicate=this->guard_predicate(event,location);
    RealSpace space=this->state_variables(location);

    RealExpression expression(0.0);
    const ExpressionInterface<Kleenean>* ptr=guard_predicate._ptr.operator->();
    const BinaryExpression<Kleenean,Gtr,Real,Real>* bptr=dynamic_cast< const BinaryExpression<Kleenean,Gtr,Real,Real>*>(ptr);
    if(bptr) {
        expression = bptr->_arg1 - bptr->_arg2;
    }
    //if(guard_predicate.operator_name()=="GTR") {
    //    List< RealExpression > subexpressions = guard_predicate.subexpressions() const;
    //    expression=subexpressions[0]-subexpressions[1];
    //} else {
    //    ARIADNE_ASSERT_MSG(false,"Current implementation can only compute dynamic of form g(x)>0; guard is "<<guard_predicate);
    //}

    return ScalarFunction(expression,space);
}





Set<RealAssignment>
HybridSystem::unordered_equations(const DiscreteValuation& state) const
{
    Set<RealVariable> variables;
    Set<RealAssignment> equations;
    for(relation_const_iterator iter=this->_algebraic_equations.begin(); iter!=this->_algebraic_equations.end(); ++iter) {
        if(evaluate(iter->loc,state)) {
            ARIADNE_ASSERT(!contains(variables,iter->lhs));
            variables.insert(iter->lhs);
            equations.insert(iter->lhs=iter->rhs);
       }
    }
    return equations;
}


List<RealAssignment>
HybridSystem::equations(const DiscreteValuation& state) const
{
    typedef Map<Identifier,Set<Identifier> >::Iterator dependencies_iterator;

    Set<Identifier> variables;
    Map<RealVariable,RealExpression> formulae;
    for(relation_const_iterator iter=this->_algebraic_equations.begin(); iter!=this->_algebraic_equations.end(); ++iter) {
        if(evaluate(iter->loc,state)) {
            ARIADNE_ASSERT(!has_key(formulae,iter->lhs));
            variables.insert(iter->lhs.name());
            formulae.insert(make_pair(iter->lhs,iter->rhs));
       }
    }

    Map<Identifier,Set<Identifier> > dependencies;
    for(Map<RealVariable,RealExpression>::ConstIterator iter=formulae.begin(); iter!=formulae.end(); ++iter) {
        dependencies.insert(iter->first.name(),iter->second.arguments().restrict(variables));
    }

    List<Identifier> ordering=flatten(dependencies);

    List<RealAssignment> equations;
    for(List<Identifier>::ConstIterator iter=ordering.begin(); iter!=ordering.end(); ++iter)
    {
        for(Map<RealVariable,RealExpression>::Iterator fiter=formulae.begin(); fiter!=formulae.end(); ++fiter) {
            if(fiter->first.name()==*iter) { equations.push_back(fiter->first=fiter->second); }
        }
    }

    return equations;
}



List<RealDynamic>
HybridSystem::dynamics(const DiscreteValuation& location) const
{
    List<RealDynamic> dynamic_equations;

    for(dynamic_const_iterator iter=this->_differential_equations.begin(); iter!=this->_differential_equations.end(); ++iter) {
        if(evaluate(iter->loc,location)) {
            dynamic_equations.push_back(dot(iter->lhs)=iter->rhs);
        }
    }

    return dynamic_equations;
}


List<StringUpdate>
HybridSystem::switching(const Event& event, const DiscreteValuation& location) const
{
    List<StringUpdate> switch_equations;

    for(switch_const_iterator iter=this->_discrete_updates.begin(); iter!=this->_discrete_updates.end(); ++iter) {
        if(iter->evnts.contains(event)) {
            if(evaluate(iter->loc,location)) {
                switch_equations.push_back(next(iter->lhs)=iter->rhs);
            }
        }
    }
    return switch_equations;
}



List<RealUpdate>
HybridSystem::resets(const Event& event, const DiscreteValuation& old_location) const
{
    List<RealUpdate> continuous_updates;

    for(jump_const_iterator iter=this->_continuous_updates.begin(); iter!=this->_continuous_updates.end(); ++iter) {
        if(iter->evnts.contains(event)) {
            if(evaluate(iter->loc,old_location)) {
                continuous_updates.push_back(next(iter->lhs)=iter->rhs);
            }
        }
    }

    return continuous_updates;
}



Map<Event,ContinuousPredicate>
HybridSystem::guards(const DiscreteValuation& location) const
{
    // Compute possibly active events, and ensure that only finitely many events can occur
    EventSet unguarded_events=EventSet::all();
    Map<Event,ContinuousPredicate> guards;

    for(guard_const_iterator iter=this->_guard_predicates.begin(); iter!=this->_guard_predicates.end(); ++iter) {
        if(evaluate(iter->loc,location)) {
            if(iter->evnts.finite()) {
                ARIADNE_ASSERT_MSG(iter->evnts.size()==1,"Guard "<<*iter<<" is active for multiple events.");
                Event e=iter->evnts.front();
                unguarded_events.remove(e);
                if(has_key(guards,e)) {
                    ContinuousPredicate& guard=guards.find(e)->second;
                    guard = guard && iter->pred;
                } else {
                    guards.insert(make_pair(e,iter->pred));
                }
            } else {
                // For an infinite set of events, the guard must be identically false
                assert(iter->pred==false);
                unguarded_events.remove(iter->evnts);
            }
        }
    }
    //ARIADNE_WARN_MSG(unguarded_events.empty(),"Events "<<unguarded_events<<" in location "<<location<<" have no guard predicate");

    return guards;
}


ContinuousPredicate
HybridSystem::guard_predicate(const Event& event, const DiscreteValuation& location) const
{
    // TODO: Make this more efficient
    const Map<Event,ContinuousPredicate> grds=this->guards(location);
    return grds[event];
}


HybridSystem parallel_composition(const HybridSystem& sys1, const HybridSystem& sys2) {
    // Parallel composition in the new model is easy; we just combine all the system rules!
    HybridSystem result;
    result._differential_equations=catenate(sys1._differential_equations,sys2._differential_equations);
    result._algebraic_equations=catenate(sys1._algebraic_equations,sys2._algebraic_equations);
    result._discrete_updates=catenate(sys1._discrete_updates,sys2._discrete_updates);
    result._continuous_updates=catenate(sys1._continuous_updates,sys2._continuous_updates);
    result._guard_predicates=catenate(sys1._guard_predicates,sys2._guard_predicates);
    result._invariant_predicates=catenate(sys1._invariant_predicates,sys2._invariant_predicates);
    result._disabled_events=catenate(sys1._disabled_events,sys2._disabled_events);
    return result;
}


OutputStream& operator<<(OutputStream& os, const HybridSystem& sys) {
    os << std::boolalpha << "HybridSystem(\n"
       << "  algebraic_equations=" << sys._algebraic_equations << ",\n"
       << "  differential_equations=" << sys._differential_equations << ",\n"
       << "  discrete_resets=" << sys._discrete_updates << ",\n"
       << "  continuous_resets=" << sys._continuous_updates << ",\n"
       << "  guard_predicates=" << sys._guard_predicates << ",\n"
       << "  invariant_predicates=" << sys._invariant_predicates << "\n)\n";
    return os;
}

OutputStream& operator<<(OutputStream& os, const HybridSystem::DifferentialEquation& de) {
    os << de.loc << " -> dot("<<de.lhs.name()<<")="<<de.rhs;
    return os;
}

OutputStream& operator<<(OutputStream& os, const HybridSystem::AlgebraicEquation& ae) {
    os << ae.loc << " -> "<<ae.lhs<<"="<<ae.rhs;
    return os;
}

OutputStream& operator<<(OutputStream& os, const HybridSystem::DiscreteUpdate& da) {
    os << da.evnts << ": " << da.loc << " -> next("<<da.lhs.name()<<")="<<da.rhs;
    return os;
}

OutputStream& operator<<(OutputStream& os, const HybridSystem::ContinuousUpdate& da) {
    os << da.evnts << ": " << da.loc << " -> next("<<da.lhs.name()<<")="<<da.rhs;
    return os;
}

OutputStream& operator<<(OutputStream& os, const HybridSystem::GuardPredicate& g) {
    os << g.evnts << ": " << g.loc << " -> "<<g.pred;
    return os;
}

OutputStream& operator<<(OutputStream& os, const HybridSystem::InvariantPredicate& inv) {
    os << inv.loc << " -> "<<inv.pred;
    return os;
}





}
