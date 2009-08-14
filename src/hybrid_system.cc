/***************************************************************************
 *            hybrid_system.cc
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

#include <map>

#include "macros.h"
#include "stlio.h"
#include "formula.h"
#include "function_interface.h"
#include "hybrid_time.h"
#include "hybrid_system.h"
#include "grid.h"

namespace Ariadne {

using std::make_pair;

std::vector<std::string> Event::_names=std::vector<std::string>();


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
        if(evaluate(iter->loc,discrete_state)) { result.insert(iter->lhs.name()); }
    }
    for(relation_const_iterator iter=this->_algebraic_equations.begin(); iter!=this->_algebraic_equations.end(); ++iter) {
        if(evaluate(iter->loc,discrete_state)) { result.insert(iter->lhs.name()); }
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
        if(evaluate(iter->loc,discrete_state)) { result.insert(iter->lhs.name()); }
    }
    return result;
}

VariableSet
HybridSystem::algebraic_variables(const DiscreteValuation& discrete_state) const
{
    VariableSet result;
    for(relation_const_iterator iter=this->_algebraic_equations.begin(); iter!=this->_algebraic_equations.end(); ++iter) {
        if(evaluate(iter->loc,discrete_state)) { result.insert(iter->lhs.name()); }
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
        for(typename Map<T,Set<T> >::iterator iter=dag.begin(); iter!=dag.end(); ++iter) {
            iter->second.erase(t);
        }
    }
    return result;
}

DiscreteValuation HybridSystem::target(const Event& event, const DiscreteValuation& source) const {
    std::cerr<<"\n";
    DiscreteValuation target;
    for(switch_const_iterator iter=_discrete_updates.begin(); iter!=this->_discrete_updates.end(); ++iter) {
        if(iter->evnts.contains(event)) {
            if(evaluate(iter->loc,source)) {
                String val=evaluate(iter->rhs,source);
                std::cerr<<iter->lhs<<":="<<val<<"="<<iter->rhs<<"\n";
                target.set(iter->lhs,val);
            }
        }
    }
    std::cerr<<"\n"<<source<<" "<<target<<"\n\n";
    return target;
}


bool HybridSystem::check_dynamic(const DiscreteValuation& location) const
{
    typedef Map<Identifier,Set<Identifier> >::iterator dependencies_iterator;

    Set<Identifier> independent_variables;
    Set<Identifier> differential_variables;
    Set<Identifier> algebraic_variables;
    Set<Identifier> auxiliary_variables;
    Set<Identifier> input_variables;
    Set<Identifier> output_variables;

    Map<Identifier,Set<Identifier> > differential_dependencies;
    Map<Identifier,Set<Identifier> > algebraic_dependencies;

    for(dynamic_const_iterator iter=this->_differential_equations.begin(); iter!=this->_differential_equations.end(); ++iter) {
        bool active;
        try {
            active=evaluate(iter->loc,location);
        } catch(...) {
            ARIADNE_ASSERT_MSG(false,"cannot determine activation of "<<*iter<<" in location "<<location);
        }
        if(active) {
            ARIADNE_ASSERT_MSG(!has_key(differential_dependencies,iter->lhs.name()),"Variable "<<iter->lhs.name()<<" is assigned in two different rules");
            differential_variables.insert(iter->lhs.name());
            differential_dependencies.insert(make_pair(iter->lhs.name(),iter->rhs.arguments()));
            adjoin(independent_variables,iter->rhs.arguments());
        }
    }

    for(relation_const_iterator iter=this->_algebraic_equations.begin(); iter!=this->_algebraic_equations.end(); ++iter) {
        if(evaluate(iter->loc,location)) {
            ARIADNE_ASSERT_MSG(!has_key(differential_dependencies,iter->lhs.name()),"");
            ARIADNE_ASSERT_MSG(!has_key(algebraic_dependencies,iter->lhs.name()),"");
            algebraic_variables.insert(iter->lhs.name());
            algebraic_dependencies.insert(make_pair(iter->lhs.name(),iter->rhs.arguments()));
            adjoin(independent_variables,iter->rhs.arguments());
       }
    }

    std::cerr<<"differential_variables="<<differential_variables<<"\n";
    std::cerr<<"algebraic_variables="<<algebraic_variables<<"\n";
    std::cerr<<"independent_variables="<<independent_variables<<"\n";

    input_variables=difference(independent_variables,join(differential_variables,algebraic_variables));
    std::cerr<<"input_variables="<<independent_variables<<"\n";
    ARIADNE_ASSERT_MSG(input_variables.empty(),"Variables "<<input_variables<<" are used, but have no defining rules");

    std::cerr<<"algebraic_dependencies="<<algebraic_dependencies<<"\n";
    for(dependencies_iterator iter=algebraic_dependencies.begin(); iter!=algebraic_dependencies.end(); ++iter) {
        restrict(iter->second,algebraic_variables);
    }
    std::vector<Identifier> ordered_algebraic_variables=flatten(algebraic_dependencies);
    std::cerr<<"ordered_algebraic_variables="<<ordered_algebraic_variables<<"\n";

    return true;
}

bool HybridSystem::check_reset(const Event& event, const DiscreteValuation& source, const DiscreteValuation& target) const
{
    Set<Identifier> source_variables=join(this->state_variables(source),this->auxiliary_variables(source));
    Set<Identifier> target_state_variables=this->state_variables(source);
    Set<Identifier> updated_variables;

    for(jump_const_iterator iter=this->_continuous_updates.begin(); iter!=this->_continuous_updates.end(); ++iter) {
        if(iter->evnts.contains(event) && evaluate(iter->loc,source)) {
            ARIADNE_ASSERT_MSG(subset(iter->rhs.arguments(),source_variables),"");
            updated_variables.insert(iter->lhs.name());
        }
    }

    std::cerr<<"source_variables="<<source_variables<<"\n";
    std::cerr<<"target_state_variables="<<target_state_variables<<"\n";
    std::cerr<<"updated_variables="<<updated_variables<<"\n";

    ARIADNE_ASSERT_MSG(subset(updated_variables,target_state_variables),"");
    ARIADNE_ASSERT_MSG(subset(target_state_variables,updated_variables),"");

    return true;
}


bool HybridSystem::check_guards(const DiscreteValuation& location) const
{
    Set<Identifier> variables=join(this->state_variables(location),this->auxiliary_variables(location));

    EventSet unguarded_events=EventSet::all();

    for(guard_const_iterator iter=this->_guard_predicates.begin(); iter!=this->_guard_predicates.end(); ++iter) {
        if(evaluate(iter->loc,location)) {
            if(iter->evnts.finite()) {
                ARIADNE_ASSERT_MSG(subset(iter->pred.arguments(),variables),"");
            } else {
                const Expression<tribool>& pred=iter->pred;
                if(!(pred==false)) { ARIADNE_ASSERT_MSG(true,"Nontrivial guard applied to an infinite set of events"); }
            }
            unguarded_events.remove(iter->evnts);
        }
    }
    ARIADNE_ASSERT_MSG(unguarded_events.size()==0,"Events "<<unguarded_events<<" have no active guard in location "<<location);
    return true;
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
    typedef Map<Identifier,Set<Identifier> >::iterator dependencies_iterator;

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
    for(Map<RealVariable,RealExpression>::const_iterator iter=formulae.begin(); iter!=formulae.end(); ++iter) {
        dependencies.insert(iter->first.name(),iter->second.arguments().restrict(variables));
    }

    List<Identifier> ordering=flatten(dependencies);

    List<RealAssignment> equations;
    for(List<Identifier>::const_iterator iter=ordering.begin(); iter!=ordering.end(); ++iter)
    {
        for(Map<RealVariable,RealExpression>::iterator fiter=formulae.begin(); fiter!=formulae.end(); ++fiter) {
            if(fiter->first.name()==*iter) { equations.push_back(fiter->first=fiter->second); }
        }
    }

    return equations;
}



List<RealDynamic>
HybridSystem::dynamic(const DiscreteValuation& location) const
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
HybridSystem::reset(const Event& event, const DiscreteValuation& old_location) const
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
    ARIADNE_ASSERT_MSG(unguarded_events.empty(),"Events "<<unguarded_events<<" in location "<<location<<" have no guard predicate");

    return guards;
}


std::ostream& operator<<(std::ostream& os, const HybridSystem& sys) {
    os << "HybridSystem(\n"
       << "  algebraic_equations=" << sys._algebraic_equations << ",\n"
       << "  differential_equations=" << sys._differential_equations << ",\n"
       << "  discrete_resets=" << sys._discrete_updates << ",\n"
       << "  continuous_resets=" << sys._continuous_updates << ",\n"
       << "  guard_predicates=" << sys._guard_predicates << ",\n"
       << "  invariant_predicates=" << sys._invariant_predicates << "\n)\n";
    return os;
}

std::ostream& operator<<(std::ostream& os, const HybridSystem::DifferentialEquation& de) {
    os << de.loc << " -> dot("<<de.lhs<<")="<<de.rhs;
    return os;
}

std::ostream& operator<<(std::ostream& os, const HybridSystem::AlgebraicEquation& ae) {
    os << ae.loc << " -> "<<ae.lhs<<"="<<ae.rhs;
    return os;
}

std::ostream& operator<<(std::ostream& os, const HybridSystem::DiscreteUpdate& da) {
    os << da.evnts << ": " << da.loc << " -> next("<<da.lhs.name()<<")="<<da.rhs;
    return os;
}

std::ostream& operator<<(std::ostream& os, const HybridSystem::ContinuousUpdate& da) {
    os << da.evnts << ": " << da.loc << " -> next("<<da.lhs.name()<<")="<<da.rhs;
    return os;
}

std::ostream& operator<<(std::ostream& os, const HybridSystem::GuardPredicate& g) {
    os << g.evnts << ": " << g.loc << " -> "<<g.pred;
    return os;
}

std::ostream& operator<<(std::ostream& os, const HybridSystem::InvariantPredicate& inv) {
    os << inv.loc << " -> "<<inv.pred;
    return os;
}





}
