/***************************************************************************
 *            hybrid_automaton-composite.cc
 *
 *  Copyright  2004-11  Alberto Casagrande, Pieter Collins
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
#include <boost/weak_ptr.hpp>

#include "config.h"

#include "macros.h"
#include "stlio.h"
#include "formula.h"
#include "expression.h"
#include "space.h"
#include "function.h"
#include "hybrid_time.h"
#include "hybrid_space.h"

#include "hybrid_automaton-composite.h"
#include "container.h"

namespace Ariadne {


typedef uint DimensionType;

DiscreteMode::DiscreteMode()
{
}

DiscreteMode::
DiscreteMode(DiscreteLocation location)
    : _location(location)
{
}

DiscreteMode::
DiscreteMode(DiscreteLocation location,
             List<RealAssignment> const& auxiliary,
             List<DottedRealAssignment> const& dynamic)
    : _location(location), _auxiliary(auxiliary), _dynamic(dynamic)
{
}



template<class K, class T1, class T2>
Map<K,Pair<T1,T2> > merge(const Map<K,T1>& m1, const Map<K,T2>& m2) {
    assert(m1.keys()==m2.keys());
    Map<K,Pair<T1,T2> > r;
    Set<K> keys=m1.keys();
    for(typename Set<K>::const_iterator iter=keys.begin(); iter!=keys.end(); ++iter) {
        r.insert(*iter,make_pair(m1[*iter],m2[*iter]));
    }
    return r;
}


template<class X> Set<Identifier> names(const Set< Variable<X> >& variables) {
    Set<Identifier> names;
    for(typename Set< Variable<X> >::const_iterator variable_iter=variables.begin(); variable_iter!=variables.end(); ++variable_iter) {
        names.insert(variable_iter->name());
    }
    return names;
}

template<class X> Set<Identifier> names(const List< Variable<X> >& variables) {
    Set<Identifier> names;
    for(typename List< Variable<X> >::const_iterator variable_iter=variables.begin(); variable_iter!=variables.end(); ++variable_iter) {
        names.insert(variable_iter->name());
    }
    return names;
}

std::ostream&
DiscreteMode::write(std::ostream& os) const
{
    const DiscreteMode& mode=*this;
    os << "DiscreteMode( "
       << "location=" << mode._location;
    if(mode._auxiliary.size()>0) {
        os << ", algebraic_equations="<<mode._auxiliary; }
    if(mode._dynamic.size()>0) {
        os << ", differential_equations="<<mode._dynamic; }
    if(mode._invariants.size()>0) {
        os << ", invariants="<<mode._invariants; }
    if(mode._guards.size()>0) {
        os << ", guards="<<mode._guards; }
    if(mode._targets.size()>0) {
        os << ", targets="<<mode._targets; }
    if(mode._targets.size()>0) {
        os << ", resets="<<mode._resets; }
    return os << " )";
}




DiscreteTransition::
DiscreteTransition(DiscreteLocation source,
                   DiscreteEvent event,
                   DiscreteLocation target,
                   const List<PrimedRealAssignment>& reset,
                   const ContinuousPredicate& guard,
                   EventKind kind)
    : _source(source), _event(event), _target(target),
      _guard(guard), _reset(reset)
{
}

std::ostream&
DiscreteTransition::write(std::ostream& os) const
{
    const DiscreteTransition& transition=*this;
    return os << "DiscreteTransition( "
              << "event=" << transition._event << ", "
              << "source=" << transition._source << ", "
              << "target=" << transition._target << ", "
              << "reset=" << transition._reset << ", "
              << "guard=" << transition._guard << " )";
}


Identifier name_composition(const List<HybridAutomaton>& components)
{
	List<HybridAutomaton>::const_iterator comp_it = components.begin();

	ARIADNE_ASSERT_MSG(comp_it != components.end(), "The components list is empty.");

	Identifier composed_name = comp_it->name();
	++comp_it;
	for (; comp_it != components.end(); ++comp_it)
		composed_name += "&" + comp_it->name();

	return composed_name;
}

EffectiveVectorFunction dynamic_function(
    Space<Real> const& space,
    List<RealAssignment> const& algebraic,
    List<DottedRealAssignment> const& differential)
{
    RealExpression default_expression;
    List<RealExpression> results(differential.size(),default_expression);
    for(uint i=0; i!=differential.size(); ++i) { results[space.index(differential[i].lhs.base())]=substitute(differential[i].rhs,algebraic); }

    return EffectiveVectorFunction(Ariadne::dimension(space),Ariadne::formula(results,algebraic,space));
}

EffectiveVectorFunction reset_function(
    Space<Real> const& space,
    List<RealAssignment> const& algebraic,
    List<PrimedRealAssignment> const& primed)
{
    RealExpression default_expression;
    List<RealExpression> results(primed.size(),default_expression);
    for(uint i=0; i!=primed.size(); ++i) { results[i]=substitute(primed[i].rhs,algebraic); }

    return EffectiveVectorFunction(Ariadne::dimension(space),Ariadne::formula(results,algebraic,space));
}

EffectiveScalarFunction constraint_function(
    Space<Real> const& space,
    List<RealAssignment> const& algebraic,
    ContinuousPredicate const& constraint,
    Sign sign)
{
    RealExpression constraint_expression=indicator(constraint,sign);
    return EffectiveScalarFunction(Ariadne::dimension(space),Ariadne::formula(constraint_expression,algebraic,space));
}

HybridAutomaton::HybridAutomaton()
    : _name("system"),_modes()
{
}


HybridAutomaton::~HybridAutomaton()
{
}

HybridAutomaton::HybridAutomaton(Identifier name)
    : _name(name),_modes()
{
}

HybridAutomaton::HybridAutomaton(
		Identifier name,
		const List<StringVariable>& discrete_variables)
    : _name(name),_modes()
{
}





void
HybridAutomaton::_new_mode(DiscreteLocation location,
                           const List<RealAssignment>& auxiliary,
                           const List<DottedRealAssignment>& dynamic)
{
    for(Map<DiscreteLocation,DiscreteMode>::const_iterator mode_iter=this->_modes.begin(); mode_iter!=this->_modes.end(); ++mode_iter) {
        if(!are_distinguishable(location,mode_iter->first)) {
            ARIADNE_THROW(IndistinguishableModeError,"HybridAutomaton::new_mode",
                          "Location "<<location<<" is indistinguishable from location "<<mode_iter->first<<
                          " in hybrid automaton with locations "<<this->locations());
        }
    }

    Set<UntypedVariable> defined_variables;
    Set<UntypedVariable> argument_variables;

    // Compute the auxiliary variables ordered by the given equations
    for(uint i=0; i!=auxiliary.size(); ++i) {
        if(defined_variables.contains(auxiliary[i].lhs)) {
            ARIADNE_THROW(SystemSpecificationError,"HybridAutomaton::new_mode",
                          "Variable "<<auxiliary[i].lhs<<" is defined twice by the algebraic equations "<<auxiliary<<" for mode "<<location);
        }
        defined_variables.insert(auxiliary[i].lhs);
        argument_variables.adjoin(auxiliary[i].rhs.arguments());
    }

    // Compute the state variables ordered by the given differential equations
    for(uint i=0; i!=dynamic.size(); ++i) {
        if(defined_variables.contains(dynamic[i].lhs.base())) {
            ARIADNE_THROW(SystemSpecificationError,"HybridAutomaton::new_mode",
                          "Variable "<<dynamic[i].lhs.base()<<" is defined by the differential equations "<<dynamic<<" for mode "<<location<<" is already defined");
        }
        defined_variables.insert(dynamic[i].lhs.base());
        argument_variables.adjoin(dynamic[i].rhs.arguments());
    }

    List<RealAssignment> sorted_auxiliary = sort(auxiliary);

    // TODO: Compute function
    DiscreteMode new_mode=DiscreteMode(location,sorted_auxiliary,dynamic);

    this->_modes.insert(location,new_mode);
}


void HybridAutomaton::_new_invariant_(DiscreteLocation location, ContinuousPredicate invariant, DiscreteEvent event)
{
    DiscreteMode& mode=this->_modes.value(location);
    mode._invariants.insert(event,invariant);
    mode._kinds.insert(event,INVARIANT);
}

void HybridAutomaton::_new_guard_(DiscreteLocation location, DiscreteEvent event, ContinuousPredicate guard, EventKind kind)
{
    DiscreteMode& mode=this->_modes.value(location);
    mode._guards.insert(event,guard);
    mode._kinds.insert(event,kind);
}


void
HybridAutomaton::_new_action(DiscreteLocation location,
                             ContinuousPredicate invariant,
                             DiscreteEvent event,
                             ContinuousPredicate guard,
                             EventKind kind)
{
    if(!this->has_location(location)) {
        ARIADNE_THROW(NonExistentModeError,"HybridAutomaton::new_action",
                      "mode "<<location<<" is not a location of the automaton with locations "<<this->locations());
    }
    DiscreteMode& mode=this->_modes.value(location);
    if(mode._kinds.has_key(event)) {
        ARIADNE_THROW(MultipleGuardError,"HybridAutomaton::new_action",
                      "Guard/invariant for event "<<event<<" is already defined in mode "<<mode);
    }
    mode._invariants.insert(event,invariant);
    mode._guards.insert(event,guard);
    mode._kinds.insert(event,kind);
}


void
HybridAutomaton::_new_update(DiscreteLocation source,
                             DiscreteEvent event,
                             DiscreteLocation target,
                             List<PrimedRealAssignment> const& reset)
{
    if(!this->has_location(source)) {
        ARIADNE_THROW(NonExistentModeError,"HybridAutomaton::new_reset",
                      "Source mode "<<source<<" is not a location of the automaton with locations "<<this->locations());
    }

    if(!this->has_location(target)) {
        ARIADNE_THROW(NonExistentModeError,"HybridAutomaton::new_reset",
                      "Target mode "<<target<<" is not a location of the automaton with locations "<<this->locations());
    }

    DiscreteMode& source_mode=this->mode(source); // Non-constant since we may wish to update input variables
    if(source_mode._targets.has_key(event)) {
        ARIADNE_THROW(MultipleTransitionError,"HybridAutomaton::new_update",
                      "Update for event "<<event<<" is already defined in mode "<<source);
    }

    Set<RealVariable> target_auxiliary_variable=make_set(this->auxiliary_variables(target));
    if(!disjoint(target_auxiliary_variable,left_hand_sides(reset))) {
        ARIADNE_THROW(OverspecifiedResetError,"HybridAutomaton::new_update",
                      "reset "<<reset<<" for event "<<event<<" in source "<<source<<
                      " specifies variables "<<intersection(target_auxiliary_variable,left_hand_sides(reset))<<
                      " which are auxiliary variables in target "<<target);
    }

    if(!unique_elements(left_hand_sides(reset))) {
        ARIADNE_THROW(OverspecifiedResetError,"HybridAutomaton::new_update",
                      "reset "<<reset<<" for event "<<event<<" in source "<<source<<
                      "overspecifies variables "<<duplicate_elements(left_hand_sides(reset)));
    }

    List<RealVariable> target_state_space=this->state_variables(target);
    Set<RealVariable> target_state_variables(target_state_space);
    Set<RealVariable> reset_variables(left_hand_sides(reset));


    source_mode._targets.insert(event,target);
    source_mode._resets.insert(event,reset);
}





Set<DiscreteLocation>
HybridAutomaton::locations() const
{
    return this->_modes.keys();
}




bool
HybridAutomaton::has_location(DiscreteLocation location) const
{
    return this->_modes.has_key(location);
}

bool
HybridAutomaton::has_mode(DiscreteLocation location) const
{
    for(Map<DiscreteLocation,DiscreteMode>::const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
    {
        if(is_restriction(mode_iter->first,location)) {
            return true;
        }
    }
    return false;
}

bool
HybridAutomaton::has_transition(DiscreteLocation source, DiscreteEvent event) const
{
   return this->has_mode(source) && this->mode(source)._targets.has_key(event);
}

bool
HybridAutomaton::has_invariant(DiscreteLocation source, DiscreteEvent event) const
{
   return this->has_mode(source) && this->mode(source)._invariants.has_key(event);
}


bool
HybridAutomaton::has_guard(DiscreteLocation source, DiscreteEvent event) const
{
   return this->has_invariant(source,event) || this->has_transition(source,event);
}



Map<DiscreteLocation,DiscreteMode> const&
HybridAutomaton::modes() const
{
    return this->_modes;
}

//Set<DiscreteMode>
//HybridAutomaton::modes() const
//{
//    return Set<DiscreteMode>(this->_modes.values());
//}

DiscreteMode&
HybridAutomaton::mode(DiscreteLocation location)
{
    const HybridAutomaton& self=*this;
    return const_cast<DiscreteMode&>(self.mode(location));

    if(!this->_modes.has_key(location)) {
        ARIADNE_THROW(SystemSpecificationError,"HybridAutomaton::mode(DiscreteLocation)",
                      location<<" is not a location of the automaton with locations "<<this->locations());
    }
    return this->_modes.find(location)->second;
}

const DiscreteMode&
HybridAutomaton::mode(DiscreteLocation location) const
{
    for(Map<DiscreteLocation,DiscreteMode>::const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
    {
        const DiscreteLocation& partial_location=mode_iter->first;
        if(is_restriction(partial_location,location)) {
            return mode_iter->second;
        }
    }

    ARIADNE_THROW(SystemSpecificationError,"HybridAutomaton::mode(DiscreteLocation)",
                  location<<" does not define a mode of the automaton with locations "<<this->locations());
}




std::ostream&
HybridAutomaton::write(std::ostream& os) const
{
    const HybridAutomaton& ha=*this;
    os << "\nHybridAutomaton( \n  modes=\n";
    Set<DiscreteMode> modes(ha.modes().values());
    for(Set<DiscreteMode>::const_iterator mode_iter=modes.begin();
            mode_iter!=modes.end(); ++mode_iter)
    {
        os << "    " <<*mode_iter<<",\n";
    }
    return os << ")\n";
}



List<RealVariable>
HybridAutomaton::state_variables(DiscreteLocation location) const {
    List<RealVariable> result;
    const DiscreteMode& mode=this->mode(location);
    for(uint i=0; i!=mode._dynamic.size(); ++i) {
        result.append(mode._dynamic[i].lhs.base());
    }
    return result;
}

List<RealVariable>
HybridAutomaton::auxiliary_variables(DiscreteLocation location) const {
    List<RealVariable> result;
    const DiscreteMode& mode=this->mode(location);
    for(uint i=0; i!=mode._auxiliary.size(); ++i) {
        result.append(mode._auxiliary[i].lhs);
    }
    return result;
}

List<RealAssignment>
HybridAutomaton::sort(const List<RealAssignment>& auxiliary) {
    List<RealVariable> lhs_array=left_hand_sides(auxiliary);
    LinkedList<RealVariable> lhs_list(lhs_array.begin(),lhs_array.end());
    Map<Identifier, Set<UntypedVariable> > dependencies;
    Set<UntypedVariable> variables(lhs_array.begin(),lhs_array.end());

    for(List<RealAssignment>::const_iterator asgn_iter=auxiliary.begin();
        asgn_iter!=auxiliary.end(); ++asgn_iter)
    {
        dependencies.insert(asgn_iter->lhs.base().name(), intersection(asgn_iter->rhs.arguments(),variables));
    }

    List<RealAssignment> sorted_auxiliary;
    sorted_auxiliary.reserve(auxiliary.size());
    if(!lhs_list.empty()) {
        bool found=false;
        for(LinkedList<RealVariable>::iterator iter=lhs_list.begin(); iter!=lhs_list.end(); ) {
            if(dependencies[iter->name()].empty()) {
                for(Map<Identifier, Set<UntypedVariable> >::iterator dep_iter=dependencies.begin(); dep_iter!=dependencies.end(); ++dep_iter) {
                    dep_iter->second.remove(static_cast<UntypedVariable>(*iter));
                }
                for(List<RealAssignment>::const_iterator asgn_iter=auxiliary.begin();
                    asgn_iter!=auxiliary.end(); ++asgn_iter)
                {
                    if(asgn_iter->lhs.name()==iter->name()) {
                        sorted_auxiliary.append(*asgn_iter);
                        break;
                    }
                }
                dependencies.erase(iter->name());
                found=true;
                LinkedList<RealVariable>::iterator next=iter;
                ++next;
                lhs_list.erase(iter);
                iter=next;
            } else {
                ++iter;
            }
        }
        if(!found) {
            ARIADNE_THROW(AlgebraicLoopError,"HybridAutomaton::sort(List<RealAssignment>)",
                          "Algebraic dependencies among variables "<<lhs_list<<" in auxiliary equations "<<auxiliary);
        }
    }
    return sorted_auxiliary;
}

Set<RealVariable>
HybridAutomaton::input_variables(DiscreteLocation location) const {
    Set<UntypedVariable> used;
    Set<RealVariable> defined;
    const DiscreteMode& mode=this->mode(location);
    for(List<RealAssignment>::const_iterator iter=mode._auxiliary.begin(); iter!=mode._auxiliary.end(); ++iter) {
        defined.insert(iter->lhs.base());
        used.adjoin(iter->rhs.arguments());
    }
    for(List<DottedRealAssignment>::const_iterator iter=mode._dynamic.begin(); iter!=mode._dynamic.end(); ++iter) {
        defined.insert(iter->lhs.base());
        used.adjoin(iter->rhs.arguments());
    }
    for(Map<DiscreteEvent,ContinuousPredicate>::const_iterator iter=mode._invariants.begin(); iter!=mode._invariants.end(); ++iter) {
        used.adjoin(iter->second.arguments());
    }
    for(Map<DiscreteEvent,ContinuousPredicate>::const_iterator iter=mode._guards.begin(); iter!=mode._guards.end(); ++iter) {
        used.adjoin(iter->second.arguments());
    }
    for(Map<DiscreteEvent,List<PrimedRealAssignment> >::const_iterator iter=mode._resets.begin(); iter!=mode._resets.end(); ++iter) {
        const List<PrimedRealAssignment>& reset=iter->second;
        for(List<PrimedRealAssignment>::const_iterator iter=reset.begin(); iter!=reset.end(); ++iter) {
            used.adjoin(iter->rhs.arguments());
        }
    }
    Set<RealVariable> result;
    for(Set<UntypedVariable>::const_iterator iter=used.begin(); iter!=used.end(); ++iter) {
        RealVariable var(iter->name());
        if(!defined.contains(var)) {
            result.insert(var);
        }
    }
    return result;
}


DiscreteLocation
HybridAutomaton::target(DiscreteLocation source, DiscreteEvent event) const {
    if(this->has_transition(source,event)) {
        return this->mode(source)._targets[event];
    } else {
        return source;
    }
}

List<RealAssignment>
HybridAutomaton::auxiliary_assignments(DiscreteLocation location) const {
    return this->mode(location)._auxiliary;
}

List<DottedRealAssignment>
HybridAutomaton::dynamic_assignments(DiscreteLocation location) const {
    return this->mode(location)._dynamic;
}

List<PrimedRealAssignment>
HybridAutomaton::reset_assignments(DiscreteLocation source, DiscreteEvent event) const {
    if(this->has_transition(source,event)) {
        return this->mode(source)._resets[event];
    } else {
        List<PrimedRealAssignment> nonjump_assignments;
        List<RealVariable> state_variables=this->state_variables(source);
        for(uint i=0; i!=state_variables.size(); ++i) {
            nonjump_assignments.append(next(state_variables[i])=state_variables[i]);
        }
        return nonjump_assignments;
    }
}

ContinuousPredicate
HybridAutomaton::invariant_predicate(DiscreteLocation location, DiscreteEvent action) const {
    const DiscreteMode& mode=this->mode(location);
    if(mode._invariants.has_key(action)) {
        return mode._invariants[action];
    } else {
        return ContinuousPredicate(true);
    }
}

ContinuousPredicate
HybridAutomaton::guard_predicate(DiscreteLocation location, DiscreteEvent event) const {
    const DiscreteMode& mode=this->mode(location);
    if(mode._guards.has_key(event)) {
        return mode._guards[event];
    } else {
        return ContinuousPredicate(true);
    }
}


HybridSpace HybridAutomaton::state_space() const {
    MonolithicHybridSpace space;
    for(Map<DiscreteLocation,DiscreteMode>::const_iterator iter=this->_modes.begin(); iter!=this->_modes.end(); ++iter) {
        const DiscreteLocation& loc=iter->first;
        space.new_location(loc,this->continuous_state_space(loc));
    }
    return space;
}

Set<DiscreteEvent> HybridAutomaton::events(DiscreteLocation location) const {
    const DiscreteMode& mode=this->mode(location);
    return join(join(mode._invariants.keys(),mode._guards.keys()),mode._targets.keys());
}

uint HybridAutomaton::dimension(DiscreteLocation location) const {
    return this->mode(location)._dynamic.size();
};

RealSpace HybridAutomaton::continuous_state_space(DiscreteLocation location) const {
    return RealSpace(left_hand_sides(this->mode(location)._dynamic));
};

EventKind HybridAutomaton::event_kind(DiscreteLocation location, DiscreteEvent event) const {
    return this->mode(location)._kinds[event];
}

EffectiveVectorFunction HybridAutomaton::auxiliary_function(DiscreteLocation location) const {
    ARIADNE_NOT_IMPLEMENTED;
}

EffectiveVectorFunction HybridAutomaton::dynamic_function(DiscreteLocation location) const {
    DiscreteMode const& mode=this->mode(location);
    return Ariadne::dynamic_function(this->continuous_state_space(location),mode._auxiliary,mode._dynamic);
}

EffectiveVectorFunction HybridAutomaton::reset_function(DiscreteLocation location, DiscreteEvent event) const {
    DiscreteMode const& mode=this->mode(location);
    return Ariadne::reset_function(this->continuous_state_space(location),mode._auxiliary,mode._resets[event]);
}

EffectiveScalarFunction HybridAutomaton::invariant_function(DiscreteLocation location, DiscreteEvent event) const {
    DiscreteMode const& mode=this->mode(location);
    return Ariadne::constraint_function(this->continuous_state_space(location),mode._auxiliary,mode._invariants[event],NEGATIVE);
}

EffectiveScalarFunction HybridAutomaton::guard_function(DiscreteLocation location, DiscreteEvent event) const {
    DiscreteMode const& mode=this->mode(location);
    return Ariadne::constraint_function(this->continuous_state_space(location),mode._auxiliary,mode._guards[event],POSITIVE);
}


void HybridAutomaton::check_mode(DiscreteLocation location) const {
    const DiscreteMode& mode=this->mode(location);

    List<RealVariable> defined_real_variables(catenate(this->state_variables(location),this->auxiliary_variables(location)));
    Set<UntypedVariable> defined_variables(make_set(defined_real_variables));

    for(List<RealAssignment>::const_iterator aux_iter=mode._auxiliary.begin(); aux_iter!=mode._auxiliary.end(); ++aux_iter) {
        if(!subset(aux_iter->rhs.arguments(),defined_variables)) {
            ARIADNE_THROW(UnderspecifiedDynamicError,"HybridAutomaton::check_mode(...)",
                          "Arguments "<<difference(aux_iter->rhs.arguments(),defined_variables)<<
                          " of "<<*aux_iter<<" are not defined in location "<<location<<" with variables "<<defined_variables);
        }
    }

    for(List<DottedRealAssignment>::const_iterator dyn_iter=mode._dynamic.begin(); dyn_iter!=mode._dynamic.end(); ++dyn_iter) {
        if(!subset(dyn_iter->rhs.arguments(),defined_variables)) {
            ARIADNE_THROW(UnderspecifiedDynamicError,"HybridAutomaton::check_mode(...)",
                          "Arguments "<<difference(dyn_iter->rhs.arguments(),defined_variables)<<
                          " of "<<*dyn_iter<<" are not defined in location "<<location<<" with variables "<<defined_variables);
        }
    }

    for(Map<DiscreteEvent,ContinuousPredicate>::const_iterator inv_iter=mode._invariants.begin(); inv_iter!=mode._invariants.end(); ++inv_iter) {
        if(!subset(inv_iter->second.arguments(),defined_variables)) {
            ARIADNE_THROW(UnderspecifiedConstraintError,"HybridAutomaton::check_mode(...)",
                          "Arguments "<<difference(inv_iter->second.arguments(),defined_variables)<<
                          " of "<<*inv_iter<<" are not defined in location "<<location<<" with variables "<<defined_variables);
        }
    }

    for(Map<DiscreteEvent,ContinuousPredicate>::const_iterator grd_iter=mode._guards.begin(); grd_iter!=mode._guards.end(); ++grd_iter) {
        if(!subset(grd_iter->second.arguments(),defined_variables)) {
            ARIADNE_THROW(UnderspecifiedConstraintError,"HybridAutomaton::check_mode(...)",
                          "Arguments "<<difference(grd_iter->second.arguments(),defined_variables)<<
                          " of "<<*grd_iter<<" are not defined in location "<<location<<" with variables "<<defined_variables);
        }
    }

    for(Map<DiscreteEvent,List<PrimedRealAssignment> >::const_iterator rst_iter=mode._resets.begin(); rst_iter!=mode._resets.end(); ++rst_iter) {
        const List<PrimedRealAssignment>& reset=rst_iter->second;
        for(List<PrimedRealAssignment>::const_iterator pra_iter=reset.begin(); pra_iter!=reset.end(); ++pra_iter) {
            if(!subset(pra_iter->rhs.arguments(),defined_variables)) {
                ARIADNE_THROW(UnderspecifiedResetError,"HybridAutomaton::check_mode(...)",
                            "Arguments "<<difference(pra_iter->rhs.arguments(),defined_variables)<<
                            " of reset"<<*pra_iter<<" are not defined in location "<<location<<" with variables "<<defined_variables);
            }
        }
    }

    for(Map<DiscreteEvent,List<PrimedRealAssignment> >::const_iterator rst_iter=mode._resets.begin(); rst_iter!=mode._resets.end(); ++rst_iter) {
        DiscreteLocation target = mode._targets[rst_iter->first];
        Set<RealVariable> reset_variables(left_hand_sides(rst_iter->second));
        Set<RealVariable> target_state_variables(this->state_variables(target));
        if(!subset(reset_variables,target_state_variables)) {
             ARIADNE_THROW(UnderspecifiedDynamicError,"HybridAutomaton::check_mode(...)",
                            "Variables "<<difference(reset_variables,target_state_variables)<<
                            " of reset"<<*rst_iter<<" are not defined in location "<<target<<" with state variables "<<target_state_variables);
        }
        if(!subset(target_state_variables,reset_variables)) {
             ARIADNE_THROW(UnderspecifiedResetError,"HybridAutomaton::check_mode(...)",
                            "Variables "<<difference(target_state_variables,reset_variables)<<
                            " in target"<<target<<" are not defined reset "<<*rst_iter);
        }
    }


}



class CompositeHybridSpace
    : public HybridSpaceInterface
{
  public:
    ~CompositeHybridSpace() { _system_ptr = 0; }
    CompositeHybridSpace(const CompositeHybridAutomaton& ha) : _system_ptr(&ha) { }
    virtual CompositeHybridSpace* clone() const { return new CompositeHybridSpace(*this); }
    virtual bool has_location(const DiscreteLocation& q) const { return this->_system_ptr->has_mode(q); }
    virtual RealSpace operator[](const DiscreteLocation& q) const { return this->_system_ptr->continuous_state_space(q); }
    virtual std::ostream& write(std::ostream& os) const { return os << "CompositeHybridSpace( " << *this->_system_ptr << " )"; }
    tribool operator==(const HybridSpaceInterface& other) const {
        const CompositeHybridSpace* chs_ptr = dynamic_cast<const CompositeHybridSpace* >(&other);
        if (!chs_ptr) return indeterminate;
        if (&*chs_ptr->_system_ptr == &*_system_ptr)
            return true;
        else return indeterminate;
    }
  private:
    const CompositeHybridAutomaton* _system_ptr;
};


CompositeHybridAutomaton::CompositeHybridAutomaton()
    : _name("system"),_components() { }

CompositeHybridAutomaton::CompositeHybridAutomaton(Identifier name)
    : _name(name),_components() { }

CompositeHybridAutomaton::CompositeHybridAutomaton(const HybridAutomaton& automaton)
    : _name(automaton.name()),_components(1u,automaton) { }

CompositeHybridAutomaton::CompositeHybridAutomaton(const List<HybridAutomaton>& components)
    : _name(name_composition(components)),_components(components) { }

CompositeHybridAutomaton::CompositeHybridAutomaton(
		Identifier name,
		const List<HybridAutomaton>& components)
    : _name(name),_components(components) { }

CompositeHybridAutomaton::~CompositeHybridAutomaton() { }

uint
CompositeHybridAutomaton::number_of_components() const
{
    return this->_components.size();
}

const HybridAutomaton&
CompositeHybridAutomaton::component(uint k) const {
    return this->_components[k];
}

namespace {

List<DottedRealVariable> dot(const List<RealVariable>& v) {
    List<DottedRealVariable> result;
    for(uint i=0; i!=v.size(); ++i) { result.append(dot(v[i])); }
    return result;
}

List<PrimedRealVariable> next(const List<RealVariable>& v) {
    List<PrimedRealVariable> result;
    for(uint i=0; i!=v.size(); ++i) { result.append(next(v[i])); }
    return result;
}

}

typedef Map<DiscreteLocation,DiscreteMode>::const_iterator ModesConstIterator;

bool has_mode(const HybridAutomaton& automaton, const DiscreteLocation& location) {
    const Map<DiscreteLocation,DiscreteMode>& modes=automaton.modes();
    for(ModesConstIterator iter=modes.begin(); iter!=modes.end(); ++iter) {
        const DiscreteMode& mode=iter->second;
        if(is_restriction(mode.location(),location)) {
            return true;
        }
    }
    return false;
}


bool
CompositeHybridAutomaton::has_mode(DiscreteLocation location) const {
    for(uint i=0; i!=this->_components.size(); ++i) {
        if(!this->_components[i].has_mode(location)) {
            return false;
        }
    }
    return true;
}

bool
CompositeHybridAutomaton::has_invariant(DiscreteLocation location, DiscreteEvent event) const
{
    for(uint i=0; i!=this->_components.size(); ++i) {
        if(this->_components[i].has_invariant(location,event)) {
            return true;
        }
    }
    return false;
}

bool
CompositeHybridAutomaton::has_guard(DiscreteLocation location, DiscreteEvent event) const
{
    for(uint i=0; i!=this->_components.size(); ++i) {
        if(this->_components[i].has_guard(location,event)) {
            return true;
        }
    }
    return false;
}

bool
CompositeHybridAutomaton::has_transition(DiscreteLocation source, DiscreteEvent event) const {
    for(uint i=0; i!=this->_components.size(); ++i) {
        if(this->_components[i].has_transition(source,event)) {
            return true;
        }
    }
    return false;
}


DiscreteMode const& CompositeHybridAutomaton::mode(DiscreteLocation location) const
{
    if(this->_cached_mode.location().values().keys()!=location.values().keys() || !(this->_cached_mode._location==location)) {
        this->_cache_mode(location);
    }
    return this->_cached_mode;
}

void CompositeHybridAutomaton::_cache_mode(DiscreteLocation location) const
{
    DiscreteMode& result=this->_cached_mode;
    result=DiscreteMode(location);

    for(List<HybridAutomaton>::const_iterator component_iter=this->_components.begin();
        component_iter!=this->_components.end(); ++component_iter)
    {
        const DiscreteMode& component_mode=component_iter->mode(location);
        result._auxiliary.append(component_mode._auxiliary);
        result._dynamic.append(component_mode._dynamic);

        for(Map<DiscreteEvent,ContinuousPredicate>::const_iterator invariant_iter=component_mode._invariants.begin();
            invariant_iter!=component_mode._invariants.end(); ++invariant_iter)
        {
            result._invariants.insert(*invariant_iter);
        }

        for(Map<DiscreteEvent,ContinuousPredicate>::const_iterator guard_iter=component_mode._guards.begin();
            guard_iter!=component_mode._guards.end(); ++guard_iter)
        {
            result._guards.insert(*guard_iter);
        }

        for(Map<DiscreteEvent,EventKind>::const_iterator kind_iter=component_mode._kinds.begin();
            kind_iter!=component_mode._kinds.end(); ++kind_iter)
        {
            result._kinds.insert(*kind_iter);
        }

        // FIXME: Need to introduce missing updates
        for(Map<DiscreteEvent,DiscreteLocation>::const_iterator target_iter=component_mode._targets.begin();
            target_iter!=component_mode._targets.end(); ++target_iter)
        {
            DiscreteLocation& target=result._targets[target_iter->first];
            target=(target,target_iter->second);
        }

        for(Map<DiscreteEvent,List<PrimedRealAssignment> >::const_iterator reset_iter=component_mode._resets.begin();
            reset_iter!=component_mode._resets.end(); ++reset_iter)
        {
            List<PrimedRealAssignment>& reset=result._resets[reset_iter->first];
            reset.append(reset_iter->second);
        }
    }
}


Set<DiscreteEvent>
CompositeHybridAutomaton::events(DiscreteLocation location) const
{
    DiscreteMode const& mode=this->mode(location);
    return join(join(mode._invariants.keys(),mode._guards.keys()),mode._targets.keys());
}

EventKind
CompositeHybridAutomaton::event_kind(DiscreteLocation location, DiscreteEvent event) const
{
    return this->mode(location)._kinds[event];
}



DiscreteLocation
CompositeHybridAutomaton::target(DiscreteLocation source, DiscreteEvent event) const {
    DiscreteLocation result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        DiscreteMode const& component_mode = this->_components[i].mode(source);
        if(component_mode._targets.has_key(event)) {
            result.adjoin(component_mode._targets[event]);
        } else {
            result.adjoin(component_mode._location);
        }
    }
    return result;
}


uint
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
    for(uint i=0; i!=this->_components.size(); ++i) {
        ARIADNE_ASSERT(this->_components[i].has_mode(location));
        // Append state space in mode i; note that this automatically checks whether a variable is already present
        result.append(this->_components[i].state_variables(location));
    }
    return result;
}

List<RealVariable>
CompositeHybridAutomaton::auxiliary_variables(DiscreteLocation location) const {
    List<RealVariable> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        ARIADNE_ASSERT(this->_components[i].has_mode(location));
        // Append state space in mode i; note that this automatically checks whether a variable is already present
        result.append(this->_components[i].auxiliary_variables(location));
    }
    return result;
}

// Find all algebraic equations valid in the location
List<RealAssignment>
CompositeHybridAutomaton::auxiliary_assignments(DiscreteLocation location) const {
    List<RealAssignment> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result.append(this->_components[i].auxiliary_assignments(location));
    }
    // TODO: Sort the result to eliminate algebraic loops
    // sort(result);
    return result;
}

List<DottedRealAssignment>
CompositeHybridAutomaton::dynamic_assignments(DiscreteLocation location) const {
    List<DottedRealAssignment> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result.append(this->_components[i].dynamic_assignments(location));
    }
    // No need to sort result since dotted variables cannot appear in the right-hand side (currently)
    return result;
}


List<PrimedRealAssignment>
CompositeHybridAutomaton::reset_assignments(DiscreteLocation location, DiscreteEvent event) const {
    List<PrimedRealAssignment> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result.append(this->_components[i].reset_assignments(location,event));
    }
    return result;
}


ContinuousPredicate
CompositeHybridAutomaton::invariant_predicate(DiscreteLocation location, DiscreteEvent event) const {
    ContinuousPredicate result(true);
    for(uint i=0; i!=this->_components.size(); ++i) {
        ContinuousPredicate guard=this->_components[i].invariant_predicate(location,event);
        if(is_constant(result,true) || is_constant(guard,false)) { result=guard; }
        else if(is_constant(guard,true) || is_constant(result,false)) { }
        else { result = result && guard; }
    }
    return result;
}


ContinuousPredicate
CompositeHybridAutomaton::guard_predicate(DiscreteLocation location, DiscreteEvent event) const {
    ContinuousPredicate result(true);
    for(uint i=0; i!=this->_components.size(); ++i) {
        ContinuousPredicate guard=this->_components[i].guard_predicate(location,event);
        if(is_constant(result,true) || is_constant(guard,false)) { result=guard; }
        else if(is_constant(guard,true) || is_constant(result,false)) { }
        else { result = result && guard; }
    }
    return result;
}





EffectiveVectorFunction
CompositeHybridAutomaton::auxiliary_function(DiscreteLocation location) const {
    RealExpression default_expression;
    Space<Real> space=this->state_variables(location);
    List<RealAssignment> algebraic=this->auxiliary_assignments(location);
    List<RealExpression> results(algebraic.size(),default_expression);
    for(uint i=0; i!=algebraic.size(); ++i) { results[i]=algebraic[i].lhs; }
    return EffectiveVectorFunction(Ariadne::dimension(space),Ariadne::formula(results,algebraic,space));
}

EffectiveVectorFunction
CompositeHybridAutomaton::dynamic_function(DiscreteLocation location) const {
    RealExpression default_expression;
    Space<Real> space=this->state_variables(location);
    List<RealAssignment> algebraic=this->auxiliary_assignments(location);
    List<DottedRealAssignment> differential=this->dynamic_assignments(location);
    List<RealExpression> results(differential.size(),default_expression);
    for(uint i=0; i!=differential.size(); ++i) { results[space.index(differential[i].lhs.base())]=substitute(differential[i].rhs,algebraic); }

    return EffectiveVectorFunction(Ariadne::dimension(space),Ariadne::formula(results,algebraic,space));
}

EffectiveVectorFunction
CompositeHybridAutomaton::reset_function(DiscreteLocation source, DiscreteEvent event) const {
    RealExpression default_expression;
    DiscreteLocation target=this->target(source,event);
    Space<Real> source_space=this->state_variables(source);
    Space<Real> target_space=this->state_variables(target);
    List<RealAssignment> algebraic=this->auxiliary_assignments(source);
    List<PrimedRealAssignment> update=this->reset_assignments(source,event);
    List<RealExpression> results(update.size(),default_expression);
    for(uint i=0; i!=update.size(); ++i) { results[target_space.index(update[i].lhs.base())]=update[i].rhs; }
    return EffectiveVectorFunction(Ariadne::dimension(source_space),Ariadne::formula(results,algebraic,source_space));
}

EffectiveScalarFunction
CompositeHybridAutomaton::invariant_function(DiscreteLocation location, DiscreteEvent event) const {
    Space<Real> space=this->state_variables(location);
    List<RealAssignment> algebraic=this->auxiliary_assignments(location);
    RealExpression invariant=indicator(invariant_predicate(location,event),NEGATIVE);
    return EffectiveScalarFunction(Ariadne::dimension(space),Ariadne::formula(invariant,algebraic,space));
}

EffectiveScalarFunction
CompositeHybridAutomaton::guard_function(DiscreteLocation location, DiscreteEvent event) const {
    Space<Real> space=this->state_variables(location);
    List<RealAssignment> algebraic=this->auxiliary_assignments(location);
    RealExpression guard=indicator(guard_predicate(location,event),POSITIVE);
    return EffectiveScalarFunction(Ariadne::dimension(space),Ariadne::formula(guard,algebraic,space));
}




std::ostream&
CompositeHybridAutomaton::write(std::ostream& os) const
{
    return os << "CompositeHybridAutomaton(\n" << this->_components << "\n)\n";
}







template<class LHS, class RHS>
List<LHS> assigned(const List< Assignment<LHS,RHS> >& assignments) {
    List<LHS> result; result.reserve(assignments.size());
    for(uint i=0; i!=assignments.size(); ++i) {
        result.append(assignments[i].lhs);
    }
    return result;
}

template<class T>
Set<RealVariable> real_arguments(const Expression<T>& expression) {
    Set<RealVariable> result;
    Set<UntypedVariable> args=expression.arguments();
    for(Set<UntypedVariable>::const_iterator iter=args.begin(); iter!=args.end(); ++iter) {
        if(iter->type()==variable_type<Real>()) {
            result.insert(RealVariable(iter->name()));
        }
    }
    return result;
}

template<class LHS, class RHS>
Set<RealVariable> real_arguments(const Assignment<LHS,RHS>& assignment) {
    return real_arguments(assignment.rhs);
}

template<class LHS, class RHS>
Set<RealVariable> real_arguments(const List< Assignment<LHS,RHS> >& assignments) {
    Set<RealVariable> result;
    for(uint i=0; i!=assignments.size(); ++i) {
        result.adjoin(real_arguments(assignments[i]));
    }
    return result;
}

template<class Tp, template<class> class Dec>
List< Variable<Tp> > base(const List< Dec<Tp> >& variables) {
    List< Variable<Tp> > result;
    for(uint i=0; i!=variables.size(); ++i) {
        result.append(variables[i].base());
    }
    return result;
}

template<class T>
Set<T> make_set(const List<T>& list) {
    return Set<T>(list.begin(),list.end());
}

template<class Var>
bool unique(const List<Var>& variables) {
    Set<Var> found;
    for(uint i=0; i!=variables.size(); ++i) {
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
    for(uint i=0; i!=variables.size(); ++i) {
        if(found.contains(variables[i])) {
            result.insert(variables[i]);
        } else {
            found.insert(variables[i]);
        }
    }
    return result;
}

namespace {

// Order the assignments so that each variable depends only on previous ones, assuming the given input variables.
List<RealAssignment> order(const List<RealAssignment>& assignments, const Set<RealVariable>& inputs)
{
    List<RealAssignment> result;

    Map< RealVariable, Set<RealVariable> > unresolved_dependencies;
    Map< RealVariable, uint > assignment_table;

    for(uint i=0; i!=assignments.size(); ++i) {
        unresolved_dependencies.insert(assignments[i].lhs,difference(real_arguments(assignments[i].rhs),inputs));
        assignment_table.insert(assignments[i].lhs,i);
    }

    while(!unresolved_dependencies.empty()) {
        // Check for an which has been resolved and add to list of ok equations
        bool found=false;
        for(Map< RealVariable, Set<RealVariable> >::iterator iter=unresolved_dependencies.begin();
            iter!=unresolved_dependencies.end(); ++iter)
        {
            if(iter->second.empty()) {
                RealVariable variable=iter->first;
                result.append(assignments[assignment_table[variable]]);
                for(Map< RealVariable, Set<RealVariable> >::iterator rem_iter=unresolved_dependencies.begin();
                        rem_iter!=unresolved_dependencies.end(); ++rem_iter) {
                    rem_iter->second.erase(variable);
                }
                unresolved_dependencies.erase(iter);
                found=true;
                break;
            }
        }

        if(!found) {
            ARIADNE_THROW(AlgebraicLoopError,"order(List<Assignment>,Set<Variable>)","Cannot order assignments "<<assignments
                            <<" with inputs "<<inputs<<" due to algebraic loops or undefined variables.");
        }
    }

    assert(result.size()==assignments.size());
    return result;
}

}


void
CompositeHybridAutomaton::check_mode(DiscreteLocation location) const {
    List<RealAssignment> equations=this->auxiliary_assignments(location);
    List<DottedRealAssignment> dynamics=this->dynamic_assignments(location);

    List<RealVariable> state_variables=base(assigned(dynamics));
    List<RealVariable> auxiliary_variables=assigned(equations);
    List<RealVariable> location_variables=catenate(state_variables,auxiliary_variables);

    if(!unique(location_variables)) {
        ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                      "Variables "<<duplicates(location_variables)<<" in location "<<location<<" with defining equations "<<equations<<" and "<<dynamics<<" have more than one defining equation.");
    }

    Set<RealVariable> result_variables=make_set(location_variables);

    Set<RealVariable> undefined_variables=difference(real_arguments(equations),result_variables);
    if(!undefined_variables.empty()) {
        ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                      "Variables "<<undefined_variables<<" are used in the definition of the algebraic equations "<<equations<<" in location "<<location<<" with state variables "<<state_variables<<", but are not defined.");
    }

    undefined_variables=difference(real_arguments(dynamics),result_variables);
    if(!undefined_variables.empty()) {
        ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                      "Variables "<<undefined_variables<<" are used in the definition of the differential equations "<<dynamics<<" in location "<<location<<" with state variables "<<state_variables<<", but are not defined.");
    }

    List<RealAssignment> ordered_equations=Ariadne::order(equations,make_set(state_variables));


/*
    Map<DiscreteEvent,ContinuousPredicate> const& invariants=this->invariant_predicates(location);
    for(Map<DiscreteEvent,ContinuousPredicate>::const_iterator invariant_iter=invariants.begin();
        invariant_iter!=invariants.end(); ++invariant_iter)
    {
        DiscreteEvent action=invariant_iter->first;
        const ContinuousPredicate& invariant=invariant_iter->second;
        if(!subset(real_arguments(invariant),result_variables)) {
            undefined_variables=difference(real_arguments(invariant),result_variables);
            ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                          "Variables "<<undefined_variables<<" are used in the invariant "<<invariant<<" with label \""<<action<<"\" in location "<<location<<" with variables "<<location_variables<<", but are not defined.");
        }
    }

    Set<DiscreteEvent> events=this->events(location);
    for(Set<DiscreteEvent>::const_iterator event_iter=events.begin(); event_iter!=events.end(); ++event_iter) {

        // Get transition information
        DiscreteEvent event=*event_iter;
        DiscreteLocation target=this->target(location,event);
        ContinuousPredicate guard=this->guard_predicate(location,event);
        List<PrimedRealAssignment> reset=this->reset_assignments(location,event);

        Set<RealVariable> target_state_variables=make_set(this->state_variables(target));
        List<RealVariable> reset_variables=base(assigned(reset));

        // Check arguments of guard
        if(!subset(real_arguments(guard),result_variables)) {
            undefined_variables=difference(real_arguments(guard),result_variables);
            ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                          "Variables "<<undefined_variables<<" are used in the guard "<<guard<<" for event \""<<event<<"\" in location "<<location<<" with variables "<<location_variables<<", but are not defined.");
        }

        // Check arguments of reset
        if(!subset(real_arguments(reset),result_variables)) {
            undefined_variables=difference(real_arguments(reset),result_variables);
            ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                          "Variables "<<undefined_variables<<" are used in the reset "<<reset<<" for event \""<<event<<"\" in location "<<location<<" with variables "<<location_variables<<", but are not defined.");
        }

        // Check duplication of reset
        if(!unique(reset_variables)) {
        ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                      "Variables "<<duplicates(reset_variables)<<" have more than one defining equation in resets "<<reset<<" for event \""<<event<<"\" in location "<<location<<".");
        }

        Set<RealVariable> updated_variables=make_set(reset_variables);

        if(!subset(updated_variables,target_state_variables)) {
            Set<RealVariable> extra_variables=difference(updated_variables,target_state_variables);
            ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                          "Variables "<<extra_variables<<" are defined in the reset "<<reset<<" for event \""<<event<<"\" in location "<<location<<", but are not state variables "<<target_state_variables<<" in the target "<<target<<".");
        }

        if(!subset(target_state_variables,updated_variables)) {
            undefined_variables=difference(target_state_variables,updated_variables);
            ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                          "Variables "<<undefined_variables<<" are state variables in location "<<target<<", but are not defined in the reset "<<reset<<" for the transition \""<<event<<"\" from location "<<location<<".");
        }

    } // Finished checking transitions
*/
    return;
}


void
CompositeHybridAutomaton::check_reachable_modes(DiscreteLocation initial_location) const
{
    Set<DiscreteLocation> initial_locations;
    initial_locations.insert(initial_location);
    this->check_reachable_modes(initial_locations);
}


void
CompositeHybridAutomaton::check_reachable_modes(const Set<DiscreteLocation>& initial_locations) const
{
    Set<DiscreteLocation> reachable_locations=this->discrete_reachability(initial_locations);
    for(Set<DiscreteLocation>::const_iterator iter=reachable_locations.begin(); iter!=reachable_locations.end(); ++iter) {
        this->check_mode(*iter);
    }
}

Set<DiscreteLocation>
CompositeHybridAutomaton::discrete_reachability(DiscreteLocation initial_location) const
{
    Set<DiscreteLocation> initial_locations;
    initial_locations.insert(initial_location);
    return this->discrete_reachability(initial_locations);
}

Set<DiscreteLocation>
CompositeHybridAutomaton::discrete_reachability(const Set<DiscreteLocation>& initial_locations) const
{
    const CompositeHybridAutomaton& automaton=*this;

    typedef uint Nat;

    Set<DiscreteLocation> reached=initial_locations;
    Set<DiscreteLocation> working=initial_locations;
    Set<DiscreteLocation> found;

    Map<DiscreteLocation,Nat> steps;
    Nat step=0u;

    for(Set<DiscreteLocation>::const_iterator initial_iter=initial_locations.begin(); initial_iter!=initial_locations.end(); ++initial_iter) {
        steps.insert(*initial_iter,step);
    }

    while(!working.empty()) {
        ++step;
        for(Set<DiscreteLocation>::const_iterator source_iter=working.begin(); source_iter!=working.end(); ++source_iter) {
            ARIADNE_LOG(5,"new_mode\n");
            DiscreteLocation location=*source_iter;
            ARIADNE_LOG(5,"  mode: "<<location<<":\n");
            ARIADNE_LOG(7,"      auxiliary="<<automaton.auxiliary_assignments(location)<<"\n");
            ARIADNE_LOG(7,"          function="<<automaton.auxiliary_function(location)<<"\n");
            ARIADNE_LOG(7,"      dynamic="<<automaton.dynamic_assignments(location)<<"\n");
            ARIADNE_LOG(7,"          function="<<automaton.dynamic_function(location)<<"\n");

/*
            Map<DiscreteEvent,ContinuousPredicate> invariants=automaton.invariant_predicates(location);
            for(Map<DiscreteEvent,ContinuousPredicate>::const_iterator invariant_iter=invariants.begin();
                    invariant_iter!=invariants.end(); ++invariant_iter) {
                ARIADNE_LOG(7,"    invariant: "<<invariant_iter->second<<"\n");
                ARIADNE_LOG(7,"      function: "<<automaton.invariant_function(location,invariant_iter->first)<<"\n");
                ARIADNE_LOG(7,"        action: "<<invariant_iter->first<<"\n");
            }
*/
            Set<DiscreteEvent> events=automaton.events(location);
            ARIADNE_LOG(5,"  events: "<<events<<"\n");
            for(Set<DiscreteEvent>::const_iterator event_iter=events.begin(); event_iter!=events.end(); ++event_iter) {
                DiscreteEvent event=*event_iter;
                ARIADNE_LOG(7,"    event:"<<event<<"\n");
                DiscreteLocation target=automaton.target(location,event);
                ARIADNE_LOG(7,"    transition: "<<event<<" -> "<<target<<"\n");
                ARIADNE_LOG(7,"        reset="<<automaton.reset_assignments(location,event)<<"\n");
                ARIADNE_LOG(7,"          function="<<automaton.reset_function(location,event)<<"\n");
                ARIADNE_LOG(7,"        guard="<<automaton.guard_predicate(location,event)<<"\n");
                ARIADNE_LOG(7,"          function="<<automaton.guard_function(location,event)<<"\n");
                if(!reached.contains(target)) {
                    found.insert(target);
                    reached.insert(target);
                    steps.insert(target,step);
                }
           }

        }
        ARIADNE_LOG(5,"\nstep "<<step<<" found: "<<found<<"\n\n");
        working.clear();
        working.swap(found);
    }

    return reached;
}

CompositeHybridAutomaton parallel_composition(const List<HybridAutomaton>& components)
{
    return CompositeHybridAutomaton(name_composition(components),components);
}

HybridAutomaton flatten(const CompositeHybridAutomaton& composite_automaton, const List<DiscreteLocation>& locations)
{
    ARIADNE_NOT_IMPLEMENTED;
}


} // namespace Ariadne
