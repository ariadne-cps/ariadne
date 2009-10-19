/***************************************************************************
 *            hybrid_automaton.cc
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

#include <map>

#include "macros.h"
#include "stlio.h"
#include "function.h"
#include "hybrid_time.h"
#include "hybrid_set.h"
#include "hybrid_automaton.h"
#include "grid_set.h"

namespace Ariadne {


typedef uint DimensionType;

class HybridSet {};


uint
AtomicDiscreteMode::
dimension() const
{
    return this->_dynamic.argument_size();
}


AtomicDiscreteMode::
AtomicDiscreteMode(DiscreteState location,
             const VectorFunction& dynamic)
    :  _location(location), _dynamic(dynamic), _invariants(), _grid(new Grid(dynamic.argument_size()))
{
}


/*
AtomicDiscreteMode::
AtomicDiscreteMode(DiscreteState location,
             const boost::shared_ptr< const VectorFunction > dynamic,
             const std::map< DiscreteEvent, boost::shared_ptr< const VectorFunction > >& invariants)
    :  _location(location), _dynamic(dynamic), _invariants(invariants), _grid(new Grid(dynamic->argument_size()))
{
    ARIADNE_ASSERT(dynamic->result_size()==dynamic->argument_size());
    for(uint i=0; i!=invariants.size(); ++i) {
        ARIADNE_ASSERT(invariants[i]->argument_size()==dynamic->argument_size());
        ARIADNE_ASSERT(invariants[i]->result_size()==1u);
    }
}
*/


std::ostream&
AtomicDiscreteMode::write(std::ostream& os) const
{
    const AtomicDiscreteMode& mode=*this;
    os << "AtomicDiscreteMode( "
       << "location=" << mode.location() << ", ";
    if(mode._algebraic_assignments.size()>0) {
        os << "algebraic_equations="<<mode._algebraic_assignments<<", "; }
    if(mode._differential_assignments.size()>0) {
        os << "differential_equations="<<mode._differential_assignments<<", "; }
    os << "dynamic=" << mode.dynamic() << ", "
       << "invariants=" << mode.invariants() << " )";
    return os;
}





AtomicDiscreteTransition::
AtomicDiscreteTransition(DiscreteEvent event,
                   const AtomicDiscreteMode& source,
                   const AtomicDiscreteMode& target,
                   const VectorFunction& reset,
                   const ScalarFunction& activation,
                   bool forced)
    : _event(event), _source(&source), _target(&target),
      _activation(activation), _reset(reset), _forced(forced)
{
    ARIADNE_ASSERT(activation.argument_size()==source.dimension());
    ARIADNE_ASSERT(reset.argument_size()==source.dimension());
    ARIADNE_ASSERT(reset.result_size()==target.dimension());
}




std::ostream&
AtomicDiscreteTransition::write(std::ostream& os) const
{
    const AtomicDiscreteTransition& transition=*this;
    return os << "AtomicDiscreteTransition( "
              << "event=" << transition.event() << ", "
              << "source=" << transition.source().location() << ", "
              << "target=" << transition.target().location() << ", "
              << "reset=" << transition.reset() << ", "
              << "activation=" << transition.activation() << ", "
              //<< "updates=" << transition._update_assignments << ", "
              //<< "predicate=" << transition._guard_predicate
              << " )";
}




AtomicHybridAutomaton::~AtomicHybridAutomaton()
{
}

AtomicHybridAutomaton::AtomicHybridAutomaton()
{
}

AtomicHybridAutomaton::AtomicHybridAutomaton(const std::string& name)
    : _name(name)
{
}





const AtomicDiscreteMode&
AtomicHybridAutomaton::new_mode(DiscreteState location,
                          const List<RealAlgebraicAssignment>& equations,
                          const List<RealDifferentialAssignment>& dynamic)
{
    if(this->has_mode(location)) {
        throw std::runtime_error("The hybrid automaton already has a mode with the given id");
    }

    RealSpace state_space;
    RealSpace auxiliary_space;
    RealSpace input_space;
    Set<String> defined_variables;
    Set<String> argument_variables;

    // Compute the auxiliary variables ordered by the given equations
    for(uint i=0; i!=equations.size(); ++i) {
        if(defined_variables.contains(equations[i].lhs.name())) {
            std::cerr<<defined_variables,equations[i].lhs.name();
            ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::new_mode",
                          "Variable "<<equations[i].lhs<<" is defined twice by the algebraic equations "<<equations<<" for mode "<<location);
        }
        auxiliary_space.append(equations[i].lhs);
        argument_variables.adjoin(equations[i].rhs.arguments());
    }

    // Compute the state variables ordered by the given differential equations
    for(uint i=0; i!=dynamic.size(); ++i) {
        if(defined_variables.contains(dynamic[i].lhs.base().name())) {
            ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::new_mode",
                          "Variable "<<dynamic[i].lhs.base()<<" is defined by the differential equations "<<dynamic<<" for mode "<<location<<" is already defined");
        }
        state_space.append(dynamic[i].lhs.base());
        argument_variables.adjoin(dynamic[i].rhs.arguments());
    }

    // Compute the input variables
    for(Set<String>::const_iterator variable_iter=argument_variables.begin();
        variable_iter!=argument_variables.end(); ++variable_iter)
    {
        if(!defined_variables.contains(*variable_iter)) {
            input_space.append(RealVariable(*variable_iter));
        }
    }

    // TODO: Compute function
    AtomicDiscreteMode new_mode=AtomicDiscreteMode(location,VectorFunction(1,1));

    new_mode._state_space=state_space;
    new_mode._auxiliary_space=auxiliary_space;
    new_mode._input_space=input_space;
    new_mode._algebraic_assignments=equations;
    new_mode._differential_assignments=dynamic;

    this->_modes.insert(new_mode);
    return this->mode(location);
}

const AtomicDiscreteMode&
AtomicHybridAutomaton::new_mode(DiscreteState location,
                          const List<RealDifferentialAssignment>& dynamic)
{
    return this->new_mode(location,List<RealAlgebraicAssignment>(),dynamic);
}

const AtomicDiscreteMode&
AtomicHybridAutomaton::new_mode(DiscreteState location,
                          const List<RealAlgebraicAssignment>& equations)
{
    return this->new_mode(location,equations,List<RealDifferentialAssignment>());
}



const AtomicDiscreteMode&
AtomicHybridAutomaton::new_invariant(DiscreteState location,
                               DiscreteEvent action,
                               const ContinuousPredicate& constraint)
{
    AtomicDiscreteMode& mode=const_cast<AtomicDiscreteMode&>(this->mode(location));
    mode._invariant_predicates.insert(make_pair(action,constraint));
    return mode;
}


const AtomicDiscreteTransition&
AtomicHybridAutomaton::new_transition(DiscreteState source,
                                DiscreteEvent event,
                                DiscreteState target,
                                const List<RealUpdateAssignment>& reset,
                                const ContinuousPredicate& guard)
{
    AtomicDiscreteMode& source_mode=const_cast<AtomicDiscreteMode&>(this->mode(source)); // Non-constant since we may wish to update input variables
    const AtomicDiscreteMode& target_mode=this->mode(target);
    RealSpace target_state_space=this->mode(target)._state_space;
    RealSpace target_auxiliary_space=this->mode(target)._auxiliary_space;
    Set<String> updated_variables;
    Set<String> argument_variables;
    for(uint i=0; i!=reset.size(); ++i) {
        if(!target_state_space.contains(reset[i].lhs.base())) {
            ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::new_transition",
                          "reset "<<reset<<" specifies variable "<<reset[i].lhs.base().name()<<
                          " which is not a state variable of the target location "<<target);
        }
        if(target_auxiliary_space.contains(reset[i].lhs.base())) {
            ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::new_transition",
                          "reset "<<reset<<" specifies variable "<<reset[i].lhs.base().name()<<
                          " which is already specified by the algebraic equations valid in "<<target);
        }
        if(updated_variables.contains(reset[i].lhs.base().name())) {
            ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::new_transition",
                          "reset "<<reset<<" specifies variable "<<reset[i].lhs.base().name()<<" twice.");
        } else {
            updated_variables.insert(reset[i].lhs.base().name());
        }
        argument_variables.adjoin(reset[i].rhs.arguments());
    }
    argument_variables.adjoin(guard.arguments());
    for(Set<String>::const_iterator variable_iter=argument_variables.begin();
        variable_iter!=argument_variables.end(); ++variable_iter)
    {
        source_mode._input_space.insert(RealVariable(*variable_iter));
    }

    VectorFunction reset_function(1,1);
    ScalarFunction guard_function(1);
    AtomicDiscreteTransition new_transition=AtomicDiscreteTransition(event,source_mode,target_mode,reset_function,guard_function,permissive);
    new_transition._update_assignments=reset;
    new_transition._guard_predicate=guard;
    this->_transitions.insert(new_transition);
    return this->transition(event,source);
}


const AtomicDiscreteTransition&
AtomicHybridAutomaton::new_transition(DiscreteState source,
                                DiscreteEvent event,
                                DiscreteState target,
                                const ContinuousPredicate& guard)
{
    return this->new_transition(source,event,target,List<RealUpdateAssignment>(),guard);
}


const AtomicDiscreteMode&
AtomicHybridAutomaton::new_mode(DiscreteState location,
                          const VectorFunction& dynamic)
{
    if(this->has_mode(location)) {
        throw std::runtime_error("The hybrid automaton already has a mode with the given id");
    }
    if(dynamic.result_size()!=dynamic.argument_size()) {
        ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::new_mode(location,dynamic)",
            "The dynamic has argument size " << dynamic.argument_size()
                << " and result size " << dynamic.result_size() << ", so does not define a vector field.");
    }
    this->_modes.insert(AtomicDiscreteMode(location,dynamic));
    return this->mode(location);
}


const AtomicDiscreteMode&
AtomicHybridAutomaton::new_invariant(DiscreteState location,
                               const ScalarFunction& invariant)
{
    ARIADNE_ASSERT(location>0);
    if(!this->has_mode(location)) {
        throw std::runtime_error("The location of the invariant must be in the automaton.");
    }
    AtomicDiscreteMode& mode=const_cast<AtomicDiscreteMode&>(this->mode(location));
    if(invariant.argument_size()!=mode.dimension()) {
        ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::new_invariant(location,invariant)",
            "The invariant has argument size " << invariant.argument_size()
                << " but the mode has state-space dimension " << mode.dimension());
    }
    DiscreteEvent invariant_event(std::string("inv")+to_string(mode._invariants.size()+1u));
    mode._invariants[invariant_event]=invariant;
    return mode;
}


const AtomicDiscreteMode&
AtomicHybridAutomaton::new_invariant(DiscreteState location,
                               const VectorFunction& invariant)
{
    if(invariant.result_size()!=1u) {
        ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::new_invariant(location,invariant)",
            "The invariant has result size " << invariant.result_size()
                << " but only scalar invariants are currently supported.");
    }
    return this->new_invariant(location,invariant[0]);
}


const AtomicDiscreteTransition&
AtomicHybridAutomaton::new_transition(const DiscreteEvent event_id,
                                const DiscreteState source_id,
                                const DiscreteState target_id,
                                const VectorFunction& reset,
                                const ScalarFunction& activation,
                                bool forced)
{
    //ARIADNE_ASSERT_MSG(event_id>0, "Transition event should be positive.");
    if(this->has_transition(event_id,source_id)) {
        throw std::runtime_error("The automaton already has a transition with the given event_id and source id.");
    }
    if(!this->has_mode(source_id)) {
        throw std::runtime_error("The source mode of the transition must be in the automaton");
    }
    if(!this->has_mode(target_id)) {
        throw std::runtime_error("The desitination mode of the transition must be in the automaton");
    }

    const AtomicDiscreteMode& this_source=this->mode(source_id);
    const AtomicDiscreteMode& this_target=this->mode(target_id);

    this->_transitions.insert(AtomicDiscreteTransition(event_id,this_source,this_target,reset,activation,forced));
    return this->transition(event_id,source_id);
}


const AtomicDiscreteTransition&
AtomicHybridAutomaton::
new_transition(DiscreteEvent event,
               DiscreteState source,
               DiscreteState target,
               const VectorFunction &reset,
               const VectorFunction &activation,
               bool forced)
{
    if(activation.result_size()!=1u) {
        ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::new_transition(...)",
            "The activation function has result size " << activation.result_size()
                << " but only scalar activations are currently supported.");
    }
    return this->new_transition(event,source,target,reset,activation[0],forced);
}


const AtomicDiscreteTransition&
AtomicHybridAutomaton::
new_forced_transition(DiscreteEvent event,
                      DiscreteState source,
                      DiscreteState target,
                      const VectorFunction &reset,
                      const VectorFunction &activation)
{
    ARIADNE_ASSERT_MSG(event>0, "Transition event should be positive.");
    if(this->has_transition(event,source)) {
        throw std::runtime_error("The automaton already has a transition with the given id and source id.");
    }
    if(!this->has_mode(source)) {
        throw std::runtime_error("The automaton does not contain a mode with ths given source id");
    }
    if(!this->has_mode(target)) {
        throw std::runtime_error("The automaton does not contain a mode with ths given desitination id");
    }
    if(activation.result_size()!=1u) {
        ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::new_transition(...)",
            "The activation function has result size " << activation.result_size()
                << " but only scalar activations are currently supported.");
    }

    const AtomicDiscreteMode& source_mode=this->mode(source);
    const AtomicDiscreteMode& target_mode=this->mode(target);
    this->_transitions.insert(AtomicDiscreteTransition(event,source_mode,target_mode,reset,activation[0],true));
    return this->transition(event,source);
}


const AtomicDiscreteTransition&
AtomicHybridAutomaton::
new_unforced_transition(DiscreteEvent event,
                        DiscreteState source,
                        DiscreteState target,
                        const VectorFunction &reset,
                        const VectorFunction &activation)
{
    ARIADNE_ASSERT_MSG(event>0, "Transition event should be positive.");
    if(this->has_transition(event,source)) {
        throw std::runtime_error("The automaton already has a transition with the given id and source id.");
    }
    if(!this->has_mode(source)) {
        throw std::runtime_error("The automaton does not contain a mode with ths given source id");
    }
    if(!this->has_mode(target)) {
        throw std::runtime_error("The automaton does not contain a mode with ths given desitination id");
    }
    if(activation.result_size()!=1u) {
        ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::new_transition(...)",
            "The activation function has result size " << activation.result_size()
                << " but only scalar activations are currently supported.");
    }

    const AtomicDiscreteMode& source_mode=this->mode(source);
    const AtomicDiscreteMode& target_mode=this->mode(target);
    this->_transitions.insert(AtomicDiscreteTransition(event,source_mode,target_mode,reset,activation[0],false));
    return this->transition(event,source);
}


void
AtomicHybridAutomaton::set_grid(DiscreteState location,
                          const Grid& grid)
{
    if(!this->has_mode(location)) {
        throw std::runtime_error("The automaton does not contain a mode with ths given location id");
    }
    AtomicDiscreteMode& mode=const_cast<AtomicDiscreteMode&>(this->mode(location));
    if(grid.dimension()!=mode.dimension()) {
        throw std::runtime_error("The mode of the automaton has a different dimension to the grid.");
    }
    mode._grid=shared_ptr<Grid>(new Grid(grid));
}

void
AtomicHybridAutomaton::set_grid(const Grid& grid)
{
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
    {
        AtomicDiscreteMode& mode=const_cast<AtomicDiscreteMode&>(*mode_iter);
        if(grid.dimension()!=mode.dimension()) {
            throw std::runtime_error("The automaton has a different dimension to the grid.");
        }
        mode._grid=shared_ptr<Grid>(new Grid(grid));
    }
}

void
AtomicHybridAutomaton::set_grid(const HybridGrid& hgrid)
{
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
    {
        AtomicDiscreteMode& mode=const_cast<AtomicDiscreteMode&>(*mode_iter);
        DiscreteState loc = mode.location();
        if(hgrid.find(loc) == hgrid.end()) {
            throw std::runtime_error("The automaton does not contain a mode with this given location id");
        }
        if(hgrid[loc].dimension()!=mode.dimension()) {
            throw std::runtime_error("The mode of the automaton has a different dimension to the grid.");
        }
        mode._grid=shared_ptr<Grid>(new Grid(hgrid[loc]));
    }
}



bool
AtomicHybridAutomaton::has_mode(DiscreteState state) const
{
    // FIXME: This is a hack since we use std::set which cannot be searched by id.
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
        {
            if(mode_iter->location()==state) {
                return true;
            }
        }
    return false;
}


bool
AtomicHybridAutomaton::has_transition(DiscreteEvent event, DiscreteState source) const
{
   return this->has_transition(source,event);
}


bool
AtomicHybridAutomaton::has_transition(DiscreteState source, DiscreteEvent event) const
{
    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
        {
            if(transition_iter->event()==event && transition_iter->source().location()==source) {
                return true;
            }
        }
    return false;
}



HybridSet
AtomicHybridAutomaton::invariant() const
{
    ARIADNE_NOT_IMPLEMENTED;
}



HybridSpace
AtomicHybridAutomaton::state_space() const
{
    HybridSpace result;
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
        {
            result[mode_iter->location()]=mode_iter->dimension();
        }
    return result;
}


const std::set< AtomicDiscreteMode >&
AtomicHybridAutomaton::modes() const
{
    return this->_modes;
}



const AtomicDiscreteMode&
AtomicHybridAutomaton::mode(DiscreteState state) const
{
    // FIXME: This is a hack; we should use a logarithmic time real search to find a mode with the given discrete state.
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
        {
            if(mode_iter->location()==state) {
                return *mode_iter;
            }
        }
    ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::mode(DiscreteState)",state<<" is not a location of the automaton with modes "<<this->modes());
}


const std::set< AtomicDiscreteTransition >&
AtomicHybridAutomaton::transitions() const
{
    return this->_transitions;
}



std::set< AtomicDiscreteTransition >
AtomicHybridAutomaton::transitions(DiscreteState source) const
{
    std::set< AtomicDiscreteTransition > result;
    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
        {
            if(transition_iter->source().location()==source) {
                result.insert(*transition_iter);
            }
        }
    return result;
}


std::map<DiscreteEvent,ScalarFunction>
AtomicHybridAutomaton::blocking_guards(DiscreteState source) const
{
    std::map<DiscreteEvent,ScalarFunction> result;
    const AtomicDiscreteMode& mode=this->mode(source);
    for(invariant_const_iterator invariant_iter=mode._invariants.begin();
        invariant_iter!=mode._invariants.end(); ++invariant_iter)
    {
        const DiscreteEvent event=invariant_iter->first;
        const ScalarFunction invariant=invariant_iter->second;
        result[event]=invariant;
    }

    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
    {
        if(transition_iter->source().location()==source && transition_iter->forced()) {
            const DiscreteEvent event=transition_iter->event();
            const ScalarFunction guard=transition_iter->activation();
            result[event]=guard;
        }
    }
    return result;
}


std::map<DiscreteEvent,ScalarFunction>
AtomicHybridAutomaton::permissive_guards(DiscreteState source) const
{
    std::map<DiscreteEvent,ScalarFunction> result;

    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
    {
        if(transition_iter->source().location()==source && !transition_iter->forced()) {
            const DiscreteEvent event=transition_iter->event();
            const ScalarFunction guard=transition_iter->activation();
            result[event]=guard;
        }
    }
    return result;
}




const AtomicDiscreteTransition&
AtomicHybridAutomaton::transition(DiscreteEvent event, DiscreteState source) const
{
    return this->transition(source,event);
}


const AtomicDiscreteTransition&
AtomicHybridAutomaton::transition(DiscreteState source, DiscreteEvent event) const
{
    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
    {
        if(transition_iter->event()==event && transition_iter->source().location()==source) {
            return *transition_iter;
        }
    }
    ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::transition(DiscreteEvent,DiscreteState)",
                  "Transition ("<<event<<","<<source<<") is not in the hybrid automaton with transitions "<<this->transitions());
}


const std::string&
AtomicHybridAutomaton::name() const
{
    return this->_name;
}

Grid
AtomicHybridAutomaton::grid(DiscreteState location) const
{
    ARIADNE_ASSERT(this->has_mode(location));
    const AtomicDiscreteMode& mode=this->mode(location);
    return Grid(mode.dimension());
}

HybridGrid
AtomicHybridAutomaton::grid() const
{
    HybridGrid result;
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
    {
        result[mode_iter->location()]=mode_iter->grid();
    }
    return result;
}

std::ostream&
AtomicHybridAutomaton::write(std::ostream& os) const
{
    const AtomicHybridAutomaton& ha=*this;
    os << "\nHybridAutomaton( \n  modes=\n";
    for(AtomicHybridAutomaton::mode_const_iterator mode_iter=ha.modes().begin();
        mode_iter!=ha.modes().end(); ++mode_iter)
    {
        os << "    " <<*mode_iter<<",\n";
    }
    os << "  transitions=\n";
    for(AtomicHybridAutomaton::transition_const_iterator transition_iter=ha.transitions().begin();
        transition_iter!=ha.transitions().end(); ++transition_iter)
    {
        os << "    " <<*transition_iter<<",\n";
    }
    return os << ")\n";
}



DiscreteState
AtomicHybridAutomaton::target(const DiscreteEvent& event, const DiscreteState& source) {
    return this->transition(event,source).target().location();
}

List<RealVariable>
AtomicHybridAutomaton::state_variables(DiscreteState location) const {
    List<RealVariable> result;
    const AtomicDiscreteMode& mode=this->mode(location);
    for(uint i=0; i!=mode._differential_assignments.size(); ++i) {
        result.append(mode._differential_assignments[i].lhs.base());
    }
    return result;
}

List<RealVariable>
AtomicHybridAutomaton::auxiliary_variables(DiscreteState location) const {
    List<RealVariable> result;
    const AtomicDiscreteMode& mode=this->mode(location);
    for(uint i=0; i!=mode._algebraic_assignments.size(); ++i) {
        result.append(mode._algebraic_assignments[i].lhs);
    }
    return result;
}


List<RealAssignment>
AtomicHybridAutomaton::algebraic_assignments(const DiscreteState& location) const {
    return this->mode(location)._algebraic_assignments;
}

List<DottedRealAssignment>
AtomicHybridAutomaton::differential_assignments(const DiscreteState& location) const {
    return this->mode(location)._differential_assignments;
}

List<PrimedRealAssignment>
AtomicHybridAutomaton::update_assignments(const DiscreteState& source, const DiscreteEvent& event) const {
    return this->transition(event,source)._update_assignments;
}

Map<DiscreteEvent,ContinuousPredicate>
AtomicHybridAutomaton::invariant_predicates(const DiscreteState& location) const {
    return this->mode(location)._invariant_predicates;
}

ContinuousPredicate
AtomicHybridAutomaton::invariant_predicate(const DiscreteState& location, const DiscreteEvent& action) const {
    return this->mode(location)._invariant_predicates[action];
}

ContinuousPredicate
AtomicHybridAutomaton::guard_predicate(const DiscreteState& source, const DiscreteEvent& event) const {
    return this->transition(event,source)._guard_predicate;
}






CompositeHybridAutomaton::CompositeHybridAutomaton(const AtomicHybridAutomaton& automaton)
    : _components(1u,automaton) { }

CompositeHybridAutomaton::CompositeHybridAutomaton(const List<AtomicHybridAutomaton>& components)
    : _components(components) { }

uint
CompositeHybridAutomaton::number_of_components() const
{
    return this->_components.size();
}

const AtomicHybridAutomaton&
CompositeHybridAutomaton::component(uint k) const {
    return this->_components[k];
}

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


VectorFunction
CompositeHybridAutomaton::output(const DiscreteLocation& location) const {
    return VectorFunction(auxiliary_variables(location),algebraic_assignments(location),state_variables(location));
}

VectorFunction
CompositeHybridAutomaton::dynamic(const DiscreteLocation& location) const {
    return VectorFunction(dot(state_variables(location)),differential_assignments(location),variables(location));
}

VectorFunction
CompositeHybridAutomaton::reset(const DiscreteLocation& source, const DiscreteEvent& event) const {
    DiscreteLocation target=this->target(source,event);
    return VectorFunction(next(state_variables(target)),update_assignments(source,event),variables(source));
}

Map<DiscreteEvent,ScalarFunction>
CompositeHybridAutomaton::invariants(const DiscreteLocation& location) const {
    ARIADNE_NOT_IMPLEMENTED;
}

ScalarFunction
CompositeHybridAutomaton::guard(const DiscreteLocation& location, const DiscreteEvent& event) const {
    return ScalarFunction(guard_predicate(location,event),variables(location));
}

bool
CompositeHybridAutomaton::has_mode(const DiscreteLocation& location) const {
    for(uint i=0; i!=this->_components.size(); ++i) {
        if(!this->_components[i].has_mode(location[i])) {
            return false;
        }
    }
    return true;
}

bool
CompositeHybridAutomaton::has_transition(const DiscreteLocation& source, const DiscreteEvent& event) const {
    for(uint i=0; i!=this->_components.size(); ++i) {
        if(this->_components[i].has_transition(source[i],event)) {
            return true;
        }
    }
    return false;
}

CompositeHybridAutomaton::DiscreteLocation
CompositeHybridAutomaton::target(const DiscreteLocation& source, const DiscreteEvent& event) const {
    DiscreteLocation result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        if(this->_components[i].has_transition(source[i],event)) {
            result.append(this->_components[i].transition(source[i],event).target().location());
        } else {
            result.append(source);
        }
    }
    return result;
}


List<RealVariable>
CompositeHybridAutomaton::variables(const DiscreteLocation& location) const {
    return catenate(this->state_variables(location),this->auxiliary_variables(location));
}

List<RealVariable>
CompositeHybridAutomaton::state_variables(const DiscreteLocation& location) const {
    List<RealVariable> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        ARIADNE_ASSERT(this->_components[i].has_mode(location[i]));
        // Append state space in mode i; note that this automatically checks whether a variable is already present
        result.append(this->_components[i].state_variables(location[i]));
    }
    return result;
}

List<RealVariable>
CompositeHybridAutomaton::auxiliary_variables(const DiscreteLocation& location) const {
    List<RealVariable> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        ARIADNE_ASSERT(this->_components[i].has_mode(location[i]));
        // Append state space in mode i; note that this automatically checks whether a variable is already present
        result.append(this->_components[i].auxiliary_variables(location[i]));
    }
    return result;
}

// Find all algebraic equations valid in the location
List<RealAssignment>
CompositeHybridAutomaton::algebraic_assignments(const DiscreteLocation& location) const {
    List<RealAssignment> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result.append(this->_components[i].algebraic_assignments(location[i]));
    }
    // TODO: Sort the result to eliminate algebraic loops
    // sort(result);
    return result;
}

List<DottedRealAssignment>
CompositeHybridAutomaton::differential_assignments(const DiscreteLocation& location) const {
    List<DottedRealAssignment> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result.append(this->_components[i].differential_assignments(location[i]));
    }
    // No need to sort result since dotted variables cannot appear in the right-hand side (currently)
    return result;
}


List<PrimedRealAssignment>
CompositeHybridAutomaton::update_assignments(const DiscreteLocation& location, const DiscreteEvent& event) const {
    List<PrimedRealAssignment> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result.append(this->_components[i].update_assignments(location[i],event));
    }
    return result;
}

ContinuousPredicate& operator&=(ContinuousPredicate& p1, const ContinuousPredicate& p2) {
    ContinuousPredicate p3 = (p1 && p2); p1=p3; return p1; }

ContinuousPredicate
CompositeHybridAutomaton::guard_predicate(const DiscreteLocation& location, const DiscreteEvent& event) const {
    ContinuousPredicate result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result &= this->_components[i].guard_predicate(location[i],event);
    }
    return result;
}


std::ostream&
CompositeHybridAutomaton::write(std::ostream& os) const
{
    return os << "CompositeHybridAutomaton(\n" << this->_components << "\n)\n";
}

}