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
DiscreteMode::
dimension() const
{
    return this->_dynamic.argument_size();
}


DiscreteMode::
DiscreteMode(DiscreteState location,
             const VectorFunction& dynamic)
    :  _location(location), _dynamic(dynamic), _invariants(), _grid(new Grid(dynamic.argument_size()))
{
}


/*
DiscreteMode::
DiscreteMode(DiscreteState location,
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
DiscreteMode::write(std::ostream& os) const
{
    const DiscreteMode& mode=*this;
    os << "DiscreteMode( "
       << "location=" << mode.location() << ", ";
    if(mode._algebraic_equations.size()>0) {
        os << "algebraic_equations="<<mode._algebraic_equations<<", "; }
    if(mode._differential_equations.size()>0) {
        os << "differential_equations="<<mode._differential_equations<<", "; }
    os << "dynamic=" << mode.dynamic() << ", "
       << "invariants=" << mode.invariants() << " )";
    return os;
}





DiscreteTransition::
DiscreteTransition(DiscreteEvent event,
                   const DiscreteMode& source,
                   const DiscreteMode& target,
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
DiscreteTransition::write(std::ostream& os) const
{
    const DiscreteTransition& transition=*this;
    return os << "DiscreteTransition( "
              << "event=" << transition.event() << ", "
              << "source=" << transition.source().location() << ", "
              << "target=" << transition.target().location() << ", "
              << "reset=" << transition.reset() << ", "
              << "activation=" << transition.activation()
              << "updates=" << transition._update_equations
              << "predicate=" << transition._guard_predicate
              << " )";
}




HybridAutomaton::~HybridAutomaton()
{
}

HybridAutomaton::HybridAutomaton()
{
}

HybridAutomaton::HybridAutomaton(const std::string& name)
    : _name(name)
{
}





const DiscreteMode&
HybridAutomaton::new_mode(DiscreteState location,
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
            ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_mode",
                          "Variable "<<equations[i].lhs<<" is defined twice by the algebraic equations "<<equations<<" for mode "<<location);
        }
        auxiliary_space.append(equations[i].lhs);
        argument_variables.adjoin(equations[i].rhs.arguments());
    }

    // Compute the state variables ordered by the given differential equations
    for(uint i=0; i!=dynamic.size(); ++i) {
        if(defined_variables.contains(dynamic[i].lhs.base().name())) {
            ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_mode",
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
    DiscreteMode new_mode=DiscreteMode(location,VectorFunction(1,1));

    new_mode._state_space=state_space;
    new_mode._auxiliary_space=auxiliary_space;
    new_mode._input_space=input_space;
    new_mode._algebraic_equations=equations;
    new_mode._differential_equations=dynamic;

    this->_modes.insert(new_mode);
    return this->mode(location);
}

const DiscreteMode&
HybridAutomaton::new_mode(DiscreteState location,
                          const List<RealDifferentialAssignment>& dynamic)
{
    return this->new_mode(location,List<RealAlgebraicAssignment>(),dynamic);
}

const DiscreteMode&
HybridAutomaton::new_mode(DiscreteState location,
                          const List<RealAlgebraicAssignment>& equations)
{
    return this->new_mode(location,equations,List<RealDifferentialAssignment>());
}



const DiscreteTransition&
HybridAutomaton::new_transition(DiscreteEvent event,
                                DiscreteState source,
                                DiscreteState target,
                                const List<RealUpdateAssignment>& reset,
                                const ContinuousPredicate& guard,
                                bool urgency)
{
    DiscreteMode& source_mode=const_cast<DiscreteMode&>(this->mode(source)); // Non-constant since we may wish to update input variables
    const DiscreteMode& target_mode=this->mode(target);
    RealSpace target_state_space=this->mode(target)._state_space;
    RealSpace target_auxiliary_space=this->mode(target)._auxiliary_space;
    Set<String> updated_variables;
    Set<String> argument_variables;
    for(uint i=0; i!=reset.size(); ++i) {
        if(!target_state_space.contains(reset[i].lhs.base())) {
            ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_transition",
                          "reset "<<reset<<" specifies variable "<<reset[i].lhs.base().name()<<
                          " which is not a state variable of the target location "<<target);
        }
        if(target_auxiliary_space.contains(reset[i].lhs.base())) {
            ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_transition",
                          "reset "<<reset<<" specifies variable "<<reset[i].lhs.base().name()<<
                          " which is already specified by the algebraic equations valid in "<<target);
        }
        if(updated_variables.contains(reset[i].lhs.base().name())) {
            ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_transition",
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
    DiscreteTransition new_transition=DiscreteTransition(event,source_mode,target_mode,reset_function,guard_function,urgency);
    new_transition._update_equations=reset;
    new_transition._guard_predicate=guard;
    this->_transitions.insert(new_transition);
    return this->transition(event,source);
}


const DiscreteMode&
HybridAutomaton::new_mode(DiscreteState location,
                          const VectorFunction& dynamic)
{
    if(this->has_mode(location)) {
        throw std::runtime_error("The hybrid automaton already has a mode with the given id");
    }
    if(dynamic.result_size()!=dynamic.argument_size()) {
        ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_mode(location,dynamic)",
            "The dynamic has argument size " << dynamic.argument_size()
                << " and result size " << dynamic.result_size() << ", so does not define a vector field.");
    }
    this->_modes.insert(DiscreteMode(location,dynamic));
    return this->mode(location);
}


const DiscreteMode&
HybridAutomaton::new_invariant(DiscreteState location,
                               const ScalarFunction& invariant)
{
    ARIADNE_ASSERT(location>0);
    if(!this->has_mode(location)) {
        throw std::runtime_error("The location of the invariant must be in the automaton.");
    }
    DiscreteMode& mode=const_cast<DiscreteMode&>(this->mode(location));
    if(invariant.argument_size()!=mode.dimension()) {
        ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_invariant(location,invariant)",
            "The invariant has argument size " << invariant.argument_size()
                << " but the mode has state-space dimension " << mode.dimension());
    }
    DiscreteEvent invariant_event(-8-mode._invariants.size());
    mode._invariants[invariant_event]=invariant;
    return mode;
}


const DiscreteMode&
HybridAutomaton::new_invariant(DiscreteState location,
                               const VectorFunction& invariant)
{
    if(invariant.result_size()!=1u) {
        ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_invariant(location,invariant)",
            "The invariant has result size " << invariant.result_size()
                << " but only scalar invariants are currently supported.");
    }
    return this->new_invariant(location,invariant[0]);
}


const DiscreteTransition&
HybridAutomaton::new_transition(const DiscreteEvent event_id,
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

    const DiscreteMode& this_source=this->mode(source_id);
    const DiscreteMode& this_target=this->mode(target_id);

    this->_transitions.insert(DiscreteTransition(event_id,this_source,this_target,reset,activation,forced));
    return this->transition(event_id,source_id);
}


const DiscreteTransition&
HybridAutomaton::
new_transition(DiscreteEvent event,
               DiscreteState source,
               DiscreteState target,
               const VectorFunction &reset,
               const VectorFunction &activation,
               bool forced)
{
    if(activation.result_size()!=1u) {
        ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_transition(...)",
            "The activation function has result size " << activation.result_size()
                << " but only scalar activations are currently supported.");
    }
    return this->new_transition(event,source,target,reset,activation[0],forced);
}


const DiscreteTransition&
HybridAutomaton::
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
        ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_transition(...)",
            "The activation function has result size " << activation.result_size()
                << " but only scalar activations are currently supported.");
    }

    const DiscreteMode& source_mode=this->mode(source);
    const DiscreteMode& target_mode=this->mode(target);
    this->_transitions.insert(DiscreteTransition(event,source_mode,target_mode,reset,activation[0],true));
    return this->transition(event,source);
}


const DiscreteTransition&
HybridAutomaton::
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
        ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_transition(...)",
            "The activation function has result size " << activation.result_size()
                << " but only scalar activations are currently supported.");
    }

    const DiscreteMode& source_mode=this->mode(source);
    const DiscreteMode& target_mode=this->mode(target);
    this->_transitions.insert(DiscreteTransition(event,source_mode,target_mode,reset,activation[0],false));
    return this->transition(event,source);
}


void
HybridAutomaton::set_grid(DiscreteState location,
                          const Grid& grid)
{
    if(!this->has_mode(location)) {
        throw std::runtime_error("The automaton does not contain a mode with ths given location id");
    }
    DiscreteMode& mode=const_cast<DiscreteMode&>(this->mode(location));
    if(grid.dimension()!=mode.dimension()) {
        throw std::runtime_error("The mode of the automaton has a different dimension to the grid.");
    }
    mode._grid=shared_ptr<Grid>(new Grid(grid));
}

void
HybridAutomaton::set_grid(const Grid& grid)
{
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
    {
        DiscreteMode& mode=const_cast<DiscreteMode&>(*mode_iter);
        if(grid.dimension()!=mode.dimension()) {
            throw std::runtime_error("The automaton has a different dimension to the grid.");
        }
        mode._grid=shared_ptr<Grid>(new Grid(grid));
    }
}

void
HybridAutomaton::set_grid(const HybridGrid& hgrid)
{
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
    {
        DiscreteMode& mode=const_cast<DiscreteMode&>(*mode_iter);
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
HybridAutomaton::has_mode(DiscreteState state) const
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
HybridAutomaton::has_transition(DiscreteEvent event, DiscreteState source) const
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
HybridAutomaton::invariant() const
{
    ARIADNE_NOT_IMPLEMENTED;
}



HybridSpace
HybridAutomaton::state_space() const
{
    HybridSpace result;
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
        {
            result[mode_iter->location()]=mode_iter->dimension();
        }
    return result;
}


const std::set< DiscreteMode >&
HybridAutomaton::modes() const
{
    return this->_modes;
}



const DiscreteMode&
HybridAutomaton::mode(DiscreteState state) const
{
    // FIXME: This is a hack; we should use a logarithmic time real search to find a mode with the given discrete state.
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
        {
            if(mode_iter->location()==state) {
                return *mode_iter;
            }
        }
    throw std::runtime_error("The hybrid automaton does not have a mode with the given id.");
}


const std::set< DiscreteTransition >&
HybridAutomaton::transitions() const
{
    return this->_transitions;
}



std::set< DiscreteTransition >
HybridAutomaton::transitions(DiscreteState source) const
{
    std::set< DiscreteTransition > result;
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
HybridAutomaton::blocking_guards(DiscreteState source) const
{
    std::map<DiscreteEvent,ScalarFunction> result;
    const DiscreteMode& mode=this->mode(source);
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
HybridAutomaton::permissive_guards(DiscreteState source) const
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




const DiscreteTransition&
HybridAutomaton::transition(DiscreteEvent event, DiscreteState source) const
{
    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
        {
            if(transition_iter->event()==event && transition_iter->source().location()==source) {
                return *transition_iter;
            }
        }
    throw std::runtime_error("The hybrid automaton does not have a transition with the given event and source.");
}


const std::string&
HybridAutomaton::name() const
{
    return this->_name;
}

Grid
HybridAutomaton::grid(DiscreteState location) const
{
    ARIADNE_ASSERT(this->has_mode(location));
    const DiscreteMode& mode=this->mode(location);
    return Grid(mode.dimension());
}

HybridGrid
HybridAutomaton::grid() const
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
operator<<(std::ostream& os, const HybridAutomaton& ha)
{
    return os << "HybridAutomaton( modes=" << ha.modes() << ", transitions=" << ha.transitions() << ")";
}





}
