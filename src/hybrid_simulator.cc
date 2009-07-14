/***************************************************************************
 *            hybrid_simulator.cc
 *
 *  Copyright  2008  Alberto Casagrande, Pieter Collins
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
 
#include "macros.h"
#include "array.h"
#include "tuple.h"
#include "stlio.h"
#include "vector.h"
#include "point.h"
#include "hybrid_set.h"
#include "expression_interface.h"
#include "function_interface.h"
#include "orbit.h"
#include "simulation_toolbox.h"
#include "evolution_parameters.h"

#include "logging.h"

#include "hybrid_time.h"
#include "hybrid_system.h"
#include "hybrid_simulator.h"



namespace Ariadne {
 

class DegenerateCrossingException { };


Simulator<HybridAutomaton>::Simulator()
    : _parameters(new EvolutionParameters()),
      _toolbox(new SimulationToolbox())
{
}


Simulator<HybridAutomaton>*
HybridSimulator<HybridAutomaton>::clone() const
{ 
    return new HybridSimulator(*this); 
}


Simulator<HybridAutomaton>::Simulator(const EvolutionParameters& p)
    : _parameters(new EvolutionParameters(p)),
      _toolbox(new SimulationToolbox())
{
}


Orbit<HybridPoint> 
Simulator<HybridAutomaton>::
orbit(const HybridAutomaton& system, 
      const HybridPoint& initial_point, 
      const HybridTime& maximum_time,
      Semantics semantics) const 
{
    const SimulationToolboxInterface& toolbox=this->_toolbox;

    Orbit<HybridPoint> result(initial_point);

    HybridTime time(0.0,0);
    
    HybridPoint current_state=initial_state;
    
    ARIADNE_ASSERT(system.has_mode(initial_state.first));

    while(time<maximum_time) {
        const DiscreteMode& current_location=current_state.first;
        const Point& current_point=current_state.second;
        std::set<DiscreteTransition>& transitions=system.transitions(current_mode);
        
        // Check for discrete events at initial time
        for(std::set<DiscreteTransition>::const_iterator transition_iter=transitions.begin();
            transition_iter!=transitions.end(); ++transition_iter) 
        {
            if(toolbox.active(*transition_iter->activation_ptr(),current_point)) {
                HybridPoint next_state(transition_iter->target(),toolbox.reset_step(transition_iter->reset(),current_point));
                orbit.insert(current_time,next_state);
                break;
            }
        }

        // Compute flow for given step size
        Point flow_point=toolbox.integration_step(automaton.mode(current_location).dynamic(),current_point,this->parameters().step_size());
        
        // Check for discrete events at final time
        std::map<DiscreteEvent,Float> active_urgent_events;
        std::map<DiscreteEvent,Float> active_nonurgent_events;
        for(std::set<DiscreteTransition>::const_iterator transition_iter=transitions.begin();
            transition_iter!=transitions.end(); ++transition_iter) 
        {
            if(toolbox.active(*transition_iter->activation_ptr(),flow_point)) {
                if(transition_iter->forced()) { active_urgent_events.insert(transition_iter->id()); }
                ele 
            }
        }
        
        
}


void 
HybridSimulator::
_evolution(EnclosureListType&,EnclosureListType&,EnclosureListType&,const SystemType& system, const EnclosureType& initial_point, const TimeType& time, Semantics semantics, bool) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


}  // namespace Ariadne

