/***************************************************************************
 *            valve-proportional-urgent.hpp
 *
 *  Copyright  2017 Luca Geretti
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

#include <cstdarg>
#include "ariadne.hpp"

using namespace Ariadne;

inline AtomicHybridAutomaton getValve()
{
    // Declare some constants.
    RealConstant K("K",2.0_decimal); // Gain of the proportional controller
    RealConstant Ref("Ref",7.0_decimal); // Reference height

    // Declare the shared system variables
    RealVariable aperture("aperture");
    RealVariable height("height");

    // Declare the events we use
    DiscreteEvent start_modulating("start_modulating");
    DiscreteEvent finished_opening("finished_opening");
    DiscreteEvent finished_closing("finished_closing");

    // Create a valve automaton
    AtomicHybridAutomaton valve_automaton("valve");

    // Declare the locations for the valve automaton
    AtomicDiscreteLocation opened("opened");
    AtomicDiscreteLocation modulated("modulated");
    AtomicDiscreteLocation closed("closed");

    // Since aperture is a known constant when the valve is open or closed,
    // specify aperture by an algebraic equation
    valve_automaton.new_mode(opened,{let(aperture)=1});
    valve_automaton.new_mode(closed,{let(aperture)=0});
    // Specify the differential equation for when the proportional control is in effect
    valve_automaton.new_mode(modulated,{let(aperture)=K*(Ref-height)});

    // Define the transitions
    valve_automaton.new_transition(modulated,finished_opening,closed,K*(Ref-height)<=0,EventKind::URGENT);
    valve_automaton.new_transition(modulated,finished_closing,opened,K*(Ref-height)>=1,EventKind::URGENT);
    valve_automaton.new_transition(opened,start_modulating,modulated,K*(Ref-height)<=1,EventKind::URGENT);
    valve_automaton.new_transition(closed,start_modulating,modulated,K*(Ref-height)>=0,EventKind::URGENT);

    return valve_automaton;
}
