/***************************************************************************
 *            valve-hysteresis-permissive.hpp
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

    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4);

    // Declare the shared system variables
    RealVariable aperture("aperture");
    RealVariable height("height");

    // Declare the events we use
    DiscreteEvent e_stop_opening("stop_opening");
    DiscreteEvent e_stop_closing("stop_closing");
    DiscreteEvent e_can_open("can_open");
    DiscreteEvent e_can_close("can_close");

    // Create the valve automaton
    AtomicHybridAutomaton valve("valve");

    // Declare the values the valve variable can have
    AtomicDiscreteLocation opening("opening");
    AtomicDiscreteLocation closed("closed");
    AtomicDiscreteLocation opened("opened");
    AtomicDiscreteLocation closing("closing");

    // Define the algebraic equations for the opened/closed locations.
    valve.new_mode(opened,{let(aperture)=1});
    valve.new_mode(closed,{let(aperture)=0});
    // Define the differential equations for the opening/closing locations.
    valve.new_mode(opening,{dot(aperture)=+1/T});
    valve.new_mode(closing,{dot(aperture)=-1/T});

    // Define the transitions: source location, event and target location;
    // then a mix of reset, guard and even kind can be present; if the event kind
    // is not specified, then also the guard can't be specified: this implicitly
    // means that the event is an input event for this automaton.
    valve.new_transition(closed,e_can_open,opening,{next(aperture)=aperture});
    valve.new_transition(opening,e_stop_opening,opened,aperture>=1,EventKind::URGENT);
    valve.new_transition(opened,e_can_close,closing,{next(aperture)=aperture});
    valve.new_transition(closing,e_stop_closing,closed,aperture<=0,EventKind::URGENT);

    return valve;
}
