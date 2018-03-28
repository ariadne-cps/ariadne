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

AtomicHybridAutomaton getValve()
{

    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4.0_decimal);

    // Declare the shared system variables
    RealVariable aperture("aperture");
    RealVariable height("height");

    // Declare the events we use
    DiscreteEvent f_opening("f_opening");
    DiscreteEvent f_closing("f_closing");
    DiscreteEvent b_opening("b_opening");
    DiscreteEvent b_closing("b_closing");

    AtomicHybridAutomaton valve("valve");

    // Declare the values the valve can variable can have
    AtomicDiscreteLocation opening("opening");
    AtomicDiscreteLocation closed("closed");
    AtomicDiscreteLocation opened("opened");
    AtomicDiscreteLocation closing("closing");

    valve.new_mode(opened,{let(aperture)=1.0_dec});
    valve.new_mode(closed,{let(aperture)=0.0_dec});
    valve.new_mode(opening,{dot(aperture)=+1/T});
    valve.new_mode(closing,{dot(aperture)=-1/T});

    valve.new_transition(closed,b_opening,opening);
    valve.new_transition(opening,f_opening,opened,aperture>=1,urgent);
    valve.new_transition(opened,b_closing,closing);
    valve.new_transition(closing,f_closing,closed,aperture<=0,urgent);

    return valve;
}
