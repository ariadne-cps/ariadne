/***************************************************************************
 *            valve-permissive.hpp
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

#include "ariadne.hpp"

using namespace Ariadne;

inline AtomicHybridAutomaton getValve1()
{

    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4.0_decimal);

    // Declare the shared system variables
    RealVariable aperture("aperture1");

    // Declare the events we use
    DiscreteEvent e_idle("idle1");
    DiscreteEvent e_can_open("can_open1");
    DiscreteEvent e_can_close("can_close1");

    AtomicHybridAutomaton valve("valve1");

    // Declare the values the valve can variable can have
    AtomicDiscreteLocation closed("closed1");
    AtomicDiscreteLocation opened("opened1");
    AtomicDiscreteLocation opening("opening1");
    AtomicDiscreteLocation closing("closing1");

    // Since aperture is a known constant when the valve is open or closed,
    // specify aperture by an algebraic equation.
    valve.new_mode(closed,{let(aperture)=0});
    valve.new_mode(opened,{let(aperture)=1});
    valve.new_mode(opening,{dot(aperture)=+1/T});
    valve.new_mode(closing,{dot(aperture)=-1/T});

    valve.new_transition(opening,e_idle,opened,aperture>=1.0_dec,EventKind::URGENT);
    valve.new_transition(closing,e_idle,closed,aperture<=0.0_dec,EventKind::URGENT);
    valve.new_transition(closed,e_can_open,opening,{next(aperture)=aperture});
    valve.new_transition(opened,e_can_close,closing,{next(aperture)=aperture});
    valve.new_transition(opening,e_can_close,closing,{next(aperture)=aperture});
    valve.new_transition(closing,e_can_open,opening,{next(aperture)=aperture});

    return valve;
}
