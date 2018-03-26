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

#include <cstdarg>
#include "ariadne.hpp"

using namespace Ariadne;

AtomicHybridAutomaton getValve3()
{

    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4.0_decimal);

    // Declare the shared system variables
    RealVariable aperture3("aperture3");

    // Declare the events we use
    DiscreteEvent e_idle("idle3");
    DiscreteEvent e_close("close3");
    DiscreteEvent e_open("open3");

    AtomicHybridAutomaton valve("valve3");

    // Declare the values the valve can variable can have
    AtomicDiscreteLocation idle3("idle3");
    AtomicDiscreteLocation opening3("opening3");
    AtomicDiscreteLocation closing3("closing3");

    // Since aperture is a known constant when the valve is open or closed,
    // specify aperture by an algebraic equation.
    valve.new_mode(idle3,{dot(aperture3)=0.0_decimal});
    valve.new_mode(opening3,{dot(aperture3)=+1/T});
    valve.new_mode(closing3,{dot(aperture3)=-1/T});

    valve.new_transition(opening3,e_idle,idle3,{next(aperture3)=1.0_dec},aperture3>=1.0_dec,urgent);
    valve.new_transition(closing3,e_idle,idle3,{next(aperture3)=0.0_dec},aperture3<=0.0_dec,urgent);
    valve.new_transition(idle3,e_open,opening3,{next(aperture3)=aperture3});
    valve.new_transition(idle3,e_close,opening3,{next(aperture3)=aperture3});

    return valve;
}
