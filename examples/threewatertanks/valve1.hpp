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

AtomicHybridAutomaton getValve1()
{

    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4.0_decimal);

    // Declare the shared system variables
    RealVariable aperture1("aperture1");

    // Declare the events we use
    DiscreteEvent e_idle("idle1");
    DiscreteEvent e_close("close1");
    DiscreteEvent e_open("open1");

    AtomicHybridAutomaton valve("valve1");

    // Declare the values the valve can variable can have
    AtomicDiscreteLocation idle1("idle1");
    AtomicDiscreteLocation opening1("opening1");
    AtomicDiscreteLocation closing1("closing1");

    // Since aperture is a known constant when the valve is open or closed,
    // specify aperture by an algebraic equation.
    valve.new_mode(idle1,{dot(aperture1)=0.0_decimal});
    valve.new_mode(opening1,{dot(aperture1)=+1/T});
    valve.new_mode(closing1,{dot(aperture1)=-1/T});

    valve.new_transition(opening1,e_idle,idle1,{next(aperture1)=1.0_dec},aperture1>=1.0_dec,urgent);
    valve.new_transition(closing1,e_idle,idle1,{next(aperture1)=0.0_dec},aperture1<=0.0_dec,urgent);
    valve.new_transition(idle1,e_open,opening1,{next(aperture1)=aperture1});
    valve.new_transition(idle1,e_close,closing1,{next(aperture1)=aperture1});

    return valve;
}
