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

AtomicHybridAutomaton getValve2()
{

    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4.0_decimal);

    // Declare the shared system variables
    RealVariable aperture2("aperture2");

    // Declare the events we use
    DiscreteEvent e_idle("idle2");
    DiscreteEvent e_close("close2");
    DiscreteEvent e_open("open2");

    AtomicHybridAutomaton valve("valve2");

    // Declare the values the valve can variable can have
    AtomicDiscreteLocation idle2("idle2");
    AtomicDiscreteLocation opening2("opening2");
    AtomicDiscreteLocation closing2("closing2");

    // Since aperture is a known constant when the valve is open or closed,
    // specify aperture by an algebraic equation.
    valve.new_mode(idle2,{dot(aperture2)=0.0_decimal});
    valve.new_mode(opening2,{dot(aperture2)=+1/T});
    valve.new_mode(closing2,{dot(aperture2)=-1/T});

    valve.new_transition(opening2,e_idle,idle2,{next(aperture2)=1.0_dec},aperture2>=1.0_dec,urgent);
    valve.new_transition(closing2,e_idle,idle2,{next(aperture2)=0.0_dec},aperture2<=0.0_dec,urgent);
    valve.new_transition(idle2,e_open,opening2,{next(aperture2)=aperture2});
    valve.new_transition(idle2,e_close,closing2,{next(aperture2)=aperture2});

    return valve;
}
