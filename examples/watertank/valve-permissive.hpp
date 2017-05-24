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

AtomicHybridAutomaton getValve()
{

    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4.0_decimal);
    RealConstant hmin("hmin",5.5_decimal);
    RealConstant hmax("hmax",8.0_decimal);
    RealConstant delta("delta",0.05_decimal);

    // Declare the shared system variables
    RealVariable aperture("aperture");
    RealVariable height("height");

    // Declare the events we use
    DiscreteEvent start_opening("start_opening");
    DiscreteEvent start_closing("start_closing");
    DiscreteEvent finished_opening("finished_opening");
    DiscreteEvent finished_closing("finished_closing");
    DiscreteEvent must_start_opening("must_start_opening");
    DiscreteEvent must_start_closing("must_start_closing");

    AtomicHybridAutomaton valve("valve");

    // Declare the values the valve can variable can have
    AtomicDiscreteLocation open("open");
    AtomicDiscreteLocation opening("opening");
    AtomicDiscreteLocation closed("closed");
    AtomicDiscreteLocation closing("closing");

    // Since aperture is a known constant when the valve is open or closed,
    // specify aperture by an algebraic equation.
    valve.new_mode(open,{let(aperture)=+1.0_decimal});
    valve.new_mode(closed,{let(aperture)=0.0_decimal});
    // Specify the differential equation for how the valve opens/closes.
    valve.new_mode(opening,{dot(aperture)=+1/T});
    valve.new_mode(closing,{dot(aperture)=-1/T});

    // Specify the invariants valid in each mode. Note that every invariant
    // must have an action label. This is used internally, for example, to
    // check non-blockingness of urgent actions.
    valve.new_invariant(open,height<=hmax+delta,must_start_closing);
    valve.new_invariant(closed,height>=hmin-delta,must_start_opening);

    valve.new_transition(closed,start_opening,opening,{next(aperture)=aperture},height<=hmin+delta,permissive);
    valve.new_transition(opening,finished_opening,open,aperture>=1,urgent);
    valve.new_transition(open,start_closing,closing,{next(aperture)=aperture},height>=hmax-delta,permissive);
    valve.new_transition(closing,finished_closing,closed,aperture<=0,urgent);

    return valve;
}
