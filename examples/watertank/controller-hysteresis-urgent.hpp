/***************************************************************************
 *            controller-urgent.hpp
 *
 *  Copyright  2018 Luca Geretti
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

AtomicHybridAutomaton getController()
{

    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant hmin("hmin",5.75_decimal);
    RealConstant hmax("hmax",7.75_decimal);

    // Declare the shared system variables
    RealVariable height("height");

    // Declare the events we use
    DiscreteEvent e_can_open("can_open");
    DiscreteEvent e_can_close("can_close");

    AtomicHybridAutomaton controller("controller");

    // Declare the values the valve can variable can have
    AtomicDiscreteLocation rising("rising");
    AtomicDiscreteLocation falling("falling");

    controller.new_mode(rising,List<RealAssignment>());
    controller.new_mode(falling,List<RealAssignment>());

    controller.new_transition(falling,e_can_open,rising,height<=hmin,EventKind::URGENT);
    controller.new_transition(rising,e_can_close,falling,height>=hmax,EventKind::URGENT);

    return controller;
}
