/***************************************************************************
 *            tank.hpp
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

inline HybridAutomaton getTank()
{
    // Declare the system constants
    RealConstant lambda("lambda",0.02_decimal);
    RealConstant rate("rate",0.3_decimal);

    // Declase the variables for the dynamics
    RealVariable aperture("aperture");
    RealVariable height("height");

    // Create the tank automaton; in this case we use the HybridAutomaton class since
    // it is more permissive than AtomicHybridAutomaton; in particular it uses
    // a DiscreteLocation instead of an AtomicDiscreteLocation; a DiscreteLocation allows
    // an empty label, which is cleaner from a logging perspective when there is only one
    // location in an automaton.
    HybridAutomaton tank("tank");

    // Declare a trivial discrete location.
    DiscreteLocation draining;

    // The water level is always given by the same dynamic.
    // The inflow is controlled by the valve aperture, the outflow depends on the
    // pressure, which is proportional to the water height.
    tank.new_mode(draining,{dot(height)=-lambda*height+rate*aperture});

    return tank;
}
