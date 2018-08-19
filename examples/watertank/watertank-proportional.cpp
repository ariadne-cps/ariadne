/***************************************************************************
 *            watertank-proportional.cpp
 *
 *  Copyright  2017  Luca Geretti
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
#include "tank.hpp"
#include "valve-proportional-urgent.hpp"

using namespace Ariadne;
using std::cout; using std::endl;

Int main(Int argc, const char* argv[])
{
    Nat evolver_verbosity = 0;
    if(argc>1) { evolver_verbosity=atoi(argv[1]); }

    // Declare the shared system variables
    RealVariable aperture("aperture");
    RealVariable height("height");

    StringVariable valve("valve");
    StringConstant opened("opened");
    StringConstant modulated("modulated");
    StringConstant closed("closed");

    HybridAutomaton tank_automaton = getTank();
    HybridAutomaton valve_automaton = getValve();
    CompositeHybridAutomaton watertank_system({tank_automaton,valve_automaton});

    cout << watertank_system << endl;
    // Compute the system evolution

    // Create a GeneralHybridEvolver object
    GeneralHybridEvolver evolver(watertank_system);
    evolver.verbosity = evolver_verbosity;

    // Set the evolution parameters
    evolver.configuration().set_maximum_enclosure_radius(3.05);
    evolver.configuration().set_maximum_step_size(0.25);

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolver::OrbitType OrbitType;

    std::cout << "Computing evolution... " << std::flush;
    Real a_max(0.0);
    
    //HybridSet initial_set({valve|closed},{0<=height<=0.0_decimal});
    HybridSet initial_set({valve|opened},{0.0_decimal<=height<=0.0_decimal});
    HybridTime evolution_time(80.0,5);
    OrbitType orbit = evolver.orbit(initial_set,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Plotting trajectory... "<<std::flush;
    Axes2d time_height_axes(0<=TimeVariable()<=80,-0.1<=height<=9.1);
    plot("watertank_proportional_t-height",time_height_axes, Colour(0.0,0.5,1.0), orbit);
    Axes2d height_aperture_axes(-0.1,height,9.1, -0.1,aperture,1.3);
    plot("watertank_proportional_height-aperture",height_aperture_axes, Colour(0.0,0.5,1.0), orbit);
    std::cout << "done." << std::endl;


}
