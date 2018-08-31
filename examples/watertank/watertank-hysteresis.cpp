/***************************************************************************
 *            watertank-hysteresis.cpp
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
#include "valve-noalgebraic-hysteresis.hpp"
#include "controller-hysteresis-permissive.hpp"

using namespace Ariadne;
using std::cout; using std::endl;

Int main(Int argc, const char* argv[])
{
    Nat evolver_verbosity=get_verbosity(argc,argv);

    // Declare the shared system variables
    RealVariable aperture("aperture");
    RealVariable height("height");

    StringVariable valve("valve");
    StringVariable controller("controller");

    StringConstant opened("opened");
    StringConstant opening("opening");
    StringConstant idle("idle");
    StringConstant rising("rising");

    HybridAutomaton tank_automaton = getTank();
    HybridAutomaton valve_automaton = getValve();
    HybridAutomaton controller_automaton = getController();
    CompositeHybridAutomaton watertank_system({tank_automaton,valve_automaton,controller_automaton});

    cout << watertank_system << endl;
    // Compute the system evolution

    // Create a GeneralHybridEvolver object
    GeneralHybridEvolver evolver(watertank_system);
    evolver.verbosity = evolver_verbosity;

    // Set the evolution parameters
    evolver.configuration().set_maximum_enclosure_radius(3.05);
    evolver.configuration().set_maximum_step_size(0.6);

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolver::OrbitType OrbitType;

    std::cout << "Computing evolution... " << std::flush;

    //HybridSet initial_set({valve|idle,controller|rising},{height==7,aperture==1});

    //HybridSet initial_set({valve|opening,controller|rising},{height==5.5_dec,aperture==0.4_dec});
    //HybridSet initial_set({valve|opened,controller|rising},{height==7});
    HybridSet initial_set({valve|opened,controller|rising},{height==7,aperture==1});
    HybridTime evolution_time(30.0,5);
    OrbitType orbit = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    std::cout << "done." << std::endl;

    std::cout << "Plotting trajectory... "<<std::flush;
    Axes2d time_height_axes(0<=TimeVariable()<=30,-0.1<=height<=9.1);
    plot("watertank-orbit",time_height_axes, Colour(0.0,0.5,1.0), orbit);
    Axes2d height_aperture_axes(-0.1,height,9.1, -0.1,aperture,1.3);
    plot("watertank-height_aperture",height_aperture_axes, Colour(0.0,0.5,1.0), orbit);
    std::cout << "done." << std::endl;

    /*
    HybridReachabilityAnalyser analyser(watertank_system,evolver);

    HybridBoxes bounding_domain;
    Box continuous_domain_2d(2,-1.0,10,0,-1.0,2.0);
    Box continuous_domain_1d(1,-1.0,10.0);
    bounding_domain.insert(valve|opening,watertank_system.state_space()[valve|opening],continuous_domain_2d);
    bounding_domain.insert(valve|closing,watertank_system.state_space()[valve|closing],continuous_domain_2d);
    bounding_domain.insert(valve|opened,watertank_system.state_space()[valve|opened],continuous_domain_1d);
    bounding_domain.insert(valve|closed,watertank_system.state_space()[valve|closed],continuous_domain_1d);
    analyser.configuration().set_bounding_domain(bounding_domain);
    */
    /*
    std::shared_ptr<HybridGrid> grid(new HybridGrid(watertank_system.state_auxiliary_space()));
    analyser.configuration().set_grid(grid);
    analyser.configuration().set_maximum_grid_depth(3);
    analyser.configuration().set_lock_to_grid_steps(2);
    analyser.configuration().set_lock_to_grid_time(80.0);
	*/
/*
    std::cout << "Computing upper reach... " << std::flush;
    HybridGridTreeSet upper_reach = analyser.upper_reach(initial_set,evolution_time);
    std::cout << "Plotting reachable sets... " << std::flush;
    plot("watertank-upper-reach", height_aperture_axes, Colour(0.0,0.5,1.0), upper_reach);
*/

    std::cout << "Discretising orbit" << std::flush;
    HybridGrid grid(watertank_system.state_auxiliary_space());
    HybridGridTreePaving hgts(grid);

    for (ListSet<HybridEnclosure>::ConstIterator it = orbit.reach().begin(); it != orbit.reach().end(); it++)
    {
        std::cout<<"."<<std::flush;
        it->state_auxiliary_set().adjoin_outer_approximation_to(hgts,4);
    }
    std::cout << "done." << std::endl;

    plot("watertank-reach", height_aperture_axes, Colour(0.0,0.5,1.0), hgts);

}
