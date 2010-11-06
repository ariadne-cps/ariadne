/***************************************************************************
 *            bouncingball.cc
 *
 *  Copyright  2008  Davide Bresolin
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
#include "ariadne.h"
#include "hybrid_evolver-working.h"

using namespace Ariadne;

int main(int argc, const char* argv[])
{
    uint evolver_verbosity=0;
    if(argc>1) { evolver_verbosity=atoi(argv[1]); }

    typedef GeneralHybridEvolver HybridEvolverType;

    /// Set the system parameters
    Real a = 0.5;  // Coefficient of restitution
    Real g = 9.8;

    /// Set the position and velocity functions.
    ScalarFunction x=ScalarFunction::coordinate(2,0);
    ScalarFunction v=ScalarFunction::coordinate(2,1);

    /// Build the Hybrid System

    /// Create a HybridAutomton object
    MonolithicHybridAutomaton ball;

    /// Create four discrete states
    DiscreteLocation freefall("freefall");
    cout << "location = " << freefall << endl << endl;

    /// Create the discrete events
    DiscreteEvent bounce("bounce");
    cout << "event = " << bounce << endl << endl;

    /// Create the dynamics
    VectorFunction dynamic((v,-g));
    cout << "dynamic = " << dynamic << endl << endl;

    /// Create the resets
    VectorFunction reset((x,-a*v));
    cout << "reset=" << reset << endl << endl;

    /// Create the guards.
    /// Guards are true when g(x) > 0
    ScalarFunction guard(-x);
    cout << "guard=" << guard << endl << endl;

    /// Build the automaton
    ball.new_mode(freefall,dynamic);
    ball.new_transition(bounce,freefall,freefall,reset,guard,impact);
    /// Finished building the automaton

    cout << "Automaton = " << ball << endl << endl;
    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolverType evolver;
    evolver.verbosity=evolver_verbosity;

    TaylorModel::set_default_maximum_degree(5);
    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = 0.05;
    evolver.parameters().maximum_step_size = 1.0/32;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolverType::EnclosureType EnclosureType;
    typedef HybridEvolverType::EnclosureListType EnclosureListType;
    typedef HybridEvolverType::OrbitType OrbitType;

    std::cout << "Computing evolution starting from location l1, x = 2.0, v = 0.0" << std::endl;

    Box initial_box(2, 1.99,2.0, 0.0,0.01);
    EnclosureType initial_enclosure(freefall,initial_box);
    Box bounding_box(2, -0.1,2.1, -10.1,10.1);

    HybridTime evolution_time(1.5,4);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(ball,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    //std::cout << "Orbit="<<orbit<<std::endl;
    std::cout << "Orbit.final()="<<orbit.final()<<std::endl;
    plot("bouncingball-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit);

    //textplot("ball-orbit.txt",orbit);

/*
    std::cout << "Computing reach set using HybridEvolver... " << std::flush;
    EnclosureListType reach = evolver.reach(ball,initial_enclosure,evolution_time);
    std::cout << "done." << std::endl;

    std::cout << "Reach="<<reach<<std::endl;
    //plot("tutorial-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit.initial());
    plot("ball-reach-evolver",bounding_box, Colour(0.0,0.5,1.0), reach);

    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.verbosity = 6;
    analyser.parameters().lock_to_grid_time = 32.0;
    analyser.parameters().maximum_grid_depth= 5;
    std::cout <<  analyser.parameters() << std::endl;

    HybridImageSet initial_set;
    initial_set[l1]=initial_box;

    HybridTime reach_time(4.0,2);

    plot("ball-initial_set1",bounding_box, Colour(0.0,0.5,1.0), initial_set);

    // Compute evolved sets (i.e. at the evolution time) and reach sets (i.e. up to the evolution time) using lower semantics.
    // These functions run a bunch of simulations with bounded approximation errors and combines the results.
    // If the desired evolution time can not be attained without exceeding the error bounds, then the run discarded (without warning)
    std::cout << "Computing lower reach set... " << std::flush;
    HybridGridTreeSet* lower_reach_set_ptr = analyser.lower_reach(ball,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("ball-lower_reach",bounding_box, Colour(0.0,0.5,1.0), *lower_reach_set_ptr);

    // Compute evolved sets and reach sets using upper semantics.
    // These functions compute over-approximations to the evolved and reachabe sets. Subdivision is used
    // as necessary to keep the local errors reasonable. The accumulated global error may be very large.
    std::cout << "Computing upper reach set... " << std::flush;
    HybridGridTreeSet* upper_reach_set_ptr = analyser.upper_reach(ball,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("ball-upper_reach",bounding_box, Colour(0.0,0.5,1.0), *upper_reach_set_ptr);
*/

}
