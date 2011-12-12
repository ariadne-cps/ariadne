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

using namespace Ariadne;

int main(int argc, const char* argv[])
{
    uint evolver_verbosity=0;
    if(argc>1) { evolver_verbosity=atoi(argv[1]); }

    typedef GeneralHybridEvolver GeneralHybridEvolverType;

    /// Set the system parameters
    Real a = 0.5;  // Coefficient of restitution
    Real g = 9.8;

    /// Set the position and velocity functions.
    RealVariable x("x");
    RealVariable v("v");

    /// Build the Hybrid System

    /// Create a HybridAutomton object
    HybridAutomaton ball;

    /// Create the discrete location
    //DiscreteLocation freefall(StringVariable("ball")|"freefall");
    DiscreteLocation freefall;
    cout << "location = " << freefall << endl << endl;

    /// Create the discrete events
    DiscreteEvent bounce("bounce");
    cout << "event = " << bounce << endl << endl;

    /// Build the automaton
    ball.new_mode(freefall,(dot(x)=v,dot(v)=-g));
    ball.new_guard(freefall,bounce,x<=0,impact);
    ball.new_update(freefall,bounce,freefall,(next(x)=x,next(v)=-a*v));
    /// Finished building the automaton

    cout << "Ball = " << ball << endl << endl;
    /// Compute the system evolution

    /// Create a GeneralHybridEvolver object
    GeneralHybridEvolverType evolver(ball);
    evolver.verbosity=evolver_verbosity;

    /// Set the evolution parameters
    evolver.configuration().set_maximum_enclosure_radius(0.05);
    evolver.configuration().set_maximum_step_size(1.0/32);
    std::cout <<  evolver.configuration() << std::endl;

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolverType::EnclosureType EnclosureType;
    typedef GeneralHybridEvolverType::EnclosureListType EnclosureListType;
    typedef GeneralHybridEvolverType::OrbitType OrbitType;

    std::cout << "Computing evolution starting from location l1, x = 2.0, v = 0.0" << std::endl;

    Box bounding_box(2, -0.1,2.1, -10.1,10.1);

    HybridExpressionSet initial_set(freefall,(2.0<=x<=2.0,v.in(0.0,0.0)));
    HybridTime evolution_time(1.5,4);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(initial_set,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    //std::cout << "Orbit="<<orbit<<std::endl;
    std::cout << "Orbit.final()="<<orbit.final()<<std::endl;
    plot("bouncingball-orbit",Axes2d(-0.1,x,2.1, -10.1,v,10.1), Colour(0.0,0.5,1.0), orbit);
    plot("bouncingball-x",Axes2d(0.0,TimeVariable(),1.5,- 0.1,x,2.1), Colour(0.0,0.5,1.0), orbit);

    //textplot("ball-orbit.txt",orbit);

/*
    std::cout << "Computing reach set using GeneralHybridEvolver... " << std::flush;
    EnclosureListType reach = evolver.reach(initial_enclosure,evolution_time);
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
