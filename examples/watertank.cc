/***************************************************************************
 *            watertank.cc
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
    uint evolver_verbosity = 0;
    if(argc>1) { evolver_verbosity=atoi(argv[1]); }

    typedef GeneralHybridEvolver GeneralHybridEvolverType;

    DRAWING_METHOD = AFFINE_DRAW;
    DRAWING_ACCURACY = 0;

    /// Set the system parameters
    Real a = -0.02;
    Real b = 0.3;
    Real T = 4.0;
    Real hmin = 5.5;
    Real Delta = 0.05;
    Real hmax = 8.0;

    /// Build the Hybrid System

    /// Create a HybridAutomton object
    MonolithicHybridAutomaton watertank_system;

    /// Create four discrete states
    AtomicDiscreteLocation l1(1);
    AtomicDiscreteLocation l2(2);
    AtomicDiscreteLocation l3(3);
    AtomicDiscreteLocation l4(4);

    /// Create the discrete events
    DiscreteEvent e12(12);
    DiscreteEvent e23(23);
    DiscreteEvent i2("i2");
    DiscreteEvent e34(34);
    DiscreteEvent e41(41);
    DiscreteEvent i4("i4");

    // Create coordinate functions in two variables.
    ScalarFunction x0=ScalarFunction::coordinate(2,0);
    ScalarFunction x1=ScalarFunction::coordinate(2,1);

    /// Create the dynamics
    VectorFunction dynamic1((a*x0+b*x1,1/T));
    cout << "dynamic1 = " << dynamic1 << endl << endl;
    VectorFunction dynamic2((a*x0+b,0));
    cout << "dynamic2 = " << dynamic2 << endl << endl;
    VectorFunction dynamic3((a*x0+b*x1,-1/T));
    cout << "dynamic3 = " << dynamic3 << endl << endl;
    VectorFunction dynamic4((a*x0,0));
    cout << "dynamic4 = " << dynamic4 << endl << endl;

    /// Create the resets
    VectorFunction reset_y_zero((x0,0));
    cout << "reset_y_zero=" << reset_y_zero << endl << endl;
    VectorFunction reset_y_one((x0,1));
    cout << "reset_y_one=" << reset_y_one << endl << endl;

    /// Create the guards.
    /// Guards are true when g(x) >= 0
    ScalarFunction guard12(x1-1);
    cout << "guard12=" << guard12 << endl << endl;
    ScalarFunction guard23(x0+(-hmax+Delta));
    cout << "guard23=" << guard23 << endl << endl;
    ScalarFunction guard34(-x1);
    cout << "guard34=" << guard34 << endl << endl;
    ScalarFunction guard41(-x0+(hmin+Delta));
    cout << "guard41=" << guard41 << endl << endl;

    /// Create the invariants.
    /// Invariants are true when c(x) <= 0
    /// Urgent transitions do not need an explicit invariant,
    /// we need only the invariants for location 2 and 4
    ScalarFunction inv2(x0+(-hmax - Delta));
    cout << "inv2=" << inv2 << endl << endl;
    ScalarFunction inv4(-x0+(hmin - Delta));
    cout << "inv4=" << inv4 << endl << endl;

    /// Build the automaton
    watertank_system.new_mode(l1,dynamic1);
    watertank_system.new_mode(l2,dynamic2);
    watertank_system.new_mode(l3,dynamic3);
    watertank_system.new_mode(l4,dynamic4);

    watertank_system.new_invariant(l2,i2,inv2);
    watertank_system.new_invariant(l4,i4,inv4);

    watertank_system.new_transition(e12,l1,l2,reset_y_one,guard12,urgent);
    watertank_system.new_transition(e23,l2,l3,reset_y_one,guard23,permissive);
    watertank_system.new_transition(e34,l3,l4,reset_y_zero,guard34,urgent);
    watertank_system.new_transition(e41,l4,l1,reset_y_zero,guard41,permissive);

    /// Finished building the automaton

    cout << "Automaton = " << watertank_system << endl << endl;


    /// Compute the system evolution

    /// Create a GeneralHybridEvolver object
    GeneralHybridEvolverType evolver;
    evolver.verbosity = evolver_verbosity;

    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = 0.25;
    evolver.parameters().maximum_step_size = 2.5;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolverType::EnclosureType EnclosureType;
    typedef GeneralHybridEvolverType::EnclosureListType EnclosureListType;
    typedef GeneralHybridEvolverType::OrbitType OrbitType;

    std::cout << "Computing evolution starting from location l1, x = 0.0, y = 0.0" << std::endl;

    Box initial_box(2, 0.0,0.001, 0.0,0.001);
    EnclosureType initial_enclosure(l1,initial_box);
    Box bounding_box(2, -0.1,9.1, -0.1,1.3);

    //HybridTime evolution_time(90.0,6);
    HybridTime evolution_time(80.0,10);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(watertank_system,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit.final_size="<<orbit.final().size()<<std::endl;
    std::cout << "Orbit.reach_size="<<orbit.reach().size()<<std::endl;
    //plot("tutorial-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit.initial());
    std::cout << "Plotting orbit... "<<std::flush;
    plot("watertank-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit);
    std::cout << "done." << std::endl;

/*
    std::cout << "Computing reach set using GeneralHybridEvolver... " << std::flush;
    EnclosureListType reach = evolver.reach(watertank_system,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Plotting reach set... "<<std::flush;
    //plot("tutorial-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit.initial());
    plot("watertank-reach-evolver",bounding_box, Colour(0.0,0.5,1.0), reach);


    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.parameters().lock_to_grid_time = 32.0;
    analyser.parameters().grid_lengths = 0.05;
    std::cout <<  analyser.parameters() << std::endl;

    HybridImageSet initial_set;
    initial_set[l2]=initial_box;

    HybridTime reach_time(64.0,2);

    plot("watertank-initial_set1",bounding_box, Colour(0.0,0.5,1.0), initial_set);

    // Compute evolved sets (i.e. at the evolution time) and reach sets (i.e. up to the evolution time) using lower semantics.
    // These functions run a bunch of simulations with bounded approximation errors and combines the results.
    // If the desired evolution time can not be attained without exceeding the error bounds, then the run discarded (without warning)
    std::cout << "Computing lower reach set... " << std::flush;
    HybridGridTreeSet* lower_reach_set_ptr = analyser.lower_reach(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-lower_reach1",bounding_box, Colour(0.0,0.5,1.0), *lower_reach_set_ptr);

    // Compute evolved sets and reach sets using upper semantics.
    // These functions compute over-approximations to the evolved and reachabe sets. Subdivision is used
    // as necessary to keep the local errors reasonable. The accumulated global error may be very large.
    std::cout << "Computing upper reach set... " << std::flush;
    HybridGridTreeSet* upper_reach_set_ptr = analyser.upper_reach(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-upper_reach1",bounding_box, Colour(0.0,0.5,1.0), *upper_reach_set_ptr);

    std::cout << "Computing evolution starting from location l1, x = 0.0, y = 0.0" << std::endl;

    Box initial_box2(2, 0.0,0.001, 0.0,0.001);
    HybridImageSet initial_set2;
    initial_set2[l1]=initial_box2;

    plot("watertank-initial_set2",bounding_box, Colour(0.0,0.5,1.0), initial_set2);

    // Compute evolved sets (i.e. at the evolution time) and reach sets (i.e. up to the evolution time) using lower semantics.
    // These functions run a bunch of simulations with bounded approximation errors and combines the results.
    // If the desired evolution time can not be attained without exceeding the error bounds, then the run discarded (without warning)
    std::cout << "Computing lower reach set... " << std::flush;
    lower_reach_set_ptr = analyser.lower_reach(watertank_system,initial_set2,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-lower_reach2",bounding_box, Colour(0.0,0.5,1.0), *lower_reach_set_ptr);

    // Compute evolved sets and reach sets using upper semantics.
    // These functions compute over-approximations to the evolved and reachabe sets. Subdivision is used
    // as necessary to keep the local errors reasonable. The accumulated global error may be very large.
    std::cout << "Computing upper reach set... " << std::flush;
    upper_reach_set_ptr = analyser.upper_reach(watertank_system,initial_set2,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-upper_reach2",bounding_box, Colour(0.0,0.5,1.0), *upper_reach_set_ptr);
*/

}
