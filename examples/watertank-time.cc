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


int main()
{

    /// Set the system parameters
    Real a = -0.02;
    Real b = 0.3;
    Real T = 1.5;
    Real hmin = 5.55;
    Real Delta = 0.05;
    Real hmax = 5.70;
    Real tmax = 20.0;
    Int dmax = 15;

    /// Build the Hybrid System

    /// Create a HybridAutomton object
    MonolithicHybridAutomaton watertank_system;

    /// Create four discrete states
    DiscreteLocation l1("q1");
    DiscreteLocation l2("q2");
    DiscreteLocation l3("q4");
    DiscreteLocation l4("q5");

    /// Create the discrete events
    DiscreteEvent e12("e12");
    DiscreteEvent e13("e13");
    DiscreteEvent e23("e23");
    DiscreteEvent e34("e34");
    DiscreteEvent e31("e31");
    DiscreteEvent e41("e41");
    DiscreteEvent i2("i2");
    DiscreteEvent i4("i4");

    /// Coordinates
    RealScalarFunction x=RealScalarFunction::coordinate(3,0);
    RealScalarFunction y=RealScalarFunction::coordinate(3,1);
    RealScalarFunction t=RealScalarFunction::coordinate(3,2);

    /// Create the dynamics
    RealVectorFunction dynamic1((a*x+b*y,1/T,1));
    RealVectorFunction dynamic2((a*x+b,0,1));
    RealVectorFunction dynamic3((a*x+b*y,-1/T,1));
    RealVectorFunction dynamic4((a*x,0,1));

    cout << "dynamic1 = " << dynamic1 << endl << endl;
    cout << "dynamic2 = " << dynamic2 << endl << endl;
    cout << "dynamic3 = " << dynamic3 << endl << endl;
    cout << "dynamic4 = " << dynamic4 << endl << endl;

    /// Create the resets
    RealVectorFunction reset_id((x,y,t));
    RealVectorFunction reset_y_zero((x,0,t));
    cout << "reset_y_zero=" << reset_y_zero << endl << endl;
    RealVectorFunction reset_y_one((x,1,t));
    cout << "reset_y_one=" << reset_y_one << endl << endl;

    /// Create the guards.
    /// Guards are true when f(x) > 0
    RealScalarFunction guard12(y-1);
    cout << "guard12=" << guard12 << endl << endl;
    RealScalarFunction guard23(x-(hmax-Delta));
    cout << "guard23=" << guard23 << endl << endl;
    RealScalarFunction guard34(-y);
    cout << "guard34=" << guard34 << endl << endl;
    RealScalarFunction guard41(-x+(hmin+Delta));
    cout << "guard41=" << guard41 << endl << endl;

    /// Create the invariants.
    /// Invariants are true when f(x) < 0
    /// forced transitions do not need an explicit invariant
    /// x < hmax + Delta
    RealScalarFunction inv2(x-(hmax+Delta));
    cout << "inv2=" << inv2 << endl << endl;
    /// x > hmin - Delta
    RealScalarFunction inv4(-x+(hmin-Delta));
    cout << "inv4=" << inv4 << endl << endl;

    /// Build the automaton
    watertank_system.new_mode(l1,dynamic1);
    watertank_system.new_mode(l2,dynamic2);
    watertank_system.new_mode(l3,dynamic3);
    watertank_system.new_mode(l4,dynamic4);

    watertank_system.new_invariant(l1,i2,inv2);
    watertank_system.new_invariant(l2,i2,inv2);
    watertank_system.new_invariant(l3,i4,inv4);
    watertank_system.new_invariant(l4,i4,inv4);

    watertank_system.new_transition(l1,e12,l2,reset_y_one,guard12,urgent);
    watertank_system.new_transition(l1,e13,l3,reset_id,guard23,permissive);
    watertank_system.new_transition(l2,e23,l3,reset_y_one,guard23,permissive);
    watertank_system.new_transition(l3,e34,l4,reset_y_zero,guard34,urgent);
    watertank_system.new_transition(l3,e31,l1,reset_id,guard41,permissive);
    watertank_system.new_transition(l4,e41,l1,reset_y_zero,guard41,permissive);


    /// Finished building the automaton

    cout << "Automaton = " << watertank_system << endl << endl;

    /// Compute the system evolution

    /// Create a GeneralHybridEvolver object
    GeneralHybridEvolver evolver;

    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = 0.1;
    evolver.parameters().maximum_step_size = 0.1;
    evolver.verbosity=1;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolver::EnclosureType HybridEnclosureType;
    typedef GeneralHybridEvolver::OrbitType OrbitType;
    typedef GeneralHybridEvolver::EnclosureListType EnclosureListType;

    std::cout << "Computing evolution starting from location l2, x = 5.0, y = 1.0, t = 0.0" << std::endl;

    RealVariableBox initial_box((5.0<=RealVariable("x")<=5.001, 1.0<=RealVariable("y")<=1.001, 0.0<=RealVariable("t")<=0.001));
    HybridSet initial_set(l2,initial_box);
    Box bounding_box(3, -0.1,9.1, -0.1,1.1, -0.1,numeric_cast<double>(tmax)+0.1);

    HybridTime evolution_time(tmax,dmax);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(watertank_system,initial_set,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit="<<orbit<<std::endl;
    Figure g;
    Box graphic_box(2, -0.1,numeric_cast<double>(tmax)+0.1, 4.1,6.1);
    g.set_bounding_box(graphic_box);
    Array<uint> p(2,2,0);
    g.set_projection_map(PlanarProjectionMap(3,2,0));

    g << fill_colour(Colour(0.0,0.5,1.0));
    g << orbit;
    g.write("watertank-time-orbit");
/*
    std::cout << "Computing reach set using GeneralHybridEvolver... " << std::flush;
    EnclosureListType reach = evolver.reach(watertank_system,initial_enclosure,evolution_time);
    std::cout << "done." << std::endl;

    std::cout << "Orbit="<<reach<<std::endl;
    //plot("tutorial-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit.initial());
    plot("watertank-reach-evolver",bounding_box, Colour(0.0,0.5,1.0), reach);


    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.parameters().lock_to_grid_time = 10.0;
    analyser.parameters().maximum_grid_depth = 5;
    std::cout <<  analyser.parameters() << std::endl;

    HybridImageSet initial_set;
    initial_set[l2]=initial_box;

    HybridTime reach_time(tmax,dmax);

    // Compute evolved sets (i.e. at the evolution time) and reach sets (i.e. up to the evolution time) using lower semantics.
    // These functions run a bunch of simulations with bounded approximation errors and combines the results.
    // If the desired evolution time can not be attained without exceeding the error bounds, then the run discarded (without warning)
    std::cout << "Computing lower reach set... " << std::flush;
    HybridGridTreeSet lower_reach_set = analyser.lower_reach(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;

    Figure g;
    Box graphic_box(2, -0.1,tmax+0.1, 4.1,6.1);
    g.set_bounding_box(graphic_box);
    Array<uint> p(2,2,0);
    g.set_projection_map(ProjectionFunction(2,3,p));

    g << fill_colour(Colour(0.0,0.5,1.0));
    g << lower_reach_set;
    g.write("watertank-time-lower");
*/

}
