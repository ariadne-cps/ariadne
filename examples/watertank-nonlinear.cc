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

//
// Definition of the dynamics
//
// the parameters of the functions are a, b and T
//
struct Opening : FunctionData<2,2,3> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
        typename R::value_type x0 = x[0];
        r[0] = - p[0] * Ariadne::sqrt(x0) + p[1] * x[1]; 
        r[1] = 1.0/p[2];        
    }
};

struct Closing : FunctionData<2,2,3> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
        typename R::value_type x0 = x[0];
        r[0] = - p[0] * Ariadne::sqrt(x0) + p[1] * x[1]; 
        r[1] = -1.0/p[2];        
    }
};


struct OpenValve : FunctionData<2,2,3> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
        typename R::value_type x0 = x[0];
        r[0] = - p[0] * Ariadne::sqrt(x0) + p[1]; 
        r[1] = 0.0;        
    }
};

struct ClosedValve : FunctionData<2,2,3> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
        typename R::value_type x0 = x[0];
        r[0] = - p[0] * Ariadne::sqrt(x0); 
        r[1] = 0;        
    }
};


int main() 
{
  
    /// Set the system parameters
    double a = 0.065;
    double b = 0.3;
    double T = 4.0;
    double hmin = 5.5;
    double Delta = 0.05;
    double hmax = 8.0;
 
    Vector<Interval> system_parameters(3);
    system_parameters[0] = a;
    system_parameters[1] = Interval(0.3,0.34);
    system_parameters[2] = T;


    /// Build the Hybrid System
  
    /// Create a HybridAutomton object
    HybridAutomaton watertank_system;
  
    /// Create four discrete states
    DiscreteState l1(1);
    DiscreteState l2(2);
    DiscreteState l3(3);
    DiscreteState l4(4);
  
    /// Create the discrete events
    DiscreteEvent e12(12);
    DiscreteEvent e23(23);
    DiscreteEvent e34(34);
    DiscreteEvent e41(41);
  
    /// Create the dynamics
    Function<Opening> dynamic1(system_parameters);
    Function<OpenValve> dynamic2(system_parameters);
    Function<Closing> dynamic3(system_parameters);
    Function<ClosedValve> dynamic4(system_parameters);
    
    cout << "dynamic1 = " << dynamic1 << endl << endl;
    cout << "dynamic2 = " << dynamic2 << endl << endl;
    cout << "dynamic3 = " << dynamic3 << endl << endl;
    cout << "dynamic4 = " << dynamic4 << endl << endl;

    /// Create the resets
    AffineFunction reset_y_zero(Matrix<Float>(2,2,1.0,0.0,0.0,0.0),Vector<Float>(2,0.0,0.0));
    cout << "reset_y_zero=" << reset_y_zero << endl << endl;
    AffineFunction reset_y_one(Matrix<Float>(2,2,1.0,0.0,0.0,0.0),Vector<Float>(2,0.0,1.0));
    cout << "reset_y_one=" << reset_y_one << endl << endl;

    /// Create the guards.
    /// Guards are true when f(x) = Ax + b > 0
    AffineFunction guard12(Matrix<Float>(1,2,0.0,1.0),Vector<Float>(1,-1.0));
    cout << "guard12=" << guard12 << endl << endl;
    AffineFunction guard23(Matrix<Float>(1,2,1.0,0.0),Vector<Float>(1, - hmax + Delta));
    cout << "guard23=" << guard23 << endl << endl;
    AffineFunction guard34(Matrix<Float>(1,2,0.0,-1.0),Vector<Float>(1,0.0));
    cout << "guard34=" << guard34 << endl << endl;
    AffineFunction guard41(Matrix<Float>(1,2,-1.0,0.0),Vector<Float>(1,hmin + Delta));
    cout << "guard41=" << guard41 << endl << endl;

    /// Create the invariants.
    /// Invariants are true when f(x) = Ax + b < 0
    /// forced transitions do not need an explicit invariant, 
    /// we need only the invariants for location 2 and 4
    AffineFunction inv2(Matrix<Float>(1,2,1.0,0.0),Vector<Float>(1, - hmax - Delta));//
    cout << "inv2=" << inv2 << endl << endl;
    AffineFunction inv4(Matrix<Float>(1,2,-1.0,0.0),Vector<Float>(1, hmin - Delta));
    cout << "inv4=" << inv4 << endl << endl;
  
    /// Build the automaton
    watertank_system.new_mode(l1,dynamic1);
    watertank_system.new_mode(l2,dynamic2);
    watertank_system.new_mode(l3,dynamic3);
    watertank_system.new_mode(l4,dynamic4);

    watertank_system.new_invariant(l2,inv2);
    watertank_system.new_invariant(l4,inv4);

    watertank_system.new_forced_transition(e12,l1,l2,reset_y_one,guard12);
    watertank_system.new_unforced_transition(e23,l2,l3,reset_y_one,guard23);
    watertank_system.new_forced_transition(e34,l3,l4,reset_y_zero,guard34);
    watertank_system.new_unforced_transition(e41,l4,l1,reset_y_zero,guard41);


    /// Finished building the automaton

    cout << "Automaton = " << watertank_system << endl << endl;

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver;

    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = 0.25;
    evolver.parameters().maximum_step_size = 0.125;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;

    std::cout << "Computing evolution starting from location l2, x = 0.0, y = 1.0" << std::endl;

    Box initial_box(2, 0.001,0.002, 1.0,1.001);
    HybridEnclosureType initial_enclosure(l2,initial_box);
    Box bounding_box(2, -0.1,9.1, -0.1,1.1);
  
    HybridTime evolution_time(64.0,6);
  
    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(watertank_system,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit="<<orbit<<std::endl;
    //plot("tutorial-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit.initial());
    plot("watertank-nonlinear-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit);

/*
    std::cout << "Computing reach set using HybridEvolver... " << std::flush;
    EnclosureListType reach = evolver.reach(watertank_system,initial_enclosure,evolution_time);
    std::cout << "done." << std::endl;

    std::cout << "Orbit="<<reach<<std::endl;
    //plot("tutorial-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit.initial());
    plot("watertank-nonlinear-reach-evolver",bounding_box, Colour(0.0,0.5,1.0), reach);


    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.parameters().lock_to_grid_time = 32.0;
    analyser.parameters().grid_lengths = 0.05;
    std::cout <<  analyser.parameters() << std::endl;

    HybridImageSet initial_set;
    initial_set[l2]=initial_box;

    HybridTime reach_time(64.0,2);

    plot("watertank-nonlinear-initial_set1",bounding_box, Colour(0.0,0.5,1.0), initial_set);

    // Compute evolved sets (i.e. at the evolution time) and reach sets (i.e. up to the evolution time) using lower semantics.
    // These functions run a bunch of simulations with bounded approximation errors and combines the results.
    // If the desired evolution time can not be attained without exceeding the error bounds, then the run discarded (without warning)
    std::cout << "Computing lower reach set... " << std::flush;
    HybridGridTreeSet* lower_reach_set_ptr = analyser.lower_reach(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-nonlinear-lower_reach1",bounding_box, Colour(0.0,0.5,1.0), *lower_reach_set_ptr);

    // Compute evolved sets and reach sets using upper semantics.
    // These functions compute over-approximations to the evolved and reachabe sets. Subdivision is used
    // as necessary to keep the local errors reasonable. The accumulated global error may be very large.
    std::cout << "Computing upper reach set... " << std::flush;
    HybridGridTreeSet* upper_reach_set_ptr = analyser.upper_reach(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-nonlinear-upper_reach1",bounding_box, Colour(0.0,0.5,1.0), *upper_reach_set_ptr);

    std::cout << "Computing evolution starting from location l1, x = 0.0, y = 0.0" << std::endl;

    Box initial_box2(2, 0.0,0.001, 0.0,0.001);
    HybridImageSet initial_set2;
    initial_set2[l1]=initial_box2;

    plot("watertank-nonlinear-initial_set2",bounding_box, Colour(0.0,0.5,1.0), initial_set2);

    // Compute evolved sets (i.e. at the evolution time) and reach sets (i.e. up to the evolution time) using lower semantics.
    // These functions run a bunch of simulations with bounded approximation errors and combines the results.
    // If the desired evolution time can not be attained without exceeding the error bounds, then the run discarded (without warning)
    std::cout << "Computing lower reach set... " << std::flush;
    lower_reach_set_ptr = analyser.lower_reach(watertank_system,initial_set2,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-nonlinear-lower_reach2",bounding_box, Colour(0.0,0.5,1.0), *lower_reach_set_ptr);

    // Compute evolved sets and reach sets using upper semantics.
    // These functions compute over-approximations to the evolved and reachabe sets. Subdivision is used
    // as necessary to keep the local errors reasonable. The accumulated global error may be very large.
    std::cout << "Computing upper reach set... " << std::flush;
    upper_reach_set_ptr = analyser.upper_reach(watertank_system,initial_set2,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-nonlinear-upper_reach2",bounding_box, Colour(0.0,0.5,1.0), *upper_reach_set_ptr);

*/

}
