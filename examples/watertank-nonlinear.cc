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
struct Opening : VectorFunctionData<3,3,3> {
    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p) {
        typename R::value_type x0 = x[0];
        r[0] = - p[0] * Ariadne::sqrt(x0) + p[1] * x[1];
        r[1] = 1.0/p[2];
        r[2] = 1.0;
    }
};

struct Closing : VectorFunctionData<3,3,3> {
    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p) {
        typename R::value_type x0 = x[0];
        r[0] = - p[0] * Ariadne::sqrt(x0) + p[1] * x[1];
        r[1] = -1.0/p[2];
        r[2] = 1.0;
    }
};


struct OpenValve : VectorFunctionData<3,3,3> {
    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p) {
        typename R::value_type x0 = x[0];
        r[0] = - p[0] * Ariadne::sqrt(x0) + p[1];
        r[1] = 0.0;
        r[2] = 1.0;
    }
};

struct ClosedValve : VectorFunctionData<3,3,3> {
    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p) {
        typename R::value_type x0 = x[0];
        r[0] = - p[0] * Ariadne::sqrt(x0);
        r[1] = 0;
        r[2] = 1.0;
    }
};

/// Function for plotting the orbit and reachability set
template<class SET> void plot(const char* filename, const int& xaxis, const int& yaxis, const int& numVariables, const Box& bbox, const Colour& fc, const SET& set, const int& MAX_GRID_DEPTH) {
    // Assigns local variables
    Figure fig;
    array<uint> xy(2,xaxis,yaxis);

    fig.set_projection_map(ProjectionFunction(xy,numVariables));
    fig.set_bounding_box(bbox);

    // If the grid must be shown
    if (MAX_GRID_DEPTH >= 0)
    {
	// The rectangle to be drawn
	Box rect = Box(numVariables);
	// Chooses the fill colour
        fig << fill_colour(Colour(1.0,1.0,1.0));

	// Gets the number of times each variable interval would be divided by 2
        int numDivisions = MAX_GRID_DEPTH / numVariables;
	// Gets the step in the x direction, by 1/2^(numDivisions+h), where h is 1 if the step is to be further divided by 2, 0 otherwise
	double step_x = 1.0/(1 << (numDivisions + ((MAX_GRID_DEPTH - numDivisions*numVariables > xaxis) ? 1 : 0)));
	// Initiates the x position to the bounding box left bound
        double pos_x = bbox[0].lower();
        // Sets the rectangle 2-nd interval to the corresponding bounding box interval (while the >2 intervals are kept at [0,0])
	rect[yaxis] = bbox[1];
        // While between the interval
        while (pos_x < bbox[0].upper())
        {
	    rect[xaxis] = Interval(pos_x,pos_x+step_x); // Sets the rectangle x coordinate
	    pos_x += step_x; // Shifts the x position
	    fig << rect; // Appends the rectangle
        }

	// Repeats for the rectangles in the y direction
	double step_y = 1.0/(1 << (numDivisions + ((MAX_GRID_DEPTH - numDivisions*numVariables > yaxis) ? 1 : 0)));
        double pos_y = bbox[1].lower();
	rect[xaxis] = bbox[0];
        while (pos_y < bbox[1].upper())
        {
	    rect[yaxis] = Interval(pos_y,pos_y+step_y);
   	    fig << rect;
	    pos_y += step_y;
        }
    }
    // Draws and creates file
    fig.set_fill_colour(fc);
    fig << set;
    fig.write(filename);
}


int main()
{

    /// Set the system parameters
    double a = 0.065;
    double b = 0.3;
    double T = 4.0;
    double hmin = 5.5;
    double Delta = 0.05;
    double hmax = 8.0;

    double tmax = 64.0;
    int jmax = 6;
    
    int max_grid_depth=3;

    Vector<Interval> system_parameters(3);
    system_parameters[0] = a;
    system_parameters[1] = b;
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
    VectorUserFunction<Opening> dynamic1(system_parameters);
    VectorUserFunction<OpenValve> dynamic2(system_parameters);
    VectorUserFunction<Closing> dynamic3(system_parameters);
    VectorUserFunction<ClosedValve> dynamic4(system_parameters);

    cout << "dynamic1 = " << dynamic1 << endl << endl;
    cout << "dynamic2 = " << dynamic2 << endl << endl;
    cout << "dynamic3 = " << dynamic3 << endl << endl;
    cout << "dynamic4 = " << dynamic4 << endl << endl;

    /// Create the resets
    VectorAffineFunction reset_y_zero(Matrix<Float>(3,3, 1.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,1.0),Vector<Float>(3, 0.0,0.0,0.0));
    cout << "reset_y_zero=" << reset_y_zero << endl << endl;
    VectorAffineFunction reset_y_one(Matrix<Float>(3,3, 1.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,1.0),Vector<Float>(3, 0.0,1.0,0.0));
    cout << "reset_y_one=" << reset_y_one << endl << endl;

    /// Create the guards.
    /// Guards are true when f(x) = Ax + b > 0
    VectorAffineFunction guard12(Matrix<Float>(1,3,0.0,1.0,0.0),Vector<Float>(1,-1.0));
    cout << "guard12=" << guard12 << endl << endl;
    VectorAffineFunction guard23(Matrix<Float>(1,3,1.0,0.0,0.0),Vector<Float>(1, - hmax + Delta));
    cout << "guard23=" << guard23 << endl << endl;
    VectorAffineFunction guard34(Matrix<Float>(1,3,0.0,-1.0,0.0),Vector<Float>(1,0.0));
    cout << "guard34=" << guard34 << endl << endl;
    VectorAffineFunction guard41(Matrix<Float>(1,3,-1.0,0.0,0.0),Vector<Float>(1,hmin + Delta));
    cout << "guard41=" << guard41 << endl << endl;

    /// Create the invariants.
    /// Invariants are true when f(x) = Ax + b < 0
    /// forced transitions do not need an explicit invariant,
    /// we need only the invariants for location 2 and 4
    VectorAffineFunction inv2(Matrix<Float>(1,3,1.0,0.0,0.0),Vector<Float>(1, - hmax - Delta));//
    cout << "inv2=" << inv2 << endl << endl;
    VectorAffineFunction inv4(Matrix<Float>(1,3,-1.0,0.0,0.0),Vector<Float>(1, hmin - Delta));
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
    evolver.parameters().maximum_enclosure_radius = 50.0;
    evolver.parameters().maximum_step_size = 1.25;
    evolver.parameters().enable_subdivisions = true;
    evolver.verbosity = 1;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;

    std::cout << "Computing evolution starting from location l2, x = 0.0, y = 1.0" << std::endl;

    Box initial_box(3, 1.0000,1.000, 1.00000,1.000, 0.0000,0.000);
    HybridEnclosureType initial_enclosure(l2,initial_box);
    Box bounding_box(3, 0.0,10.0, -0.1,1.1, 0.0,tmax);

    HybridTime evolution_time(tmax,jmax);

    std::cout << "Computing orbit with lower semantics... " << std::flush;
    OrbitType orbit = evolver.orbit(watertank_system,initial_enclosure,evolution_time,LOWER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit.reach.size()="<<orbit.reach().size()<<std::endl;
    std::cout << "Orbit.final.size()="<<orbit.final().size()<<std::endl;
    //plot("tutorial-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit.initial());
    plot("watertank-nonlinear-lower-orbit-tx", 2,0, 3, bounding_box, Colour(0.0,0.5,1.0), orbit, -1);
    plot("watertank-nonlinear-lower-orbit-ty", 2,1, 3, bounding_box, Colour(0.0,0.5,1.0), orbit, -1);
    plot("watertank-nonlinear-lower-orbit-xy", 0,1, 3, bounding_box, Colour(0.0,0.5,1.0), orbit, -1);
    HybridGridTreeSet hgr = outer_approximation(orbit.reach(),watertank_system.grid(),max_grid_depth);
    plot("watertank-nonlinear-lower-orbit-xy-grid", 0,1, 3, bounding_box, Colour(0.0,0.5,1.0), hgr, -1);

    std::cout << "Computing orbit with upper semantics... " << std::flush;
    orbit = evolver.orbit(watertank_system,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit.reach.size()="<<orbit.reach().size()<<std::endl;
    std::cout << "Orbit.final.size()="<<orbit.final().size()<<std::endl;
    //plot("tutorial-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit.initial());
    plot("watertank-nonlinear-upper-orbit-tx", 2,0, 3, bounding_box, Colour(0.0,0.5,1.0), orbit, -1);
    plot("watertank-nonlinear-upper-orbit-ty", 2,1, 3, bounding_box, Colour(0.0,0.5,1.0), orbit, -1);
    plot("watertank-nonlinear-upper-orbit-xy", 0,1, 3, bounding_box, Colour(0.0,0.5,1.0), orbit, -1);
    hgr = outer_approximation(orbit.reach(),watertank_system.grid(),max_grid_depth);
    plot("watertank-nonlinear-upper-orbit-xy-grid", 0,1, 3, bounding_box, Colour(0.0,0.5,1.0), hgr, -1);

    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.parameters().lock_to_grid_time = 100.0;
    analyser.parameters().lock_to_grid_steps = 10;
    analyser.parameters().maximum_grid_depth = max_grid_depth;
    analyser.verbosity=3;
    std::cout <<  analyser.parameters() << std::endl;

    HybridImageSet initial_set;
    initial_set[l2]=initial_box;

    HybridTime reach_time(tmax,jmax);

    // Compute evolved sets (i.e. at the evolution time) and reach sets (i.e. up to the evolution time) using lower semantics.
    // These functions run a bunch of simulations with bounded approximation errors and combines the results.
    // If the desired evolution time can not be attained without exceeding the error bounds, then the run discarded (without warning)
    std::cout << "Computing lower evolve set... " << std::flush;
    HybridGridTreeSet evolve_set = analyser.lower_evolve(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-nonlinear-lower_evolve", 2,0, 3, bounding_box, Colour(0.0,0.5,1.0), evolve_set, -1);

    std::cout << "Computing lower reach set... " << std::flush;
    HybridGridTreeSet reach_set = analyser.lower_reach(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-nonlinear-lower_reach", 2,0, 3, bounding_box, Colour(0.0,0.5,1.0), reach_set, -1);

    std::cout << "Computing lower reach and evolve set... " << std::flush;
    make_lpair(reach_set,evolve_set) = analyser.lower_reach_evolve(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-nonlinear-lower_reach_evolve1", 2,0, 3, bounding_box, Colour(0.0,0.5,1.0), reach_set, -1);
    plot("watertank-nonlinear-lower_reach_evolve2", 2,0, 3, bounding_box, Colour(0.0,0.5,1.0), evolve_set, -1);

    std::cout << "Computing upper evolve set... " << std::flush;
    evolve_set = analyser.upper_evolve(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-nonlinear-upper_evolve", 2,0, 3, bounding_box, Colour(0.0,0.5,1.0), evolve_set, -1);

    std::cout << "Computing upper reach set... " << std::flush;
    reach_set = analyser.upper_reach(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-nonlinear-upper_reach", 2,0, 3, bounding_box, Colour(0.0,0.5,1.0), reach_set, -1);

    std::cout << "Computing upper reach and evolve set... " << std::flush;
    make_lpair(reach_set,evolve_set) = analyser.upper_reach_evolve(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-nonlinear-upper_reach_evolve1", 2,0, 3, bounding_box, Colour(0.0,0.5,1.0), reach_set, -1);
    plot("watertank-nonlinear-upper_reach_evolve2", 2,0, 3, bounding_box, Colour(0.0,0.5,1.0), evolve_set, -1);

    analyser.parameters().lock_to_grid_time = 100.0;
    analyser.parameters().lock_to_grid_steps = 1;


    // bould bounding box to ensure termination fo chain_reach
    HybridBoxes hbox = bounding_boxes(watertank_system.state_space(), bounding_box);

    std::cout << "Computing chain reach set.... " << std::flush;
    reach_set = analyser.chain_reach(watertank_system,initial_set,hbox);
    std::cout << "done." << std::endl;
    plot("watertank-nonlinear-chainreach", 0,1, 3, bounding_box, Colour(0.0,0.5,1.0), reach_set, -1);


}
