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
	Float step_x = 1.0/(1 << (numDivisions + ((MAX_GRID_DEPTH - numDivisions*numVariables > xaxis) ? 1 : 0)));
	// Initiates the x position to the bounding box left bound
        Float pos_x = bbox[0].lower();
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
	Float step_y = 1.0/(1 << (numDivisions + ((MAX_GRID_DEPTH - numDivisions*numVariables > yaxis) ? 1 : 0)));
        Float pos_y = bbox[1].lower();
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

    double tmax = 80.0;
    int jmax = 6;


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
    DiscreteEvent e34(34);
    DiscreteEvent e41(41);

    /// Coordinates
    ScalarFunction x=ScalarFunction::coordinate(3,0);
    ScalarFunction y=ScalarFunction::coordinate(3,1);
    ScalarFunction t=ScalarFunction::coordinate(3,2);

    /// Create the dynamics

    ScalarFunction zero=ScalarFunction::constant(3,0.0);
    ScalarFunction one=ScalarFunction::constant(3,1.0);

    /// The dynamic for an opening valve
    VectorFunction opening_d((-a*Ariadne::sqrt(x)+b*y,one/T,one));
    /// The dynamic for a closing valve
    VectorFunction closing_d((-a*Ariadne::sqrt(x)+b*y,-one/T,one));
    /// The dynamic for an opened valve
    VectorFunction opened_d((-a*Ariadne::sqrt(x)+b,zero,one));
    /// The dynamic for an opened valve
    VectorFunction closed_d((-a*Ariadne::sqrt(x),zero,one));

    cout << "opening dynamic = " << opening_d << endl << endl;
    cout << "closing dynamic = " << closing_d << endl << endl;
    cout << "opened dynamic = " << opened_d << endl << endl;
    cout << "closed dynamic = " << closed_d << endl << endl;

    /// Create the resets
    VectorAffineFunction reset_y_zero(Matrix<Real>(3,3, 1.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,1.0),Vector<Real>(3, 0.0,0.0,0.0));
    cout << "reset_y_zero=" << reset_y_zero << endl << endl;
    VectorAffineFunction reset_y_one(Matrix<Real>(3,3, 1.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,1.0),Vector<Real>(3, 0.0,1.0,0.0));
    cout << "reset_y_one=" << reset_y_one << endl << endl;

    /// Create the guards.
    /// Guards are true when f(x) = Ax + b > 0
    VectorAffineFunction guard12(Matrix<Real>(1,3,0.0,1.0,0.0),Vector<Real>(1,-1.0));
    cout << "guard12=" << guard12 << endl << endl;
    VectorAffineFunction guard23(Matrix<Real>(1,3,1.0,0.0,0.0),Vector<Real>(1, - hmax + Delta));
    cout << "guard23=" << guard23 << endl << endl;
    VectorAffineFunction guard34(Matrix<Real>(1,3,0.0,-1.0,0.0),Vector<Real>(1,0.0));
    cout << "guard34=" << guard34 << endl << endl;
    VectorAffineFunction guard41(Matrix<Real>(1,3,-1.0,0.0,0.0),Vector<Real>(1,hmin + Delta));
    cout << "guard41=" << guard41 << endl << endl;

    /// Create the invariants.
    /// Invariants are true when f(x) = Ax + b < 0
    /// forced transitions do not need an explicit invariant,
    /// we need only the invariants for location 2 and 4
    VectorAffineFunction inv2(Matrix<Real>(1,3,1.0,0.0,0.0),Vector<Real>(1, - hmax - Delta));//
    cout << "inv2=" << inv2 << endl << endl;
    VectorAffineFunction inv4(Matrix<Real>(1,3,-1.0,0.0,0.0),Vector<Real>(1, hmin - Delta));
    cout << "inv4=" << inv4 << endl << endl;

    /// Build the automaton
    watertank_system.new_mode(l1,opening_d);
    watertank_system.new_mode(l2,opened_d);
    watertank_system.new_mode(l3,closing_d);
    watertank_system.new_mode(l4,closed_d);

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
    evolver.parameters().maximum_enclosure_radius = 0.5;
    evolver.parameters().maximum_step_size = 1.25;
    evolver.parameters().enable_subdivisions = true;
    evolver.verbosity = 1;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;

    std::cout << "Computing evolution starting from location l2, x = 0.0, y = 1.0" << std::endl;

    Box initial_box(3, 6.0,6.0, 1.0,1.0, 0.0000,0.0000);
    HybridEnclosureType initial_enclosure(l2,initial_box);
    Box bounding_box(3, 0.0,10.0, -0.1,1.1, 0.0,tmax);

    HybridTime evolution_time(tmax,jmax);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(watertank_system,initial_enclosure,evolution_time,LOWER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit.reach.size()="<<orbit.reach().size()<<std::endl;
    std::cout << "Orbit.final.size()="<<orbit.final().size()<<std::endl;
    //plot("tutorial-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit.initial());
    plot("watertank-nonlinear-orbit-tx", 2,0, 3, bounding_box, Colour(0.0,0.5,1.0), orbit, -1);
    plot("watertank-nonlinear-orbit-ty", 2,1, 3, bounding_box, Colour(0.0,0.5,1.0), orbit, -1);
    plot("watertank-nonlinear-orbit-xy", 0,1, 3, bounding_box, Colour(0.0,0.5,1.0), orbit, -1);
    // textplot("watertank-nonlinear-orbit", orbit);
    // textplot("watertank-nonlinear-final", orbit.final());


    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.parameters().lock_to_grid_time = 32.0;
    analyser.verbosity=5;
    std::cout <<  analyser.parameters() << std::endl;

    HybridImageSet initial_set;
    initial_set[l2]=initial_box;

    HybridTime reach_time(64.0,6);

    // Compute evolved sets (i.e. at the evolution time) and reach sets (i.e. up to the evolution time) using lower semantics.
    // These functions run a bunch of simulations with bounded approximation errors and combines the results.
    // If the desired evolution time can not be attained without exceeding the error bounds, then the run discarded (without warning)
    std::cout << "Computing lower reach set... " << std::flush;
    HybridGridTreeSet lower_reach_set = analyser.lower_reach(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-nonlinear-lower_reach1", 2,0, 3, bounding_box, Colour(0.0,0.5,1.0), lower_reach_set, -1);
/*
    // Compute evolved sets and reach sets using upper semantics.
    // These functions compute over-approximations to the evolved and reachabe sets. Subdivision is used
    // as necessary to keep the local errors reasonable. The accumulated global error may be very large.
    std::cout << "Computing upper reach set... " << std::flush;
    HybridGridTreeSet upper_reach_set = analyser.upper_reach(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-nonlinear-upper_reach1",bounding_box, Colour(0.0,0.5,1.0), upper_reach_set);

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
