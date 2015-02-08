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

Figure& operator<<(Figure& fig, const ListSet<HybridEnclosure>& set) {
    for(ListSet<HybridEnclosure>::ConstIterator set_iter=set.begin(); set_iter!=set.end(); ++set_iter) {
        fig << set_iter->continuous_set(); } return fig; }

Figure& operator<<(Figure& fig, const HybridGridTreeSet& set) {
    for(Map<DiscreteLocation,GridTreeSet>::ConstIterator loc_iter=set.locations_begin(); loc_iter!=set.locations_end(); ++loc_iter) {
        fig << loc_iter->second; } return fig; }

/// Function for plotting the orbit and reachability set
template<class SET> Void plot(const char* filename, const Int& xaxis, const Int& yaxis, const Int& numVariables, const ExactBox& bbox, const Colour& fc, const SET& set, const Int& MAX_GRID_DEPTH) {
    // Assigns local variables
    Figure fig;

    fig.set_projection_map(PlanarProjectionMap(numVariables,xaxis,yaxis));
    fig.set_bounding_box(bbox);

    // If the grid must be shown
    if (MAX_GRID_DEPTH >= 0)
    {
        // The rectangle to be drawn
        ExactBox rect = ExactBox(numVariables);
        // Chooses the fill colour
        fig << fill_colour(Colour(1.0,1.0,1.0));

        // Gets the number of times each variable interval would be divided by 2
        Int numDivisions = MAX_GRID_DEPTH / numVariables;
        // Gets the step in the x direction, by 1/2^(numDivisions+h), where h is 1 if the step is to be further divided by 2, 0 otherwise
        Float64 step_x = 1.0/(1 << (numDivisions + ((MAX_GRID_DEPTH - numDivisions*numVariables > xaxis) ? 1 : 0)));
        // Initiates the x position to the bounding box left bound
        Float64 pos_x = bbox[0].lower().raw();
        // Sets the rectangle 2-nd interval to the corresponding bounding box interval (while the >2 intervals are kept at [0,0])
        rect[yaxis] = bbox[1];
        // While between the interval
        while (pos_x < bbox[0].upper().raw())
        {
            rect[xaxis] = ExactInterval(pos_x,pos_x+step_x); // Sets the rectangle x coordinate
            pos_x += step_x; // Shifts the x position
            fig << rect; // Appends the rectangle
        }

        // Repeats for the rectangles in the y direction
        Float64 step_y = 1.0/(1 << (numDivisions + ((MAX_GRID_DEPTH - numDivisions*numVariables > yaxis) ? 1 : 0)));
        Float64 pos_y = bbox[1].lower().raw();
        rect[xaxis] = bbox[0];
        while (pos_y < bbox[1].upper().raw())
        {
            rect[yaxis] = ExactInterval(pos_y,pos_y+step_y);
            fig << rect;
            pos_y += step_y;
        }
    }
    // Draws and creates file
    fig.set_fill_colour(fc);
    fig << set;
    fig.write(filename);
}


Int main()
{

    /// Set the system parameters
    Real a = 0.065_dec;
    Real b = 0.3_dec;
    Real T = 4.0_dec;
    Real hmin = 5.5_dec;
    Real Delta = 0.05_dec;
    Real hmax = 8.0_dec;

    Real zero = 0;
    Real one = 1;

    double tmax = 80.0;
    Int jmax = 6;


    /// Build the Hybrid System

    /// Create a HybridAutomton object
    HybridAutomaton watertank_system;

    /// Create four discrete states
    StringVariable valve("valve");
    DiscreteLocation opening(valve|"opening");
    DiscreteLocation open(valve|"open");
    DiscreteLocation closing(valve|"closing");
    DiscreteLocation closed(valve|"closed");

    /// Create the discrete events
    DiscreteEvent finish_opening("finish_opening");
    DiscreteEvent start_closing("start_closing");
    DiscreteEvent must_start_closing("must_start_closing");
    DiscreteEvent finish_closing("finish_closing");
    DiscreteEvent start_opening("start_opening");
    DiscreteEvent must_start_opening("must_start_opening");

    /// Coordinates
    RealVariable height("height");
    RealVariable aperture("aperture");
    TimeVariable time;
    /// Create the dynamics

    /// The dynamic for an opening valve
    DottedRealAssignments opening_dynamic = { dot(height)=-a*sqrt(height)+b*aperture, dot(aperture)=1/T };
    /// The dynamic for a closing valve
    DottedRealAssignments closing_dynamic = { dot(height)=-a*sqrt(height)+b*aperture, dot(aperture)=-1/T };
    /// The dynamic for an opened valve
    DottedRealAssignments open_dynamic = { dot(height)=-a*sqrt(height)+b, dot(aperture)=zero };
    /// The dynamic for an opened valve
    DottedRealAssignments closed_dynamic = { dot(height)=-a*sqrt(height), dot(aperture)=zero };

    cout << "opening dynamic = " << opening_dynamic << endl << endl;
    cout << "closing dynamic = " << closing_dynamic << endl << endl;
    cout << "open dynamic = " << open_dynamic << endl << endl;
    cout << "closed dynamic = " << closed_dynamic << endl << endl;

    /// Create the resets
    PrimedRealAssignments reset_aperture_zero = { next(height)=height, next(aperture)=0 };
    cout << "reset_aperture_zero=" << reset_aperture_zero << endl << endl;
    PrimedRealAssignments reset_aperture_one = { next(height)=height, next(aperture)=1 };
    cout << "reset_aperture_one=" << reset_aperture_one << endl << endl;

    /// Create the guards.
    /// Guards are true when f(height) >= 0
    ContinuousPredicate finish_opening_guard(aperture-1>=0);
    cout << "finish_opening_guard=" << finish_opening_guard << endl << endl;
    ContinuousPredicate start_closing_guard(height-(hmax-Delta)>=0);
    cout << "start_closing_guard=" << start_closing_guard << endl << endl;
    ContinuousPredicate finish_closing_guard(-aperture>=0);
    cout << "finish_closing_guard=" << finish_closing_guard << endl << endl;
    ContinuousPredicate start_opening_guard(-height+(hmin+Delta)>=0);
    cout << "start_opening_guard=" << start_opening_guard << endl << endl;

    /// Create the invariants.
    /// Invariants are true when f(height) <= 0
    /// forced transitions do not need an explicit invariant,
    /// we need only the invariants for location 2 and 4
    ContinuousPredicate start_closing_invariant(height-(hmax+Delta)<=0);
    cout << "start_closing_invariant=" << start_closing_invariant << endl << endl;
    ContinuousPredicate start_opening_invariant(-height+(hmin-Delta)<=0);
    cout << "start_opening_invariant=" << start_opening_invariant << endl << endl;

    /// Build the automaton
    watertank_system.new_mode(opening,opening_dynamic);
    watertank_system.new_mode(open,open_dynamic);
    watertank_system.new_mode(closing,closing_dynamic);
    watertank_system.new_mode(closed,closed_dynamic);

    watertank_system.new_invariant(open,start_closing_invariant,must_start_closing);
    watertank_system.new_invariant(closed,start_opening_invariant,must_start_opening);

    watertank_system.new_transition(opening,finish_opening,open,reset_aperture_one,finish_opening_guard,urgent);
    watertank_system.new_transition(open,start_closing,closing,reset_aperture_one,start_closing_guard,permissive);
    watertank_system.new_transition(closing,finish_closing,closed,reset_aperture_zero,finish_closing_guard,urgent);
    watertank_system.new_transition(closed,start_opening,opening,reset_aperture_zero,start_opening_guard,permissive);


    /// Finished building the automaton

    cout << "Automaton = " << watertank_system << endl << endl;

    /// Compute the system evolution

    /// Create a GeneralHybridEvolver object
    GeneralHybridEvolver evolver(watertank_system);

    /// Set the evolution parameters
    evolver.configuration().set_maximum_enclosure_radius(0.5);
    evolver.configuration().set_maximum_step_size(1.25);
    evolver.configuration().set_enable_subdivisions(true);
    evolver.verbosity = 1;
    std::cout <<  evolver.configuration() << std::endl;

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolver::EnclosureType HybridEnclosureType;
    typedef GeneralHybridEvolver::OrbitType OrbitType;
    typedef GeneralHybridEvolver::EnclosureListType EnclosureListType;

    std::cout << "Computing evolution starting from location l2, x = 6.0, y = 1.0" << std::endl;

    //RealVariablesBox initial_box((height==0.5, aperture==0.0));
    Decimal h0l(0.5),h0u(0.5001), a0l(0.0),a0u(0.0001);
    RealVariablesBox initial_box((h0l<=height<=h0u, a0l<=aperture<=a0u));
    HybridSet initial_set(opening,initial_box);

    HybridTime evolution_time(tmax,jmax);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(initial_set,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Plotting orbit... "<<std::flush;
    Axes2d height_aperture_axes(-0.1,height,9.1, -0.1,aperture,1.3);
    Axes2d time_height_axes(0.0,time,tmax, -0.1,height,9.1);

    std::cout << "Orbit.reach.size()="<<orbit.reach().size()<<std::endl;
    std::cout << "Orbit.final.size()="<<orbit.final().size()<<std::endl;
    //plot("tutorial-orbit",bounding_box, Colour(0.0,0.5,1.0), orbit.initial());
    plot("watertank-nonlinear-orbit", height_aperture_axes, Colour(0.0,0.5,1.0), orbit, Colour(0.0,1.0,1.0), orbit.final());
    plot("watertank-nonlinear-height", time_height_axes, Colour(0.0,0.5,1.0), orbit, Colour(0.0,1.0,1.0), orbit.final());
    //plot("watertank-nonlinear-orbit-time-height", 2,0, 3, bounding_box, Colour(0.0,0.5,1.0), orbit.reach(), -1);
    //plot("watertank-nonlinear-orbit-time-aperture", 2,1, 3, bounding_box, Colour(0.0,0.5,1.0), orbit.reach(), -1);
    //plot("watertank-nonlinear-orbit-height-aperture", 0,1, 3, bounding_box, Colour(0.0,0.5,1.0), orbit.reach(), -1);
    // textplot("watertank-nonlinear-orbit", orbit);
    // textplot("watertank-nonlinear-final", orbit.final());

    std::cout << "Discretising orbit" << std::flush;
    Int depth = 3;
    HybridScaling scaling( {height|1.0, aperture|0.125} );
    HybridGrid grid(watertank_system.state_space(),scaling);
    HybridGridTreeSet hgts(grid);
    for (ListSet<HybridEnclosure>::ConstIterator it = orbit.reach().begin(); it != orbit.reach().end(); it++)
    {
        std::cout<<"."<<std::flush;
        it->adjoin_outer_approximation_to(hgts,depth);
    }
    std::cout << "done." << std::endl;

    plot("watertank-nonlinear-reach", height_aperture_axes, Colour(0.0,0.5,1.0), hgts);

    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(watertank_system,GeneralHybridEvolver(watertank_system));
    analyser.configuration().set_lock_to_grid_time(32.0);
    analyser.verbosity=5;
    std::cout <<  analyser.configuration() << std::endl;

    HybridTime reach_time(64.0,6);

    std::cout << "Omitting computation of global upper and lower reach set." << std::endl;

/*
    // Compute evolved sets (i.e. at the evolution time) and reach sets (i.e. up to the evolution time) using lower semantics.
    // These functions run a bunch of simulations with bounded approximation errors and combines the results.
    // If the desired evolution time can not be attained without exceeding the error bounds, then the run discarded (without warning)

    std::cout << "Computing lower reach set... " << std::flush;
    HybridGridTreeSet lower_reach_set = analyser.lower_reach(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-nonlinear-lower_reach1", 2,0, 3, bounding_box, Colour(0.0,0.5,1.0), lower_reach_set, -1);
    // Compute evolved sets and reach sets using upper semantics.
    // These functions compute over-approximations to the evolved and reachabe sets. Subdivision is used
    // as necessary to keep the local errors reasonable. The accumulated global error may be very large.
    std::cout << "Computing upper reach set... " << std::flush;
    HybridGridTreeSet upper_reach_set = analyser.upper_reach(watertank_system,initial_set,reach_time);
    std::cout << "done." << std::endl;
    plot("watertank-nonlinear-upper_reach1",bounding_box, Colour(0.0,0.5,1.0), upper_reach_set);

    std::cout << "Computing evolution starting from location l1, x = 0.0, y = 0.0" << std::endl;

    ExactBox initial_box2(2, 0.0,0.001, 0.0,0.001);
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
