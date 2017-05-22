/***************************************************************************
 *            watertank-monolithic.cpp
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
#include "ariadne.hpp"

#include "expression/expression_set.hpp"

using namespace Ariadne;

Int main(Int argc, const char* argv[])
{
    Nat evolver_verbosity = 0;
    if(argc>1) { evolver_verbosity=atoi(argv[1]); }

    typedef GeneralHybridEvolver GeneralHybridEvolverType;


    DRAWING_METHOD = AFFINE_DRAW;
    DRAWING_ACCURACY = 1;

    /// Set the system parameters
    RealConstant a("a", Decimal(-0.02));
    RealConstant b("b", 0.3_dec);
    RealConstant T("T", 4.0_dec);
    RealConstant hmin("hmin", 5.5_dec);
    RealConstant Delta("Delta", 0.05_dec);
    RealConstant hmax("hmax", 8.0_dec);

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

    // Create coordinate functions in two variables.
    RealVariable height("height");
    RealVariable aperture("aperture");

    /// Create the dynamics
    DottedRealAssignment tank_dynamic(dot(height)=a*height+b*aperture);
    DottedRealAssignment valve_opening_dynamic(dot(aperture)=1/T);
    DottedRealAssignment valve_closing_dynamic(dot(aperture)=-1/T);
    DottedRealAssignment valve_constant_dynamic(dot(aperture)=0.0_dec);
//    RealAssignment valve_open_dynamic(let(aperture)=1.0_dec);
    RealAssignment valve_open_dynamic(let(aperture)=1.01_dec);
    RealAssignment valve_closed_dynamic(let(aperture)=0.0_dec);

    /// Create the resets
    PrimedRealAssignment tank_reset(next(height)=height);
    PrimedRealAssignment valve_reset(next(aperture)=aperture);
    PrimedRealAssignment valve_open_reset(next(aperture)=1.0_dec);
    PrimedRealAssignment valve_closed_reset(next(aperture)=0.0_dec);
    cout << "tank_reset=" << tank_reset << endl << endl;
    cout << "valve_reset=" << valve_reset << endl << endl;

    /// Create the guards.
    /// Guards are true when g(x) >= 0
    ContinuousPredicate finish_opening_guard(aperture>=1.01_dec);
    cout << "finish_opening_guard=" << finish_opening_guard << endl << endl;
    ContinuousPredicate start_closing_guard(height>=hmax-Delta);
    cout << "start_closing_guard=" << start_closing_guard << endl << endl;
    ContinuousPredicate finish_closing_guard(aperture<=0);
    cout << "finish_closing_guard=" << finish_closing_guard << endl << endl;
    ContinuousPredicate start_opening_guard(height<=hmin+Delta);
    cout << "start_opening_guard=" << start_opening_guard << endl << endl;

    /// Create the invariants.
    /// Invariants are true when c(x) <= 0
    /// Urgent transitions do not need an explicit invariant,
    /// we need only the invariants for location 2 and 4
    ContinuousPredicate start_closing_invariant(height<=hmax + Delta);
    cout << "start_closing_invariant=" << start_closing_invariant << endl << endl;
    ContinuousPredicate start_opening_invariant(height>=hmin - Delta);
    cout << "start_opening_invariant=" << start_opening_invariant << endl << endl;

    /// Build the automaton
    watertank_system.new_mode(opening,{tank_dynamic,valve_opening_dynamic});
//    watertank_system.new_mode(open,{tank_dynamic,valve_constant_dynamic});
    watertank_system.new_mode(open,{valve_open_dynamic},{tank_dynamic});
    watertank_system.new_mode(closing,{tank_dynamic,valve_closing_dynamic});
//    watertank_system.new_mode(closed,{tank_dynamic,valve_constant_dynamic});
    watertank_system.new_mode(closed,{valve_closed_dynamic},{tank_dynamic});

    watertank_system.new_invariant(open,start_closing_invariant,must_start_closing);
    watertank_system.new_invariant(closed,start_opening_invariant,must_start_opening);

//    watertank_system.new_transition(opening,finish_opening,open,{tank_reset,valve_open_reset},finish_opening_guard,urgent);
    watertank_system.new_transition(opening,finish_opening,open,{tank_reset},finish_opening_guard,urgent);
    watertank_system.new_transition(open,start_closing,closing,{tank_reset,valve_reset},start_closing_guard,permissive);
//    watertank_system.new_transition(closing,finish_closing,closed,{tank_reset,valve_closed_reset},finish_closing_guard,urgent);
    watertank_system.new_transition(closing,finish_closing,closed,{tank_reset},finish_closing_guard,urgent);
    watertank_system.new_transition(closed,start_opening,opening,{tank_reset,valve_reset},start_opening_guard,permissive);

    /// Finished building the automaton

    cout << "Automaton = " << watertank_system << endl << endl;

    /// Compute the system evolution

    /// Create a GeneralHybridEvolver object
    GeneralHybridEvolverType evolver(watertank_system);
    evolver.verbosity = evolver_verbosity;

    /// Set the evolution parameters
    evolver.configuration().set_maximum_enclosure_radius(0.5);
    evolver.configuration().set_maximum_step_size(2.5);
    evolver.configuration().set_maximum_spacial_error(1e-3);
    std::cout <<  evolver.configuration() << std::endl;

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolverType::EnclosureType EnclosureType;
    typedef GeneralHybridEvolverType::EnclosureListType EnclosureListType;
    typedef GeneralHybridEvolverType::OrbitType OrbitType;

    //HybridSet initial_set(opening,(height==0.0, aperture==0.0));
    Dyadic eps(1.0/256);
    HybridSet initial_set(opening,{0<=height<=eps, 0<=aperture<=eps});
    std::cout << "Initial set = " << initial_set << "\n" ;
    HybridTime evolution_time(80.0,10);
    std::cout << "Evolution time = "  << evolution_time << "\n" ;

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(initial_set,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit.final_size="<<orbit.final().size()<<std::endl;
    std::cout << "Orbit.reach_size="<<orbit.reach().size()<<std::endl;

    std::cout << "Plotting orbit... "<<std::flush;
    Axes2d height_aperture_axes(-0.1,height,9.1, -0.1,aperture,1.3);
    plot("watertank-orbit", height_aperture_axes, Colour(0.0,0.5,1.0), orbit, Colour(0.0,1.0,1.0), orbit.final());//, Colour(1.0,0.0,0.0),orbit.final()[9]
    Axes2d time_height_axes(0,TimeVariable(),80, -0.1,height,9.1);
    plot("watertank-height", time_height_axes, Colour(0.0,0.5,1.0), orbit, Colour(0.0,1.0,1.0), orbit.final());
    Axes2d time_aperture_axes(0,TimeVariable(),80, -0.1,aperture,1.31);
    plot("watertank-aperture", time_aperture_axes, Colour(0.0,0.5,1.0), orbit, Colour(0.0,1.0,1.0), orbit.final());
    std::cout << "done." << std::endl;


//    HybridGrid grid(watertank_system.state_space());
    HybridGrid grid(watertank_system.state_auxiliary_space());
    std::cerr << "grid=" << grid << std::endl;

    std::cout << "Discretising orbit" << std::flush;
    HybridGridTreeSet hgts(grid);
    for (ListSet<EnclosureType>::ConstIterator it = orbit.reach().begin(); it != orbit.reach().end(); it++)
    {
        std::cout<<"."<<std::flush;
        it->state_auxiliary_set().adjoin_outer_approximation_to(hgts,4);
    }
    std::cout << "done." << std::endl;

    plot("watertank-reach", height_aperture_axes, Colour(0.0,0.5,1.0), hgts);


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
    // These functions run a bunch of simulations with singleton approximation errors and combines the results.
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

    ExactBoxType initial_box2(2, 0.0,0.001, 0.0,0.001);
    HybridImageSet initial_set2;
    initial_set2[l1]=initial_box2;

    plot("watertank-initial_set2",bounding_box, Colour(0.0,0.5,1.0), initial_set2);

    // Compute evolved sets (i.e. at the evolution time) and reach sets (i.e. up to the evolution time) using lower semantics.
    // These functions run a bunch of simulations with singleton approximation errors and combines the results.
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
