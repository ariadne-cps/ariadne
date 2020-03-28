/***************************************************************************
 *            hybrid_evolution_tutorial.cpp
 *
 *  Copyright  2020  Luca Geretti
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

#include <ariadne.hpp>

using namespace Ariadne;

HybridAutomaton get_tank()
{
    // Declare the system constants
    RealConstant lambda("lambda",0.02_decimal);
    RealConstant rate("rate",0.3_decimal);

    // Declare the variables for the dynamics
    RealVariable aperture("aperture");
    RealVariable height("height");

    // Declare the variable for the automaton name
    StringVariable tank("tank");

    // Create the tank automaton
    HybridAutomaton automaton(tank.name());

    // Declare a trivial discrete location (we use an empty label since there is only one location).
    DiscreteLocation draining;

    // The water level is always given by the same dynamic.
    // The inflow is controlled by the valve aperture, the outflow depends on the
    // pressure, which is proportional to the water height.
    automaton.new_mode(draining,{dot(height)=-lambda*height+rate*aperture});

    return automaton;
}

HybridAutomaton get_valve()
{
    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4);

    // Declare the shared system variables
    RealVariable aperture("aperture");
    RealVariable height("height");

    // Declare the events we use
    DiscreteEvent e_stop_opening("stop_opening");
    DiscreteEvent e_stop_closing("stop_closing");
    DiscreteEvent e_can_open("can_open");
    DiscreteEvent e_can_close("can_close");

    // Declare the variable for the automaton name
    StringVariable valve("valve");

    // Create the valve automaton
    HybridAutomaton automaton(valve.name());

    // Declare the values the valve variable can have
    DiscreteLocation opening(valve|"opening");
    DiscreteLocation closed(valve|"closed");
    DiscreteLocation opened(valve|"opened");
    DiscreteLocation closing(valve|"closing");

    // Define the algebraic equations for the opened/closed locations.
    automaton.new_mode(opened,{let(aperture)=1});
    automaton.new_mode(closed,{let(aperture)=0});
    // Define the differential equations for the opening/closing locations.
    automaton.new_mode(opening,{dot(aperture)=+1/T});
    automaton.new_mode(closing,{dot(aperture)=-1/T});

    // Define the transitions: source location, event and target location;
    // then a mix of reset, guard and event kind can be present; if the event kind
    // is not specified, then also the guard can't be specified: this implicitly
    // means that the event is an input event for this automaton.
    automaton.new_transition(closed,e_can_open,opening,{next(aperture)=aperture});
    automaton.new_transition(opening,e_stop_opening,opened,aperture>=1,EventKind::URGENT);
    automaton.new_transition(opened,e_can_close,closing,{next(aperture)=aperture});
    automaton.new_transition(closing,e_stop_closing,closed,aperture<=0,EventKind::URGENT);

    return automaton;
}

HybridAutomaton get_controller()
{
    // Declare some constants
    RealConstant hmin("hmin",5.75_decimal);
    RealConstant hmax("hmax",7.75_decimal);
    RealConstant delta("delta",0.02_decimal);

    // Declare the shared system variables
    RealVariable height("height");

    // Declare the events we use
    DiscreteEvent e_can_open("can_open");
    DiscreteEvent e_can_close("can_close");
    DiscreteEvent e_must_open("must_open");
    DiscreteEvent e_must_close("must_close");

    // Declare the variable for the automaton name
    StringVariable controller("controller");

    // Create the controller automaton
    HybridAutomaton automaton(controller.name());

    // Declare the locations for the controller
    DiscreteLocation rising(controller|"rising");
    DiscreteLocation falling(controller|"falling");

    // Instantiate modes for each location; since no dynamics is present, we create an empty list
    automaton.new_mode(rising,List<RealAssignment>());
    automaton.new_mode(falling,List<RealAssignment>());

    // Specify the invariants valid in each mode. Note that every invariant
    // must have an action label. This is used internally, for example, to
    // check non-blockingness of urgent actions.
    automaton.new_invariant(falling,height>=hmin-delta,e_must_open);
    automaton.new_invariant(rising,height<=hmax+delta,e_must_close);

    // Specify the transitions, starting from the source location, according to an event, to a target location;
    // Following those arguments you specify a guard and whether the event is permissive or urgent.
    automaton.new_transition(falling,e_can_open,rising,height<=hmin+delta,EventKind::PERMISSIVE);
    automaton.new_transition(rising,e_can_close,falling,height>=hmax-delta,EventKind::PERMISSIVE);

    return automaton;
}

Void simulate_evolution(const CompositeHybridAutomaton& system, const Nat& log_verbosity)
{
    // Re-introduce the shared system variables required for the initial set
    RealVariable aperture("aperture");
    RealVariable height("height");

    StringVariable valve("valve");
    StringVariable controller("controller");

    StringConstant opened("opened");
    StringConstant rising("rising");

    // Create a simulator object
    HybridSimulator simulator;
    simulator.set_step_size(0.01);
    simulator.verbosity = log_verbosity;

    // Set an initial point for the simulation
    HybridRealPoint initial_point({valve|opened,controller|rising},{height=7});

    // Set the maximum simulation time
    HybridTime termination_time(30.0,5);

    // Compute a simulation trajectory
    std::cout << "Computing simulation trajectory... \n" << std::flush;
    Orbit<HybridApproximatePoint> orbit = simulator.orbit(system,initial_point,termination_time);
    std::cout << "done." << std::endl;

    // Plot the simulation trajectory using all different projections
    std::cout << "Plotting simulation trajectory... " << std::flush;
    Axes2d time_height_axes(0<=TimeVariable()<=30,5<=height<=9);
    plot("simulation_t-height.png",time_height_axes, Colour(1.0,0.75,0.5), orbit);
    Axes2d time_aperture_axes(0<=TimeVariable()<=30,-0.1<=aperture<=1.1);
    plot("simulation_t-aperture.png",time_aperture_axes, Colour(1.0,0.75,0.5), orbit);
    Axes2d height_aperture_axes(5<=height<=9,-0.1<=aperture<=1.1);
    plot("simulation_height-aperture",height_aperture_axes, Colour(1.0,0.75,0.5), orbit);
    std::cout << "done." << std::endl;
}

GeneralHybridEvolver create_evolver(const CompositeHybridAutomaton& system, const Nat& log_verbosity)
{
    // Create a GeneralHybridEvolver object
    GeneralHybridEvolver evolver(system);

    // Set the evolver configuration
    evolver.configuration().set_maximum_enclosure_radius(3.0);
    evolver.configuration().set_maximum_step_size(0.25);
    evolver.verbosity=log_verbosity;

    std::cout << evolver.configuration() << std::endl;

    return evolver;
}

Void compute_evolution(const GeneralHybridEvolver& evolver) {

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolver::OrbitType OrbitType;

    // Re-introduce the shared system variables required for the initial set
    RealVariable aperture("aperture");
    RealVariable height("height");

    StringVariable valve("valve");
    StringVariable controller("controller");

    StringConstant opened("opened");
    StringConstant rising("rising");

    // Define the initial set, by supplying the location as a list of locations for each composed automata, and
    // the continuous set as a list of variable assignments for each variable controlled on that location
    // (the assignment can be either a singleton value using the == symbol or an interval using the <= symbols)
    HybridSet initial_set({valve|opened,controller|rising},{6.9_decimal<=height<=7});
    // Define the evolution time: continuous time and maximum number of transitions
    HybridTime termination_time(30.0,5);
    // Compute the orbit using upper semantics
    std::cout << "Computing evolution... \n" << std::flush;
    OrbitType orbit = evolver.orbit(initial_set,termination_time,Semantics::UPPER);
    std::cout << "done." << std::endl;

    // Plot the trajectory using two different projections
    std::cout << "Plotting trajectory... " << std::flush;
    Axes2d time_height_axes(0<=TimeVariable()<=30,5<=height<=9);
    plot("finite_evolution_t-height",time_height_axes, Colour(1.0,0.75,0.5), orbit);
    Axes2d time_aperture_axes(0<=TimeVariable()<=30,-0.1<=aperture<=1.1);
    plot("finite_evolution_t-aperture.png",time_aperture_axes, Colour(1.0,0.75,0.5), orbit);
    Axes2d height_aperture_axes(5<=height<=9,-0.1<=aperture<=1.1);
    plot("finite_evolution_height-aperture",height_aperture_axes, Colour(1.0,0.75,0.5), orbit);
    std::cout << "done." << std::endl;
}

HybridReachabilityAnalyser create_analyser(const GeneralHybridEvolver& evolver, const Nat& log_verbosity)
{
    // Silence the evolver
    evolver.verbosity=0;
    // Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);

    //  Set the analyser configuration
    analyser.configuration().set_maximum_grid_depth(6);
    analyser.configuration().set_lock_to_grid_time(5);
    analyser.verbosity=log_verbosity;

    std::cout << analyser.configuration() << std::endl;

    return analyser;
}

Void compute_reachability(const HybridReachabilityAnalyser& analyser) {

    // Re-introduce the shared system variables required for the initial set
    RealVariable aperture("aperture");
    RealVariable height("height");

    StringVariable valve("valve");
    StringVariable controller("controller");

    StringConstant opened("opened");
    StringConstant rising("rising");

    // Define the initial set
    HybridSet initial_set({valve|opened,controller|rising},{6.9_decimal<=height<=7});

    // Compute over-approximation to finite-time reachable set using upper semantics.
    std::cout << "Computing outer chain reach set... \n" << std::flush;
    HybridGridTreePaving outer_chain_reach = analyser.outer_chain_reach(initial_set);
    std::cout << "done." << std::endl;

    std::cout << "Plotting trajectory... " << std::flush;
    Axes2d height_aperture_axes(5<=height<=9,-0.1<=aperture<=1.1);
    plot("outer_chain_reach.png", height_aperture_axes, Colour(0.0,0.5,1.0), outer_chain_reach);
    std::cout << "done." << std::endl;
}

Int main(Int argc, const char* argv[])
{
    // Acquire the verbosity value from the command line
    Nat log_verbosity = get_verbosity(argc,argv);

    // Create the composed automaton
    CompositeHybridAutomaton watertank_system({get_tank(),get_valve(),get_controller()});

    // Choose a compact output representation for systems
    CompositeHybridAutomaton::set_default_writer(new CompactCompositeHybridAutomatonWriter());

    // Print the system description on the command line
    std::cout << "System: " << watertank_system << std::endl;

    // Compute an approximate simulation of the system evolution
    simulate_evolution(watertank_system,log_verbosity);

    // Create an evolver object
    GeneralHybridEvolver evolver = create_evolver(watertank_system,log_verbosity);

    // Compute the system evolution
    compute_evolution(evolver);

    // Create an analyser object
    HybridReachabilityAnalyser analyser = create_analyser(evolver,log_verbosity);

    // Compute the system reachability
    compute_reachability(analyser);
}
