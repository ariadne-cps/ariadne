/***************************************************************************
 *            tutorial.cc
 *
 *  Copyright  2008  Pieter Collins
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


//! \file tutorial.cc

#include "ariadne.h"

#include "hybrid_automaton-composite.h"
#include "hybrid_set.h"
#include "hybrid_evolver.h"
#include "hybrid_simulator.h"
#include "hybrid_graphics.h"

template<class T> void write(const char* filename, const T& t) {
    std::ofstream ofs(filename); ofs << t; ofs.close();
}

using namespace Ariadne;

typedef GeneralHybridEvolver HybridEvolverType;

// The automaton has two modes, both
// with two variables, room_temperature and time_of_day
// The dynamics is given by the differential equations
//   dot(T) = K(T_ext(t) - T(t)) + P delta(q(t) = On) \f$
// where
//   T_ext(t) = T_av - T_amp \cos(2 pi t)
// is the external temperature and
//   q(t) in {On,Off}
// is the discrete state of the heater.
//
// The parameters are the insulation coefficient K, the average
// external temperature T_av, the amplitude of the
// temperature fluctuations T_amp and the power of the heater P.
//
// The heater is controlled by a thermostat. The heater turns off whenever
// the temperature rises above a fixed threshold T_Off and
// turns on nondeterministically for a temperature \T between
// T^+_On and \f$T^-_On.
//
// The clock time tau is reset to zero whenever tau becomes
// equal to one.

// System variables
//   time-of-day t (d)
//   room temperature T (C)

// System paramters
//   Heating power P
//   Thermal coefficient K
//   Average external temperature Te
//   Amplitude of external temperature fluctuations Ta
//   Temperature at which the heater is turned off Toff
//   Temperature below which the heater may be turned on Tonact
//   Temperature below which the heater must be turned on Toninv



CompositeHybridAutomaton create_heating_system()
{
    // Set the system dynamic parameters
    RealConstant P("P",4.0);
    RealConstant K("K",1.0);
    RealConstant Tav("Tav",16.0);
    RealConstant Tamp("Tamp",8.0);

    // Set the system control parameters
    RealConstant Tmax("Tmax",23.0);
    RealConstant Tmin("Tmin",14.0);
    RealConstant Toff("Toff",21.0);
    RealConstant Ton_upper("Ton_upper",15.125);
    RealConstant Ton_lower("Ton_lower",14.875);

    // Create the discrete states
    StringVariable heating("heating");
    StringConstant on("on");
    StringConstant off("off");

    // Create the discrete events
    DiscreteEvent switch_on("switch_on");
    DiscreteEvent switch_off("switch_off");
    DiscreteEvent midnight("midnight");

    // Declare the system variables.
    RealVariable T("T");
    RealVariable t("t");

    // Create the heater subsystem
    HybridAutomaton heater;
    heater.new_mode( heating|on, (dot(T)=P+K*(Tav-Tamp*Ariadne::cos(2.0*pi*t)-T)) );
    heater.new_mode( heating|off, (dot(T)=K*(Tav-Tamp*Ariadne::cos(2.0*pi*t)-T)) );
    heater.new_invariant( heating|off, T<=Ton_lower, switch_on );
    heater.new_transition( (heating|off), switch_on, (heating|on), (next(T)=T), T<=Ton_upper, permissive );
    heater.new_transition( (heating|on), switch_off, heating|off, (next(T)=T), T>=Toff, urgent );

    // Create the clock subsystem
    HybridAutomaton clock;
    clock.new_mode( (dot(t)=1.0) );
    clock.new_transition( midnight, next(t)=0.0, t>=1.0, urgent );

    CompositeHybridAutomaton heating_system((clock,heater));
    std::cout << "heating_system=" << heating_system << "\n" << "\n";

    return heating_system;
}

HybridEvolverType create_evolver()
{
    // Create a GeneralHybridEvolver object
    HybridEvolverType evolver;

    // Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = 0.25;
    //evolver.parameters().maximum_step_size = 0.125;
    evolver.parameters().maximum_step_size = 0.5;
    evolver.verbosity=1;
    cout <<  evolver.parameters() << endl << endl;

    return evolver;
}


void compute_evolution(const CompositeHybridAutomaton& heating_system, const GeneralHybridEvolver& evolver)
{

    // Redefine the two discrete states
    AtomicDiscreteVariable clock("clock");
    StringVariable heater("heater");
    DiscreteLocation heating_off(StringVariable("heating")|"off");
    DiscreteLocation heating_on(StringVariable("heating")|"on");
    RealVariable T("T");
    RealVariable t("t");


    // Create a simulator object.
    HybridSimulator simulator;
    simulator.set_step_size(0.03125);

    // Set an initial point for the simulation
    HybridPoint initial_point(heating_off, (t|0.0,T|18.0) );
    cout << "initial_point=" << initial_point << endl;
    // Set the maximum simulation time
    HybridTime simulation_time(6.5,9);
    cout << "simulation_time=" << simulation_time << endl;

    // Compute a simulation trajectory
    cout << "Computing simulation trajectory... \n" << flush;
    Orbit<HybridPoint> trajectory = simulator.orbit(heating_system,initial_point,simulation_time);
    cout << "    done." << endl;
    // Write the simulation trajectory to standard output and plot.
    cout << "Writing simulation trajectory... " << flush;
    write("tutorial-trajectory.txt",trajectory);
    cout << "done." << endl;
    cout << "Plotting simulation trajectory... " << flush;
    hplot("tutorial-trajectory.png",(0.0<=t<=1.0,14.0<=T<=23.0), Colour(0.0,0.5,1.0), trajectory);
    cout << "done." << endl << endl;


    // Set the initial set.
    HybridSet initial_set(heating_off, (16.0<=T<=16.0625,0.0<=t<=0.015625) );
    cout << "initial_set=" << initial_set << endl;
    // Compute the initial set as a validated enclosure.
    HybridEnclosure initial_enclosure = evolver.enclosure(heating_system,initial_set);
    cout << "initial_enclosure="<<initial_enclosure << endl << endl;

    // Set the maximum evolution time
    HybridTime evolution_time(1.5,4);
    cout << "evolution_time=" << evolution_time << endl;

    // Compute a validated orbit.
    cout << "Computing orbit... \n" << flush;
    Orbit<HybridEnclosure> orbit = evolver.orbit(heating_system,initial_set,evolution_time,UPPER_SEMANTICS);
    cout << "    done." << endl;
    // Write the validated orbit to standard output and plot.
    cout << "Writing orbit... " << flush;
    write("tutorial-orbit.txt",orbit);
    cout << "done." << endl;
    cout << "Plotting orbit... " << flush;
    plot("tutorial-orbit.png",Box(2, 0.0,1.0, 14.0,23.0), Colour(0.0,0.5,1.0), orbit);
    cout << "done." << endl << endl;


    // Compute reachable and evolved sets
    cout << "Computing reach and evolve sets... \n" << flush;
    ListSet<HybridEnclosure> reach,evolve;
    make_lpair(reach,evolve) = evolver.reach_evolve(heating_system,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    cout << "    done." << endl;
    // Write the orbit to standard output and plot.
    cout << "Plotting reach and evolve sets... " << flush;
    plot("tutorial-reach_evolve.png",Box(2, 0.0,1.0, 14.0,23.0),
         Colour(0.0,0.5,1.0), reach, Colour(0.0,0.25,0.5), initial_enclosure, Colour(0.25,0.0,0.5), evolve);
    plot("tutorial-reach_evolve-off.png",Box(2, 0.0,1.0, 14.0,23.0),
         Colour(0.0,0.5,1.0), reach[heating_off], Colour(0.0,0.25,0.5), initial_enclosure, Colour(0.25,0.0,0.5), evolve[heating_on]);
    plot("tutorial-reach_evolve-on.png",Box(2, 0.0,1.0, 14.0,23.0),
         Colour(0.0,0.5,1.0), reach[heating_on], Colour(0.0,0.25,0.5), initial_enclosure, Colour(0.25,0.0,0.5), evolve[heating_on]);
    cout << "done." << endl;


}


void compute_reachable_sets(const CompositeHybridAutomaton& heating_system, const GeneralHybridEvolver& evolver)
{
/*
    // Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.parameters().initial_grid_density=5;
    analyser.parameters().initial_grid_depth=6;
    analyser.parameters().maximum_grid_depth=6;


    // Define the initial set
    HybridImageSet initial_set;
    AtomicDiscreteLocation heater_off(2);
    Box initial_box(2, 0.0,0.015625/4, 16.0,16.0+0.0625/16);
    initial_set[heater_off]=initial_box;

    // Set the maximum evolution time
    HybridTime reach_time(1.5,4);


    // Compute lower-approximation to finite-time evolved set using lower-semantics.
    std::cout << "Computing lower evolve set... " << std::flush;
    HybridGridTreeSet lower_evolve_set = analyser.lower_evolve(heating_system,initial_set,reach_time);
    std::cout << "done." << std::endl;

    // Compute lower-approximation to finite-time reachable set using lower-semantics.
    std::cout << "Computing lower reach set... " << std::flush;
    HybridGridTreeSet lower_reach_set = analyser.lower_reach(heating_system,initial_set,reach_time);
    std::cout << "done." << std::endl;

    plot("tutorial-lower_reach_evolve.png",Box(2, 0.0,1.0, 14.0,21.0),
         Colour(0.0,0.5,1.0), lower_reach_set,
         Colour(0.0,0.25,0.5), initial_set,
         Colour(0.25,0.0,0.5), lower_evolve_set);

    // Compute over-approximation to finite-time evolved set using upper semantics.
    // Subdivision is used as necessary to keep the local errors reasonable.
    // The accumulated global error may be very large.
    std::cout << "Computing upper evolve set... " << std::flush;
    HybridGridTreeSet upper_evolve_set = analyser.upper_evolve(heating_system,initial_set,reach_time);
    std::cout << "done." << std::endl;

    // Compute over-approximation to finite-time reachable set using upper semantics.
    std::cout << "Computing upper reach set... " << std::flush;
    HybridGridTreeSet upper_reach_set = analyser.upper_reach(heating_system,initial_set,reach_time);
    std::cout << "done." << std::endl;

    plot("tutorial-upper_reach_evolve.png",Box(2, 0.0,1.0, 14.0,21.0),
         Colour(0.0,0.5,1.0), upper_reach_set,
         Colour(0.0,0.25,0.5), initial_set,
         Colour(0.25,0.0,0.5), upper_evolve_set);

    // Compute over-approximation to infinite-time chain-reachable set using upper semantics.
    std::cout << "Computing chain reach set... " << std::flush;
    HybridGridTreeSet chain_reach_set = analyser.chain_reach(heating_system,initial_set);
    std::cout << "done." << std::endl;
    plot("tutorial-chain_reach.png",Box(2, 0.0,1.0, 14.0,21.0), Colour(0.0,0.5,1.0), chain_reach_set);
*/
}



void compute_reachable_sets_with_serialisation(const CompositeHybridAutomaton& heating_system, const HybridReachabilityAnalyser& analyser)
{
/*
    // Define the initial set
    HybridImageSet initial_set;
    AtomicDiscreteLocation heater_off(2);
    Box initial_box(2, 0.0,0.015625, 16.0,16.0625);
    initial_set[heater_off]=initial_box;


    // Compute the reach set for times between tlower and tupper.
    // The intermediate set is stored to an archive file and used to build the initial set for the reach step
    // Note that because of peculiarities in the Boost serialization library,
    // the object to be serialized must be declared const.
    Float tlower=0.25; Float tupper=0.75;
    HybridTime transient_time(tlower,4);
    HybridTime recurrent_time(tupper-tlower,16);

    const HybridGridTreeSet upper_intermediate_set = analyser.upper_evolve(heating_system,initial_set,transient_time);
    plot("tutorial-upper_intermediate.png",Box(2, 0.0,1.0, 14.0,18.0), Colour(0.0,0.5,1.0), upper_intermediate_set);

    std::ofstream output_file_stream("tutorial-transient.txt");
    text_oarchive output_archive(output_file_stream);
    output_archive << upper_intermediate_set;
    output_file_stream.close();

    HybridGridTreeSet rebuilt_upper_intermediate_set;

    std::ifstream input_file_stream("tutorial-transient.txt");
    text_iarchive input_archive(input_file_stream);
    input_archive >> rebuilt_upper_intermediate_set;
    input_file_stream.close();

    HybridGridTreeSet upper_recurrent_set = analyser.upper_reach(heating_system,initial_set,recurrent_time);
    plot("tutorial-upper_recurrent.png",Box(2, 0.0,1.0, 14.0,18.0), Colour(0.0,0.5,1.0), upper_recurrent_set);
*/
}




int main()
{
    // Create the system
    CompositeHybridAutomaton heating_system=create_heating_system();
    std::cerr<<heating_system<<"\n";

    // Create the analyser classes
    HybridEvolverType evolver=create_evolver();
    std::cerr<<evolver<<"\n";

    // Compute the system evolution
    compute_evolution(heating_system,evolver);
    compute_reachable_sets(heating_system,evolver);
}
