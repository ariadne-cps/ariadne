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

template<class T> void write(const char* filename, const T& t) {
    std::ofstream ofs(filename); ofs << t; ofs.close();
}

using namespace Ariadne;


const Float pi_flt = Ariadne::pi<Float>(); 
const Interval pi_ivl = Ariadne::pi<Interval>(); 

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
//   x[0]: time-of-day tau (d)
//   x[1]: room temperature T (C)

// System paramters
//   p[0]: Heating power P 
//   p[1]: Thermal coefficient K
//   p[2]: Average external temperature Te
//   p[3]: Amplitude of external temperature fluctuations Ta 
//   p[4]: Temperature at which the heater is turned off Toff
//   p[5]: Temperature below which the heater may be turned on Tonact
//   p[6]: Temperature below which the heater must be turned on Toninv

// System dynamic when the heater is on.
struct HeaterOn : FunctionData<2,2,4> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
        r[0] = 1.0; 
        r[1] = p[0] + p[1] * ( p[2] - p[3] * Ariadne::cos(2*pi_flt*x[0]) - x[1] ); // Need explicit Ariadne::cos due to bug in g++.
    }
};

// System dynamic when the heater is off.
struct HeaterOff : FunctionData<2,2,4> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
        typename R::value_type t0=2*pi_flt*x[0];
        r[0] = 1.0; 
        r[1] = p[1] * ( p[2] - p[3] * Ariadne::cos(t0) - x[1] ); // Need explicit Ariadne::cos due to bug in g++.
    }  
};


HybridAutomaton create_heating_system()
{
    // Create a HybridAutomton object
    HybridAutomaton heating_system;
  
    // Set the system dynamic parameters
    Float P=1.0;
    Float k=1.0;
    Float Tav=16.0;
    Float Tamp=8.0;
    
    // Set the system control parameters
    Float Toff=22.0;
    Float Ton_upper=16.0;
    Float Ton_lower=15.0;
  
    // Create the two discrete state 
    DiscreteState heater_on(1);
    DiscreteState heater_off(2);
  
    // Create the discrete events
    DiscreteEvent switch_on(1);
    DiscreteEvent switch_off(2);
    DiscreteEvent midnight(3);
  
    // Create the dynamics
    Function<HeaterOn> heater_on_dynamic(Vector<Float>(4, P,k,Tav,Tamp));
    Function<HeaterOff> heater_off_dynamic(Vector<Float>(4, P,k,Tav,Tamp));

    // Create the resets
    IdentityFunction heater_turn_off_reset(2);
    IdentityFunction heater_turn_on_reset(2);
    AffineFunction midnight_reset(Matrix<Float>(2,2,0.0,0.0,0.0,1.0),Vector<Float>(2,0.0,0.0));

    // Create the guards.
    AffineFunction heater_turn_off_guard(Matrix<Float>(1,2,0.0,1.0),Vector<Float>(1,-Toff));
    AffineFunction heater_turn_on_activation(Matrix<Float>(1,2,0.0,-1.0),Vector<Float>(1,Ton_upper));
    AffineFunction heater_turn_on_invariant(Matrix<Float>(1,2,0.0,-1.0),Vector<Float>(1,Ton_lower));
    AffineFunction midnight_guard(Matrix<Float>(1,2,1.0,0.0),Vector<Float>(1,-1.0));
  
    // Create the system modes
    heating_system.new_mode(heater_on,heater_on_dynamic);
    heating_system.new_mode(heater_off,heater_off_dynamic);
  
    // Create the system transitions for switching the heater
    heating_system.new_invariant(heater_off,heater_turn_on_invariant);
    heating_system.new_unforced_transition(switch_on,heater_off,heater_on,heater_turn_on_reset,heater_turn_on_activation);
    heating_system.new_forced_transition(switch_off,heater_on,heater_off,heater_turn_off_reset,heater_turn_off_guard);

    // Create the system transitions for resetting the clock at midnight
    heating_system.new_forced_transition(midnight,heater_on,heater_on,midnight_reset,midnight_guard);
    heating_system.new_forced_transition(midnight,heater_off,heater_off,midnight_reset,midnight_guard);

    return heating_system;
}

HybridEvolver create_evolver()
{
    // Create a HybridEvolver object
    HybridEvolver evolver;

    // Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = 0.25;
    evolver.parameters().maximum_step_size = 0.0625;
    cout <<  evolver.parameters() << endl;

    return evolver;
}


void compute_evolution(const HybridAutomaton& heating_system, const HybridEvolver& evolver)
{
    // Redefine the two discrete states
    DiscreteState heater_on(1);
    DiscreteState heater_off(2);

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::EnclosureListType HybridEnclosureListType;
    typedef HybridEvolver::OrbitType OrbitType;

    // Define the initial set
    Box initial_box(2, 0.0,0.015625, 16.0,16.0625);
    HybridEnclosureType initial(heater_off,initial_box);
  
    // Set the maximum evolution time
    HybridTime evolution_time(1.5,4);
  
    // Compute reachable and evolved sets
    cout << "Computing evolved sets... " << flush;
    HybridEnclosureListType reach,evolve;
    make_lpair(reach,evolve)=evolver.reach_evolve(heating_system,initial,evolution_time,UPPER_SEMANTICS);
    plot("tutorial-reach_evolve.png",Box(2, 0.0,1.0, 14.0,19.0),
         Colour(0.0,0.5,1.0), reach, Colour(0.0,0.25,0.5), initial, Colour(0.25,0.0,0.5), evolve);
    cout << "done." << endl;

    // Compute the orbit.
    cout << "Computing orbit... " << flush;
    OrbitType orbit = evolver.orbit(heating_system,initial,evolution_time,UPPER_SEMANTICS);
    cout << "done." << endl;

    // Write the orbit to standard output and plot.
    write("tutorial-orbit.txt",orbit);
    plot("tutorial-orbit.png",Box(2, 0.0,1.0, 14.0,19.0), Colour(0.0,0.5,1.0), orbit);
}


void compute_reachable_sets(const HybridAutomaton& heating_system, const HybridEvolver& evolver)
{
    // Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.parameters().initial_grid_density=10;
    analyser.parameters().initial_grid_depth=12;
    analyser.parameters().maximum_grid_depth=12;


    // Define the initial set
    HybridImageSet initial_set;
    DiscreteState heater_off(2);
    Box initial_box(2, 0.0,0.015625/4, 16.0,16.0+0.0625/16);
    initial_set[heater_off]=initial_box;

    // Set the maximum evolution time
    HybridTime reach_time(1.5,4);


    // Compute lower-approximation to finite-time evolved set using lower-semantics.
    std::cout << "Computing lower evolve set... " << std::flush;
    HybridGridTreeSet* lower_evolve_set_ptr = analyser.lower_evolve(heating_system,initial_set,reach_time);
    std::cout << "done." << std::endl;

    // Compute lower-approximation to finite-time reachable set using lower-semantics.
    std::cout << "Computing lower reach set... " << std::flush;
    HybridGridTreeSet* lower_reach_set_ptr = analyser.lower_reach(heating_system,initial_set,reach_time);
    std::cout << "done." << std::endl;

    plot("tutorial-lower_reach_evolve.png",Box(2, 0.0,1.0, 14.0,19.0),
         Colour(0.0,0.5,1.0), *lower_reach_set_ptr,
         Colour(0.0,0.25,0.5), initial_set,
         Colour(0.25,0.0,0.5), *lower_evolve_set_ptr);

    // Compute over-approximation to finite-time evolved set using upper semantics.
    // Subdivision is used as necessary to keep the local errors reasonable. 
    // The accumulated global error may be very large.
    std::cout << "Computing upper evolve set... " << std::flush;
    HybridGridTreeSet* upper_evolve_set_ptr = analyser.upper_evolve(heating_system,initial_set,reach_time);
    std::cout << "done." << std::endl;

    // Compute over-approximation to finite-time reachable set using upper semantics.
    std::cout << "Computing upper reach set... " << std::flush;
    HybridGridTreeSet* upper_reach_set_ptr = analyser.upper_reach(heating_system,initial_set,reach_time);
    std::cout << "done." << std::endl;

    plot("tutorial-upper_reach_evolve.png",Box(2, 0.0,1.0, 14.0,19.0),
         Colour(0.0,0.5,1.0), *upper_reach_set_ptr,
         Colour(0.0,0.25,0.5), initial_set,
         Colour(0.25,0.0,0.5), *upper_evolve_set_ptr);

    // Compute over-approximation to infinite-time chain-reachable set using upper semantics.
    HybridGridTreeSet* chain_reach_set_ptr = 0;
    std::cout << "Computing chain reach set... " << std::flush;
    chain_reach_set_ptr = analyser.chain_reach(heating_system,initial_set);
    std::cout << "done." << std::endl;
    plot("tutorial-chain_reach.png",Box(2, 0.0,1.0, 14.0,19.0), Colour(0.0,0.5,1.0), *chain_reach_set_ptr);
}



void compute_reachable_sets_with_serialisation(const HybridAutomaton& heating_system, const HybridReachabilityAnalyser& analyser)
{
    // Define the initial set
    HybridImageSet initial_set;
    DiscreteState heater_off(2);
    Box initial_box(2, 0.0,0.015625, 16.0,16.0625);
    initial_set[heater_off]=initial_box;


    // Compute the reach set for times between tlower and tupper.
    // The intermediate set is stored to an archive file and used to build the initial set for the reach step
    // Note that because of peculiarities in the Boost serialization library, 
    // the object to be serialized must be declared const.
    Float tlower=0.25; Float tupper=0.75;
    HybridTime transient_time(tlower,4);
    HybridTime recurrent_time(tupper-tlower,16);

    const HybridGridTreeSet* upper_intermediate_set_ptr = analyser.upper_evolve(heating_system,initial_set,transient_time);
    const HybridGridTreeSet upper_intermediate_set = *upper_intermediate_set_ptr;
    plot("tutorial-upper_intermediate.png",Box(2, 0.0,1.0, 14.0,18.0), Colour(0.0,0.5,1.0), *upper_intermediate_set_ptr);

    std::ofstream output_file_stream("tutorial-transient.txt");
    text_oarchive output_archive(output_file_stream);
    //output_archive << *upper_intermediate_set_ptr;
    output_archive << upper_intermediate_set;
    output_file_stream.close();

    HybridGridTreeSet rebuilt_upper_intermediate_set;

    std::ifstream input_file_stream("tutorial-transient.txt");
    text_iarchive input_archive(input_file_stream);
    input_archive >> rebuilt_upper_intermediate_set;
    input_file_stream.close();

    HybridGridTreeSet* upper_recurrent_set_ptr = analyser.upper_reach(heating_system,initial_set,recurrent_time);
    plot("tutorial-upper_recurrent.png",Box(2, 0.0,1.0, 14.0,18.0), Colour(0.0,0.5,1.0), *upper_recurrent_set_ptr);
}




int main() 
{
    // Create the system
    HybridAutomaton heating_system=create_heating_system();
    
    // Create the analyser classes
    HybridEvolver evolver=create_evolver();
    HybridReachabilityAnalyser reachability_analysier(evolver);
    
    // Compute the system evolution
    compute_evolution(heating_system,evolver);
    compute_reachable_sets(heating_system,evolver);
}
