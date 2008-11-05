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

#include <cstdarg>
#include "ariadne.h"

using namespace Ariadne;
using Ariadne::cos;



Float pi_=Ariadne::pi<Float>();



//! Declare nonlinear function objects

struct HeaterOn : FunctionData<2,2,4> {
  template<class R, class A, class P> static void 
  compute(R& r, const A& x, const P& p) {
    //    r[0] = p[0] + p[1] * ( p[2] - p[3] * Ariadne::cos(2*pi_*x[1]) - x[0] );
    typename R::value_type t0=2*pi_*x[0];
    r[0] = 1.0; 
    r[1] = p[0] + p[1] * ( p[2] - p[3] * (1-sqr(t0)/2) - x[1] );
  }
};

struct HeaterOff : FunctionData<2,2,4> {
  template<class R, class A, class P> static void 
  compute(R& r, const A& x, const P& p) {
    //    r[0] = p[1] * ( p[2] - p[3] * Ariadne::cos(2*pi_*x[1]) - x[0] );
    typename R::value_type t0=2*pi_*x[0];
    r[0] = 1.0; 
    r[1] = p[1] * ( p[2] - p[3] * (1-sqr(t0)/2) - x[1] );
  }  
};



int main() 
{
  //! The automaton has two modes, both 
  //! with two variables, room_temperature and time_of_day
  //! The dynamics is given by the differential equations
  //! \f$ \dot{T} = K(T_\mathrm{ext}(t) - T(t)) + P \delta(q(t) = \mathsc{on}) \f$
  //! where \f$T_\mathrm{ext}(t) = T_\mathrm{av} - A \cos(2\pi t)\f$ is the
  //! external temperature and \f$ q(t)\in\{\mathsc{on},\mathsc{off}\}\f$ is the
  //! discrete state of the heater.
  //! The parameters are the insulation coefficient \f$k\f$, the average 
  //! external temperature \f$T_\mathrm{av}\f$, the amplitude of the 
  //! temperature fluctuations \f$A\f$ and the power of the heater \f$P\f$.
  //!
  //! The heater is controlled by a thermostat. The heater turns off whenever
  //! the temperature rises above a fixed threshold \f$T_\mathrm{off}\f$ and
  //! turns on nondeterministically for a temperature \f$T\f$ between 
  //! \f$T^+_\mathrm{on}\f$ and \f$T^-_\mathrm{on}\f$.
  //! The clock time \f$\tau\f$ is reset to zero whenever \f$\tau\f$ becomes 
  //! equal to one.
 
  //! The system variables are
  //!  x[0] = \tau
  //!  x[1] = T

  //! The system parameters are
  //!  p[0] = P
  //!  p[1] = K
  //!  p[2] = T_{av}
  //!  p[3] = A
  
  //! The hard-coded parameters are
  //!  p[4] = T_{off}
  //!  p[5] = T^+_{on}
  //!  p[6] = T^-_{on}

  
  
  //! Set the system parameters
  Vector<Interval> system_parameters(4);
  system_parameters[0] = 1.0;
  system_parameters[1] = 1.0;
  system_parameters[2] = 16.0;
  system_parameters[3] = 8.0;
  
  Float Toff=22.0;
  Float Ton_upper=16.0;
  Float Ton_lower=15.0;

  //! Build the Hybrid System
  
  //! Create a HybridAutomton object
  HybridAutomaton heating_system;
  
  //! Create the two discrete state 
  DiscreteState heater_on(1);
  DiscreteState heater_off(2);
  
  //! Create the discrete events
  DiscreteEvent switch_on(1);
  DiscreteEvent switch_off(2);
  DiscreteEvent midnight(3);
  
  //! Create the dynamics
  Function<HeaterOn> heater_on_dynamic(system_parameters);
  Function<HeaterOff> heater_off_dynamic(system_parameters);

  //! Create the resets
  IdentityFunction heater_turn_off_reset(2);
  IdentityFunction heater_turn_on_reset(2);
  AffineFunction midnight_reset(Matrix<Float>(2,2,0.0,0.0,0.0,1.0),Vector<Float>(2,0.0,0.0));

  //! Create the guards.
  //! Note that an affine predicate is of the form \f$ Ax+b \gtrless 0 \f$, 
  //! and a guard is active (i.e. a discrete event is allowed/forced) when the value is positive. 
  //! Hence the condition \f$x_i>c\f$ is given by \f$e_i\cdot x + (-c) > 0\f$, and the condition
  //! \f$ x_i<c\f$ by \f$ -e_i\cdot x + c > 0\f$.
  AffineFunction heater_turn_off_guard(Matrix<Float>(1,2,0.0,1.0),Vector<Float>(1,-Toff));
  AffineFunction heater_turn_on_activation(Matrix<Float>(1,2,0.0,-1.0),Vector<Float>(1,Ton_upper));
  AffineFunction heater_turn_on_invariant(Matrix<Float>(1,2,0.0,-1.0),Vector<Float>(1,Ton_lower));
  AffineFunction midnight_guard(Matrix<Float>(1,2,1.0,0.0),Vector<Float>(1,-1.0));
  
  //! Build the automaton
  heating_system.new_mode(heater_on,heater_on_dynamic);
  heating_system.new_mode(heater_off,heater_off_dynamic);
  
  //heating_system.new_invariant(heater_off,heater_turn_on_invariant);
  //heating_system.new_unforced_transition(switch_on,heater_off,heater_on,heater_turn_on_reset,heater_turn_on_activation);
  //heating_system.new_forced_transition(switch_off,heater_on,heater_off,heater_turn_off_reset,heater_turn_off_guard);

  //heating_system.new_forced_transition(midnight,heater_on,heater_on,midnight_reset,midnight_guard);
  //heating_system.new_forced_transition(midnight,heater_off,heater_off,midnight_reset,midnight_guard);

  //! Compute the system evolution


  //! Create a HybridEvolver object
  HybridEvolver evolver;

  //! Set the evolution parameters
  evolver.parameters().maximum_enclosure_radius = 0.25;
  evolver.parameters().maximum_step_size = 0.0625;
  std::cout <<  evolver.parameters() << std::endl;

  // Declare the type to be used for the system evolution
  typedef HybridEvolver::EnclosureType HybridEnclosureType;
  typedef HybridEvolver::OrbitType OrbitType;

  Box initial_box(2, 0.0,0.015625, 16.0,16.0625);
  HybridEnclosureType initial_enclosure(heater_off,initial_box);
  
  HybridTime evolution_time(0.25,4);
  
  std::cout << "Computing orbit... " << std::flush;
  OrbitType orbit = evolver.orbit(heating_system,initial_enclosure,evolution_time,UPPER_SEMANTICS);
  std::cout << "done." << std::endl;

  std::cout << "Orbit="<<orbit<<std::endl;
  std::cout << "Orbit="<<orbit<<std::endl;
  //plot("tutorial-orbit.png",Box(2, 0.0,1.0, 14.0,18.0), Colour(0.0,0.5,1.0), orbit.initial());
  plot("tutorial-orbit.png",Box(2, 0.0,1.0, 14.0,18.0), Colour(0.0,0.5,1.0), orbit);



  //! Create a ReachabilityAnalyser object
  HybridReachabilityAnalyser analyser(evolver);

  HybridImageSet initial_set;
  initial_set[heater_off]=initial_box;

  HybridTime reach_time(0.25,4);

  HybridGridTreeSet* upper_reach_set_ptr = analyser.upper_reach(heating_system,initial_set,reach_time);
  plot("tutorial-initial_set.png",Box(2, 0.0,1.0, 14.0,18.0), Colour(0.0,0.5,1.0), initial_set);
  plot("tutorial-upper_reach.png",Box(2, 0.0,1.0, 14.0,18.0), Colour(0.0,0.5,1.0), *upper_reach_set_ptr);

}
