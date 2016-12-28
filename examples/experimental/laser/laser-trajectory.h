/*****************************************************************************************************
 *            laser-trajectory.h
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the trajectory of the laser.
 *
 *****************************************************************************************************/

#include "ariadne.h"

#ifndef LASER_TRAJECTORY_H_
#define LASER_TRAJECTORY_H_

namespace Ariadne {

AtomicHybridAutomaton getLaserTrajectory()
{
    /// Build the Hybrid System

    /// Create a HybridAutomaton object
	AtomicHybridAutomaton automaton("trajectory");

    // Parameters
    RealConstant velocity("velocity",0.092_dec); // Velocity in modulus
    RealConstant width("width",0.0046_dec); // Width of the cut

    /// Modes

    AtomicDiscreteLocation scanning("scanning");

    // Variables

    RealVariable x("x"); // X position
    RealVariable vx("vx"); // X velocity

    DiscreteEvent switch_left("switch_left");
    DiscreteEvent switch_right("switch_right");

	automaton.new_mode(scanning, {dot(x)=vx,dot(vx)=0});

	automaton.new_guard(DiscreteLocation(automaton|scanning),switch_right,x<=0,impact);
	automaton.new_update(DiscreteLocation(automaton|scanning),switch_right,DiscreteLocation(automaton|scanning),{next(x)=0,next(vx)=velocity});
	automaton.new_transition(scanning,switch_left,scanning,{next(x)=width,next(vx)=-velocity},x>=width);
	//automaton.new_transition(scanning,switch_right,scanning,{next(x)=0,next(vx)=velocity},x<=0);

	return automaton;
}

}

#endif
