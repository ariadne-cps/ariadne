/*****************************************************************************************************
 *            laser-trajectory.hpp
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the trajectory of the laser.
 *
 *****************************************************************************************************/

#include "ariadne.hpp"

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

	automaton.new_guard(scanning,switch_right,x<=0,EventKind::IMPACT);
	automaton.new_update(scanning,switch_right,scanning,{next(x)=0,next(vx)=velocity});

	automaton.new_guard(scanning,switch_left,x>=width,EventKind::IMPACT);
	automaton.new_update(scanning,switch_left,scanning,{next(x)=width,next(vx)=-velocity});

	return automaton;
}

}

#endif
