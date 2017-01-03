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

HybridAutomaton getLaserTrajectory()
{
    /// Build the Hybrid System

    /// Create a HybridAutomaton object
	HybridAutomaton automaton("trajectory");

    // Parameters
    RealConstant velocity("velocity",0.092_dec); // Velocity in modulus
    RealConstant width("width",0.0046_dec); // Width of the cut

    /// Modes

    DiscreteLocation scanning("scanning");

    // Variables

    RealVariable x("x"); // X position
    RealVariable vx("vx"); // X velocity

    DiscreteEvent switch_left("switch_left");
    DiscreteEvent switch_right("switch_right");

	automaton.new_mode(scanning, {dot(x)=vx,dot(vx)=0});

	automaton.new_guard(scanning,switch_right,x<=0,impact);
	automaton.new_update(scanning,switch_right,scanning,{next(x)=0,next(vx)=velocity});

	automaton.new_guard(scanning,switch_left,x>=width,impact);
	automaton.new_update(scanning,switch_left,scanning,{next(x)=width,next(vx)=-velocity});

	return automaton;
}

}

#endif
