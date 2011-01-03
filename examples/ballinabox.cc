/***************************************************************************
 *            ballinabox.cc
 *
 *  An example for checking the behavior when more than
 *  one guard is active.
 *
 *  Copyright  2010  Luca Geretti
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

int main(int argc, const char* argv[])
{
    uint evolver_verbosity=0;
    if(argc>1) { evolver_verbosity=atoi(argv[1]); }

    typedef GeneralHybridEvolver GeneralHybridEvolverType;

    /// Set the system parameters
    Real ax = 0.75;  // Coefficient of restitution for bounces with vertical walls
    Real ay = 0.5;  // Coefficient of restitution for bounces with horizontal walls

    /// Set the position and velocity functions.
    RealScalarFunction x=RealScalarFunction::coordinate(4,0); // X position
    RealScalarFunction vx=RealScalarFunction::coordinate(4,1); // X velocity
    RealScalarFunction y=RealScalarFunction::coordinate(4,2); // Y position
    RealScalarFunction vy=RealScalarFunction::coordinate(4,3); // Y velocity

    /// Build the Hybrid System

    /// Create a HybridAutomton object
    MonolithicHybridAutomaton ball;

    /// Create the discrete states
    DiscreteLocation free("free");

    /// Create the discrete events
    DiscreteEvent b_yu("b_yu"); // When the ball bounced on the upper boundary for y
    DiscreteEvent b_yl("b_yl"); // When the ball bounced on the lower boundary for y
    DiscreteEvent b_xu("b_xu"); // When the ball bounced on the upper boundary for x
    DiscreteEvent b_xl("b_xl"); // When the ball bounced on the lower boundary for x

    /// Create the dynamics
    RealVectorFunction free_d((vx,0,vy,0));

    /// Create the resets
    RealVectorFunction x_r((x,-ax*vx,y,vy)); // Bounces on a boundary for x
    RealVectorFunction y_r((x,vx,y,-ay*vy)); // Bounces on a boundary for y

    /// Create the guards
    /// Guards are true when g(x) > 0
    RealScalarFunction xu_g(x-1); // For when the upper boundary for x is reached
    RealScalarFunction xl_g(-1-x); // For when the lower boundary for x is reached
    RealScalarFunction yu_g(y-1); // For when the upper boundary for y is reached
    RealScalarFunction yl_g(-1-y); // For when the lower boundary for y is reached

    ball.new_mode(free,free_d);

    ball.new_transition(free,b_yu,free,y_r,yu_g,impact); // Bounces on the upper boundary for y
    ball.new_transition(free,b_xu,free,x_r,xu_g,impact); // Bounces on the upper boundary for x
	ball.new_transition(free,b_yl,free,y_r,yl_g,impact); // Bounces on the lower boundary for y
    ball.new_transition(free,b_xl,free,x_r,xl_g,impact); // Bounces on the lower boundary for x

    /// Finished building the automaton

    cout << "Automaton = " << ball << endl << endl;
    /// Compute the system evolution

    /// Create a GeneralHybridEvolver object
    GeneralHybridEvolverType evolver;
    evolver.verbosity=evolver_verbosity;

    TaylorModelAccuracy::set_default_maximum_degree(5);
    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = 1.0;
    evolver.parameters().maximum_step_size = 16.0;
    //Uncomment out to use small time-steps
    //evolver.parameters().maximum_step_size = 1.0/16;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolverType::EnclosureType EnclosureType;
    typedef GeneralHybridEvolverType::EnclosureListType EnclosureListType;
    typedef GeneralHybridEvolverType::OrbitType OrbitType;

    std::cout << "Computing evolution..." << std::endl;

    Box initial_box(4, 0.48,0.52, 2.5,2.5, 0.78,0.81, 1.0,1.0);
    EnclosureType initial_enclosure(free,initial_box);
    Box bounding_box(4, -1.0,1.0, -10.0,10.0, -1.0,1.0, -10.0,10.0);

    HybridTime evolution_time(10.0,3);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(ball,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit.final()="<<orbit.final()<<std::endl;
    Array<uint> xy(4,0,2);
    plot("ballinabox-orbit",PlanarProjectionMap(4,0,2),bounding_box, Colour(0.0,0.5,1.0), orbit);

}
