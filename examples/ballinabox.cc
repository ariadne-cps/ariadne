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
#include "hybrid/hybrid_graphics.h"

using namespace Ariadne;

int main(int argc, const char* argv[])
{
    uint evolver_verbosity=0;
    if(argc>1) { evolver_verbosity=atoi(argv[1]); }

    typedef GeneralHybridEvolver GeneralHybridEvolverType;

    /// Set the system parameters
    Real ax = 0.75_bin;  // Coefficient of restitution for bounces with vertical walls
    Real ay = 0.5_bin;  // Coefficient of restitution for bounces with horizontal walls

    /// Set the position and velocity functions.
    RealVariable x("x"); // X position
    RealVariable vx("vx"); // X velocity
    RealVariable y("y"); // Y position
    RealVariable vy("vy"); // Y velocity

    /// Build the Hybrid System

    /// Create a HybridAutomton object
    HybridAutomaton ball;

    /// Create the discrete states
    DiscreteLocation free;

    /// Create the discrete events
    DiscreteEvent b_yu("b_yu"); // When the ball bounced on the upper boundary for y
    DiscreteEvent b_yl("b_yl"); // When the ball bounced on the lower boundary for y
    DiscreteEvent b_xu("b_xu"); // When the ball bounced on the upper boundary for x
    DiscreteEvent b_xl("b_xl"); // When the ball bounced on the lower boundary for x

    /// Create the dynamics
    DottedRealAssignments free_d((dot(x)=vx,dot(vx)=0,dot(y)=vy,dot(vy)=0));

    /// Create the resets
    PrimedRealAssignments x_r((next(x)=x,next(vx)=-ax*vx,next(y)=y,next(vy)=vy)); // Bounces on a boundary for x
    PrimedRealAssignments y_r((next(x)=x,next(vx)=vx,next(y)=y,next(vy)=-ay*vy)); // Bounces on a boundary for y

    /// Create the guards
    /// Guards are true when g(x) > 0
    ContinuousPredicate xu_g(x>=1); // For when the upper boundary for x is reached
    ContinuousPredicate xl_g(x<=-1); // For when the lower boundary for x is reached
    ContinuousPredicate yu_g(y>=1);// For when the upper boundary for y is reached
    ContinuousPredicate yl_g(y<=-1); // For when the lower boundary for y is reached

    ball.new_mode(free_d);

    ball.new_transition(b_yu,y_r,yu_g,impact); // Bounces on the upper boundary for y
    ball.new_transition(b_xu,x_r,xu_g,impact); // Bounces on the upper boundary for x
	ball.new_transition(b_yl,y_r,yl_g,impact); // Bounces on the lower boundary for y
    ball.new_transition(b_xl,x_r,xl_g,impact); // Bounces on the lower boundary for x

    /// Finished building the automaton

    cout << "Automaton = " << ball << endl << endl;
    /// Compute the system evolution

    /// Create a GeneralHybridEvolver object
    GeneralHybridEvolverType evolver(ball);
    typedef GeneralHybridEvolver::OrbitType OrbitType;
    evolver.verbosity=evolver_verbosity;

    /// Set the evolution parameters
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(16.0);
    //Uncomment out to use small time-steps
    //evolver.parameters().maximum_step_size = 1.0/16;
    std::cout <<  evolver.configuration() << std::endl;

    std::cout << "Computing evolution..." << std::endl;

    Decimal x0l(0.48), x0u(0.52), vx0(2.5), y0l(0.78), y0u(0.81), vy0(1.0);
    RealVariablesBox initial_box((x0l<=x<=x0u, vx0<=vx<=vx0, y0l<=y<=y0u, vy0<=vy<=vy0));
    HybridSet initial_set(free,initial_box);

    HybridTime evolution_time(10.0,3);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(initial_set,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit.final()="<<orbit.final()<<std::endl;
    Axes2d axes(-1.0<=x<=1.0,-1.0<=y<=1.0);
    plot("ballinabox-orbit",axes, Colour(0.0,0.5,1.0), orbit);

}
