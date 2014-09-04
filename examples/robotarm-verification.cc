/***************************************************************************
 *            robotarm.cc
 *
 *  Copyright  2010  Davide Bresolin
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
#include "valuation.h"
#include "operators.h"

using namespace Ariadne;

// Declare the type to be used for the system evolution
typedef GeneralHybridEvolver::EnclosureType HybridEnclosureType;

int main(int argc, char** argv)
{
    uint verbosity = 0;
    if(argc>1) { verbosity=atoi(argv[1]); }

    // Set up the evolution parameters and grid
    double time(5.0);
    int steps = 3;
    Real step_size(1e-1);


    // Definition of the reference trajectory
    // Constants
    Real per = Real(time);
    Real omega = 2*pi/per;
    // Variables
    RealVariable x("x");
    RealVariable vx("vx");
    RealVariable t("t");

    RealExpression dot_t=RealExpression::constant(1);

    // Algebraic expression for the reference trajectory
    RealExpression ddot_xd = omega/per*sin(omega*t);
    std::cout << "Expression for dot dot xd = " << ddot_xd << std::endl;

//   return 1;

    // First and second derivatives of the reference trajectory
//    RealExpression dot_xd = 5.0*Ax*t*t*t*t + 4.0*Bx*t*t*t + 3.0*Cx*t*t;
    RealExpression dot_xd = -1/per*cos(omega*t) + 1/per;
    std::cout << "Expression for dot xd = " << dot_xd << std::endl;

//    RealExpression ddot_xd = 20.0*Ax*t*t*t + 12.0*Bx*t*t + 6.0*Cx*t;
    RealExpression xd = -1/(2*pi)*sin(omega*t) + 1/per*t;
    std::cout << "Algebraic expression for xd = " << xd << std::endl;


    // Definition of the arm system
    // Constants for the dynamics
    Real h_scaling = 2;
    Real m = 100;            Real b = h_scaling*500;            Real k = h_scaling*h_scaling*2500;
    Real ke = 1000;          Real kFx = h_scaling*h_scaling*10/m;         Real kFz = (h_scaling*h_scaling/4)*10/m;

    // Constants for the contact point
    Real xc = 0.95_dec;
    Real delta = 0.03_dec;

    // Dynamics for mode free
    Real fifth=0.2_dec;
    RealExpression dot_x = vx;
    RealExpression dot_vx = ddot_xd + fifth*b/m * (dot_xd - vx) + fifth*k/m * (xd - x);
    std::cout << "Expression for dot vx = " << dot_vx << std::endl;

    // Dynamics for mode contact
    RealExpression cdot_vx = ddot_xd + fifth*b/m * (dot_xd - vx) + fifth*k/m * (xd - x) + kFx * ke * (xc - x); // * cos(phi);

    //
    // Definition of the Hybrid Automaton for the Robot Arm
    //
    HybridAutomaton robotarm_automaton;

    // First mode: free run
    StringVariable robotarm("robot_arm");
    DiscreteLocation free(robotarm|"free");
    DottedRealAssignments free_dyn(( dot(x)=vx, dot(vx) = ddot_xd + fifth*b/m * (dot_xd - vx) + fifth*k/m * (xd - x), dot(t)=1.0));
    std::cout << "free_dyn = " << free_dyn << std::endl;
    robotarm_automaton.new_mode(free, free_dyn);
    DiscreteEvent finv("finv");
    robotarm_automaton.new_invariant(free, (x <= xc + delta), finv);

    // Second mode: contact
    DiscreteLocation contact(robotarm|"contact");
    robotarm_automaton.new_mode(contact, (dot(x)=vx,dot(vx)=ddot_xd + fifth*b/m * (dot_xd - vx) + fifth*k/m * (xd - x) + kFx * ke * (xc - x),dot(t)=1));

    // transition from free to contact
    DiscreteEvent f2c("f2c");
    PrimedRealAssignments reset_id( (next(x)=x,next(vx)=vx,next(t)=t) );
    robotarm_automaton.new_transition(free, f2c, contact, reset_id, (x >= xc - delta), permissive);

    std::cout << "RobotArm = " << std::endl;
    std::cout << robotarm << std::endl << std::endl;


    //
    // COMPUTE THE EVOLUTION OF THE SYSTEM
    //

    // Set up the evaluators
    GeneralHybridEvolver evolver(robotarm_automaton);
    evolver.verbosity = verbosity;
    evolver.configuration().set_maximum_enclosure_radius(10.0);
    evolver.configuration().set_maximum_step_size(2.5);
    evolver.configuration().set_enable_reconditioning(true);
    evolver.configuration().set_maximum_spacial_error(1e-3);

    std::cout << "Evolution parameters:" << evolver.configuration() << std::endl;

    // Define the initial box
    RealVariablesBox initial_box((x==0,vx==0,t==0));

    cout << "initial_box=" << initial_box << endl;


    HybridSet initial_set(free, initial_box);
    cout << "initial_set=" << initial_set << endl << endl;

    Semantics semantics=UPPER_SEMANTICS;

    // Compute the reachable sets
    Orbit<HybridEnclosure> orbit = evolver.orbit(initial_set,HybridTime(time,steps),semantics);
    cout << std::endl;

    cout << "\norbit.final=\n" << orbit.final() << endl << endl;
    // cout << "\norbit=\n" << orbit << endl << endl;

    // Print the intial, evolve and reach sets
    // cout << "Plotting sets" << endl;
    // cout << "evolve_set=" << hybrid_evolve_set << endl;
    // cout << "reach_set=" << hybrid_reach_set << endl;
    std::cout << "Plotting..." << std::flush;
    HybridFigure fig;
    // Plotting x
    fig.set_axes(Axes2d(0.0<=t<=time,0.0<=x<=1.0));
    fig.set_bounds(t,0.0,time);
    fig << line_style(true) << fill_colour(cyan);
    fig << orbit;
    fig << fill_colour(red) << orbit.final();
//    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-verification-orbit-x");
    fig.clear();

    // Plotting vx
    fig.set_axes(Axes2d(0.0<=t<=time,-1.0<=vx<=1.0));
    fig.set_bounds(t,0.0,time);
    fig << line_style(true) << fill_colour(cyan);
    fig << orbit;
    fig << fill_colour(red) << orbit.final();
//    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-verification-orbit-vx");


    std::cout << " done." << endl;

    return 0;

}
