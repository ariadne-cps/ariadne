/***************************************************************************
 *            robotarm.cc
 *
 *  Copyright  2013  Davide Bresolin
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
#include "expression/valuation.h"
#include "numeric/operators.h"

using namespace Ariadne;

// Declare the type to be used for the system evolution
typedef GeneralHybridEvolver::EnclosureType HybridEnclosureType;

// Auxiliary function for plotting
Void plot_results(const Orbit<HybridEnclosure>& orbit,
  double delta, double time, const RealVariable& t,
  const RealVariable& x, const RealVariable& vx,
  const RealVariable& z, const RealVariable& vz,
  const RealVariable& phi, const RealVariable& vphi) {
    // Print the intial, evolve and reach sets
    // cout << "Plotting sets" << endl;
    // cout << "evolve_set=" << hybrid_evolve_set << endl;
    // cout << "reach_set=" << hybrid_reach_set << endl;
    std::cout << "Plotting..." << std::flush;
    HybridFigure fig;
    char filename[50];
    // Plotting x
    fig.set_axes(Axes2d(0.0<=t<=time,0.0<=x<=1.0));
    fig.set_bounds(t,0.0,time);
    fig << line_style(true) << fill_colour(cyan);
    fig << orbit;
    fig << fill_colour(red) << orbit.final();
//    fig << fill_colour(blue) << initial_set;
    sprintf(filename,"robotarm-verification-%f-x.png",delta);
    fig.write(filename);

    // Plotting vx
    fig.clear();
    fig.set_axes(Axes2d(0.0<=t<=time,-0.2<=vx<=0.5));
    fig.set_bounds(t,0.0,time);
    fig << line_style(true) << fill_colour(cyan);
    fig << orbit;
    fig << fill_colour(red) << orbit.final();
//    fig << fill_colour(blue) << initial_set;
    sprintf(filename,"robotarm-verification-%f-vx.png",delta);
    fig.write(filename);


    // Plotting z
    fig.clear();
    fig.set_axes(Axes2d(0.0<=t<=time,0.0<=z<=0.2));
    fig.set_bounds(t,0.0,time);
    fig << line_style(true) << fill_colour(cyan);
    fig << orbit;
    fig << fill_colour(red) << orbit.final();
//    fig << fill_colour(blue) << initial_set;
    sprintf(filename,"robotarm-verification-%f-z.png",delta);
    fig.write(filename);

    // Plotting vz
    fig.clear();
    fig.set_axes(Axes2d(0.0<=t<=time,-0.02<=vz<=0.1));
    fig.set_bounds(t,0.0,time);
    fig << line_style(true) << fill_colour(cyan);
    fig << orbit;
    fig << fill_colour(red) << orbit.final();
//    fig << fill_colour(blue) << initial_set;
    sprintf(filename,"robotarm-verification-%f-vz.png",delta);
    fig.write(filename);

    // Plotting phi
    fig.clear();
    fig.set_axes(Axes2d(0.0<=t<=time,0.0<=phi<=pi/2));
    fig.set_bounds(t,0.0,time);
    fig << line_style(true) << fill_colour(cyan);
    fig << orbit;
    fig << fill_colour(red) << orbit.final();
//    fig << fill_colour(blue) << initial_set;
    sprintf(filename,"robotarm-verification-%f-phi.png",delta);
    fig.write(filename);

    // Plotting vphi
    fig.clear();
    fig.set_axes(Axes2d(0.0<=t<=time,-0.7<=vphi<=0.1));
    fig.set_bounds(t,0.0,time);
    fig << line_style(true) << fill_colour(cyan);
    fig << orbit;
    fig << fill_colour(red) << orbit.final();
//    fig << fill_colour(blue) << initial_set;
    sprintf(filename,"robotarm-verification-%f-vphi.png",delta);
    fig.write(filename);

    std::cout << " done." << endl;
}

Int main(Int argc, char** argv)
{
    Nat verbosity = 0;
    if(argc>1) { verbosity=atoi(argv[1]); }

    // Redirect std::cerr to file
    std::ofstream errlog("robotarm-verification-errors.log");
    std::cerr.rdbuf(errlog.rdbuf());

    DRAWING_METHOD = AFFINE_DRAW;
    DRAWING_ACCURACY = 0;

    // Set up the evolution parameters and grid
    double time(5.0);
    Int steps = 3;
    Real step_size(1e-1);

    // Definition of the reference trajectory
    // Constants
    Real period ( time );
    Real omega = 2*pi/period;

    // Variables
    RealVariable x("x");
    RealVariable vx("vx");
    RealVariable z("z");
    RealVariable vz("vz");
    RealVariable phi("phi");
    RealVariable vphi("vphi");
    RealVariable t("t");
    RealVariable xc("xc");
    RealVariable zc("zc");

    RealExpression dot_t=RealExpression::constant(1);

    // Algebraic expression for the reference trajectory
    RealExpression xd = -1/(2*pi)*sin(omega*t) + 1/period*t;
    RealExpression zd = -1/(2*pi*period)*sin(omega*t) + 1/(period*period)*t;
    RealExpression phid = sin(omega*t)/4 - pi/(2*period)*t+pi/2;
    std::cout << "Algebraic expression for xd = " << xd << std::endl;
    std::cout << "Algebraic expression for zd = " << zd << std::endl;
    std::cout << "Algebraic expression for phid = " << phid << std::endl;

//   return 1;

    // First and second derivatives of the reference trajectory
    RealExpression dot_xd = -1/period*cos(omega*t) + 1/period;
    RealExpression dot_zd = -1/(period*period)*cos(omega*t) + 1/(period*period);
    RealExpression dot_phid = pi/(2*period)*cos(omega*t) - pi/(2*period);
    std::cout << "Expression for dot xd = " << dot_xd << std::endl;
    std::cout << "Expression for dot zd = " << dot_zd << std::endl;
    std::cout << "Expression for dot phid = " << dot_phid << std::endl;

    RealExpression ddot_xd = omega/period*sin(omega*t);
    RealExpression ddot_zd = omega/(period*period)*sin(omega*t);
    RealExpression ddot_phid = (pi*pi)/(period*period)*sin(omega*t);
    std::cout << "Expression for dot dot xd = " << ddot_xd << std::endl;
    std::cout << "Expression for dot dot zd = " << ddot_zd << std::endl;
    std::cout << "Expression for dot dot phid = " << ddot_phid << std::endl;

    // Definition of the arm system
    // Constants for the dynamics
    Real h_scaling = 2;
    Real m = 100;            Real b = h_scaling*500;            Real k = h_scaling*h_scaling*2500;
    Real ke = 1000;          Real kFx = h_scaling*h_scaling*10/m;         Real kFz = (h_scaling*h_scaling/4)*10/m;
    Real Fp = 25;

    // Constants for the contact point
    Real c(0.95_dec); Rational delta(0.05_dec);
    //
    // Definition of the Hybrid Automaton for the Robot Arm
    //
    HybridAutomaton robotarm_automaton;

    // First mode: free run
    StringVariable robotarm("robot_arm");
    DiscreteLocation free(robotarm|"free");
    Real fifth=0.2_dec; Real half=0.5_bin;
    DottedRealAssignments free_dyn(
       (dot(x)   = vx,   dot(vx)   = ddot_xd + fifth*b/m * (dot_xd - vx) + fifth*k/m * (xd - x),
        dot(z)   = vz,   dot(vz)   = ddot_zd + b/m * (dot_zd - vz) + k/m * (zd - z),
        dot(phi) = vphi, dot(vphi) = ddot_phid + b/m * (dot_phid - vphi) + k/m * (phid - phi),
        dot(t)   = 1,  dot(xc)   = 0,   dot(zc)  = 0));
    std::cout << "free_dyn = " << free_dyn << std::endl;
    robotarm_automaton.new_mode(free, free_dyn);
    DiscreteEvent finv("finv");
    robotarm_automaton.new_invariant(free, (x <= c + Real(delta)), finv);

    // Second mode: contact
    DiscreteLocation contact(robotarm|"contact");
    DottedRealAssignments contact_dyn(
       (dot(x)   = vx,   dot(vx)   = ddot_xd + fifth*b/m * (dot_xd - vx) + fifth*k/m * (xd - x) + kFx * ke * (xc - x),
        dot(z)   = vz,   dot(vz)   = ddot_zd + b/m * (dot_zd - vz) + k/m * (zd - z) + kFz * ke * (zc - z),
        dot(phi) = vphi, dot(vphi) = ddot_phid + b/m * (dot_phid - vphi) + k/m * (phid - phi),
        dot(t)   = 1,  dot(xc)   = 0,   dot(zc)  = 0));
    std::cout << "contact_dyn = " << contact_dyn << std::endl;
    robotarm_automaton.new_mode(contact, contact_dyn);

    // Third mode: puncturing
    DiscreteLocation punct(robotarm|"punct");
    DottedRealAssignments punct_dyn(
       (dot(x)   = vx,   dot(vx)   = fifth*b/m * (half - vx) + fifth*k/m * (xc + half*t - half*period - x) + kFx * ke * (xc - x),
        dot(z)   = vz,   dot(vz)   = b/m * (- vz) + k/m * (zc - z) + kFz * ke * (zc - z),
        dot(phi) = vphi, dot(vphi) = b/m * (- vphi) + k/m * (- phi),
        dot(t)   = 1,  dot(xc)   = 0,   dot(zc)  = 0));
    std::cout << "punct_dyn = " << punct_dyn << std::endl;
    robotarm_automaton.new_mode(punct, punct_dyn);
    DiscreteEvent pinv("pinv");
    robotarm_automaton.new_invariant(punct, (x <= Fp/ke + xc), pinv);

    // transition from free to contact
    DiscreteEvent f2c("f2c");
    PrimedRealAssignments reset_id(
        (next(x)=x, next(vx)=vx, next(z)=z, next(vz)=vz, next(phi)=phi, next(vphi)=vphi,
         next(t)=t, next(xc)=x,  next(zc)=z) );
    robotarm_automaton.new_transition(free, f2c, contact, reset_id, (x >= c - Real(delta)), permissive);

    // transition contact to punct
    DiscreteEvent c2p("c2p");
    robotarm_automaton.new_transition(contact, c2p, punct, reset_id, (t >= period), urgent);


    std::cout << "RobotArm = " << std::endl;
    std::cout << robotarm_automaton << std::endl << std::endl;

    //
    // COMPUTE THE EVOLUTION OF THE SYSTEM
    //

    // Define the initial box
    RealVariablesBox initial_box({x==0, vx==0, z==0, vz==0,
        phi==pi/2, vphi==0, t==0, xc==0, zc==0});

    cout << "initial_box=" << initial_box << endl;

    HybridSet initial_set(free, initial_box);
    cout << "initial_set=" << initial_set << endl << endl;

    Semantics semantics=UPPER_SEMANTICS;

    // safe set
    ExactIntervalType zsafe(0.19,0.20);
    // Inital range for delta
    Rational deltamin ( 0.0_dec );
    Rational deltamax ( 0.1_dec );
    // accuracy
    Float64Value eps = 1e-2_exact;
    Int iter=0; // iteration count

    Orbit<HybridEnclosure> orbit;
    while(deltamax - deltamin >= eps) {
        // Set up the evaluators
        GeneralHybridEvolver evolver(robotarm_automaton);
        evolver.verbosity = verbosity;
        evolver.configuration().set_maximum_enclosure_radius(10.0);
        evolver.configuration().set_maximum_step_size(2.5);
        evolver.configuration().set_enable_reconditioning(true);
        evolver.configuration().set_maximum_spacial_error(1e-3);

        if(iter==0) { std::cout << "Evolution parameters:" << evolver.configuration() << std::endl; }

        std::cout << "Computing evolution for delta = " << delta << " .... " <<  std::flush;
        // Compute the reachable sets
        double etime = 10.0;
        orbit = evolver.orbit(initial_set,HybridTime(etime,steps),semantics);
        cout << "done" << std::endl;

//    cout << "\norbit.final=\n" << orbit.final() << endl << endl;
//    cout << "\norbit=\n" << orbit << endl << endl;
        // Show the value of the bounding box of the final set for x and z
        UpperBoxType bbox = orbit.reach()[punct].bounding_box();
        std::cout << "Bounding box of the final set:" << std::endl;
        std::cout << "      x = " << bbox[0] << std::endl;
        std::cout << "      z = " << bbox[2] << std::endl;
        std::cout << "    phi = " << bbox[4] << std::endl;

        std::cout << "Safe set for z: " << zsafe << std::endl;

        // if(iter==0) {
        //plot_results(orbit, delta, time, t, x, vx, z, vz, phi, vphi); // }


        if(definitely(subset(bbox[2],zsafe))) {
            std::cout << "      system is safe, increasing delta." << std::endl  << std::endl;
            deltamin = delta;
        } else {
            std::cout << "      system is possibly unsafe, decreasing delta" << std::endl  << std::endl;
            deltamax = delta;
        }
        iter++;
        delta = (deltamax+deltamin)/2;
        // change the invariant of location free
        robotarm_automaton.set_invariant(free, (x <= c + Real(delta)), finv);
        robotarm_automaton.set_guard(free, f2c, (x >= c - Real(delta)), permissive);
    }

    std::cout << "Final value for delta is " << deltamin << ", found after " << iter << " iterations" << std::endl << std::endl;

    return 0;
}


