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
#include "valuation.h"
#include "operators.h"

using namespace Ariadne;

// Declare the type to be used for the system evolution
typedef GeneralHybridEvolver::EnclosureType HybridEnclosureType;

// Auxiliary function for plotting
void plot_results(const Orbit<HybridEnclosure>& orbit, 
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

int main(int argc, char** argv)
{
    uint verbosity = 0;
    if(argc>1) { verbosity=atoi(argv[1]); }

    // Redirect std::cerr to file
    std::ofstream errlog("robotarm-verification-errors.log");
    std::cerr.rdbuf(errlog.rdbuf());

    DRAWING_METHOD = AFFINE_DRAW;
    DRAWING_ACCURACY = 0;

    // Set up the evolution parameters and grid
    double time(5.0);
    int steps = 3;
    Real step_size(1e-1);

    // Definition of the reference trajectory
    // Constants
    Real Pi = 3.1415926535897931;
    // Constants
    Real per = time;
    Real omega = 2*pi/per;

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

    RealExpression dot_t((1.0));

    // Algebraic expression for the reference trajectory
    RealExpression xd = -1.0/(2*pi)*sin(omega*t) + 1/per*t;
    RealExpression zd = -1.0/(2*pi*per)*sin(omega*t) + 1/(per*per)*t;
    RealExpression phid = 1.0/4*sin(omega*t) - pi/(2*per)*t+pi/2;
    std::cout << "Algebraic expression for xd = " << xd << std::endl;
    std::cout << "Algebraic expression for zd = " << zd << std::endl;
    std::cout << "Algebraic expression for phid = " << phid << std::endl;
   
//   return 1;
   
    // First and second derivatives of the reference trajectory
    RealExpression dot_xd = -1/per*cos(omega*t) + 1/per;
    RealExpression dot_zd = -1/(per*per)*cos(omega*t) + 1/(per*per);
    RealExpression dot_phid = pi/(2*per)*cos(omega*t) - pi/(2*per);       
    std::cout << "Expression for dot xd = " << dot_xd << std::endl;
    std::cout << "Expression for dot zd = " << dot_zd << std::endl;
    std::cout << "Expression for dot phid = " << dot_phid << std::endl;

    RealExpression ddot_xd = omega/per*sin(omega*t);
    RealExpression ddot_zd = omega/(per*per)*sin(omega*t);
    RealExpression ddot_phid = (pi*pi)/(per*per)*sin(omega*t);       
    std::cout << "Expression for dot dot xd = " << ddot_xd << std::endl;
    std::cout << "Expression for dot dot zd = " << ddot_zd << std::endl;
    std::cout << "Expression for dot dot phid = " << ddot_phid << std::endl;

    // Definition of the arm system
    // Constants for the dynamics
    Real h_scaling = 2.0;
    Real m = 100.0;            Real b = h_scaling*500.0;            Real k = h_scaling*h_scaling*2500.0;
    Real ke = 1000.0;          Real kFx = h_scaling*h_scaling*10.0/m;         Real kFz = (h_scaling*h_scaling/4.0)*10.0/m;
    Real Fp = 25.0;

    // Constants for the contact point
    double c = 0.95;             double delta=0.05;
    Real rc = c;                
    Real rdelta = delta;

    //
    // Definition of the Hybrid Automaton for the Robot Arm
    //
    HybridAutomaton robotarm_automaton;

    // First mode: free run
    StringVariable robotarm("robot_arm");
    DiscreteLocation free(robotarm|"free");
    DottedRealAssignments free_dyn(
       (dot(x)   = vx,   dot(vx)   = ddot_xd + 0.2*b/m * (dot_xd - vx) + 0.2*k/m * (xd - x), 
        dot(z)   = vz,   dot(vz)   = ddot_zd + b/m * (dot_zd - vz) + k/m * (zd - z),
        dot(phi) = vphi, dot(vphi) = ddot_phid + b/m * (dot_phid - vphi) + k/m * (phid - phi),
        dot(t)   = 1.0,  dot(xc)   = 0.0,   dot(zc)  = 0.0));
    std::cout << "free_dyn = " << free_dyn << std::endl;
    robotarm_automaton.new_mode(free, free_dyn);
    DiscreteEvent finv("finv");
    robotarm_automaton.new_invariant(free, (x <= rc + rdelta), finv);

    // Second mode: contact
    DiscreteLocation contact(robotarm|"contact");
    DottedRealAssignments contact_dyn(
       (dot(x)   = vx,   dot(vx)   = ddot_xd + 0.2*b/m * (dot_xd - vx) + 0.2*k/m * (xd - x) + kFx * ke * (xc - x), 
        dot(z)   = vz,   dot(vz)   = ddot_zd + b/m * (dot_zd - vz) + k/m * (zd - z) + kFz * ke * (zc - z),
        dot(phi) = vphi, dot(vphi) = ddot_phid + b/m * (dot_phid - vphi) + k/m * (phid - phi),
        dot(t)   = 1.0,  dot(xc)   = 0.0,   dot(zc)  = 0.0));
    std::cout << "contact_dyn = " << contact_dyn << std::endl;
    robotarm_automaton.new_mode(contact, contact_dyn);

    // Third mode: puncturing
    DiscreteLocation punct(robotarm|"punct");
    DottedRealAssignments punct_dyn(
       (dot(x)   = vx,   dot(vx)   = 0.2*b/m * (0.5 - vx) + 0.2*k/m * (xc + 0.5*t - 0.5*per - x) + kFx * ke * (xc - x), 
        dot(z)   = vz,   dot(vz)   = b/m * (- vz) + k/m * (zc - z) + kFz * ke * (zc - z),
        dot(phi) = vphi, dot(vphi) = b/m * (- vphi) + k/m * (- phi),
        dot(t)   = 1.0,  dot(xc)   = 0.0,   dot(zc)  = 0.0));
    std::cout << "punct_dyn = " << punct_dyn << std::endl;
    robotarm_automaton.new_mode(punct, punct_dyn);
    DiscreteEvent pinv("pinv");
    robotarm_automaton.new_invariant(punct, (x <= Fp/ke + xc), pinv);
    
    // transition from free to contact
    DiscreteEvent f2c("f2c");
    PrimedRealAssignments reset_id( 
        (next(x)=x, next(vx)=vx, next(z)=z, next(vz)=vz, next(phi)=phi, next(vphi)=vphi,
         next(t)=t, next(xc)=x,  next(zc)=z) );
    robotarm_automaton.new_transition(free, f2c, contact, reset_id, (x >= rc - rdelta), permissive);

    // transition contact to punct
    DiscreteEvent c2p("c2p");
    robotarm_automaton.new_transition(contact, c2p, punct, reset_id, (t >= per), urgent);


    std::cout << "RobotArm = " << std::endl;
    std::cout << robotarm_automaton << std::endl << std::endl;

    //
    // COMPUTE THE EVOLUTION OF THE SYSTEM
    //

    // Define the initial box
    RealVariablesBox initial_box((x==0.0, vx==0.0, z==0.0, vz==0.0, 
        phi==pi.get_d()/2, vphi==0.0, t==0.0, xc==0.0, zc==0.0));

    cout << "initial_box=" << initial_box << endl;

    HybridSet initial_set(free, initial_box);
    cout << "initial_set=" << initial_set << endl << endl;

    Semantics semantics=UPPER_SEMANTICS;

    // safe set
    Interval zsafe(0.19,0.20);
    // Inital range for delta
    double deltamin = 0.0;
    double deltamax = 0.1;
    // accuracy
    double eps = 1e-2;
    int iter=0; // iteration count
         
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
        
        std::cout << "Computing evolution for delta = " << rdelta << " .... " <<  std::flush;
        // Compute the reachable sets
        double etime = 10.0;
        orbit = evolver.orbit(initial_set,HybridTime(etime,steps),semantics);
        cout << "done" << std::endl;

//    cout << "\norbit.final=\n" << orbit.final() << endl << endl;
//    cout << "\norbit=\n" << orbit << endl << endl;
        // Show the value of the bounding box of the final set for x and z
        Box bbox = orbit.reach()[punct].bounding_box();
        std::cout << "Bounding box of the final set:" << std::endl;
        std::cout << "      x = " << bbox[0] << std::endl;
        std::cout << "      z = " << bbox[2] << std::endl;
        std::cout << "    phi = " << bbox[4] << std::endl;
        
        std::cout << "Safe set for z: " << zsafe << std::endl;

        // if(iter==0) { 
        //plot_results(orbit, delta, time, t, x, vx, z, vz, phi, vphi); // }


        if(subset(bbox[2],zsafe)) {
            std::cout << "      system is safe, increasing delta." << std::endl  << std::endl;
            deltamin = delta;
        } else {
            std::cout << "      system is possibly unsafe, decreasing delta" << std::endl  << std::endl;
            deltamax = delta;
        }
        iter++;
        delta = (deltamax+deltamin)/2;
        rdelta = delta;
        // change the invariant of location free
        robotarm_automaton.set_invariant(free, (x <= rc + rdelta), finv);
        robotarm_automaton.set_guard(free, f2c, (x >= rc - rdelta), permissive);
    }
    
    std::cout << "Final value for delta is " << deltamin << ", found after " << iter << " iterations" << std::endl << std::endl;

    return 0;
}


