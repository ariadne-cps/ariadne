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

    // Constants for the contact point
    double c = 0.95;             double delta=0.03;
    Real rc = c;                
    RealExpression rdelta=Real(delta);

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

    // transition from free to contact
    DiscreteEvent f2c("f2c");
    PrimedRealAssignments reset_id( 
        (next(x)=x, next(vx)=vx, next(z)=z, next(vz)=vz, next(phi)=phi, next(vphi)=vphi,
         next(t)=t, next(xc)=x,  next(zc)=z) );
    robotarm_automaton.new_transition(free, f2c, contact, reset_id, (x >= rc - rdelta), permissive);

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
    RealVariablesBox initial_box((x==0.0, vx==0.0, z==0.0, vz==0.0, 
        phi==pi.get_d()/2, vphi==0.0, t==0.0, xc==0.0, zc==0.0));

    cout << "initial_box=" << initial_box << endl;


    HybridSet initial_set(free, initial_box);
    cout << "initial_set=" << initial_set << endl << endl;

    Semantics semantics=UPPER_SEMANTICS;

    // Compute the reachable sets
    cout << "Computing evolution ... " << std::flush;
    Orbit<HybridEnclosure> orbit = evolver.orbit(initial_set,HybridTime(time,steps),semantics);
    cout << "done" << std::endl;

//    cout << "\norbit.final=\n" << orbit.final() << endl << endl;
//    cout << "\norbit=\n" << orbit << endl << endl;
    // Show the value of the bounding box of the final set for x and z
    Box bbox = orbit.reach()[free].bounding_box();
    std::cout << "Bounding box for location free:" << std::endl;
    std::cout << "      x = " << bbox[0] << std::endl;
    std::cout << "      z = " << bbox[2] << std::endl;
    std::cout << "    phi = " << bbox[4] << std::endl << std::endl;
    bbox = orbit.reach()[contact].bounding_box();
    std::cout << "Bounding box for location contact:" << std::endl;
    std::cout << "      x = " << bbox[0] << std::endl;
    std::cout << "      z = " << bbox[2] << std::endl;
    std::cout << "    phi = " << bbox[4] << std::endl << std::endl;

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

    // Plotting vx
    fig.clear();
    fig.set_axes(Axes2d(0.0<=t<=time,-0.2<=vx<=0.5));
    fig.set_bounds(t,0.0,time);
    fig << line_style(true) << fill_colour(cyan);
    fig << orbit;
    fig << fill_colour(red) << orbit.final();
//    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-verification-orbit-vx");

    
    // Plotting z
    fig.clear();
    fig.set_axes(Axes2d(0.0<=t<=time,0.0<=z<=0.2));
    fig.set_bounds(t,0.0,time);
    fig << line_style(true) << fill_colour(cyan);
    fig << orbit;
    fig << fill_colour(red) << orbit.final();
//    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-verification-orbit-z");

    // Plotting vz
    fig.clear();
    fig.set_axes(Axes2d(0.0<=t<=time,-0.02<=vz<=0.1));
    fig.set_bounds(t,0.0,time);
    fig << line_style(true) << fill_colour(cyan);
    fig << orbit;
    fig << fill_colour(red) << orbit.final();
//    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-verification-orbit-vz");

    // Plotting phi
    fig.clear();
    fig.set_axes(Axes2d(0.0<=t<=time,0.0<=phi<=pi/2));
    fig.set_bounds(t,0.0,time);
    fig << line_style(true) << fill_colour(cyan);
    fig << orbit;
    fig << fill_colour(red) << orbit.final();
//    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-verification-orbit-phi");

    // Plotting vphi
    fig.clear();
    fig.set_axes(Axes2d(0.0<=t<=time,-0.7<=vphi<=0.1));
    fig.set_bounds(t,0.0,time);
    fig << line_style(true) << fill_colour(cyan);
    fig << orbit;
    fig << fill_colour(red) << orbit.final();
//    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-verification-orbit-vphi");


    std::cout << " done." << endl;

    return 0;

    // Set up the evaluators
//     GeneralHybridEvolver evolver;
//     evolver.verbosity = 0;
//     evolver.parameters().maximum_enclosure_radius = 10.0;
//     evolver.parameters().maximum_step_size = 2.5;
//     evolver.parameters().enable_reconditioning = true;
//     evolver.parameters().maximum_spacial_error = 1e-3;
// 
//     std::cout << "Evolution parameters:" << evolver.parameters() << std::endl;
// 
//     // Define the initial box
//     Box initial_box(numvar);
//     initial_box[0] = Interval(0.0,0.0);
//     initial_box[1] = Interval(0.0,0.0);
//     initial_box[2] = Interval(10.0,10.0);
//     initial_box[3] = Interval(0.0,0.0);
//     initial_box[4] = Interval(pi<Float>()/2.0,pi<Float>()/2.0);
//     initial_box[5] = Interval(0.0,0.0);
//     initial_box[6] = Interval(0.0,0.0);
//  
//     cout << "initial_box=" << initial_box << endl;
// 
// 
//     HybridEnclosureType initial_set(AtomicDiscreteLocation("free"), initial_box);
//     cout << "initial_set=" << initial_set << endl << endl;
// 
//     Semantics semantics=UPPER_SEMANTICS;
// 
//     Interval zsafe(11.85,11.95);
//     Interval zrange(0.0,12.0);
//     int i=0;
//     int imax=2;
//     
//     while(!subset(zrange,zsafe) && i <= imax) {
//         std::cout << "Computing evolution for delta = " << delta << " .... " <<  std::flush;
//         // Compute the reachable sets
//         Orbit<HybridEnclosureType> orbit = evolver.orbit(robotarm,initial_set,HybridTime(time,steps),semantics);
//         cout << "done." << std::endl;
//     
//         // cout << "\norbit=\n" << orbit << endl << endl;
//     
//         // Show the value of the bounding box of the final set for x and z
//         Box bbox = orbit.final()[contact].bounding_box();
//         zrange = bbox[2];
//         std::cout << "Bounding box of the final set:" << std::endl;
//         std::cout << "      x = " << bbox[0] << std::endl;
//         std::cout << "      z = " << zrange << std::endl << std::endl;
//         
//         std::cout << "Safe set for z: " << zsafe << std::endl;
//         
//         if(subset(zrange,zsafe)) {
//             std::cout << "      system is safe, exiting the loop." << std::endl;
//             std::cout << "      the current value of delta is " << delta << std::endl;
//         } else {
//             std::cout << "      system is possibly unsafe, decreasing delta" << std::endl;
//             delta = delta - 0.05;
//             rdelta = RealExpression::constant(numvar,Real(delta));
//             robotarm.set_transition(free, f2c, contact, reset_f2c, RealExpression(x - rc + rdelta), permissive);
//             robotarm.set_invariant(free, ifree, RealExpression(x - rc - rdelta)); 
//             i++;
//         }        
//     }
/*    

    // Print the intial, evolve and reach sets
    // cout << "Plotting sets" << endl;
    // cout << "evolve_set=" << hybrid_evolve_set << endl;
    // cout << "reach_set=" << hybrid_reach_set << endl;
    std::cout << "Plotting..." << std::flush;
    Figure fig;
    // Plotting x
    fig.set_projection_map(PlanarProjectionMap(numvar,6,0));
    fig.set_bounding_box(Box(numvar, 0.0,10.0, 0.0,1.0, 10.0,12.0, 0.0,1.0, 0.0,1.6, 0.0,1.0, 0.0,time, 0.0,0.0, 0.0,0.0));
    fig << line_style(false) << fill_colour(yellow);
    fig << Box(numvar, c-delta,c+delta, 0.0,1.0, 10.0,12.0, 0.0,1.0, 0.0,1.6, 0.0,1.0, -0.1,time+0.1, 0.0,0.0, 0.0,0.0);
    fig << line_style(true) << fill_colour(cyan);
    fig << orbit;
    fig.write("robotarm-verification-orbit-x");

    // Plotting z
    fig.clear();
    fig.set_projection_map(PlanarProjectionMap(numvar,6,2));
    fig << line_style(true) << fill_colour(cyan) << orbit;
    fig.write("robotarm-verification-orbit-z");    

    // Plotting phi
    fig.clear();
    fig.set_projection_map(PlanarProjectionMap(numvar,6,4));
    fig << line_style(true) << fill_colour(cyan) << orbit;
    fig.write("robotarm-verification-orbit-phi");    


    return 0;
*/

}
