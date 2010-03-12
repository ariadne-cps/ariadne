/***************************************************************************
 *            watertank.cc
 *
 *  Copyright  2008  Davide Bresolin
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

// Declare the type to be used for the system evolution
typedef HybridEvolver::EnclosureType HybridEnclosureType;


template<class X> Matrix<X> inverse2x2(const Matrix<X>& m) {
    ARIADNE_ASSERT_MSG(m.row_size() == 2 && m.column_size() == 2,
        "ERROR in inverse2x2: input matrix should be of size 2x2.");
        
    Matrix<X> invm(2,2);
    // Compute the determinant
    X det = m[0][0]*m[1][1] - m[0][1]*m[1][0];    
    // Compute transposed matrix of cofactors
    invm[0][0] = m[1][1];      invm[0][1] = -m[0][1];
    invm[1][0] = -m[1][0];     invm[1][1] = m[0][0];
    // Inverse = 1/Det(m) * tcof(m);
    invm = invm / det;
    
    return invm;    
}

void simplify(Matrix<RealExpression>& m) {
    for(uint i = 0; i < m.row_size(); i++) {
        for(uint j = 0; j < m.column_size(); j++) {
            // std::cout << "Simplyfing m[" << i << "][" << j << "]" << std::endl;
            m[i][j].simplify();
            // std::cout << std::endl << "RESULT = " << m[i][j] << std::endl << std::endl;
        }
    }
}

ListSet<HybridEnclosureType> apply(const VectorFunction& f, const ListSet<HybridEnclosureType>& ls) {
    ListSet<HybridEnclosureType> ret;
    for(ListSet<HybridEnclosureType>::const_iterator iter=ls.begin(); iter != ls.end() ; iter++) {
        ret.adjoin(iter->first, apply(f, TaylorSet((iter->second).bounding_box())));
    }
    
    return ret;
}

Float max_radius(const ListSet<HybridEnclosureType>& ls) {
    Float rad=0.0;
    for(ListSet<HybridEnclosureType>::const_iterator iter=ls.begin(); iter != ls.end() ; iter++) {
        if(iter->second.radius() > rad)
            rad = iter->second.radius();
    }
    
    return rad;
}

ListSet<Box> bounding_boxes(const ListSet<TaylorSet>& ls) {
    ListSet<Box> ret;
    for(ListSet<TaylorSet>::const_iterator iter=ls.begin(); iter != ls.end() ; iter++) {
        ret.push_back(iter->bounding_box());
    }
    
    return ret;
}

int main(int argc, char** argv) 
{
    // Definition of the arm system
    // Constants
    Float r1 = 0.25;        Float r2  = 0.15;       // m
    Float J1 = 0.0980;      Float J2 = 0.0115;      // kg m^2
    Float m1 = 1.90;        Float m2 = 0.93;        // kg
    Float Jm = 3.3e-6;      // kg m^2
    Float Jg1 = 0.0002;     Float Jg2 = 0.0005;     // kg m^2
    Float Jm1 = Jm + Jg1;   Float Jm2 = Jm + Jg2;   // kg m^2
    Float N1 = 90.0;        Float N2 = 220.0;       // coupling ratio
    
    // State variables
    RealVariable x1("x1");  // theta1: shoulder angle
    RealVariable v1("v1");  // dot theta1
    RealVariable x2("x2");  // theta2: elbow angle
    RealVariable v2("v2");  // dot theta2
    // Input variables
    RealVariable u1("u1");  // Input torque for the shoulder
    RealVariable u2("u2");  // Input torque for the elbow
    // Time
    RealVariable t("t");

    /* 
     * Dynamics is expressed as a function of ddot theta1 and ddot theta2 (angular accelerations)
     *
     *  / u1 \   / m11 m12 \ / a1 \   / g1 \
     *  |    | = |         | |    | + |    |
     *  \ u2 /   \ m21 m22 / \ a2 /   \ g2 /
     */
    
    // Matrix of the coefficients for ddot theta1 and ddot theta2
    Matrix<RealExpression> m(2,2);   
    m[0][0] = (Jm1 + (J1 + J2)/(N1*N1) + (m1*r1*r1 + m2*r2*r2 + 4.0*m2*r1*r1)/(4.0*N1*N1)) +
              (m2*r1*r2)/(N1*N1)*cos(x2/N2);
    m[0][1] = J2/(N1*N2) + (m2*r2*r2)/(4.0*N1*N2) + (m2*r1*r2)/(2.0*N1*N2) * cos(x2/N2);              
    m[1][0] = m[0][1];
    m[1][1] = Jm2 + J2/(N2*N2) + (m2*r2*r2)/(4.0*N2*N2);
    std::cout << "m11 = " << m[0][0] << std::endl;
    std::cout << "m12 = " << m[0][1] << std::endl;
    std::cout << "m21 = " << m[1][0] << std::endl;
    std::cout << "m22 = " << m[1][1] << std::endl;
    
    // "Constant" terms in theta1, theta2, dot theta1, dot theta2
    RealExpression g1 = (m2*r1*r2)/(N1*N2) * ((v1*v2)/N1 + (v2*v2)/(2.0*N2)) * sin(-x2/N2);
    RealExpression g2 = (m2*r1*r2)/(2.0*N1*N2) * (v1*v1)/(N1) * sin(x2/N2);
    std::cout << "g1 = " << g1 << std::endl;
    std::cout << "g2 = " << g2 << std::endl;
    
    /*
     * Invert M to obtain dynamics for dot dot theta1, dot dot theta2 in function of 
     * u1, u2, theta1, theta2, dot theta1, dot theta2.
     *
     * dot dot theta = M^-1 (u - g)
     */
    Matrix<RealExpression> invm = inverse2x2(m);
    simplify(invm);
    std::cout << "invm11 = " << invm[0][0] << std::endl;
    std::cout << "invm12 = " << invm[0][1] << std::endl;
    std::cout << "invm21 = " << invm[1][0] << std::endl;
    std::cout << "invm22 = " << invm[1][1] << std::endl;

    /*
     * Define differential equations for dot dot theta1, dot dot theta2, dot theta1 and dot theta2.
     */
    RealExpression dv1 = invm[0][0] * (u1 - g1) + invm[0][1] * (u2 - g2);
    RealExpression dv2 = invm[1][0] * (u1 - g1) + invm[1][1] * (u2 - g2);
    RealExpression dx1 = v1;
    RealExpression dx2 = v2;
    std::cout << "dot x1 = " << dx1 << std::endl;
    std::cout << "dot x2 = " << dx2 << std::endl;
    dv1.simplify();
    dv2.simplify();
    std::cout << "dot v1 = " << dv1 << std::endl;
    std::cout << "dot v2 = " << dv2 << std::endl;
    /*
     * Inputs u1 and u2 are constant
     */
    RealExpression du1 = 0.0;
    RealExpression du2 = 0.0;
    std::cout << "dot u1 = " << du1 << std::endl;
    std::cout << "dot u2 = " << du2 << std::endl;
    // t is a clock
    RealExpression dt= 1.0;
    std::cout << "dot t = " << dt << std::endl << std::endl;
        
    /*
     * Set up the VectorFunction representing the dynamics;
     */
    List<RealExpression> exprlist;
    exprlist.append(dx1);       exprlist.append(dx2);
    exprlist.append(dv1);       exprlist.append(dv2);
    exprlist.append(du1);       exprlist.append(du2);
    exprlist.append(dt);
    std::cout << "Expression list = " << exprlist << std::endl;
    List<RealVariable> varlist;
    varlist.append(x1);         varlist.append(x2);
    varlist.append(v1);         varlist.append(v2);
    varlist.append(u1);         varlist.append(u2);
    varlist.append(t);
    std::cout << "Variable list = " << varlist << std::endl << std::endl;
  
    VectorFunction dyn(exprlist,varlist);
    std::cout << "Dynamics = " << std::endl << dyn << std::endl << std::endl;

    // Transitions that keeps theta1 and theta2 in [0, k*2pi]
    int k = 4;
    RealExpression x1_leq_zero = -x1;               // x1 <= 0.0
    RealExpression x1_geq_k2pi = x1 - k*2.0*pi<Float>();      // x1 >= 2pi
    RealExpression x2_leq_zero = -x2;               // x1 <= 0.0
    RealExpression x2_geq_k2pi = x2 - k*2.0*pi<Float>();      // x1 >= 2pi
    std::cout << "Guards:" << std::endl;
    std::cout << "  x1 <= 0.0 : " << x1_leq_zero << std::endl;
    std::cout << "  x1 >= 2pi : " << x1_geq_k2pi << std::endl;
    std::cout << "  x2 <= 0.0 : " << x2_leq_zero << std::endl;
    std::cout << "  x2 >= 2pi : " << x2_geq_k2pi << std::endl;
        
    // Reset Functions
    RealExpression zero = 0.0;
    RealExpression ktwopi = k*2*pi<Float>();
    RealExpression idx1 = x1;
    RealExpression idx2 = x2;   
    RealExpression idv1 = v1;
    RealExpression idv2 = v2;
    RealExpression idu1 = u1;
    RealExpression idu2 = u2;
    RealExpression idt = t;
    
    // VectorFunctions for the resets
    exprlist[0] = zero;     exprlist[1] = idx2;
    exprlist[2] = idv1;     exprlist[3] = idv2;
    exprlist[4] = idu1;     exprlist[5] = idu2;
    exprlist[6] = idt;
    VectorFunction reset_x1_zero(exprlist, varlist);
    exprlist[0] = ktwopi;     exprlist[1] = idx2;
    VectorFunction reset_x1_k2pi(exprlist, varlist);
    exprlist[0] = idx1;      exprlist[1] = zero;
    VectorFunction reset_x2_zero(exprlist, varlist);
    exprlist[0] = idx1;      exprlist[1] = ktwopi;
    VectorFunction reset_x2_k2pi(exprlist, varlist);
    
    std::cout << "Resets:" << std::endl;
    std::cout << "    x1 := 0 : " << reset_x1_zero << std::endl;
    std::cout << "  x1 := 2pi : " << reset_x1_k2pi << std::endl;
    std::cout << "    x2 := 0 : " << reset_x2_zero << std::endl;
    std::cout << "  x2 := 2pi : " << reset_x2_k2pi << std::endl;
    
    //
    // Definition of the Hybrid Automaton for the Robot Arm
    //
    HybridAutomaton robotarm("RobotArm");
    
    // Two modes: initial speed-up and free run
    DiscreteState speedup("speedup");
    robotarm.new_mode(speedup, dyn);
    DiscreteState free("free");
    robotarm.new_mode(free, dyn);
    
    // Transitions
    DiscreteEvent x1l0("x1l0");
    DiscreteEvent x1g2pi("x1g2pi");
    DiscreteEvent x2l0("x2l0");
    DiscreteEvent x2g2pi("x2g2pi");
/*
    robotarm.new_transition(x1l0, free, free, reset_x1_k2pi, ScalarFunction(x1_leq_zero, varlist), true);
    robotarm.new_transition(x1g2pi, free, free, reset_x1_zero, ScalarFunction(x1_geq_k2pi, varlist), true);
    robotarm.new_transition(x2l0, free, free, reset_x2_k2pi, ScalarFunction(x2_leq_zero, varlist), true);
    robotarm.new_transition(x2g2pi, free, free, reset_x2_zero, ScalarFunction(x2_geq_k2pi, varlist), true);
*/
    // Transitions that puts u1 and u2 to zero after t0 seconds
    exprlist[0] = idx1;         exprlist[1] = idx2;
    exprlist[4] = zero;         exprlist[5] = zero;
    VectorFunction set_u1u2_zero(exprlist, varlist);
    Float t0 = 5.0;
    RealExpression  t_geq_t0 = t - t0;
    DiscreteEvent stop("stop");
    robotarm.new_transition(stop, speedup, free, set_u1u2_zero, ScalarFunction(t_geq_t0, varlist), true);

    std::cout << "RobotArm = " << std::endl;
    std::cout << robotarm << std::endl << std::endl;
        
    /*
     * COMPUTE THE EVOLUTION OF THE SYSTEM
     */
      
    // Set up the evolution parameters and grid
    Float time(10.0);
    if(argc >= 4)   // read initial value for time from the arguments
        time = atof(argv[3]);
    int steps = 2;
    if(argc >= 5)   // read initial value for steps from the arguments
        steps = atoi(argv[4]);
    Float step_size(1.25e-2);
    if(argc >= 6)
        step_size = atof(argv[5]);  // read step size from the arguments        
    Vector<Float> enclosure_cell(7,0.5);

    EvolutionParameters parameters;
    parameters.maximum_enclosure_cell=enclosure_cell;
    parameters.maximum_step_size=step_size;
    parameters.enable_set_model_reduction=true;

    std::cout << "Evolution parameters:" << parameters << std::endl;

    // Set up the evaluators
    HybridEvolver evolver(parameters);
    evolver.verbosity = 1;
    
    // Define the initial box
    Float iu1 = 5e-5;
    if(argc >= 2)   // read initial value for u1 from the arguments
        iu1 = atof(argv[1]);
    Float iu2 = 5e-5;
    if(argc >= 3)   // read initial value for u2 from the arguments
        iu2 = atof(argv[2]);
    
    Box initial_box(7, 0.5*pi<Float>(),0.5*pi<Float>(), 0.0,0.0, 0.0,0.0, 0.0,0.0, iu1,iu1, iu2,iu2, 0.0,0.0);
    cout << "initial_box=" << initial_box << endl;

    // Over-approximate the initial set by a grid cell
    HybridEnclosureType initial_set(speedup, initial_box);
    cout << "initial_set=" << initial_set << endl << endl;

    Semantics semantics=UPPER_SEMANTICS;

    // Compute the reachable sets
    Orbit<HybridEnclosureType> orbit = evolver.orbit(robotarm,initial_set,HybridTime(time,steps),semantics);
    cout << std::endl;
    
    Box bbox = orbit.reach().bounding_box();
    
    // cout << "\norbit=\n" << orbit << endl << endl;

    // Print the intial, evolve and reach sets
    // cout << "Plotting sets" << endl;
    // cout << "evolve_set=" << hybrid_evolve_set << endl;
    // cout << "reach_set=" << hybrid_reach_set << endl;
    std::cout << "Plotting..." << std::flush;
    Figure fig;
    // Plotting x1 and x2
    fig.set_projection_map(PlanarProjectionMap(7,6,0));
    fig.set_bounding_box(Box(2, 0.0,time, bbox[0].lower(),bbox[0].upper()));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-orbit-x1");
    fig.clear();
    fig.set_projection_map(PlanarProjectionMap(7,6,1));
    fig.set_bounding_box(Box(2, 0.0,time, bbox[1].lower(),bbox[1].upper()));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-orbit-x2");

    // Plotting v1 and v2
    fig.clear();
    fig.set_projection_map(PlanarProjectionMap(7,6,2));
    fig.set_bounding_box(Box(2, 0.0,time, bbox[2].lower(),bbox[2].upper()));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-orbit-v1");
    fig.clear();
    fig.set_projection_map(PlanarProjectionMap(7,6,3));
    fig.set_bounding_box(Box(2, 0.0,time, bbox[3].lower(),bbox[3].upper()));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-orbit-v2");

    // Plotting u1 and u2
    fig.clear();
    fig.set_projection_map(PlanarProjectionMap(7,6,4));
    fig.set_bounding_box(Box(2, 0.0,time, bbox[4].lower(),bbox[4].upper()));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-orbit-u1");
    fig.clear();
    fig.set_projection_map(PlanarProjectionMap(7,6,5));
    fig.set_bounding_box(Box(2, 0.0,time, bbox[5].lower(),bbox[5].upper()));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-orbit-u2");

    /*
     *  Convert polar coordinates theta1, theta2 into cartesian coordinates
     *  for plotting the trajectory
     */
    
    // Definition of the conversion function
    RealExpression px1 = r1 * cos(x1) + r2 * cos(x1 + x2);
    RealExpression px2 = r1 * sin(x1) + r2 * sin(x1 + x2);
    RealExpression pt = t;
    
    List<RealExpression> projlist;
    projlist.append(px1);       projlist.append(px2);
    projlist.append(pt);
    VectorFunction proj(projlist, varlist);
    
    HybridEnclosureType cartesian_initial(orbit.initial().location(), apply(proj, orbit.initial().continuous_state_set()));
    ListSet<HybridEnclosureType> cartesian_final = apply(proj, orbit.final());
    ListSet<HybridEnclosureType> cartesian_reach = apply(proj, orbit.reach());

    // Plotting cartesian trajectory
    fig.clear();
    fig.set_bounding_box(Box(2, -0.5,0.5, -0.5,0.5));
    fig.set_projection_map(PlanarProjectionMap(3,0,1));
    fig << line_style(true) << fill_colour(cyan) << cartesian_reach;
    fig << fill_colour(red) << cartesian_final;
    fig << fill_colour(blue) << cartesian_initial;
    fig.write("robotarm-orbit-xy");
    fig.clear();
    fig.set_bounding_box(Box(2, -0.0,time, -0.5,0.5));
    fig.set_projection_map(PlanarProjectionMap(3,2,0));
    fig << line_style(true) << fill_colour(cyan) << cartesian_reach;
    fig << fill_colour(red) << cartesian_final;
    fig << fill_colour(blue) << cartesian_initial;
    fig.write("robotarm-orbit-xt");
    fig.clear();
    fig.set_bounding_box(Box(2, -0.0,time, -0.5,0.5));
    fig.set_projection_map(PlanarProjectionMap(3,2,1));
    fig << line_style(true) << fill_colour(cyan) << cartesian_reach;
    fig << fill_colour(red) << cartesian_final;
    fig << fill_colour(blue) << cartesian_initial;
    fig.write("robotarm-orbit-yt");


    std::cout << "done!" << std::endl;

/*    
    std::cout << "Cartesian final = " << cartesian_final << std::endl;
    std::cout << "Polar final = " << orbit.final() << std::endl;
    
    
      
    std::cout << "Projection function = " << proj << std::endl;
    std::cout << "           evaluate = " << proj.evaluate(orbit.initial().range()) << std::endl;
    
    std::cout << "Initial set = " << orbit.initial() << std::endl;
    std::cout << "     domain = " << orbit.initial().domain() << std::endl;
    std::cout << "      range = " << orbit.initial().range() << std::endl;
    std::cout << "     models = " << orbit.initial().models() << std::endl;

    TaylorSet cartesian_initial = apply(proj, orbit.initial());
    std::cout << "Projected set = " << cartesian_initial << std::endl;
    std::cout << "       domain = " << cartesian_initial.domain() << std::endl;
    std::cout << "        range = " << cartesian_initial.range() << std::endl;
    std::cout << "       models = " << cartesian_initial.models() << std::endl;
*/
    
    
    
    return 0;
     
}
