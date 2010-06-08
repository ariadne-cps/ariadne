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
    
    Float pmin = 0.5;                                 // Sampling period
    Float pmax = 1.25;                                 // Sampling period
    Float lper = 30.0;                              // Switching period
    
    // State variables
    RealVariable x1("x1");  // theta1: shoulder angle
    RealVariable v1("v1");  // dot theta1
    RealVariable x2("x2");  // theta2: elbow angle
    RealVariable v2("v2");  // dot theta2
    // Input variables
    RealVariable u1("u1");      // Input torque for the shoulder
    RealVariable u2("u2");      // Input torque for the elbow
    RealVariable uc1("u1");    // Input torque from the controller
    RealVariable uc2("u2");    // Input torque from the controller
  
    // Time and clocks
    RealVariable t("t");
    RealVariable l("l");
    RealVariable c("c");        // Sampling clock for the network

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
    
    //
    // Invert M to obtain dynamics for dot dot theta1, dot dot theta2 in function of 
    // u1, u2, theta1, theta2, dot theta1, dot theta2.
    //
    // dot dot theta = M^-1 (u - g)
    //
    Matrix<RealExpression> invm = inverse2x2(m);
    simplify(invm);
    std::cout << "invm11 = " << invm[0][0] << std::endl;
    std::cout << "invm12 = " << invm[0][1] << std::endl;
    std::cout << "invm21 = " << invm[1][0] << std::endl;
    std::cout << "invm22 = " << invm[1][1] << std::endl;

    //
    // Define differential equations for dot dot theta1, dot dot theta2, dot theta1 and dot theta2.
    //
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
    //
    // Inputs uc1 and uc2 are obtained from a given trajectory
    //
    // Algebraic expression for the trajectory
    RealExpression alg_x1 = pi<Float>() * cos(2* pi<Float>() / 120.0 * t);             // theta1
//    RealExpression alg_x1 = 2*pi<Float>()/sqr(60.0) * sqr(t);             // theta1
    RealExpression alg_x2 = pi<Float>() * cos(2* pi<Float>() / 30.0 * t);             // theta1
    RealExpression alg_v1 = derivative(alg_x1, t);     // dot theta1      
    RealExpression alg_v2 = derivative(alg_x2, t);     // dot theta2
    RealExpression alg_a1 = derivative(alg_v1, t);                    // dot dot theta1 
    RealExpression alg_a2 = derivative(alg_v2, t);                    // dot dot theta2
    std::cout << "Algebraic expression for a1 = " << alg_a1 << std::endl;
    std::cout << "Algebraic expression for a2 = " << alg_a2 << std::endl;
    std::cout << "Algebraic expression for v1 = " << alg_v1 << std::endl;
    std::cout << "Algebraic expression for v2 = " << alg_v2 << std::endl;
    std::cout << "Algebraic expression for x1 = " << alg_x1 << std::endl;
    std::cout << "Algebraic expression for x2 = " << alg_x2 << std::endl << std::endl;
    RealExpression alg_u1= simplify(m[0][0]*alg_a1 + m[1][0]*alg_a2 + g1);
    std::cout << "Algebraic expression for uc1 = " << alg_u1 << std::endl;
    alg_u1 = alg_u1.substitute(v1, alg_v1);
    std::cout << "                             = " << alg_u1 << std::endl;
    alg_u1 = alg_u1.substitute(v2, alg_v2);
    std::cout << "                             = " << alg_u1 << std::endl;
    alg_u1 = alg_u1.substitute(x1, alg_x1);
    std::cout << "                             = " << alg_u1 << std::endl;
    alg_u1 = simplify(alg_u1.substitute(x2, alg_x2));
    std::cout << "                             = " << alg_u1 << std::endl;
    RealExpression alg_u2= simplify(m[1][0]*alg_a1 + m[1][1]*alg_a2 + g2);
    alg_u2 = alg_u2.substitute(v1, alg_v1);
    alg_u2 = alg_u2.substitute(v2, alg_v2);
    alg_u2 = alg_u2.substitute(x1, alg_x1);
    alg_u2 = simplify(alg_u2.substitute(x2, alg_x2));
    std::cout << "Algebraic expression for uc2 = " << alg_u2 << std::endl;
       
    RealExpression duc1 = simplify(derivative(alg_u1,t));
    RealExpression duc2 = simplify(derivative(alg_u2,t));

    std::cout << "dot u1 = " << duc1 << std::endl;
    std::cout << "dot u2 = " << duc2 << std::endl;

//    return 0;

    // t is a clock
    RealExpression dt= 1.0;
    std::cout << "dot t = " << dt << std::endl << std::endl;

//    return 0;
 
    //
    // Definition of the Hybrid Automaton for the Robot Arm
    //
    HybridIOAutomaton robotarm("RobotArm");
    robotarm.add_input_var(u1);
    robotarm.add_input_var(u2);
    robotarm.add_output_var(x1);
    robotarm.add_output_var(v1);
    robotarm.add_output_var(x2);
    robotarm.add_output_var(v2);   
    // events to keep x1 and x2 between -pi and pi
    DiscreteEvent x1l("x1l");    DiscreteEvent x1u("x1u");
    DiscreteEvent x2l("x2l");    DiscreteEvent x2u("x2u");
    robotarm.add_internal_event(x1l);   robotarm.add_internal_event(x1u);
    robotarm.add_internal_event(x2l);   robotarm.add_internal_event(x2u);
    
    // One mode: free run
    DiscreteState free("free");
    robotarm.new_mode(free);
    robotarm.set_dynamics(free, x1, dx1);
    robotarm.set_dynamics(free, v1, dv1);
    robotarm.set_dynamics(free, x2, dx2);
    robotarm.set_dynamics(free, v2, dv2);
    // transitions that resets x1
    robotarm.new_forced_transition(x1l, free, free, indicator((x1 <= -pi<Float>()) && (v1 <= 0.0), positive));
    robotarm.set_reset(x1l, free, x1, pi<Float>());
    robotarm.set_reset(x1l, free, x2, x2);
    robotarm.set_reset(x1l, free, v1, v1);
    robotarm.set_reset(x1l, free, v2, v2);
    robotarm.new_forced_transition(x1u, free, free, indicator((x1 >= pi<Float>()) && (v1 >= 0.0), positive));
    robotarm.set_reset(x1u, free, x1, -pi<Float>());
    robotarm.set_reset(x1u, free, x2, x2);
    robotarm.set_reset(x1u, free, v1, v1);
    robotarm.set_reset(x1u, free, v2, v2);
    // transitions that resets x2
    robotarm.new_forced_transition(x2l, free, free, indicator((x2 <= -pi<Float>()) && (v2 <= 0.0), positive));
    robotarm.set_reset(x2l, free, x2, pi<Float>());
    robotarm.set_reset(x2l, free, x1, x1);
    robotarm.set_reset(x2l, free, v1, v1);
    robotarm.set_reset(x2l, free, v2, v2);
    robotarm.new_forced_transition(x2u, free, free, indicator((x2 >= pi<Float>()) && (v2 >= 0.0), positive));
    robotarm.set_reset(x2u, free, x2, -pi<Float>());
    robotarm.set_reset(x2u, free, x1, x1);
    robotarm.set_reset(x2u, free, v1, v1);
    robotarm.set_reset(x2u, free, v2, v2);
    
    std::cout << "RobotArm = " << std::endl;
    std::cout << robotarm << std::endl << std::endl;

    //
    // Definition of the Hybrid Automaton for the controller
    //
    HybridIOAutomaton controller("Controller");
    controller.add_output_var(uc1);
    controller.add_output_var(uc2);
    controller.add_internal_var(t);     // global time
    controller.add_internal_var(l);     // local clock
    DiscreteEvent lreset("lreset");
    controller.add_internal_event(lreset);     // local clock
    
    // Two modes: uc1 and uc2 increase
    DiscreteState incr("incr");
    controller.new_mode(incr);
    controller.set_dynamics(incr, uc1, duc1);
    controller.set_dynamics(incr, uc2, duc2);
    controller.set_dynamics(incr, t, dt);
    controller.set_dynamics(incr, l, dt);

//     DiscreteState decr("decr");
//     controller.new_mode(decr);
//     controller.set_dynamics(decr, uc1, -5e-7);
//     controller.set_dynamics(decr, uc2, -2.5e-6);
//     controller.set_dynamics(decr, t, dt);
//     controller.set_dynamics(decr, l, dt);

    // Two transitions that resets l and switch mode
    RealExpression l_geq_lper = l - lper;     // l >= lper
    RealExpression l_res_zero = 0.0;    // l := 0
    RealExpression uc1id = uc1;     // uc1 := uc1
    RealExpression uc2id = uc2;     // uc2 := uc2
    RealExpression tid = t;         // t := t
//     controller.new_forced_transition(lreset, incr, decr, l_geq_lper);
//     controller.set_reset(lreset, incr, l, l_res_zero);
//     controller.set_reset(lreset, incr, uc1, uc1id);
//     controller.set_reset(lreset, incr, uc2, uc2id);
//     controller.set_reset(lreset, incr, t, tid);
//     controller.new_forced_transition(lreset, decr, incr, l_geq_lper);
//     controller.set_reset(lreset, decr, l, l_res_zero);
//     controller.set_reset(lreset, decr, uc1, uc1id);
//     controller.set_reset(lreset, decr, uc2, uc2id);
//     controller.set_reset(lreset, decr, t, tid);
        
    std::cout << "Controller = " << std::endl;
    std::cout << controller << std::endl << std::endl;
    
    //
    // Hybrid Automaton for the network
    //
//     HybridIOAutomaton network("Network");
//     network.add_input_var(uc1);
//     network.add_input_var(uc2);
//     network.add_output_var(u1);
//     network.add_output_var(u2);
//     network.add_internal_var(c);
//     DiscreteEvent sample_e("sample_e");
//     network.add_internal_event(sample_e);
//     
//     // One mode: sample uc1 and uc2
//     DiscreteState sample("sample");
//     network.new_mode(sample);
//     network.set_dynamics(sample, u1, 0.0);
//     network.set_dynamics(sample, u2, 0.0);
//     network.set_dynamics(sample, c, dt);
// 
//     // One transition that resets c and update the value of u1 and u2
//     RealExpression c_res_zero = 0.0;    // c := 0
//     RealExpression u1_eq_uc1 = uc1;     // u1 := uc1
//     RealExpression u2_eq_uc2 = uc2;     // u2 := uc2
//     network.new_forced_transition(sample_e, sample, sample, indicator(c >= pmin, positive));
//     // network.new_invariant(sample, indicator(c <= pmax, negative));
//     network.set_reset(sample_e, sample, c, c_res_zero);
//     network.set_reset(sample_e, sample, u1, u1_eq_uc1);
//     network.set_reset(sample_e, sample, u2, u2_eq_uc2);
// 
//     std::cout << "Network = " << std::endl;
//     std::cout << network << std::endl << std::endl;
//         
    //
    // COMPOSE THE SYSTEM
    //
//    HybridIOAutomaton contrnet = compose("contrnet", controller, network, incr, sample);
//    HybridIOAutomaton ncs = compose("NCS", contrnet, robotarm, DiscreteState("incr,sample"), free);
    HybridIOAutomaton ncs = compose("NCS", controller, robotarm, incr, free);
    RealSpace spc;
    HybridAutomaton system;
    make_lpair<HybridAutomaton, RealSpace>(system, spc) = make_monolithic_automaton(ncs);

    std::cout << "Complete System = " << std::endl;
    std::cout << system << std::endl << std::endl;

        
    //
    // COMPUTE THE EVOLUTION OF THE SYSTEM
    //
      
    // Set up the evolution parameters and grid
    Float time(60.0);
    if(argc >= 2)   // read initial value for time from the arguments
        time = atof(argv[1]);
    int steps = 10;
    if(argc >= 3)   // read initial value for steps from the arguments
        steps = atoi(argv[2]);
    Float step_size(1.25e-2);
    if(argc >= 4)
        step_size = atof(argv[3]);  // read step size from the arguments        
    Vector<Float> enclosure_cell(spc.dimension(),1.0);

    EvolutionParameters parameters;
    parameters.maximum_enclosure_cell=enclosure_cell;
    parameters.maximum_step_size=step_size;
    parameters.enable_set_model_reduction=true;
    parameters.enable_premature_termination=true;

    std::cout << "Evolution parameters:" << parameters << std::endl;

    // Set up the evaluators
    HybridEvolver evolver(parameters);
    evolver.verbosity = 1;
    
    // Define the initial box
    //
    // The default initial values for u1 and u2 are given by alg_u1(0.0) and alg_u2(0.0)
    //
    ContinuousValuation<Real> val;
    val.set(t, 0.0);
    std::cout << "Initial valuation: ";    
    Real iv1 = 0.0;
    std::cout << ", v1 = " << iv1;
    Real iv2 = 0.0;
    std::cout << ", v2 = " << iv2;
    Real ix1 = evaluate(alg_x1, val);;
    std::cout << ", x1 = " << ix1;
    Real ix2 = evaluate(alg_x2, val);;
    std::cout << ", x2 = " << ix2;
    Real iuc1 = evaluate(alg_u1, val);
    std::cout << "uc1 = " << iuc1;
    Real iuc2 = evaluate(alg_u2, val);
    std::cout << ", uc2 = " << iuc2 << endl; 
    Real iu1 = iuc1;
    std::cout << "u1 = " << iu1;
    Real iu2 = iuc2;
    std::cout << ", u2 = " << iu2;
    
    Box initial_box(spc.dimension());
    initial_box[spc.index(x1)] = Interval(ix1.lower(),ix1.upper());
    initial_box[spc.index(v1)] = Interval(iv1.lower(),iv1.upper());
    initial_box[spc.index(x2)] = Interval(ix2.lower(),ix2.upper());
    initial_box[spc.index(v2)] = Interval(iv2.lower(),iv2.upper());
    initial_box[spc.index(u1)] = Interval(iu1.lower(),iu1.upper());
    initial_box[spc.index(u2)] = Interval(iu2.lower(),iu2.upper());
    initial_box[spc.index(uc1)] = Interval(iuc1.lower(),iuc1.upper());
    initial_box[spc.index(uc2)] = Interval(iuc2.lower(),iuc2.upper());
    initial_box[spc.index(t)] = Interval(0.0,0.0);
//    initial_box[spc.index(c)] = Interval(0.0,0.0);
    initial_box[spc.index(l)] = Interval(0.5*lper,0.5*lper);
    
    cout << "initial_box=" << initial_box << endl;


    // Over-approximate the initial set by a grid cell
//    HybridEnclosureType initial_set(DiscreteState("incr,sample,free"), initial_box);
    HybridEnclosureType initial_set(DiscreteState("incr,free"), initial_box);
    cout << "initial_set=" << initial_set << endl << endl;

    Semantics semantics=UPPER_SEMANTICS;

    // Compute the reachable sets
    Orbit<HybridEnclosureType> orbit = evolver.orbit(system,initial_set,HybridTime(time,steps),semantics);
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
    fig.set_projection_map(PlanarProjectionMap(spc.dimension(),spc.index(t),spc.index(x1)));
    fig.set_bounding_box(Box(2, 0.0,time, -pi<Float>(), pi<Float>()));
    // Plot the reference trajectory for x1
    fig << line_style(true) << fill_colour(yellow);
    Float tstep = time/100.0;
    Box bx(spc.dimension());
    for(Float tp = 0.0; tp <= time; tp += tstep) {
        Interval ti(tp-0.1*tstep, tp+0.1*tstep);
        val.set(t, Real(ti));
        Interval x1i = evaluate(alg_x1, val);
        bx[spc.index(x1)] = x1i;
        bx[spc.index(t)] = ti;
        fig << bx;
    }
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-orbit-x1");
    fig.clear();
    fig.set_projection_map(PlanarProjectionMap(spc.dimension(),spc.index(t),spc.index(x2)));
    fig.set_bounding_box(Box(2, 0.0,time, -pi<Float>(), pi<Float>()));
    // Plot the reference trajectory for x2
    fig << line_style(true) << fill_colour(yellow);
    for(Float tp = 0.0; tp <= time; tp += tstep) {
        Interval ti(tp-0.1*tstep, tp+0.1*tstep);
        val.set(t, Real(ti));
        Interval x2i = evaluate(alg_x2, val);
        bx[spc.index(x2)] = x2i;
        bx[spc.index(t)] = ti;
        fig << bx;
    }
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-orbit-x2");

    // Plotting v1 and v2
    fig.clear();
    fig.set_projection_map(PlanarProjectionMap(spc.dimension(),spc.index(t),spc.index(v1)));
    fig.set_bounding_box(Box(2, 0.0,time, bbox[spc.index(v1)].lower(),bbox[spc.index(v1)].upper()));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-orbit-v1");
    fig.clear();
    fig.set_projection_map(PlanarProjectionMap(spc.dimension(),spc.index(t),spc.index(v2)));
    fig.set_bounding_box(Box(2, 0.0,time, bbox[spc.index(v2)].lower(),bbox[spc.index(v2)].upper()));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-orbit-v2");

    // Plotting u1 and u2
    fig.clear();
    fig.set_projection_map(PlanarProjectionMap(spc.dimension(),spc.index(t),spc.index(u1)));
    fig.set_bounding_box(Box(2, 0.0,time, bbox[spc.index(u1)].lower(),bbox[spc.index(u1)].upper()));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-orbit-u1");
    fig.clear();
    fig.set_projection_map(PlanarProjectionMap(spc.dimension(),spc.index(t),spc.index(u2)));
    fig.set_bounding_box(Box(2, 0.0,time, bbox[spc.index(u2)].lower(),bbox[spc.index(u2)].upper()));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-orbit-u2");

    // Plotting uc1 and uc2
    fig.clear();
    fig.set_projection_map(PlanarProjectionMap(spc.dimension(),spc.index(t),spc.index(uc1)));
    fig.set_bounding_box(Box(2, 0.0,time, bbox[spc.index(uc1)].lower(),bbox[spc.index(uc1)].upper()));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-orbit-uc1");
    fig.clear();
    fig.set_projection_map(PlanarProjectionMap(spc.dimension(),spc.index(t),spc.index(uc2)));
    fig.set_bounding_box(Box(2, 0.0,time, bbox[spc.index(uc2)].lower(),bbox[spc.index(uc2)].upper()));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-orbit-uc2");

    // Plotting c
//     fig.clear();
//     fig.set_projection_map(PlanarProjectionMap(spc.dimension(),spc.index(t),spc.index(c)));
//     fig.set_bounding_box(Box(2, 0.0,time, bbox[spc.index(c)].lower(),bbox[spc.index(c)].upper()));
//     fig << line_style(true) << fill_colour(cyan) << orbit.reach();
//     fig << fill_colour(magenta) << orbit.intermediate();
//     fig << fill_colour(red) << orbit.final();
//     fig << fill_colour(blue) << initial_set;
//     fig.write("robotarm-orbit-c");

    // Plotting l
    fig.clear();
    fig.set_projection_map(PlanarProjectionMap(spc.dimension(),spc.index(t),spc.index(l)));
    fig.set_bounding_box(Box(2, 0.0,time, bbox[spc.index(l)].lower(),bbox[spc.index(l)].upper()));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("robotarm-orbit-l");


    //
    //  Convert polar coordinates theta1, theta2 into cartesian coordinates
    //  for plotting the trajectory
    //
    
    // Definition of the conversion function
    RealExpression px1 = r1 * cos(x1) + r2 * cos(x1 + x2);
    RealExpression px2 = r1 * sin(x1) + r2 * sin(x1 + x2);
    RealExpression pt = t;
    
    List<RealExpression> projlist;
    projlist.append(px1);       projlist.append(px2);
    projlist.append(pt);
    VectorFunction proj(projlist, spc);
    
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

    
//     std::cout << "Cartesian final = " << cartesian_final << std::endl;
//     std::cout << "Polar final = " << orbit.final() << std::endl;
//       
//     std::cout << "Projection function = " << proj << std::endl;
//     std::cout << "           evaluate = " << proj.evaluate(orbit.initial().range()) << std::endl;
//     
//     std::cout << "Initial set = " << orbit.initial() << std::endl;
//     std::cout << "     domain = " << orbit.initial().domain() << std::endl;
//     std::cout << "      range = " << orbit.initial().range() << std::endl;
//     std::cout << "     models = " << orbit.initial().models() << std::endl;
// 
//     TaylorSet cartesian_initial = apply(proj, orbit.initial());
//     std::cout << "Projected set = " << cartesian_initial << std::endl;
//     std::cout << "       domain = " << cartesian_initial.domain() << std::endl;
//     std::cout << "        range = " << cartesian_initial.range() << std::endl;
//     std::cout << "       models = " << cartesian_initial.models() << std::endl;

    return 0;

}
