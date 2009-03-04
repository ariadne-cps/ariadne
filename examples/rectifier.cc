/***************************************************************************
 *            rectifier.cc
 ****************************************************************************/

#include <cstdarg>
#include "ariadne.h"

using namespace Ariadne;

/// Control variables:
/// t: time
/// vi: input (sinusoidal) voltage
/// vo: output (rectified) voltage

/// Non-affine dynamics

/// Dynamics for the case of both diodes being off
/// t'=1, vi'= A*cos(2*pi*f*t), vo'=-vo/(Rl*Cl)
struct offoff_df : FunctionData<3,3,5> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
	r[0] = 1.0;
    r[1] = p[0]*2*pi<Float>()*p[1]*Ariadne::cos(2*pi<Float>()*p[1]*x[0]);
	r[2] = -x[2]/(p[4]*p[3]);
    }
};

/// Dynamics for the case of the first diode being on, the second being off
/// t'=1, vi'= A*cos(2*pi*f*t), vo'=-vo/(Rl*Cl)+(vi-vo)/(Ron*Cl)
struct onoff_df : FunctionData<3,3,5> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
	r[0] = 1.0;
    r[1] = p[0]*2*pi<Float>()*p[1]*Ariadne::cos(2*pi<Float>()*p[1]*x[0]);
	r[2] = -x[2]/(p[4]*p[3]) + (x[1]-x[2])/(p[2]*p[3]);
    }
};

/// Dynamics for the case of the first diode being off, the second being on
/// t'=1, vi'= A*cos(2*pi*f*t), vo'=-vo/(Rl*Cl)-(vi+vo)/(Ron*Cl)
struct offon_df : FunctionData<3,3,5> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
	r[0] = 1.0;
    r[1] = p[0]*2*pi<Float>()*p[1]*Ariadne::cos(2*pi<Float>()*p[1]*x[0]);
	r[2] = -x[2]/(p[4]*p[3]) - (x[1]+x[2])/(p[2]*p[3]);
    }
};

/// Dynamics for the case of both diodes being on
/// t'=1, vi'= A*cos(2*pi*f*t), vo'=-vo/(Rl*Cl)-2*vo/(Ron*Cl)
struct onon_df : FunctionData<3,3,5> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
	r[0] = 1.0;
    r[1] = p[0]*2*pi<Float>()*p[1]*Ariadne::cos(2*pi<Float>()*p[1]*x[0]);
	r[2] = -x[2]/(p[4]*p[3]) -2*x[2]/(p[2]*p[3]);
    }
};

int main() 
{    
    /// Introduces the parameters
    Vector<Float> dynamics_parameters(5);
    dynamics_parameters[0] = 4.0; /// Amplitude of the input voltage, Vi
    dynamics_parameters[1] = 50.0; /// Sinusoid frequency, f
    dynamics_parameters[2] = 0.1; /// Diode resistance when on, Ron
    dynamics_parameters[3] = 0.0001; /// Load capacitance, Cl
    dynamics_parameters[4] = 1000.0; /// Load resistance, Rl
//    dynamics_parameters[3] = 0.05; /// Load capacitance, Cl
//    dynamics_parameters[4] = 1.0; /// Load resistance, Rl

    /// Build the Hybrid System
  
    /// Create a HybridAutomton object
    HybridAutomaton rectifier;
  
    /// Create the discrete states
    DiscreteState offoff(1);
    DiscreteState onoff(2);
    DiscreteState offon(3);
    DiscreteState onon(4);

    /// Create the discrete events
    DiscreteEvent resettime(1);
    DiscreteEvent jump1(2), jump2(3), jump3(4);
    
    /// Create the resets

    /// Reset the time (t^=0,vi^=vi,vo^=vo)
    AffineFunction resettime_r(Matrix<Float>(3,3,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0),
                               Vector<Float>(3,0.0,0.0,0.0));
    /// Do nothing (t^=t,vi^=vi,vo^=vo)
    IdentityFunction noop_r(3); 
//    AffineFunction noop_r(Matrix<Float>(3,3, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,1.0,0.0),
//                               Vector<Float>(3, 0.0,0.0,0.0));

    /// Create the guards
    
    /// Guard for the reset of time (t>=1/f)
    AffineFunction resettime_g(Matrix<Float>(1,3,1.0,0.0,0.0),Vector<Float>(1,-1/dynamics_parameters[1]));
    /// Guard for the jump from onoff to offoff (vi-vo<=0)
    AffineFunction onoff_offoff_g(Matrix<Float>(1,3,0.0,-1.0,1.0),Vector<Float>(1,0.0));
    /// Guard for the jump from offon to offoff (-vi-vo<=0)
    AffineFunction offon_offoff_g(Matrix<Float>(1,3,0.0,1.0,1.0),Vector<Float>(1,0.0));
    /// Guard for the jump from offoff to onoff (vi-vo>=0)
    AffineFunction offoff_onoff_g(Matrix<Float>(1,3,0.0,1.0,-1.0),Vector<Float>(1,0.0));
    /// Guard for the jump from onon to onoff (-vi-vo<=0)
    AffineFunction onon_onoff_g(Matrix<Float>(1,3,0.0,1.0,1.0),Vector<Float>(1,0.0));
    /// Guard for the jump from offoff to offon (-vi-vo>=0)
    AffineFunction offoff_offon_g(Matrix<Float>(1,3,0.0,-1.0,-1.0),Vector<Float>(1,0.0));
    /// Guard for the jump from onon to offon (vi-vo<=0)
    AffineFunction onon_offon_g(Matrix<Float>(1,3,0.0,-1.0,-1.0),Vector<Float>(1,0.0));
    /// Guard for the jump from offon to onon (vi-vo>=0) 
    AffineFunction offon_onon_g(Matrix<Float>(1,3,0.0,1.0,-1.0),Vector<Float>(1,0.0));
    /// Guard for the jump from onoff to onon (-vi-vo>=0) 
    AffineFunction onoff_onon_g(Matrix<Float>(1,3,0.0,-1.0,-1.0),Vector<Float>(1,0.0));
 
    /// Create the dynamics
    Function<offoff_df> offoff_d(dynamics_parameters);
    Function<onoff_df> onoff_d(dynamics_parameters);
    Function<offon_df> offon_d(dynamics_parameters);
    Function<onon_df> onon_d(dynamics_parameters);
  
    /// Build the automaton
    
    /// Locations
    rectifier.new_mode(offoff,offoff_d);
    rectifier.new_mode(onoff,onoff_d);
    rectifier.new_mode(offon,offon_d);
    rectifier.new_mode(onon,onon_d);
    /// OffOff events
    rectifier.new_forced_transition(resettime,offoff,offoff,resettime_r,resettime_g);
    rectifier.new_forced_transition(jump1,offoff,onoff,noop_r,offoff_onoff_g);
    rectifier.new_forced_transition(jump2,offoff,offon,noop_r,offoff_offon_g);
//    rectifier.new_forced_transition(jump3,offoff,onon,noop_r,to_onon_g);
    /// OnOff events
    rectifier.new_forced_transition(resettime,onoff,onoff,resettime_r,resettime_g);
    rectifier.new_forced_transition(jump1,onoff,offoff,noop_r,onoff_offoff_g);
//    rectifier.new_forced_transition(jump2,onoff,offon,noop_r,to_offon_g);
//    rectifier.new_forced_transition(jump3,onoff,onon,noop_r,to_onon_g);
    /// OffOn events
    rectifier.new_forced_transition(resettime,offon,offon,resettime_r,resettime_g);
    rectifier.new_forced_transition(jump1,offon,offoff,noop_r,offon_offoff_g);
/*    rectifier.new_forced_transition(jump2,offon,onoff,noop_r,to_onoff_g);
    rectifier.new_forced_transition(jump3,offon,onon,noop_r,to_onon_g);
    /// OnOn events
    rectifier.new_forced_transition(resettime,onon,onon,resettime_r,resettime_g);
    rectifier.new_forced_transition(jump1,onon,offoff,noop_r,to_offoff_g);
    rectifier.new_forced_transition(jump2,onon,onoff,noop_r,to_onoff_g);
    rectifier.new_forced_transition(jump3,onon,offon,noop_r,to_offon_g);
*/

    /// Finished building the automaton

    cout << "Automaton = " << rectifier << endl << endl;

    /// Compute the system evolution

    global_verbosity = 1;

    /// Create a HybridEvolver object
    HybridEvolver evolver;

    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = 0.001/dynamics_parameters[1];
    evolver.parameters().maximum_step_size = 0.001/dynamics_parameters[1];
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;

    std::cout << "Computing evolution starting from location offoff, t = 0.0, vi = 0.0, vo = 1.0" << std::endl;

    Box initial_box(3, 0.0,0.0, 0.0,0.0, dynamics_parameters[0],dynamics_parameters[0]);
    HybridEnclosureType initial_enclosure(offoff,initial_box);

//    Box initial_box(3, 0.01,0.01, 3.80423,3.80423, 3.36842,3.36842);
//    HybridEnclosureType initial_enclosure(onoff,initial_box);

//    Box initial_box(3, 0.04,0.04, -3.80422,-3.80422, 3.34775,3.34775);
//    HybridEnclosureType initial_enclosure(offon,initial_box);

    
    std::cout << "Initial set=" << initial_enclosure << std::endl;
  
    HybridTime evolution_time(2.0/dynamics_parameters[1],8);
//    HybridTime evolution_time(0.01678,3); 

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(rectifier,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit.initial="<<orbit.initial()<<std::endl;
    std::cout << "Orbit.final="<<orbit.final()<<std::endl;

///    Box graphic_box(2,-dynamics_parameters[0]-0.1,dynamics_parameters[0]+0.1,-dynamics_parameters[0]-0.1,dynamics_parameters[0]+0.1);
    Box graphic_box(2,-0.1/dynamics_parameters[1],1.0/dynamics_parameters[1]*1.1,-dynamics_parameters[0]*1.1,dynamics_parameters[0]*1.1);
 //   Box graphic_box2(2,-dynamics_parameters[0]*1.1,dynamics_parameters[0]*1.1,-dynamics_parameters[0]*1.1,dynamics_parameters[0]*1.1);
    Box graphic_box2(2,-dynamics_parameters[0],dynamics_parameters[0],-3.5,dynamics_parameters[0]);

    Figure g;
    array<uint> tvin(2,0,1);
    array<uint> tvout(2,0,2);
    array<uint> vinvout(2,1,2);


    g.set_bounding_box(graphic_box);
    g.set_projection_map(ProjectionFunction(tvin,3));

    g << fill_colour(Colour(0.0,0.5,1.0));
    g << orbit;
    g.write("rectifier_orbit_t_vin");
    
    g.clear();
    g.set_bounding_box(graphic_box);
    g.set_projection_map(ProjectionFunction(tvout,3));

    g << fill_colour(Colour(0.0,0.5,1.0));
    g << orbit;
    g.write("rectifier_orbit_t_vout");
    
    g.clear();
    g.set_bounding_box(graphic_box2);
    g.set_projection_map(ProjectionFunction(vinvout,3));
/*    
    g << fill_colour(Colour(1.0,1.0,1.0));
    for(int i = -40; i < 40; i++) {
        g << Box(3, 0.0,0.0, -4.0,4.0, 0.1*i,0.1*(i+1));
        g << Box(3, 0.0,0.0, 0.1*i,0.1*(i+1), -4.0,4.0);
    }
*/
    g << fill_colour(Colour(0.0,0.5,1.0));
    g << orbit;
    g.write("rectifier_orbit_vin_vout");
/*
    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.parameters().lock_to_grid_time = 1.0/(dynamics_parameters[1]);
    analyser.parameters().maximum_grid_depth= 10 + 12 * round(std::log(dynamics_parameters[1]));
    std::cout <<  analyser.parameters() << std::endl;

    HybridImageSet initial_set;
    initial_set[offoff]=initial_box;
    HybridTime reach_time(2.0/dynamics_parameters[1],8);
//    HybridTime reach_time(0.001,1);

    std::cout << "Computing upper reach set... " << std::flush;
    HybridGridTreeSet* upper_reach_set_ptr = analyser.upper_reach(rectifier,initial_set,reach_time);
    //HybridGridTreeSet* upper_reach_set_ptr = analyser.chain_reach(rectifier,initial_set);
    std::cout << "done." << std::endl;

    g.clear();
    g.set_bounding_box(graphic_box);
    g.set_projection_map(ProjectionFunction(tvin,3));

    g << fill_colour(Colour(0.0,0.5,1.0));
    g << *upper_reach_set_ptr;
    g.write("rectifier_reach_t_vin");
    
    g.clear();
    g.set_bounding_box(graphic_box);
    g.set_projection_map(ProjectionFunction(tvout,3));

    g << fill_colour(Colour(0.0,0.5,1.0));
    g << *upper_reach_set_ptr;
    g.write("rectifier_reach_t_vout");
    
    g.clear();
    g.set_bounding_box(graphic_box2);
    g.set_projection_map(ProjectionFunction(vinvout,3));
    
    g << fill_colour(Colour(1.0,1.0,1.0));
    for(int i = -40; i < 40; i++) {
        g << Box(3, 0.0,0.0, -4.0,4.0, 0.1*i,0.1*(i+1));
        g << Box(3, 0.0,0.0, 0.1*i,0.1*(i+1), -4.0,4.0);
    }
    g << fill_colour(Colour(0.0,0.5,1.0));
    g << *upper_reach_set_ptr;
    g.write("rectifier_reach_vin_vout");
*/

}
