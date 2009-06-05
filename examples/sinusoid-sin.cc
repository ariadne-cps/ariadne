/***************************************************************************
 *            sinusoid_cos_timed_infinite.cc
 ****************************************************************************/

#include <cstdarg>
#include "ariadne.h"

using namespace Ariadne;

/// Function for the generation of a sinusoid
struct GenSinusoid : FunctionData<2,2,0> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
	      r[0] = 1.0;
        r[1] = Ariadne::sin(x[0]+pi<Float>()/2);
    }
};

int main() 
{    
    Vector<Interval> system_parameters(0);

    /// Build the Hybrid System
  
    /// Create a HybridAutomton object
    HybridAutomaton sinusoid;
  
    /// Create the discrete state
    DiscreteState l1(1);

    /// Create the discrete event
    DiscreteEvent e1(1);
 
    /// Create the resets
    AffineFunction reset(Matrix<Float>(2,2,0.0,0.0,0.0,1.0),Vector<Float>(2,0.0,0.0));
    cout << "reset=" << reset << endl << endl;

    /// Create the guards.
    /// Guards are true when f(x) = Ax + b > 0
    AffineFunction guard(Matrix<Float>(1,2,1.0,0.0),Vector<Float>(1,-2*pi<Float>()));
    cout << "guard=" << guard << endl << endl;
 
    /// Create the dynamics
    Function<GenSinusoid> sinusoid_dynamic(system_parameters);
    
    cout << "dynamic = " << sinusoid_dynamic << endl << endl;
  
    /// Build the automaton
    sinusoid.new_mode(l1,sinusoid_dynamic);
    sinusoid.new_forced_transition(e1,l1,l1,reset,guard);

    /// Finished building the automaton

    cout << "Automaton = " << sinusoid << endl << endl;

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver;

    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = 0.5;
    evolver.parameters().maximum_step_size = 0.05;
    std::cout <<  evolver.parameters() << std::endl;
    evolver.verbosity=1;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;

    std::cout << "Computing evolution starting from location l1, t = 0.0, x =0.0" << std::endl;

    Box initial_box(2, 0.0,0.0,0.0,0.0);
    HybridEnclosureType initial_enclosure(l1,initial_box);
  
    HybridTime evolution_time(6*pi<Float>(),4);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(sinusoid,initial_enclosure,evolution_time,LOWER_SEMANTICS);
    std::cout << "done." << std::endl;

    //std::cout << "Orbit="<<orbit<<std::endl;
    
    Figure g;
    Box graphic_box(2, -0.01,2*pi<Float>()+0.01, -1.1,1.6);
    g.set_bounding_box(graphic_box);
    array<uint> p(2,1,0);
    g.set_projection_map(ProjectionFunction(p,3));

    g << fill_colour(Colour(0.0,0.5,1.0));
    g << orbit;
    g.write("sinusoid_sin_orbit");
/*
    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.verbosity = 3;
    analyser.parameters().lock_to_grid_time = 2*pi<Float>();
    analyser.parameters().maximum_grid_depth+=2;
    std::cout <<  analyser.parameters() << std::endl;

    HybridImageSet initial_set;
    initial_set[l1]=initial_box;

    // HybridTime reach_time(3*pi<Float>(),1);

    std::cout << "Computing lower reach set... " << std::flush;
    HybridGridTreeSet reach_set;
    HybridGridTreeSet evolve_set;
    make_lpair(reach_set,evolve_set) = analyser.lower_reach_evolve(sinusoid,initial_set,evolution_time);
    std::cout << "done." << std::endl;
    
    std::cout << "reach_set.size()=" << reach_set.size() << std::endl;
    std::cout << "evolve_set.size()=" << evolve_set.size() << std::endl;
   
    g.clear();
    Box graphic_box2(2, -0.1,2*pi<Float>()+0.1, -1.1,1.6);
    g.set_bounding_box(graphic_box2);
    array<uint> p2(2,1,0);
    g.set_projection_map(ProjectionFunction(p2,3));

    g << fill_colour(Colour(0.0,0.5,1.0));
    g << reach_set;
    g << fill_colour(Colour(1.0,0.0,0.0));
    g << evolve_set;
    g.write("sinusoid_sin_lower");

    std::cout << "Computing upper reach set... " << std::flush;
    make_lpair(reach_set,evolve_set) = analyser.upper_reach_evolve(sinusoid,initial_set,evolution_time);
    std::cout << "done." << std::endl;
    
    std::cout << "reach_set.size()=" << reach_set.size() << std::endl;
    std::cout << "evolve_set.size()=" << evolve_set.size() << std::endl;
   
    g.clear();
    g.set_bounding_box(graphic_box2);
    g.set_projection_map(ProjectionFunction(p2,3));

    g << fill_colour(Colour(0.0,0.5,1.0));
    g << reach_set;
    g << fill_colour(Colour(1.0,0.0,0.0));
    g << evolve_set;
    g.write("sinusoid_sin_upper");
 
    std::cout << "Computing chainreach set... " << std::flush;
    reach_set = analyser.chain_reach(sinusoid,initial_set);
    std::cout << "done." << std::endl;
    
    std::cout << "reach_set.size()=" << reach_set.size() << std::endl;

    g.clear();
    g.set_bounding_box(graphic_box2);
    g.set_projection_map(ProjectionFunction(p2,3));

    g << fill_colour(Colour(0.0,0.5,1.0));
    g << reach_set;
    g.write("sinusoid_sin_chain");
*/ 
}
