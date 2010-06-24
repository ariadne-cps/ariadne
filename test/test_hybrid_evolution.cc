/***************************************************************************
 *            test_hybrid_evolution.cc
 *
 *  Copyright  2006-8  Pieter Collins
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

#include <fstream>
#include <iostream>

#include "tuple.h"
#include "vector.h"
#include "matrix.h"
#include "function.h"
#include "taylor_set.h"
#include "taylor_function.h"
#include "box.h"
#include "zonotope.h"
#include "list_set.h"
#include "evolution_parameters.h"
#include "orbit.h"
#include "hybrid_time.h"
#include "hybrid_set.h"
#include "hybrid_evolver.h"
#include "graphics_interface.h"
#include "graphics.h"
#include "logging.h"

#include "models.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

int evolver_verbosity=0;

class TestHybridEvolution
{
    typedef Vector<Float> FVector;
    typedef Matrix<Float> FMatrix;

    static const bool non_urgent=false;
    static const bool urgent=true;
  private:
    static HybridAutomaton system();
  public:
    void test() const;
    void test_constant_derivative_system() const;
    void test_bouncing_ball() const;
    void test_affine_system() const;
};

HybridAutomaton
TestHybridEvolution::system()
{
    const DiscreteState location1(1);
    const DiscreteState location2(2);
    const DiscreteEvent event3(3);
    const DiscreteEvent event4(4);

    HybridAutomaton automaton("Affine Hysteresis System");
    double adata[]={-0.5,-1.0,1.0,-0.5};
    double bdata[]={1.0,0.0};
    Matrix<Float> A(2,2,adata);
    Vector<Float> b(2,bdata);
    VectorAffineFunction dynamic1(A,3*b);
    VectorAffineFunction dynamic2(A,-b);
    IdentityFunction reset(2);

    Matrix<Float> c(1,2,bdata);
    Vector<Float> d(1,Float(1.0));
    VectorAffineFunction guard3(c,-d);
    VectorAffineFunction guard4(-c,-d);
    VectorAffineFunction activation4(-c,-d);
    VectorAffineFunction invariant2(-c,-1.125*d);

    automaton.new_mode(location1,dynamic1);
    automaton.new_mode(location2,dynamic2);
    //automaton.new_invariant(location2,invariant2);
    automaton.new_transition(event3,location1,location2,reset,guard3,urgent);
    automaton.new_transition(event4,location2,location1,reset,guard4,urgent);
    //automaton.new_unforced_transition(event4,location2,location1,reset,activation4);

    cout << "Finished creating hybrid automaton." << endl;

    return automaton;
}

void TestHybridEvolution::test_constant_derivative_system() const
{
    // Test the system (d(x),d(y))=(1,0) with reset (x',y')=(x-2,y) when x+y>0
    // Starting in a small box near the origin, the system should return to
    // the initial condition after time 2
    DiscreteState q1(1); DiscreteState q2(2); DiscreteEvent e(1);
    VectorAffineFunction d(FMatrix(2,2, 0.,0.,0.,0.),FVector(2, 1.0,0.));
    VectorAffineFunction r(FMatrix(2,2, 1.,0.,0.,1.),FVector(2, -2.,0.));
    VectorAffineFunction g(FMatrix(1,2, 1.,0.,0.,0.),FVector(1, -1.25));

    HybridAutomaton automaton("Constant Derivative System");
    automaton.new_mode(q1,d);
    automaton.new_mode(q2,d);
    automaton.new_transition(e,q1,q2,r,g,urgent);

    TaylorSet initial_enclosure(Box(2, -0.0625,0.0625, -0.0625,+0.0625));
    HybridTaylorSet initial_set(q1,initial_enclosure);

    HybridEvolver evolver;
    evolver.verbosity=evolver_verbosity;
	evolver.parameters().maximum_enclosure_cell=Vector<Float>(2,0.5);

    ARIADNE_TEST_PRINT(automaton);
    ARIADNE_TEST_PRINT(initial_set);

    {
        // Test continuous evolution without any jumps
        HybridTime evolution_time(0.5,1);
        ARIADNE_TEST_PRINT(evolution_time);
        Orbit<HybridTaylorSet> orbit=evolver.orbit(automaton,initial_set,evolution_time);
        ARIADNE_TEST_PRINT(orbit);
        ListSet<HybridTaylorSet> final_set=evolver.evolve(automaton,initial_set,evolution_time);
        ARIADNE_TEST_PRINT(final_set);
        HybridTaylorSet expected_final_set(q1,Box(2, +0.4375,+0.5625, -0.0625,+0.0625));
        ARIADNE_TEST_PRINT(expected_final_set);
        ARIADNE_TEST_COMPARE(norm(final_set[q1][0].models()-expected_final_set.second.models()),<,1e-15);
    }

    {
        // Test continuous evolution with a single transverse jump
        HybridTime evolution_time(2.0,2);
        ARIADNE_TEST_PRINT(evolution_time);

        Orbit<HybridTaylorSet> orbit=evolver.orbit(automaton,initial_set,evolution_time);
        ARIADNE_TEST_PRINT(orbit);

        ListSet<HybridTaylorSet> final_set=evolver.evolve(automaton,initial_set,evolution_time);
        ARIADNE_TEST_PRINT(final_set);
        HybridTaylorSet expected_final_set(q2,Box(2, -0.0625,+0.0625, -0.0625,+0.0625));
        ARIADNE_TEST_PRINT(expected_final_set);

        ARIADNE_TEST_COMPARE(norm(final_set[q2][0].models()-expected_final_set.second.models()),<,1e-14);
    }

}

void TestHybridEvolution::test_bouncing_ball() const
{
    /// Set the system parameters
    double a = 0.5; // Coefficient of restitution
    double g = 9.8; // Constant of gravity
    double x0 = 5.0; // Initial height
    double r0 = 1.0/16; // Initial box radius

    /// Create the system functions
    DiscreteState q1(1);
    DiscreteState q2(2);
    DiscreteEvent e12(12);
    //VectorAffineFunction dynamic(FMatrix(3,3, 0.,1.,0., 0.,0.,0., 0.,0.,0.), FVector(3, 0.0, -g, 1.0));
    //VectorAffineFunction reset(FMatrix(3,3, 1.0,0.0 ,0.0,-a,0.0, 0.0,0.0,1.0), FVector(3, 0.0,0.0,0.0));
    //VectorAffineFunction guard(FMatrix(1,3, -1.0,0.0,0.0), FVector(1, 0.0));
    VectorAffineFunction dynamic(FMatrix(2,2, 0.,1., 0.,0.), FVector(2, 0.0, -g));
    VectorAffineFunction reset(FMatrix(2,2, 1.0,0.0 ,0.0,-a), FVector(2, 0.0,0.0));
    VectorAffineFunction guard(FMatrix(1,2, -1.0,0.0), FVector(1, 0.0));

    /// Build the automaton
    HybridAutomaton automaton;
    automaton.new_mode(q1,dynamic);
    automaton.new_mode(q2,dynamic);
    automaton.new_transition(e12,q1,q2,reset,guard,urgent);

    //TaylorSet initial_enclosure(Box(3, x0-r0,x0+r0, -r0,+r0, 0.0,0.0));
    TaylorSet initial_enclosure(Box(2, x0-r0,x0+r0, -r0,+r0));
    HybridTaylorSet initial_set(q1,initial_enclosure);

    HybridEvolver evolver;
    evolver.verbosity=evolver_verbosity;
    evolver.parameters().hybrid_maximum_step_size[1]=0.125;
    evolver.parameters().hybrid_maximum_step_size[2]=0.125;
	evolver.parameters().maximum_enclosure_cell=Vector<Float>(2,0.5);

    ARIADNE_TEST_PRINT(automaton);
    ARIADNE_TEST_PRINT(initial_set);

    // Test continuous evolution without any jumps
    HybridTime evolution_time(1.5,2);
    ARIADNE_TEST_PRINT(evolution_time);
    Orbit<HybridTaylorSet> orbit=evolver.orbit(automaton,initial_set,evolution_time);
    ARIADNE_TEST_PRINT(orbit);
    ARIADNE_TEST_EVALUATE(orbit.intermediate()[q2]);
    ARIADNE_TEST_EVALUATE(orbit.reach()[q2]);
    ARIADNE_TEST_EVALUATE(orbit.final()[q2]);

    plot("test_hybrid_evolution-bouncing_ball-orbit",Box(2,-1.,6.,-12.,8.),
         Colour(.99,.99,.99),Box(2,-2.,0.,-20.,+20.),Colour(0,0.5,1),orbit.reach()[q1],Colour(0,1,1),orbit.reach()[q2]);
}



void TestHybridEvolution::test_affine_system() const
{
    cout << __PRETTY_FUNCTION__ << endl;

    const DiscreteState location1(1);
    const DiscreteState location2(2);
    const DiscreteEvent event3(3);
    const DiscreteEvent event4(4);

    typedef TaylorSet EnclosureType;
    typedef HybridBasicSet<TaylorSet> HybridEnclosureType;

    // Set up the evolution parameters and grid
    Float step_size(0.5);
    Float max_enclosure_width(0.5);

    EvolutionParameters parameters;
    parameters.maximum_enclosure_cell=Vector<Float>(2,max_enclosure_width);
    parameters.hybrid_maximum_step_size[1]=step_size;
    parameters.hybrid_maximum_step_size[2]=step_size;
    parameters.hybrid_maximum_step_size[3]=step_size;
    parameters.hybrid_maximum_step_size[4]=step_size;

    // Set up the evaluators
    HybridEvolver evolver(parameters);
    evolver.verbosity = evolver_verbosity;
	evolver.parameters().maximum_enclosure_cell=Vector<Float>(2,0.5);

    // Make a hybrid automaton for the Van der Pol equation
    HybridAutomaton automaton=system();
    ARIADNE_TEST_PRINT(automaton);

    // Define the initial box
    Box initial_box(2, -0.01,0.01, 0.49,0.51);
    cout << "initial_box=" << initial_box << endl;
    EnclosureType initial_set=TaylorSet(initial_box);
    cout << "initial_set=" << initial_set << endl << endl;
    HybridEnclosureType initial_hybrid_set(location1,initial_set);
    HybridTime hybrid_evolution_time(0.25,1);


    // Compute the reachable sets
    cout << "Computing orbit... "<<std::flush;
    Orbit<HybridEnclosureType> orbit=evolver.orbit(automaton,initial_hybrid_set,hybrid_evolution_time);
    cout << "done"<<std::endl;
    ListSet<HybridEnclosureType> hybrid_evolve_set,hybrid_intermediate_set,hybrid_reach_set;
    hybrid_evolve_set = orbit.final();
    hybrid_intermediate_set = orbit.intermediate();
    hybrid_reach_set = orbit.reach();

    ARIADNE_TEST_PRINT(hybrid_evolve_set);
    ARIADNE_TEST_PRINT(hybrid_reach_set);
    ARIADNE_TEST_PRINT(hybrid_intermediate_set);

    cout << "Plotting sets... " << flush;
    Figure fig;
    fig.set_bounding_box(Box(2, -0.25, 0.75, 0.0, 1.0));
    fig << line_style(true);
    fig << hybrid_reach_set;
    fig << fill_colour(cyan) << hybrid_reach_set;
    fig << fill_colour(magenta) << hybrid_intermediate_set;
    fig << fill_colour(blue) << hybrid_evolve_set;
    fig << fill_colour(red) << initial_set;
    fig.write("test_hybrid_evolution-affine");
    cout << "done" << endl;

    cout << "Plotting orbit... " << flush;
    fig.clear();
    fig << orbit;
    fig.write("test_hybrid_evolution-orbit");
    cout << "done" << endl;

}


void TestHybridEvolution::test() const
{
    //ARIADNE_TEST_CALL(test_constant_derivative_system());
    ARIADNE_TEST_CALL(test_bouncing_ball());
    //ARIADNE_TEST_CALL(test_affine_system());
}





// Test the HybridEvolver class on two-dimensional examples with simple flows and resets
class TestHybridEvolver
{
  private:
    DiscreteState q1,q2;
    DiscreteEvent e;
    HybridEvolver evolver;
    ScalarFunction z,o,x,y;
    ScalarFunction x0,y0,t;
  public:
    TestHybridEvolver();
    HybridAutomaton make_hybrid_automaton(const ScalarFunction& guard);

    void test();
    void test_transverse_linear_crossing();
    void test_transverse_cubic_crossing();
    void test_transverse_cube_root_crossing();

};

TestHybridEvolver::TestHybridEvolver()
    : evolver()
{
    // Set up convenience variables
    q1=DiscreteState(1);
    q2=DiscreteState(2);
    e=DiscreteEvent(3);

    z=ScalarFunction::constant(2,0.0);
    o=ScalarFunction::constant(2,1.0);
    x=ScalarFunction::variable(2,0);
    y=ScalarFunction::variable(2,1);
    x0=ScalarFunction::variable(3,0);
    y0=ScalarFunction::variable(3,1);
    t=ScalarFunction::variable(3,2);
}

HybridAutomaton TestHybridEvolver::make_hybrid_automaton(const ScalarFunction& guard)
{
    HybridAutomaton system;
    system.new_mode(q1,VectorFunction(join(o,z)));
    system.new_mode(q2,VectorFunction(join(z,o)));
    system.new_transition(e,q1,q2,IdentityFunction(2),VectorFunction(1u,guard),true);
    return system;
}

void TestHybridEvolver::test_transverse_linear_crossing()
{
    Float r=1.0/8;
    Float tol=1e-5;
    ScalarFunction guard=x+y/2-1;
    HybridAutomaton system=make_hybrid_automaton(guard);
    Box initial_box(2, -r,+r, -r,+r);
    HybridTaylorSet initial_set(q1,initial_box);
    HybridTime evolution_time(2.0,3);

	evolver.parameters().maximum_enclosure_cell=Vector<Float>(2,0.5);
	HybridSpace hspace = system.state_space();
	for (HybridSpace::locations_const_iterator loc_it = hspace.locations_begin(); loc_it != hspace.locations_end(); loc_it++)
		evolver.parameters().hybrid_maximum_step_size[loc_it->first] = 1.0;

    ListSet<HybridTaylorSet> evolved_set=evolver.evolve(system,initial_set,evolution_time,UPPER_SEMANTICS);

    ScalarFunction ct=-guard; // Crossing time
    VectorFunction f=join(x+ct,y+2-ct);
    // VectorFunction f=join(x+ct,y+ct);
    Vector<Interval> tolerance(2,Interval(-tol,+tol));
    TaylorSet expected_evolved_set(f,initial_box);
    ARIADNE_TEST_ASSERT(!evolved_set.empty());
    if(!evolved_set.empty()){
        ARIADNE_TEST_BINARY_PREDICATE(refines,expected_evolved_set.models(),evolved_set[q2][0].models());
        ARIADNE_TEST_BINARY_PREDICATE(refines,evolved_set[q2][0].models(),expected_evolved_set.models()+tolerance);
        plot("test_hybrid_evolution-transverse_linear_crossing",Box(2, -1.0,3.0, -1.0,3.0),
             Colour(0,0,1),evolved_set[q2][0]);
    }
}

void TestHybridEvolver::test_transverse_cubic_crossing()
{
    Float r=1.0/8;
    Float tol=1e-5;
    ScalarFunction guard=x-(1+y/2+y*y*y);
    HybridAutomaton system=make_hybrid_automaton(guard);
    Box initial_box(2, -r,+r, -r,+r);
    HybridTaylorSet initial_set(q1,initial_box);
    HybridTime evolution_time(2.0,3);

	evolver.parameters().maximum_enclosure_cell=Vector<Float>(2,0.5);

    ListSet<HybridTaylorSet> evolved_set=evolver.evolve(system,initial_set,evolution_time,UPPER_SEMANTICS);

    ScalarFunction ct=-guard; // Crossing time

    VectorFunction f=join(x+ct,y+2-ct);
    Vector<Interval> tolerance(2,Interval(-tol,+tol));
    TaylorSet expected_evolved_set(f,initial_box);
    ARIADNE_TEST_ASSERT(!evolved_set.empty());
    if(!evolved_set.empty()){
        ARIADNE_TEST_BINARY_PREDICATE(refines,expected_evolved_set.models(),evolved_set[q2][0].models());
        ARIADNE_TEST_BINARY_PREDICATE(refines,evolved_set[q2][0].models(),expected_evolved_set.models()+tolerance);
        plot("test_hybrid_evolution-transverse_cubic_crossing",Box(2, -1.0,3.0, -1.0,3.0),
             Colour(0,0,1),evolved_set[q2][0]);
    }
}

void TestHybridEvolver::test_transverse_cube_root_crossing()
{
    Float r=1.0/32;
    Float tol=1e-5;
    ScalarFunction guard=((x-1)*(x-1)+1.0)*(x-1)-y-1./64;
    HybridAutomaton system=make_hybrid_automaton(guard);
    Box initial_box(2, -r,+r, -r,+r);
    HybridTaylorSet initial_set(q1,initial_box);
    HybridTime evolution_time(2.0,3);

	evolver.parameters().maximum_enclosure_cell=Vector<Float>(2,0.5);

    ScalarFunction ct=y-pow(y,3)+3*pow(y,5)-12*pow(y,7)+55*pow(y,9)-273*pow(y,11)+1-x;
    VectorFunction f=join(x+ct,y+2-ct);
    Vector<Interval> tolerance(2,Interval(-tol,+tol));
    TaylorSet expected_evolved_set(f,initial_box);
    
    ListSet<HybridTaylorSet> evolved_set=evolver.evolve(system,initial_set,evolution_time,UPPER_SEMANTICS);
    
    ARIADNE_TEST_ASSERT(!evolved_set.empty());
    
    if(!evolved_set.empty()) {
        ARIADNE_TEST_BINARY_PREDICATE(refines,evolved_set[q2][0].models(),expected_evolved_set.models()+tolerance);
    }
}

void TestHybridEvolver::test() {
    ARIADNE_TEST_CALL(test_transverse_linear_crossing());
    ARIADNE_TEST_CALL(test_transverse_cubic_crossing());
    ARIADNE_TEST_CALL(test_transverse_cube_root_crossing());
}


int main(int argc, const char* argv[])
{
    if(argc>=1) { evolver_verbosity=atoi(argv[0]); }

    //std::cerr<<"SKIPPED "; return 1;
    //TestHybridEvolution().test();
    TestHybridEvolver().test();
    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

