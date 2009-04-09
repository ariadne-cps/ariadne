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
#include "graphics.h"
#include "logging.h"

#include "models.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

int verbosity=0;

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
    void test_affine_system() const;
};

int main()
{
    //std::cerr<<"SKIPPED "; return 1;
    TestHybridEvolution().test();
    std::cerr<<"INOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

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
    AffineFunction dynamic1(A,3*b);
    AffineFunction dynamic2(A,-b);
    IdentityFunction reset(2);

    Matrix<Float> c(1,2,bdata);
    Vector<Float> d(1,Float(1.0));
    AffineFunction guard3(c,-d);
    AffineFunction guard4(-c,-d);
    AffineFunction activation4(-c,-d);
    AffineFunction invariant2(-c,-1.125*d);

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
    AffineFunction d(FMatrix(2,2, 0.,0.,0.,0.),FVector(2, 1.0,0.));
    AffineFunction r(FMatrix(2,2, 1.,0.,0.,1.),FVector(2, -2.,0.));
    AffineFunction g(FMatrix(1,2, 1.,0.,0.,0.),FVector(1, -1.25));

    HybridAutomaton automaton("Constant Derivative System");
    automaton.new_mode(q1,d);
    automaton.new_mode(q2,d);
    automaton.new_transition(e,q1,q2,r,g,urgent);

    TaylorSet initial_enclosure(Box(2, -0.0625,0.0625, -0.0625,+0.0625));
    HybridTaylorSet initial_set(q1,initial_enclosure);

    HybridEvolver evolver;

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


void TestHybridEvolution::test_affine_system() const
{
    cout << __PRETTY_FUNCTION__ << endl;

    const DiscreteState location1(1);
    const DiscreteState location2(2);
    const DiscreteEvent event3(3);
    const DiscreteEvent event4(4);

    typedef TaylorSet EnclosureType;
    typedef pair<DiscreteState,TaylorSet> HybridEnclosureType;

    // Set up the evolution parameters and grid
    Float step_size(0.125);
    Float enclosure_radius(0.25);

    EvolutionParameters parameters;
    parameters.maximum_enclosure_radius=enclosure_radius;
    parameters.maximum_step_size=step_size;

    // Set up the evaluators
    HybridEvolver evolver(parameters);
    evolver.verbosity = verbosity;


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
    ARIADNE_TEST_CALL(test_constant_derivative_system());
    ARIADNE_TEST_CALL(test_affine_system());
}