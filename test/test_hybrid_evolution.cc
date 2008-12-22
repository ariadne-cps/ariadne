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
#include "sparse_differential.h"
#include "differential_vector.h"
#include "function.h"
#include "approximate_taylor_model.h"
#include "box.h"
#include "zonotope.h"
#include "list_set.h"
#include "evolution_parameters.h"
#include "orbit.h"
#include "hybrid_set.h"
#include "hybrid_evolver.h"
#include "graphics.h"
#include "logging.h"

#include "models.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

class TestHybridEvolution
{
  private:
    static HybridAutomaton system();
  public:
    void test() const;
};

int main() 
{
    TestHybridEvolution().test();
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
    automaton.new_forced_transition(event3,location1,location2,reset,guard3);
    automaton.new_forced_transition(event4,location2,location1,reset,guard4);
    //automaton.new_unforced_transition(event4,location2,location1,reset,activation4);

    cout << "Finished creating hybrid automaton." << endl;

    return automaton;
}

void TestHybridEvolution::test() const
{
    cout << __PRETTY_FUNCTION__ << endl;

    const DiscreteState location1(1);
    const DiscreteState location2(2);
    const DiscreteEvent event3(3);
    const DiscreteEvent event4(4);

    typedef ApproximateTaylorModel EnclosureType;
    typedef std::pair<DiscreteState,ApproximateTaylorModel> HybridEnclosureType;

    // Set up the evolution parameters and grid
    Float time(3.0);
    Float step_size(0.125);
    Float enclosure_radius(0.25);
    
    EvolutionParameters parameters;
    parameters.maximum_enclosure_radius=enclosure_radius;
    parameters.maximum_step_size=step_size;

    // Set up the evaluators
    HybridEvolver evolver(parameters);
  
    // Define the initial box
    Box initial_box(2); 
    initial_box[0]=Interval(-0.01,0.01);
    initial_box[1]=Interval(0.49,0.51);
    //initial_box[0]=Interval(-0.02,0.02);
    //initial_box[1]=Interval(0.48,0.52);

    // cout << "initial_box=" << initial_box << endl;


    // Make a hybrid automaton for the Van der Pol equation

    HybridAutomaton automaton=system();
    ARIADNE_TEST_PRINT(automaton);

    // Over-approximate the initial set by a grid cell
    Vector<Interval> unit_box(2,Interval(-1,1));
    EnclosureType initial_set=ApproximateTaylorModel(unit_box,ScalingFunction(initial_box),4,1);
    // cout << "initial_set=" << initial_set << endl << endl;
    HybridEnclosureType initial_hybrid_set(location1,initial_set);
    HybridTime hybrid_time(time,5);

  
    // Compute the reachable sets
    Orbit<HybridEnclosureType> orbit=evolver.orbit(automaton,initial_hybrid_set,hybrid_time);
    ListSet<HybridEnclosureType> hybrid_evolve_set,hybrid_intermediate_set,hybrid_reach_set;
    hybrid_evolve_set = orbit.final();
    hybrid_intermediate_set = orbit.intermediate();
    hybrid_reach_set = orbit.reach();

    ARIADNE_TEST_PRINT(hybrid_evolve_set);
    ARIADNE_TEST_PRINT(hybrid_reach_set);
    ARIADNE_TEST_PRINT(hybrid_intermediate_set);
  
    cout << "Plotting sets... " << flush;
    Figure fig;
    fig << line_style(true);
    fig << fill_colour(cyan) << hybrid_reach_set;
    fig << fill_colour(magenta) << hybrid_intermediate_set;
    fig << fill_colour(blue) << hybrid_evolve_set;
    fig << fill_colour(blue) << initial_set;
    fig.write("test_hybrid_evolution-affine");
    cout << "done" << endl;

}
