/***************************************************************************
 *            test_hybrid_simulator.cc
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
#include "evolution_parameters.h"
#include "orbit.h"
#include "hybrid_time.h"
#include "hybrid_set.h"
#include "hybrid_automaton.h"
#include "hybrid_simulator.h"
#include "graphics.h"
#include "logging.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

int verbosity=0;

class TestHybridSimulator
{
  private:
    static MonolithicHybridAutomaton system();
  public:
    void test() const;
};

int main()
{
    TestHybridSimulator().test();
    std::cerr<<"INOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

MonolithicHybridAutomaton
TestHybridSimulator::system()
{
    const AtomicDiscreteLocation location1(1);
    const AtomicDiscreteLocation location2(2);
    const DiscreteEvent event3(3);
    const DiscreteEvent event4(4);

    MonolithicHybridAutomaton automaton("Affine Hysteresis System");
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
    automaton.new_forced_transition(event3,location1,location2,reset,guard3);
    automaton.new_forced_transition(event4,location2,location1,reset,guard4);
    //automaton.new_unforced_transition(event4,location2,location1,reset,activation4);

    cout << "Finished creating hybrid automaton." << endl;

    return automaton;
}

void TestHybridSimulator::test() const
{
    cout << __PRETTY_FUNCTION__ << endl;

    const AtomicDiscreteLocation location1(1);
    const AtomicDiscreteLocation location2(2);
    const DiscreteEvent event3(3);
    const DiscreteEvent event4(4);

    // Set up the simulator parameters and grid
    Float step_size(0.125);
    Float enclosure_radius(0.25);

    EvolutionParameters parameters;
    parameters.maximum_step_size=step_size;

    // Set up the evaluators
    HybridSimulator simulator(parameters);
    simulator.verbosity = verbosity;


    // Make a hybrid automaton for the Van der Pol equation
    MonolithicHybridAutomaton automaton=system();
    ARIADNE_TEST_PRINT(automaton);

    // Define the initial box
    Point initial_point(2, -0.00, 0.50);
    cout << "initial_point=" << initial_point << endl;
    HybridPoint initial_hybrid_point(location1,initial_point);
    HybridTime simulation_time(0.25,1);


    // Compute the reachable sets
    cout << "Computing orbit... "<<std::flush;
    Orbit<HybridPoint> hybrid_orbit=simulator.orbit(automaton,initial_hybrid_point,simulation_time);
    cout << "done"<<std::endl;

    ARIADNE_TEST_PRINT(hybrid_orbit);


    cout << "Plotting orbit... " << flush;
    Figure fig;
    fig << hybrid_orbit;
    fig.write("test_hybrid_simulator-orbit");
    cout << "done" << endl;

}
