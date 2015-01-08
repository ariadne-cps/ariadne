/***************************************************************************
 *            test_hybrid_simulator.cc
 *
 *  Copyright  2006-11  Pieter Collins
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

#include "config.h"
#include "expression/expression.h"
#include "expression/space.h"
#include "hybrid/hybrid_set.h"
#include "hybrid/hybrid_orbit.h"
#include "hybrid/hybrid_time.h"
#include "hybrid/hybrid_automaton.h"
#include "hybrid/hybrid_simulator.h"
#include "output/graphics.h"
#include "utility/logging.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

Int verbosity=0;

class TestHybridSimulator
{
  private:
    static HybridAutomaton system();
  public:
    Void test() const;
};

Int main()
{
    TestHybridSimulator().test();
    std::cerr<<"INOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

HybridAutomaton
TestHybridSimulator::system()
{
    const DiscreteLocation location1(1);
    const DiscreteLocation location2(2);
    const DiscreteEvent event3(3);
    const DiscreteEvent event4(4);
    const RealVariable x("x");
    const RealVariable y("y");
    const RealVariable z("z");

    HybridAutomaton automaton;
    Dyadic a(-0.5); Dyadic b(1.0); Dyadic cx1(3.0); Dyadic cx2(-1.0); Dyadic cz(1.0);
    automaton.new_mode(location1,(dot(x)= a*x-b*y+cx1,dot(y)=b*x+a*y) );
    automaton.new_mode(location2,(dot(x)= a*x-b*y+cx2,dot(y)=b*x+a*y,dot(z)=cz));
    automaton.new_transition(location1,event3,location2,(next(x)=x,next(y)=y,next(z)=y),x>=1,urgent);
    automaton.new_transition(location2,event4,location1,(next(x)=x,next(y)=y),x<=-1,urgent);

    cout << "Finished creating hybrid automaton." << endl;

    return automaton;
}

Void TestHybridSimulator::test() const
{
    cout << __PRETTY_FUNCTION__ << endl;

    const DiscreteLocation location1(1);
    const DiscreteLocation location2(2);
    const DiscreteEvent event3(3);
    const DiscreteEvent event4(4);
    const RealVariable x("x");
    const RealVariable y("y");

    // Set up the simulator parameters and grid
    Float step_size(0.125);
    Float enclosure_radius(0.25);

    // Set up the evaluators
    HybridSimulator simulator;
    simulator.set_step_size(0.0625);
    simulator.verbosity = verbosity;


    // Make a hybrid automaton for the Van der Pol equation
    HybridAutomaton automaton=system();
    ARIADNE_TEST_PRINT(automaton);

    // Define the initial box
    RealSpace space((x,y));
    ExactPoint initial_point = ExactPoint{-0.00, 0.50};
    cout << "initial_point=" << initial_point << endl;
    HybridPoint initial_hybrid_point(location1,space,initial_point);
    HybridTime simulation_time(2.25,3);


    // Compute the reachable sets
    cout << "Computing orbit... "<<std::flush;
    Orbit<HybridPoint> hybrid_orbit=simulator.orbit(automaton,initial_hybrid_point,simulation_time);
    cout << "done"<<std::endl;

    ARIADNE_TEST_PRINT(hybrid_orbit);


/*
    cout << "Plotting orbit... " << flush;
    Figure fig;
    fig << hybrid_orbit;
    fig.write("test_hybrid_simulator-orbit");
    cout << "done" << endl;
*/

}
