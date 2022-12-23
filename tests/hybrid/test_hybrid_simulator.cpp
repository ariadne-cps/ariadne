/***************************************************************************
 *            test_hybrid_simulator.cpp
 *
 *  Copyright  2006-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <iostream>

#include "config.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/space.hpp"
#include "hybrid/hybrid_set.hpp"
#include "hybrid/hybrid_paving.hpp"
#include "hybrid/hybrid_orbit.hpp"
#include "hybrid/hybrid_time.hpp"
#include "hybrid/hybrid_automata.hpp"
#include "hybrid/hybrid_simulator.hpp"
#include "io/figure.hpp"
#include "conclog/logging.hpp"

#include "../test.hpp"

using namespace ConcLog;

using namespace Ariadne;
using namespace std;

class TestHybridSimulator
{
  private:
    HybridAutomaton system() const
    {
        const DiscreteLocation location1(1);
        const DiscreteLocation location2(2);
        const DiscreteEvent event3("e3");
        const DiscreteEvent event4("e4");
        const RealVariable x("x");
        const RealVariable y("y");
        const RealVariable z("z");

        HybridAutomaton automaton;
        Dyadic a(-0.5_x); Dyadic b(1.0_x); Dyadic cx1(3.0_x); Dyadic cx2(-1.0_x); Dyadic cz(1.0_x);
        automaton.new_mode(location1,{dot(x)= a*x-b*y+cx1,dot(y)=b*x+a*y} );
        automaton.new_mode(location2,{dot(x)= a*x-b*y+cx2,dot(y)=b*x+a*y,dot(z)=cz});
        automaton.new_transition(location1,event3,location2,{next(x)=x,next(y)=y,next(z)=y},x>=1,EventKind::URGENT);
        automaton.new_transition(location2,event4,location1,{next(x)=x,next(y)=y},x<=-1,EventKind::URGENT);

        return automaton;
    }

  public:
    Void test() const {
        ARIADNE_TEST_CALL(test_run());
    }

    Void test_run() const {
        typedef HybridSimulator::ApproximatePointType ApproximatePointType;
        typedef HybridSimulator::HybridApproximatePointType HybridApproximatePointType;

        const DiscreteLocation location1(1);
        const DiscreteLocation location2(2);
        const DiscreteEvent event3("e3");
        const DiscreteEvent event4("e4");
        const RealVariable x("x");
        const RealVariable y("y");

        // Make a hybrid automaton for the Van der Pol equation
        HybridAutomaton automaton=system();
        ARIADNE_TEST_PRINT(automaton);

        // Set up the evaluators
        HybridSimulator simulator(automaton);
        simulator.configuration().set_step_size(0.0625);

        cout << "configuration=" << simulator.configuration() << endl;

        // Define the initial box
        RealSpace space={x,y};
        RealPoint initial_point = Point<Real>{-0.00_x, 0.50_x};
        cout << "initial_point=" << initial_point << endl;
        ApproximatePointType approximate_initial_point = ApproximatePointType(initial_point,dp);
        HybridApproximatePointType initial_hybrid_point(location1,space,approximate_initial_point);
        HybridTime simulation_time(2.25_x,3);

        // Compute the reachable sets
        cout << "Computing orbit... "<<std::flush;
        auto orbit1=simulator.orbit(initial_hybrid_point,simulation_time);
        auto orbit2=simulator.orbit(HybridBoundedConstraintSet(location1,{x==0,y==0.5_dec}),simulation_time);
        cout << "done"<<std::endl;

        ARIADNE_TEST_EQUALS(orbit1.curve(0).size(),orbit2.curve(0).size());
    }
};

int main()
{
    TestHybridSimulator().test();
    return ARIADNE_TEST_FAILURES;
}