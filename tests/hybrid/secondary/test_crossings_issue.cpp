/***************************************************************************
 *            test_crossings_issue.cpp
 *
 *  Copyright  2008-18 Luca Geretti, Pieter Collins
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

#include "ariadne.hpp"
#include "../../test.hpp"

using namespace std;
using namespace Ariadne;

class TestCrossingsIssue
{
  public:
	TestCrossingsIssue(unsigned int verbosity);
    void test();
  private:
    unsigned int _verbosity;
    void test_crossings();
};

TestCrossingsIssue::TestCrossingsIssue(unsigned int verbosity) : _verbosity(verbosity) { }

void TestCrossingsIssue::test()
{
    // Set the accuracy of numerical output
    FloatDPApproximation::set_output_places(7);

    ARIADNE_TEST_CALL(test_crossings());
}

void TestCrossingsIssue::test_crossings()
{
    /// Build the Hybrid System
    RealConstant v("v",0.092_dec);
    RealConstant L("L",0.00025_dec);
    RealConstant M("M",0.0023_dec);

    /// Create a HybridAutomaton object
    AtomicHybridAutomaton automaton("compute_crossings_issue");

    /// Modes

    AtomicDiscreteLocation far("far");
    AtomicDiscreteLocation close("close");

    // Variables

    RealVariable x("x"); // X position
    TimeVariable t;

    automaton.new_mode(far, {dot(x)=v});
    automaton.new_mode(close, {dot(x)=v});

    DiscreteEvent comes("comes");
    DiscreteEvent leaves("leaves");

    automaton.new_transition(far,comes,close,{next(x)=x+L/4},sqr(x-M)<=sqr(L),EventKind::URGENT);
	automaton.new_transition(close,leaves,far,{next(x)=x+L/2},sqr(x-M)>=sqr(L),EventKind::URGENT);

    // Create a GeneralHybridEvolver object
    GeneralHybridEvolver evolver(automaton);
    evolver.verbosity = _verbosity;

    // Set the evolution parameters
    evolver.configuration().set_maximum_enclosure_radius(0.5);
    evolver.configuration().set_maximum_step_size(0.10);

    typedef GeneralHybridEvolver::OrbitType OrbitType;

    Real initial_x(0);
    HybridSet initial_set({automaton|far},{x==initial_x});

    Nat max_n = 3; Real max_t = 0.046_dec;
    HybridTime evolution_time(max_t,max_n);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    std::cout << "done." << std::endl;

    Real max_x = initial_x + max_t * v + 2*L;
    HybridBoxSet guard_set(automaton|far,{t.in(0,2*max_t),x.in(M-L,M+L)});
    plot("crossings_issue",Axes2d(0.0,t,1.1*max_t.get_d(), 0.0,x,1.1*max_x.get_d()), Colour(1.0,0.5,0.0), guard_set, Colour(0.0,0.5,1.0), orbit);

    for (auto enclosure : orbit.final()) {
        enclosure.reduce();
//        ARIADNE_TEST_ASSERT( definitely(enclosure.is_empty()) or enclosure.previous_events().size() == 2);
        if (not (enclosure.previous_events().size()==2 or definitely(enclosure.is_empty()))) {
            ListSet<HybridEnclosure> split_enclosures = enclosure.split();
            for (auto split_enclosure : split_enclosures) {
                split_enclosure.reduce();
                ARIADNE_TEST_ASSERT( split_enclosure.is_empty() or split_enclosure.previous_events().size() == 2);
            }
        }
    }
}

int main(Int argc, const char* argv[]) {

    auto verbosity=get_verbosity(argc,argv);

    TestCrossingsIssue(verbosity).test();
    return ARIADNE_TEST_FAILURES;
}

