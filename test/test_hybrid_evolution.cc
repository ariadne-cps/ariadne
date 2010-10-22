/***************************************************************************
 *            test_hybrid_evolution.cc
 *
 *  Copyright  2006-9  Pieter Collins
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
#include "hybrid_evolver-constrained.h"
#include "hybrid_evolver-working.h"
#include "graphics_interface.h"
#include "graphics.h"
#include "logging.h"

#include "models.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

int evolver_verbosity=0;


ScalarFunction c=ScalarFunction::constant(2,1.0);
ScalarFunction x0=ScalarFunction::coordinate(2,0);
ScalarFunction x1=ScalarFunction::coordinate(2,1);
DiscreteLocation q("q");
DiscreteEvent e("e");

Colour reach_set_colour(0.25,0.25,0.50);
Colour intermediate_set_colour(0.50,0.50,0.75);
Colour final_set_colour(0.75,0.75,1.00);
Colour initial_set_colour(0.75,0.75,1.00);
Colour guard_set_colour(0.75,0.75,0.75);

// Test evolution of realistic hybrid systems
class TestHybridEvolution
{
    HybridEvolverInterface& evolver;
  public:
    TestHybridEvolution(const HybridEvolverInterface&);
  public:
    void test() const;
    void test_bouncing_ball() const;
    void test_water_tank() const;
};

TestHybridEvolution::TestHybridEvolution(const HybridEvolverInterface& _evolver)
    : evolver(*_evolver.clone()) { }

void TestHybridEvolution::test() const {
    ARIADNE_TEST_CALL(test_bouncing_ball());
    ARIADNE_TEST_CALL(test_water_tank());
};

void TestHybridEvolution::test_bouncing_ball() const {
    MonolithicHybridAutomaton bouncing_ball;
    ScalarFunction c=ScalarFunction::constant(2,1.0);
    ScalarFunction x=ScalarFunction::coordinate(2,0);
    ScalarFunction v=ScalarFunction::coordinate(2,1);
    double lambda=0.5;
    bouncing_ball.new_mode(q,(v,-c));
    ARIADNE_TEST_WARN("IMPACT case not tested properly; instead adding small constant");
    bouncing_ball.new_transition(e,q,q,(x+0.000001,-lambda*v),-x,impact);

    double height=2.0;
    double radius=1.0/64;
    HybridBox initial(q,Box(2, height-radius,height+radius, -radius,+radius));
    HybridTime time(4.5,3);

    Orbit<HybridEnclosure> orbit=evolver.orbit(bouncing_ball,HybridEnclosure(initial),time,UPPER_SEMANTICS);

    ARIADNE_TEST_PRINT(orbit);
    if(orbit.final().size()!=1u) {
        ARIADNE_TEST_WARN("orbit.final().size()="<<orbit.final().size()<<"; expected 1. "
                          "This may indicate over-zealous splitting, and/or errors in detecting the end conditions.");
    }

    HybridEnclosure orbit_final_set=*orbit.final().begin();
    HybridBox expected_orbit_final_bounding_box=HybridBox(q,Box(2, 0.12,0.13, -0.04,0.04));
    ARIADNE_TEST_BINARY_PREDICATE(subset,orbit_final_set,expected_orbit_final_bounding_box);
    ARIADNE_TEST_PRINT(orbit.final().size());
    ARIADNE_TEST_PRINT(orbit_final_set.bounding_box());
    ARIADNE_TEST_PRINT(expected_orbit_final_bounding_box);

    Box bounding_box(2, -0.5,+2.5, -4.0, +4.0);
    plot("test_hybrid_evolution-bouncing_ball",bounding_box,
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}


void TestHybridEvolution::test_water_tank() const {
    // Declare some constants
    Real T(4.0);
    Real hmin(5.5);
    Real hmax(8.0);
    Real delta(0.05);
    Real lambda(0.02);
    Real rate(0.3);

    ScalarFunction height=ScalarFunction::coordinate(2,0);
    ScalarFunction aperture=ScalarFunction::coordinate(2,1);
    ScalarFunction zero=ScalarFunction::constant(2,0.0);
    ScalarFunction one=ScalarFunction::constant(2,1.0);

    DiscreteLocation open("open");
    DiscreteLocation opening("opening");
    DiscreteLocation closed("closed");
    DiscreteLocation closing("closing");

    DiscreteEvent start_opening("start_opening");
    DiscreteEvent start_closing("start_closing");
    DiscreteEvent finished_opening("finished_opening");
    DiscreteEvent finished_closing("finished_closing");

    // Create the tank object
    MonolithicHybridAutomaton watertank("Watertank");

    watertank.new_mode(open,(-lambda*height+rate*aperture,zero));
    watertank.new_mode(closed,(-lambda*height+rate*aperture,zero));
    watertank.new_mode(opening,(-lambda*height+rate*aperture,one/T));
    watertank.new_mode(closing,(-lambda*height+rate*aperture,-one/T));

    watertank.new_transition(start_opening,closed,opening,(height,aperture),hmin-height,urgent);
    watertank.new_transition(start_closing,open,closing,(height,aperture),height-hmax,urgent);
    watertank.new_transition(finished_opening,opening,open,(height,aperture),aperture-1,urgent);
    watertank.new_transition(finished_closing,closing,closed,(height,aperture),-aperture,urgent);

    DiscreteLocation initial_location=opening;
    Box initial_box(2, 0.0,0.05, 0.0,0.01);
    HybridEnclosure initial_enclosure(initial_location,initial_box);

    //HybridTime evolution_time(80.0,5);
    HybridTime evolution_time(80.0,8);

    dynamic_cast<HybridEvolverBase&>(evolver).parameters().maximum_step_size=1.0;

    Orbit<HybridEnclosure> orbit = evolver.orbit(watertank,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    if(orbit.final().size()!=1u) {
        ARIADNE_TEST_WARN("orbit.final().size()="<<orbit.final().size()<<"; expected 1. "
                          "This may indicate over-zealous splitting, and/or errors in detecting the end conditions.");
    }
    HybridEnclosure final_enclosure=HybridEnclosure(*orbit.final().begin());
    ARIADNE_TEST_PRINT(final_enclosure.bounding_box());
    ARIADNE_TEST_BINARY_PREDICATE(subset,final_enclosure,HybridBox(open,Box(2, 7.7,8.0, 0.999,1.001)));

    Box bounding_box(2, -0.1,9.1, -0.3,1.3);
    plot("test_hybrid_evolution-water_tank",bounding_box,
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());

}


int main(int argc, const char* argv[])
{
    if(argc>1) { evolver_verbosity=atoi(argv[1]); }
    TaylorModel::set_default_sweep_threshold(1e-12);
    GeneralHybridEvolver evolver;
    evolver.verbosity=evolver_verbosity;
    evolver.parameters().maximum_step_size=0.125;
    TestHybridEvolution(evolver).test();
    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

