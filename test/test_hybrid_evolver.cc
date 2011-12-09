/***************************************************************************
 *            test_hybrid_evolver.cc
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

#include "config.h"
#include "tuple.h"
#include "vector.h"
#include "matrix.h"
#include "function.h"
#include "box.h"
#include "list_set.h"
#include "evolution_parameters.h"
#include "integrator.h"
#include "orbit.h"
#include "graphics_interface.h"
#include "graphics.h"
#include "hybrid_automaton.h"
#include "hybrid_time.h"
#include "hybrid_set.h"
#include "hybrid_evolver.h"
#include "hybrid_graphics.h"
#include "logging.h"


#include "hybrid_automaton-composite.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

int evolver_verbosity=0;


RealVariable v0=RealVariable("x0");
RealVariable v1=RealVariable("x1");
RealScalarFunction z=RealScalarFunction::constant(2,0.0);
RealScalarFunction c=RealScalarFunction::constant(2,1.0);
RealScalarFunction x0=RealScalarFunction::coordinate(2,0);
RealScalarFunction x1=RealScalarFunction::coordinate(2,1);
DiscreteLocation q("q");
DiscreteEvent e("e");

Colour reach_set_colour(0.25,0.25,0.50);
Colour intermediate_set_colour(0.50,0.50,0.75);
Colour final_set_colour(0.75,0.75,1.00);
Colour initial_set_colour(0.75,0.75,1.00);
Colour guard_set_colour(0.75,0.75,0.75);

inline std::string operator+(const char* cstr, const std::string& str) { return std::string(cstr)+str; }
inline const char* cstr(const std::string& str) { return str.c_str(); }

class TestHybridEvolver
{
  private:
    shared_ptr<HybridEvolverBase> evolver_ptr;
    std::string evolver_name;
  public:
    TestHybridEvolver(const HybridEvolverInterface& evolver, const String& name);
    void test_all() const;
    void test_flow() const;
    void test_affine_flow() const;
    void test_exact_final_time() const;
    void test_maximum_steps() const;
    void test_urgent_event() const;
    void test_empty_interior() const;
    void test_partial_event() const;
    void test_step_size_event() const;
    void test_initially_active_event() const;
    void test_initially_active_attracting_event() const;
    void test_initially_active_repelling_event() const;
    void test_impact() const;
    void test_tangency() const;
    void test_simultaneous_events() const;
    void test_creep() const;
    void test_unwind() const;
    void test_permissive() const;

    void test_splitting_on_urgent_event() const;
    void test_affine_flow_system() const;
    void test_transverse_linear_crossing() const;
    void test_transverse_cubic_crossing() const;
    void test_transverse_cube_root_crossing() const;
    void test_constant_derivative_system() const;

};

TestHybridEvolver::TestHybridEvolver(const HybridEvolverInterface& evolver, const String& name)
    : evolver_ptr(dynamic_cast<HybridEvolverBase*>(evolver.clone()))
    , evolver_name(name)
{
    DRAWING_METHOD = AFFINE_DRAW;
    DRAWING_ACCURACY = 1;
}

void TestHybridEvolver::test_all() const {
    ARIADNE_TEST_PRINT(*evolver_ptr);
    ARIADNE_TEST_CALL(test_flow());
    ARIADNE_TEST_CALL(test_affine_flow());
    ARIADNE_TEST_CALL(test_exact_final_time());
    ARIADNE_TEST_CALL(test_maximum_steps());
    ARIADNE_TEST_CALL(test_urgent_event());
    ARIADNE_TEST_CALL(test_step_size_event());
    ARIADNE_TEST_CALL(test_empty_interior());
    ARIADNE_TEST_CALL(test_partial_event());
    ARIADNE_TEST_CALL(test_initially_active_event());
    ARIADNE_TEST_CALL(test_initially_active_attracting_event());
    ARIADNE_TEST_CALL(test_initially_active_repelling_event());
    ARIADNE_TEST_CALL(test_simultaneous_events());
    ARIADNE_TEST_CALL(test_impact());
    ARIADNE_TEST_CALL(test_creep());
    ARIADNE_TEST_CALL(test_unwind());
    ARIADNE_TEST_CALL(test_tangency());
    ARIADNE_TEST_CALL(test_permissive());

    ARIADNE_TEST_CALL(test_splitting_on_urgent_event());
    ARIADNE_TEST_CALL(test_affine_flow_system());
    ARIADNE_TEST_CALL(test_transverse_linear_crossing());
    ARIADNE_TEST_CALL(test_transverse_cubic_crossing());
    ARIADNE_TEST_CALL(test_transverse_cube_root_crossing());
    ARIADNE_TEST_CALL(test_constant_derivative_system());
}

void TestHybridEvolver::test_flow() const {
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q,(c,c/2));
    RealSpace space=automaton.continuous_state_space(q);
    HybridBox initial(q,space,Box(2, -0.125,0.125, -0.125,0.125));
    HybridTime time(2.5,3);
    evolver_ptr->parameters().maximum_step_size=1.0;

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);
    ARIADNE_TEST_CHECK_WARN(orbit.final().size(),1u);
    ARIADNE_TEST_CHECK_WARN(orbit.reach().size(),3u);
    ARIADNE_TEST_CHECK_WARN(orbit.intermediate().size(),2u);

    ARIADNE_TEST_PRINT(orbit);

    plot(cstr("test_hybrid_evolver-"+evolver_name+"-flow"),Axes2d(-0.5,v0,+3.5, -1.0,v1, +3.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}


void TestHybridEvolver::test_exact_final_time() const {
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q,(c,c/2));
    RealSpace space=automaton.continuous_state_space(q);

    HybridBox initial(q,space,Box(2, -0.125,0.125, -0.125,0.125));
    HybridTime time(2.0,1);

    evolver_ptr->parameters().maximum_step_size=1.0;

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);

    ARIADNE_TEST_CHECK_WARN(orbit.final().size(),1u);
    ARIADNE_TEST_CHECK_WARN(orbit.reach().size(),2u);

    plot(cstr("test_hybrid_evolver-"+evolver_name+"-exact_final_time"),Axes2d(-0.5,v0,+3.5, -1.0,v1,+3.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}


//! A simple test for the maximum number of steps
void TestHybridEvolver::test_maximum_steps() const {
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q,(c,c/2));
    automaton.new_transition(q,e,q,(x0-1.5,x1),x0+x1/16-0.5,urgent);
    //automaton.new_transition(e,q,q,(x0-1.5,x1),x0-0.5,urgent);
    // FIXME: Change so that hitting coordinate guard is not an error.

    RealSpace space=automaton.continuous_state_space(q);
    HybridBox initial(q,space,Box(2, -0.125,0.125, -0.125,0.125));
    HybridTime time(1.0,1);

    evolver_ptr->parameters().maximum_step_size=1.0;

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);
    ARIADNE_TEST_PRINT(orbit);

    ARIADNE_TEST_CHECK_WARN(orbit.final().size(),1u);
    ARIADNE_TEST_CHECK_WARN(orbit.reach().size(),1u);

    Axes2d axes(-1.5,space[0],+2.5, -1.0,space[1],+3.0);
    plot(cstr("test_hybrid_evolver-"+evolver_name+"-maximum_steps"),axes,
         //guard_set_colour,RealBoundedConstraintSet(bounding_box,x0+x1/16-0.5>=0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}

//! A simple test for an urgent event
void TestHybridEvolver::test_urgent_event() const {
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q,(c,c/2));
    automaton.new_transition(q,e,q,(x0-1.5,x1),x0+x1/16-0.5,urgent);
    //automaton.new_transition(e,q,q,(x0-1.5,x1),x0-0.5,urgent);
    // FIXME: Change so that hitting coordinate guard is not an error.

    RealSpace space=automaton.continuous_state_space(q);
    HybridBox initial(q,space,Box(2, -0.125,0.125, -0.125,0.125));
    HybridTime time(1.0,2);

    evolver_ptr->parameters().maximum_step_size=1.0;

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);
    ARIADNE_TEST_PRINT(orbit);

    ARIADNE_TEST_CHECK_WARN(orbit.final().size(),1u);
    ARIADNE_TEST_CHECK_WARN(orbit.reach().size(),2u);

    Axes2d axes(-1.5,space[0],+2.5, -1.0,space[1],+3.0);
    plot(cstr("test_hybrid_evolver-"+evolver_name+"-urgent_event"),axes,
         //guard_set_colour,RealBoundedConstraintSet(bounding_box,x0+x1/16-0.5>=0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}


// A test to ensure that an event which is active at the final time does actually occur
void TestHybridEvolver::test_empty_interior() const {
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q,(c,c/2));
    automaton.new_transition(q,e,q,(x0-2,x1),x0-1,urgent);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBox initial(q,space,Box(2, -0.125,0.25, -0.125,0.125));
    HybridTime time(2.0,3);

    evolver_ptr->parameters().maximum_step_size=2.0;

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);

    ARIADNE_TEST_CHECK_WARN(orbit.final().size(),1u);
    ARIADNE_TEST_CHECK_WARN(orbit.reach().size(),2u);

    plot(cstr("test_hybrid_evolver-"+evolver_name+"-empty_interior"),Axes2d(-1.5,space[0],+2.5, -1.0,space[1],+3.0),
         //guard_set_colour,Box(2,1.0,8.0,-8.0,+8.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}

// A test to ensure that an event which is active at the final time does actually occur
void TestHybridEvolver::test_partial_event() const {
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q,(c,c/2));
    automaton.new_transition(q,e,q,(x0-2,x1),x0-x1/16-2,urgent);
    //automaton.new_transition(q,e,q,(x0-2,x1),x0-2,urgent);
    //FIXME: Need to allow domain of TaylorFunction to have empty interior

    RealSpace space=automaton.continuous_state_space(q);
    HybridBox initial(q,space,Box(2, -0.125,0.25, -0.125,0.125));
    HybridTime time(2.0,3);

    evolver_ptr->parameters().maximum_step_size=2.0;

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);

    ARIADNE_TEST_CHECK_WARN(orbit.final().size(),2u);
    ARIADNE_TEST_CHECK_WARN(orbit.reach().size(),2u);

    Axes2d axes(-1.5,v0,+2.5, -1.0,v1,+3.0);
    plot(cstr("test_hybrid_evolver-"+evolver_name+"-partial_event"),axes,
         //guard_set_colour,RealBoundedConstraintSet(bounding_box,x0-x1/16-2>=0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}


// A test to ensure that an event which would be active after a time step
// but is avoided because the evolution is completed before the event would occur,
// does not occur
void TestHybridEvolver::test_step_size_event() const {
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q,(c,c/2));
    automaton.new_transition(q,e,q,(x0-2,x1),x0-2.0,urgent);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBox initial(q,space,Box(2, -0.125,0.125, -0.125,0.125));
    HybridTime time(1.0,3);

    evolver_ptr->parameters().maximum_step_size=2.0;

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);

    ARIADNE_TEST_CHECK_WARN(orbit.final().size(),1u);
    ARIADNE_TEST_CHECK_WARN(orbit.reach().size(),1u);

    plot(cstr("test_hybrid_evolver-"+evolver_name+"-step_size_event"),Axes2d(-0.5,v0,+2.5, -1.0,v1,+3.0),
         //guard_set_colour,Box(2,2.0,8.0,-8.0,+8.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}


void TestHybridEvolver::test_initially_active_event() const {
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q,(c,c));
    automaton.new_transition(q,e,q,(x0+1,x1),-x0,urgent);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBox initial(q,space,Box(2, -1.625,-1.375, -0.125,0.125));
    HybridTime time(1.0,4);

    evolver_ptr->parameters().maximum_step_size=2.0;

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);
    ARIADNE_TEST_CHECK_WARN(orbit.final().size(),1u);

    // There should be two components of the reachable set which come from initially active events, and one from flowing
    ARIADNE_TEST_CHECK_WARN(orbit.reach().size(),3u);

    plot(cstr("test_hybrid_evolver-"+evolver_name+"-initially_active"),Axes2d(-2.0,v0,+2.0, -1.0,v1,+2.0),
         //guard_set_colour,Box(2,-8.0,0.0,-8.0,+8.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());


}

void TestHybridEvolver::test_initially_active_attracting_event() const {
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q,(-0.5*c,c));
    automaton.new_transition(q,e,q,(x0+1.0,x1),-x0-x1*1.0/256,urgent);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBox initial(q,space,Box(2, -0.125,0.25, -0.125,0.125));
    HybridTime time(1.0,4);

    evolver_ptr->parameters().maximum_step_size=2.0;

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);

    plot(cstr("test_hybrid_evolver-"+evolver_name+"-initially_active_attracting"),Axes2d(-1.0,v0,+2.0, -1.0,v1,+2.0),
         //guard_set_colour,Box(2,-1.0,0.0,-8.0,+8.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());


}

void TestHybridEvolver::test_initially_active_repelling_event() const {
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q,(+0.5*c,c));
    automaton.new_transition(q,e,q,(x0+1,x1),-x0,urgent);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBox initial(q,space,Box(2, -0.125,0.25, -0.125,0.125));
    HybridTime time(1.0,2);

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);

    plot(cstr("test_hybrid_evolver-"+evolver_name+"-initially_active_repelling"),Axes2d(-1.0,v0,+2.0, -1.0,v1,+2.0),
         guard_set_colour,HybridBox(q,space,Box(2,-1.0,0.0,-8.0,+8.0)),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}



void TestHybridEvolver::test_impact() const {
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q,(x1,Real(0)*c));
    automaton.new_transition(q,e,q,(x0,x1-2),x0-1,impact);
    //automaton.new_transition(q,e,q,(x0+0.001*x1-0.0004,x1-2),x0-1,urgent);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBox initial(q,space,Box(2, 0.4375,0.5625, 0.9375,1.0625));
    HybridTime time(2.0,3);

    evolver_ptr->parameters().maximum_step_size=2.0;

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);
    //ARIADNE_TEST_CHECK(orbit.final().size(),2u);

    plot(cstr("test_hybrid_evolver-"+evolver_name+"-impact"),Axes2d(-3.0,v0,+2.0, -4.0,v1,+2.0),
         //guard_set_colour,Box(2,1.0,8.0,-8.0,+8.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}

void TestHybridEvolver::test_tangency() const {
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q,(c,z));
    automaton.new_transition(q,e,q,(x0,x1-1),x1-sqr(x0),urgent);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBox initial(q,space,Box(2, -1.125,-0.875, -0.25,0.25));
    HybridTime time(2.0,3);

    evolver_ptr->parameters().maximum_step_size=2.0;

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);
    ARIADNE_TEST_CHECK_WARN(orbit.final().size(),2u);
    ARIADNE_TEST_CHECK_WARN(orbit.reach().size(),3u);
    ARIADNE_TEST_PRINT(orbit);

    Axes2d axes(-2.0,v0,+2.0, -2.0,v1,+1.0);
    plot(cstr("test_hybrid_evolver-"+evolver_name+"-tangency"),axes,
         //guard_set_colour,RealBoundedConstraintSet(bounding_box,x1-x0*x0>=0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}

#include "affine_set.h"
#include <include/hybrid_set.h>

void TestHybridEvolver::test_simultaneous_events() const {
    MonolithicHybridAutomaton automaton;
    DiscreteEvent e1("e1");
    DiscreteEvent e2("e2");
    automaton.new_mode(q,(c,c));
    automaton.new_transition(q,e1,q,(x0-1,x1-2),x0-1.0,urgent);
    automaton.new_transition(q,e2,q,(x0-2,x1-1),x1-1.0,urgent);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBox initial(q,space,Box(2, -0.25,0.125, -0.125,0.25));
    HybridTime time(2.5,4);

    evolver_ptr->parameters().maximum_step_size=4.0;

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);
    ARIADNE_TEST_CHECK_WARN(orbit.final().size(),2u);

    plot(cstr("test_hybrid_evolver-"+evolver_name+"-simultaneous_events"),Axes2d(-3.0,v0,+2.0, -3.0,v1,+2.0),
         //guard_set_colour,Box(2,1.0,8.0,-8.0,+8.0),
         //guard_set_colour,Box(2,-8.0,+8.0,1.0,8.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}


void TestHybridEvolver::test_creep() const {
    MonolithicHybridAutomaton automaton;
    DiscreteEvent e("e");
    automaton.new_mode(q,(c,c));
    automaton.new_transition(q,e,q,(x0-1,x1),x0-1.0,urgent);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBox initial(q,space,Box(2, -0.25,0.125, -0.125,0.25));
    HybridTime time(1.5,4);

    evolver_ptr->parameters().maximum_step_size=1.0;

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);
    ARIADNE_TEST_PRINT(orbit);
    ARIADNE_TEST_PRINT(*orbit.reach().begin());
    ARIADNE_TEST_PRINT(HybridEnclosure(*orbit.final().begin()).bounding_box());
    ARIADNE_TEST_CHECK_WARN(orbit.final().size(),1u);
    ARIADNE_TEST_CHECK_WARN(orbit.reach().size(),3u);
    ARIADNE_TEST_BINARY_PREDICATE(inside,HybridEnclosure(*orbit.final().begin()),HybridBox(q,space,Box(2, 0.24,0.635, 1.365,1.76)));

    plot(cstr("test_hybrid_evolver-"+evolver_name+"-creep"),Axes2d(-1.5,v0,+1.5, -0.5,v1,+3.5),
         //guard_set_colour,Box(2,1.0,8.0,-8.0,+8.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}

void TestHybridEvolver::test_unwind() const {
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q,(c,c));
    //automaton.new_transition(q,e,q,(x0-3,x1-1),x0-1,urgent);
    automaton.new_transition(q,e,q,(x0-3,x1-1),x0-x1/16-1,urgent);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBox initial(q,space,Box(2, -0.25,0.125, -0.125,0.25));
    HybridTime time(3.0,4);

    evolver_ptr->parameters().maximum_step_size=2.0;

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);
    ARIADNE_TEST_CHECK_WARN(orbit.final().size(),1u);

    Axes2d axes(-2.5,v0,+1.5, -0.5,v1,+2.5);
    plot(cstr("test_hybrid_evolver-"+evolver_name+"-unwind"),axes,
         //guard_set_colour,RealBoundedConstraintSet(bounding_box,x0-x1/16-1>=0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}

void TestHybridEvolver::test_permissive() const {
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q,(c,c/2));
    DiscreteEvent i("i");
    // Invariant has form f(x)<=0
    automaton.new_invariant(q,i,(x0-2.125));
    automaton.new_transition(q,e,q,(x0-3,x1),x0-1.875,permissive);
    ARIADNE_TEST_PRINT(automaton);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBox initial(q,space,Box(2, -0.125,0.125, -0.125,0.125));
    HybridTime time(3.0,4);

    evolver_ptr->parameters().maximum_step_size=1.0;

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);
    ARIADNE_TEST_PRINT(orbit);
    ARIADNE_TEST_CHECK_WARN(orbit.final().size(),1u);

    HybridBox guard_set(q,space,Box(2, 2-0.125,2+0.125, -0.5,+2.5));
    Axes2d axes(-4.5,v0,+2.5, -0.5,v1,+2.5);
    plot(cstr("test_hybrid_evolver-"+evolver_name+"-permissive"),axes,
         guard_set_colour,guard_set,
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}







void
TestHybridEvolver::test_affine_flow() const
{
    evolver_ptr->parameters().maximum_step_size=0.5;
    RealVariable x0("x");
    RealVariable x1("x1");
    DiscreteLocation q("q");

    HybridAutomaton system;
    system.new_mode(q,(dot(x0)=1.0-x0-2*x1,dot(x1)=+2*x0-x1));

    const double r=0.125;
    HybridExpressionSet initial_set(q,(2-r<=x0<=2+r,-r<=x1<=r));
    HybridEnclosure initial_enclosure=this->evolver_ptr->enclosure(system,initial_set);

    HybridTime evolution_time=HybridTime(1.5,2);
    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(system,initial_enclosure,evolution_time);
    ARIADNE_TEST_PRINT(orbit);
    ARIADNE_TEST_EQUAL(orbit.final().size(),1);
    ARIADNE_TEST_COMPARE(orbit.reach().size(),>=,3);
    if(orbit.reach().size()>12) { ARIADNE_TEST_WARN("Continuous evolution may use very short step size."); }
    else if(orbit.reach().size()!=3u) { ARIADNE_TEST_WARN("Continuous evolution may use shorter step size than necessary."); }

    Axes2d axes(-1.0,x0,3.0, -2.0,x1,2.0);
    plot(cstr("test_hybrid_evolver-"+evolver_name+"-affine_flow"), axes,
         reach_set_colour,orbit.reach(), intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(), initial_set_colour,orbit.initial());
}


void
TestHybridEvolver::test_splitting_on_urgent_event() const
{
    evolver_ptr->parameters().maximum_step_size=2.0;

    RealVariable x("x");
    RealVariable y("y");
    StringVariable q("");
    DiscreteLocation upwards(q|"upwards");
    DiscreteLocation downwards(q|"downwards");
    DiscreteEvent changeup("changeup");
    DiscreteEvent changedown("changedown");
    DiscreteEvent block("block");

    HybridAutomaton system;
    system.new_mode(upwards,(dot(x)=1.0,dot(y)=1.0));
    system.new_mode(downwards,(dot(x)=1.0,dot(y)=-1.0));

    system.new_guard(upwards,changedown,y>=2.0,urgent);
    system.new_transition(upwards,changedown,downwards,(next(x)=x+1,next(y)=y));

    ARIADNE_TEST_PRINT(system);

    DiscreteLocation initial_location(upwards);
    const double r=0.125;
    RealVariableBox initial_box((-r<=x<=+r, -r<=y<=+r));
    HybridEnclosure initial_enclosure=this->evolver_ptr->enclosure(system,HybridExpressionSet(initial_location,initial_box));
    HybridTime evolution_time=HybridTime(4.0,2);
    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(system,initial_enclosure,evolution_time);
    ARIADNE_TEST_PRINT(orbit);
    ARIADNE_TEST_CHECK_WARN(orbit.final().size(),1u);
    ARIADNE_TEST_CHECK_WARN(orbit.reach().size(),2u);

    Axes2d axes(-1.0,x,11.0, -6.0,y,6.0);
    HybridExpressionSet guard_set(upwards,(-1.0<=x<=11.0,2.0<=y<=6.0));
    plot(cstr("test_hybrid_evolver-"+evolver_name+"-splitting_on_urgent_event"), axes, guard_set_colour,guard_set,
         reach_set_colour,orbit.reach(), intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(), initial_set_colour,orbit.initial());
}


void
TestHybridEvolver::test_affine_flow_system() const
{
    evolver_ptr->parameters().maximum_step_size=8.0;

    // A hybrid automaton with a single component and two locations,
    // with an affine flow in each, and affine guards and resets.
    // This should be very easy to analyse numerically, and is there to test
    // switching logic

    RealVariable x("x");
    RealVariable y("y");
    StringVariable q("");
    DiscreteLocation upwards(q|"upwards");
    DiscreteLocation downwards(q|"downwards");
    DiscreteEvent changeup("changeup");
    DiscreteEvent changedown("changedown");
    DiscreteEvent block("block");

    HybridAutomaton system;
    system.new_mode(upwards,(dot(x)=1.0,dot(y)=1.0));
    system.new_mode(downwards,(dot(x)=1.0,dot(y)=-1.0));

    system.new_guard(upwards,changedown,y>=2.0,urgent);
    system.new_invariant(downwards,y>=-2.0,block);
    system.new_guard(downwards,changeup,y<=-1.5,permissive);
    system.new_transition(upwards,changedown,downwards,(next(x)=x+1,next(y)=y));
    system.new_transition(downwards,changeup,upwards,(next(x)=x+1,next(y)=y));

    ARIADNE_TEST_PRINT(system);


    DiscreteLocation initial_location(upwards);
    Real r=0.125;
    RealVariableBox initial_box((-r<=x<=+r, -r<=y<=+r));
    HybridEnclosure initial_enclosure=this->evolver_ptr->enclosure(system,HybridExpressionSet(initial_location,initial_box));
    HybridTime evolution_time(0.0,0u);

    Axes2d axes(-1.0,x,3.0, -1.0,y,3.0);
    HybridExpressionSet guard_set_changeup(downwards,(-1.0<=x<=11.0,-2.0<=y<=-1.5));

    evolution_time=HybridTime(1.0,3);
    ARIADNE_TEST_PRINT(evolution_time);
    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(system,initial_enclosure,evolution_time);
    ARIADNE_TEST_PRINT(orbit);
    ARIADNE_TEST_EQUAL(orbit.final().size(),1u);
    ARIADNE_TEST_EQUAL(orbit.reach().size(),1u);

    evolution_time=HybridTime(2.0,3);
    ARIADNE_TEST_PRINT(evolution_time);
    orbit=evolver_ptr->orbit(system,initial_enclosure,evolution_time);
    ARIADNE_TEST_PRINT(orbit);
    ARIADNE_TEST_EQUAL(orbit.reach().size(),2u);
    ARIADNE_TEST_EQUAL(orbit.final().size(),2u);
    plot(cstr("test_hybrid_evolver-"+evolver_name+"-hysterestis-t=2"), axes, guard_set_colour,guard_set_changeup,
         reach_set_colour,orbit.reach(), intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(), initial_set_colour,orbit.initial());

    evolution_time=HybridTime(4.0,3);
    ARIADNE_TEST_PRINT(evolution_time);
    orbit=evolver_ptr->orbit(system,initial_enclosure,evolution_time);
    ARIADNE_TEST_PRINT(orbit);
    ARIADNE_TEST_EQUAL(orbit.reach().size(),2u);
    ARIADNE_TEST_EQUAL(orbit.final().size(),1u);

    evolution_time=HybridTime(6.0,3);
    ARIADNE_TEST_PRINT(evolution_time);
    orbit=evolver_ptr->orbit(system,initial_enclosure,evolution_time);
    ARIADNE_TEST_PRINT(orbit);
    ARIADNE_TEST_EQUAL(orbit.reach().size(),3u);
    ARIADNE_TEST_EQUAL(orbit.final().size(),2u);

    plot(cstr("test_hybrid_evolver-"+evolver_name+"-hysterestis-t=6"), axes, guard_set_colour,guard_set_changeup,
         reach_set_colour,orbit.reach(), intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(), initial_set_colour,orbit.initial());

    evolution_time=HybridTime(8.0,3);
    ARIADNE_TEST_PRINT(evolution_time);
    orbit=evolver_ptr->orbit(system,initial_enclosure,evolution_time);
    ARIADNE_TEST_PRINT(orbit);
    ARIADNE_TEST_EQUAL(orbit.reach().size(),3u);
    ARIADNE_TEST_EQUAL(orbit.final().size(),1u);
    plot(cstr("test_hybrid_evolver-"+evolver_name+"-hysterestis-t=8"), axes, guard_set_colour,guard_set_changeup,
         reach_set_colour,orbit.reach(), intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(), initial_set_colour,orbit.initial());
}




void TestHybridEvolver::test_constant_derivative_system() const
{
    // Test the system (d(x),d(y))=(1,0) with reset (x',y')=(x-2,y) when x+y>0
    // Starting in a small box near the origin, the system should return to
    // the initial condition after time 2
    DiscreteLocation q1("q1"), q2("q2");
    DiscreteEvent e("e");
    RealVariable x("x"), y("y");

    HybridAutomaton system;
    system.new_mode(q1,(dot(x)=1,dot(y)=0));
    system.new_mode(q2,(dot(x)=1,dot(y)=0));
    system.new_transition(q1,e,q2,(prime(x)=x-2,prime(y)=y),x+y>=0,urgent);

    Real r=0.0625;
    RealVariableBox initial_box((-r<=x<=+r,-r<=y<=+r));
    HybridExpressionSet initial_set(q1,initial_box);


    ARIADNE_TEST_PRINT(system);
    ARIADNE_TEST_PRINT(initial_set);

    // Test continuous evolution with a single transverse jump
    HybridTime evolution_time(2.0,2);
    ARIADNE_TEST_PRINT(evolution_time);

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(system,initial_set,evolution_time,UPPER_SEMANTICS);
    ARIADNE_TEST_PRINT(orbit);

    ListSet<HybridEnclosure> final_set=orbit.final();
    ARIADNE_TEST_PRINT(final_set);
    HybridEnclosure expected_final_set=evolver_ptr->enclosure(system,HybridExpressionSet(q2,initial_box));
    ARIADNE_TEST_PRINT(expected_final_set);

    ARIADNE_TEST_COMPARE(norm(final_set[0].space_function()-expected_final_set.space_function()),<,1e-14);

}



void TestHybridEvolver::test_transverse_linear_crossing() const
{
    DiscreteLocation q1("q1"), q2("q2");
    RealVariable x("x"), y("y");
    Real r=1.0/8;
    Float tol=1e-5;
    RealExpression guard=x+y/2-1;

    HybridAutomaton system;
    system.new_mode(q1,(dot(x)=1,dot(y)=0));
    system.new_mode(q2,(dot(x)=0,dot(y)=1));
    system.new_transition(q1,e,q2,(prime(x)=x,prime(y)=y),guard>=0,urgent);
    HybridExpressionSet initial_set(q1,(-r<=x<=r,-r<=y<=r));
    HybridTime evolution_time(2.0,3);

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(system,initial_set,evolution_time,UPPER_SEMANTICS);
    Enclosure final_enclosure=orbit.final()[0].continuous_state_set();

    RealSpace final_space=system.continuous_state_space(q2);
    RealExpression ct=-guard; // Crossing time
    //RealVectorFunction function=make_function((x+ct,y+2-ct),final_space);
    RealVectorFunction function=join(make_function(x+ct,final_space),make_function(y+2-ct,final_space));
    Enclosure expected_final_enclosure(initial_set.continuous_state_set(final_space),evolver_ptr->function_factory());
    expected_final_enclosure.apply_map(function);

    Vector<Interval> tolerance(2,Interval(-tol,+tol));
    ARIADNE_TEST_BINARY_PREDICATE(refines,expected_final_enclosure.space_function()[0],final_enclosure.space_function()[0]);
    ARIADNE_TEST_BINARY_PREDICATE(refines,final_enclosure.space_function()[0],expected_final_enclosure.space_function()[0]+tolerance[0]);

    Axes2d axes(-1.0<=x<=3.0, -1.0<=y<=3.0);
    HybridExpressionSet guard_set(q1,(-1.0<=x<=3.0,-1.0<=y<=3.0),guard>=0);
    plot(cstr("test_hybrid_evolver-"+evolver_name+"-transverse_linear_crossing"), axes, guard_set_colour,guard_set,
         reach_set_colour,orbit.reach(), intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(), initial_set_colour,orbit.initial());
}

void TestHybridEvolver::test_transverse_cubic_crossing() const
{
    DiscreteLocation q1("q1"), q2("q2");
    RealVariable x("x"), y("y");
    Real r=1.0/8;
    Float tol=1e-5;
    RealExpression guard=x-(1+y/2+y*y*y);
    HybridAutomaton system;
    system.new_mode(q1,(dot(x)=1,dot(y)=0));
    system.new_mode(q2,(dot(x)=0,dot(y)=1));
    system.new_transition(q1,e,q2,(prime(x)=x,prime(y)=y),guard>=0,urgent);
    HybridExpressionSet initial_set(q1,(-r<=x<=r,-r<=y<=r));
    HybridTime evolution_time(2.0,3);

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(system,initial_set,evolution_time,UPPER_SEMANTICS);

    Axes2d axes(-1.0<=x<=3.0, -1.0<=y<=3.0);
    HybridExpressionSet guard_set(q1,(-1.0<=x<=3.0,-1.0<=y<=3.0),guard>=0);
    HybridEnclosure guard_enclosure=evolver_ptr->enclosure(system,guard_set);
    plot(cstr("test_hybrid_evolver-"+evolver_name+"-transverse_cubic_crossing"), axes, guard_set_colour,guard_enclosure,
         reach_set_colour,orbit.reach(), intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(), initial_set_colour,orbit.initial());
}

void TestHybridEvolver::test_transverse_cube_root_crossing() const
{
    DiscreteLocation q1("q1"), q2("q2");
    RealVariable x("x"), y("y");
    Real r=1.0/32;
    Float tol=1e-5;
    RealExpression guard=((x-1)*(x-1)+1.0)*(x-1)-y-1./64;
    HybridAutomaton system;
    system.new_mode(q1,(dot(x)=1,dot(y)=0));
    system.new_mode(q2,(dot(x)=0,dot(y)=1));
    system.new_transition(q1,e,q2,(prime(x)=x,prime(y)=y),guard>=0,urgent);
    HybridExpressionSet initial_set(q1,(-r<=x<=r,-r<=y<=r));
    HybridTime evolution_time(2.0,3);

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(system,initial_set,evolution_time,UPPER_SEMANTICS);

    Axes2d axes(-1.0<=x<=3.0, -1.0<=y<=3.0);
    HybridExpressionSet guard_set(q1,(-1.0<=x<=3.0,-1.0<=y<=3.0),guard>=0);
    HybridEnclosure guard_enclosure=evolver_ptr->enclosure(system,guard_set);
    plot(cstr("test_hybrid_evolver-"+evolver_name+"-transverse_cube_root_crossing"), axes, guard_set_colour,guard_enclosure,
         reach_set_colour,orbit.reach(), intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(), initial_set_colour,orbit.initial());

}




int main(int argc, const char* argv[])
{
    int evolver_verbosity=get_verbosity(argc,argv);

    DRAWING_METHOD = AFFINE_DRAW; DRAWING_ACCURACY = 2u;

    GeneralHybridEvolver evolver;
    evolver.set_integrator(TaylorSeriesIntegrator(1e-3));
    evolver.verbosity=evolver_verbosity;
    //TestHybridEvolver(evolver,"general_hybrid_evolver").test_unwind();
    //TestHybridEvolver(evolver,"general_hybrid_evolver").test_transverse_only();
    ARIADNE_TEST_CALL(TestHybridEvolver(evolver,"general").test_all());

    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

