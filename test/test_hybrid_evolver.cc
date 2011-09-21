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

#include "tuple.h"
#include "vector.h"
#include "matrix.h"
#include "function.h"
#include "taylor_set.h"
#include "taylor_function.h"
#include "box.h"
#include "list_set.h"
#include "evolution_parameters.h"
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

class TestSimpleHybridEvolver
{
  private:
    shared_ptr<HybridEvolverBase> evolver_ptr;
    std::string evolver_name;
  public:
    TestSimpleHybridEvolver(const HybridEvolverInterface& evolver, const String& name);
    void test_all() const;
    void test_flow() const;
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
};

TestSimpleHybridEvolver::TestSimpleHybridEvolver(const HybridEvolverInterface& evolver, const String& name)
    : evolver_ptr(dynamic_cast<HybridEvolverBase*>(evolver.clone()))
    , evolver_name(name)
{
    DRAWING_METHOD = AFFINE_DRAW;
    DRAWING_ACCURACY = 1;
}

void TestSimpleHybridEvolver::test_all() const {
    ARIADNE_TEST_PRINT(*evolver_ptr);
    ARIADNE_TEST_CALL(test_flow());
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
}

void TestSimpleHybridEvolver::test_flow() const {
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

    plot(cstr("test_"+evolver_name+"-flow"),Axes2d(-0.5,v0,+3.5, -1.0,v1, +3.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}


void TestSimpleHybridEvolver::test_exact_final_time() const {
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q,(c,c/2));
    RealSpace space=automaton.continuous_state_space(q);

    HybridBox initial(q,space,Box(2, -0.125,0.125, -0.125,0.125));
    HybridTime time(2.0,1);

    evolver_ptr->parameters().maximum_step_size=1.0;

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);

    ARIADNE_TEST_CHECK_WARN(orbit.final().size(),1u);
    ARIADNE_TEST_CHECK_WARN(orbit.reach().size(),2u);

    plot(cstr("test_"+evolver_name+"-exact_final_time"),Axes2d(-0.5,v0,+3.5, -1.0,v1,+3.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}


//! A simple test for the maximum number of steps
void TestSimpleHybridEvolver::test_maximum_steps() const {
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
    plot(cstr("test_"+evolver_name+"-maximum_steps"),axes,
         //guard_set_colour,BoundedConstraintSet(bounding_box,x0+x1/16-0.5>=0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}

//! A simple test for an urgent event
void TestSimpleHybridEvolver::test_urgent_event() const {
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
    plot(cstr("test_"+evolver_name+"-urgent_event"),axes,
         //guard_set_colour,BoundedConstraintSet(bounding_box,x0+x1/16-0.5>=0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}


// A test to ensure that an event which is active at the final time does actually occur
void TestSimpleHybridEvolver::test_empty_interior() const {
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

    plot(cstr("test_"+evolver_name+"-empty_interior"),Axes2d(-1.5,space[0],+2.5, -1.0,space[1],+3.0),
         //guard_set_colour,Box(2,1.0,8.0,-8.0,+8.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}

// A test to ensure that an event which is active at the final time does actually occur
void TestSimpleHybridEvolver::test_partial_event() const {
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
    plot(cstr("test_"+evolver_name+"-partial_event"),axes,
         //guard_set_colour,BoundedConstraintSet(bounding_box,x0-x1/16-2>=0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}


// A test to ensure that an event which would be active after a time step
// but is avoided because the evolution is completed before the event would occur,
// does not occur
void TestSimpleHybridEvolver::test_step_size_event() const {
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

    plot(cstr("test_"+evolver_name+"-step_size_event"),Axes2d(-0.5,v0,+2.5, -1.0,v1,+3.0),
         //guard_set_colour,Box(2,2.0,8.0,-8.0,+8.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}


void TestSimpleHybridEvolver::test_initially_active_event() const {
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

    plot(cstr("test_"+evolver_name+"-initially_active"),Axes2d(-2.0,v0,+2.0, -1.0,v1,+2.0),
         //guard_set_colour,Box(2,-8.0,0.0,-8.0,+8.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());


}

void TestSimpleHybridEvolver::test_initially_active_attracting_event() const {
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q,(-0.5*c,c));
    automaton.new_transition(q,e,q,(x0+1.0,x1),-x0-x1*1.0/256,urgent);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBox initial(q,space,Box(2, -0.125,0.25, -0.125,0.125));
    HybridTime time(1.0,4);

    evolver_ptr->parameters().maximum_step_size=2.0;

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);

    plot(cstr("test_"+evolver_name+"-initially_active_attracting"),Axes2d(-1.0,v0,+2.0, -1.0,v1,+2.0),
         //guard_set_colour,Box(2,-1.0,0.0,-8.0,+8.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());


}

void TestSimpleHybridEvolver::test_initially_active_repelling_event() const {
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q,(+0.5*c,c));
    automaton.new_transition(q,e,q,(x0+1,x1),-x0,urgent);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBox initial(q,space,Box(2, -0.125,0.25, -0.125,0.125));
    HybridTime time(1.0,2);

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(automaton,initial,time,UPPER_SEMANTICS);

    plot(cstr("test_"+evolver_name+"-initially_active_repelling"),Axes2d(-1.0,v0,+2.0, -1.0,v1,+2.0),
         //guard_set_colour,Box(2,-1.0,0.0,-8.0,+8.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}



void TestSimpleHybridEvolver::test_impact() const {
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

    plot(cstr("test_"+evolver_name+"-impact"),Axes2d(-3.0,v0,+2.0, -4.0,v1,+2.0),
         //guard_set_colour,Box(2,1.0,8.0,-8.0,+8.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}

void TestSimpleHybridEvolver::test_tangency() const {
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
    plot(cstr("test_"+evolver_name+"-tangency"),axes,
         //guard_set_colour,BoundedConstraintSet(bounding_box,x1-x0*x0>=0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}

#include "affine_set.h"

void TestSimpleHybridEvolver::test_simultaneous_events() const {
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

    plot(cstr("test_"+evolver_name+"-simultaneous_events"),Axes2d(-3.0,v0,+2.0, -3.0,v1,+2.0),
         //guard_set_colour,Box(2,1.0,8.0,-8.0,+8.0),
         //guard_set_colour,Box(2,-8.0,+8.0,1.0,8.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}


void TestSimpleHybridEvolver::test_creep() const {
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
    ARIADNE_TEST_BINARY_PREDICATE(subset,HybridEnclosure(*orbit.final().begin()),HybridBox(q,space,Box(2, 0.24,0.635, 1.365,1.76)));

    plot(cstr("test_"+evolver_name+"-creep"),Axes2d(-1.5,v0,+1.5, -0.5,v1,+3.5),
         //guard_set_colour,Box(2,1.0,8.0,-8.0,+8.0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}

void TestSimpleHybridEvolver::test_unwind() const {
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
    plot(cstr("test_"+evolver_name+"-unwind"),axes,
         //guard_set_colour,BoundedConstraintSet(bounding_box,x0-x1/16-1>=0),
         reach_set_colour,orbit.reach(),
         intermediate_set_colour,orbit.intermediate(),
         final_set_colour,orbit.final(),
         initial_set_colour,orbit.initial());
}

/*

class TestContraintHybridEvolver
{
  private:
    shared_ptr<HybridEvolverBase> evolver_ptr;
    std::string evolver_name;
  private:
    static CompositeHybridAutomaton affine_flow_system();
  public:
    TestContraintHybridEvolver(const HybridEvolverInterface& evolver, const String& name);
  public:
    void test() const;
    void test_flow_only() const;
    void test_affine_flow_system() const;
    void test_splitting_on_urgent_event() const;
};


TestContraintHybridEvolver::TestContraintHybridEvolver(const HybridEvolverInterface& evolver, const String& name)
    : evolver_ptr(dynamic_cast<HybridEvolverBase*>(evolver.clone()))
    , evolver_name(name)
{
}

CompositeHybridAutomaton
TestContraintHybridEvolver::affine_flow_system() {
    // A hybrid automaton with a single component and two locations,
    // with an affine flow in each, and affine guards and resets.
    // This should be very easy to analyse numerically, and is there to test
    // switching logic

    RealVariable x("x");
    RealVariable y("y");
    StringConstant upwards("upwards");
    StringConstant downwards("downwards");
    DiscreteEvent changeup("changeup");
    DiscreteEvent changedown("changedown");
    DiscreteEvent block("block");

    AtomicHybridAutomaton affine("Affine Flow Automaton");
    affine.new_mode(upwards,(dot(x)=1.0,dot(y)=1.0));
    affine.new_mode(downwards,(dot(x)=1.0,dot(y)=-1.0));

    affine.new_urgent_guard(upwards,changedown,y>=2.0);
    affine.new_invariant(downwards,y>=-2.0,block);
    affine.new_permissive_guard(downwards,changeup,y<=-1.5);
    affine.new_transition(upwards,changedown,downwards,(next(x)=x+1,next(y)=y));
    affine.new_transition(downwards,changeup,upwards,(next(x)=x+1,next(y)=y));

    return CompositeHybridAutomaton(affine);
}



void
TestContraintHybridEvolver::test() const
{
    test_flow_only();
    //test_splitting_on_urgent_event();
    //test_affine_flow_system();
}


void
TestContraintHybridEvolver::test_flow_only() const
{
    evolver_ptr->parameters().maximum_step_size=0.5;
    CompositeHybridAutomaton system=affine_flow_system();

    DiscreteLocation upwards("upwards");
    DiscreteLocation initial_location(upwards);
    const double r=0.125;
    Box initial_box(2,-r,+r, -r,+r);
    HybridEnclosure initial_enclosure=this->evolver_ptr->enclosure(HybridBox(initial_location,initial_box));

    HybridTime evolution_time=HybridTime(1.5,2);
    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(system,initial_enclosure,evolution_time);
    ARIADNE_TEST_PRINT(orbit);
    ARIADNE_TEST_EQUAL(orbit.final().size(),1u);
    ARIADNE_TEST_EQUAL(orbit.reach().size(),3u);

    Box bounding_box(2, -1.0,11.0, -6.0, 6.0);
    plot("test_constraint_hybrid_evolver-flow",bounding_box, Colour(0.75,0.75,0.75),Box(2,-1.0,11.0,-2.0,-1.5),Colour(0.25,0.25,0.5),orbit.reach(),Colour(0.5,0.5,0.75),orbit.final(),Colour(0.5,0.5,0.75),orbit.initial());
}


void
TestContraintHybridEvolver::test_splitting_on_urgent_event() const
{
    evolver_ptr->parameters().maximum_step_size=2.0;
    CompositeHybridAutomaton system=affine_flow_system();

    DiscreteLocation upwards("upwards");
    DiscreteLocation initial_location(upwards);
    const double r=0.125;
    Box initial_box(2,-r,+r, -r,+r);
    HybridEnclosure initial_enclosure=this->evolver_ptr->enclosure(HybridBox(initial_location,initial_box));
    HybridTime evolution_time=HybridTime(4.0,2);
    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(system,initial_enclosure,evolution_time);
    ARIADNE_TEST_PRINT(orbit);
    ARIADNE_TEST_EQUAL(orbit.final().size(),1u);
    ARIADNE_TEST_EQUAL(orbit.reach().size(),2u);

    Box bounding_box(2, -1.0,11.0, -6.0, 6.0);
    plot("test_constraint_hybrid_evolver-urgent",bounding_box, Colour(0.75,0.75,0.75),Box(2,-1.0,11.0,-2.0,-1.5),Colour(0.25,0.25,0.5),orbit.reach(),Colour(0.5,0.5,0.75),orbit.final(),Colour(0.5,0.5,0.75),orbit.initial());
}


void
TestContraintHybridEvolver::test_affine_flow_system() const
{
    evolver_ptr->parameters().maximum_step_size=8.0;

    CompositeHybridAutomaton system=affine_flow_system();
    ARIADNE_TEST_PRINT(system);


    DiscreteLocation upwards("upwards");
    DiscreteLocation initial_location(upwards);
    double r=0.125;
    Box initial_box(2,-r,+r, -r,+r);
    HybridEnclosure initial_enclosure=this->evolver_ptr->enclosure(HybridBox(initial_location,initial_box));
    ARIADNE_TEST_PRINT(initial_enclosure);

    HybridTime evolution_time(0.0,0u);

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
    plot("test_constraint_hybrid_evolver-affine-t=2",Box(2, -1.0,11.0, -6.0, 6.0), Colour(0.75,0.75,0.75),Box(2,-1.0,11.0,-2.0,-1.5),Colour(0.25,0.25,0.5),orbit.reach(),Colour(0.5,0.5,0.75),orbit.final(),Colour(0.5,0.5,0.75),orbit.initial());

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

    plot("test_constraint_hybrid_evolver-affine-t=6",Box(2, -1.0,11.0, -6.0, 6.0), Colour(0.75,0.75,0.75),Box(2,-1.0,11.0,-2.0,-1.5),Colour(0.25,0.25,0.5),orbit.reach(),Colour(0.5,0.5,0.75),orbit.final(),Colour(0.5,0.5,0.75),orbit.initial());

    evolution_time=HybridTime(8.0,3);
    ARIADNE_TEST_PRINT(evolution_time);
    orbit=evolver_ptr->orbit(system,initial_enclosure,evolution_time);
    ARIADNE_TEST_PRINT(orbit);
    ARIADNE_TEST_EQUAL(orbit.reach().size(),3u);
    ARIADNE_TEST_EQUAL(orbit.final().size(),1u);
    plot("test_constraint_hybrid_evolver-affine-t=8",Box(2, -1.0,11.0, -6.0, 6.0), Colour(0.75,0.75,0.75),Box(2,-1.0,11.0,-2.0,-1.5),Colour(0.25,0.25,0.5),orbit.reach(),Colour(0.5,0.5,0.75),orbit.final(),Colour(0.5,0.5,0.75),orbit.initial());
}

*/

/*
class TestHybridEvolution
{
    typedef Vector<Float> FVector;
    typedef Matrix<Float> FMatrix;

    static const bool non_urgent=false;
    static const bool urgent=true;
  private:
    static MonolithicHybridAutomaton system();
  public:
    void test() const;
    void test_constant_derivative_system() const;
    void test_bouncing_ball() const;
    void test_affine_system() const;
};



class TestHybridEvolution
{
    typedef Vector<Float> FVector;
    typedef Matrix<Float> FMatrix;

    static const bool non_urgent=false;
    static const bool urgent=true;
  private:
    static MonolithicHybridAutomaton system();
  public:
    void test() const;
    void test_constant_derivative_system() const;
    void test_bouncing_ball() const;
    void test_affine_system() const;
};

MonolithicHybridAutomaton
TestHybridEvolution::system()
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
    AtomicDiscreteLocation q1(1); AtomicDiscreteLocation q2(2); DiscreteEvent e(1);
    VectorAffineFunction d(FMatrix(2,2, 0.,0.,0.,0.),FVector(2, 1.0,0.));
    VectorAffineFunction r(FMatrix(2,2, 1.,0.,0.,1.),FVector(2, -2.,0.));
    VectorAffineFunction g(FMatrix(1,2, 1.,0.,0.,0.),FVector(1, -1.25));

    MonolithicHybridAutomaton automaton("Constant Derivative System");
    automaton.new_mode(q1,d);
    automaton.new_mode(q2,d);
    automaton.new_transition(e,q1,q2,r,g,urgent);

    TaylorConstrainedImageSet initial_enclosure(Box(2, -0.0625,0.0625, -0.0625,+0.0625));
    HybridTaylorConstrainedImageSet initial_set(q1,initial_enclosure);

    HybridEvolver evolver;
    evolver.verbosity=evolver_verbosity;

    ARIADNE_TEST_PRINT(automaton);
    ARIADNE_TEST_PRINT(initial_set);

    {
        // Test continuous evolution without any jumps
        HybridTime evolution_time(0.5,1);
        ARIADNE_TEST_PRINT(evolution_time);
        Orbit<HybridTaylorConstrainedImageSet> orbit=evolver.orbit(automaton,initial_set,evolution_time);
        ARIADNE_TEST_PRINT(orbit);
        ListSet<HybridTaylorConstrainedImageSet> final_set=evolver.evolve(automaton,initial_set,evolution_time);
        ARIADNE_TEST_PRINT(final_set);
        HybridTaylorConstrainedImageSet expected_final_set(q1,Box(2, +0.4375,+0.5625, -0.0625,+0.0625));
        ARIADNE_TEST_PRINT(expected_final_set);
        ARIADNE_TEST_COMPARE(norm(final_set[q1][0].models()-expected_final_set.second.models()),<,1e-15);
    }

    {
        // Test continuous evolution with a single transverse jump
        HybridTime evolution_time(2.0,2);
        ARIADNE_TEST_PRINT(evolution_time);

        Orbit<HybridTaylorConstrainedImageSet> orbit=evolver.orbit(automaton,initial_set,evolution_time);
        ARIADNE_TEST_PRINT(orbit);

        ListSet<HybridTaylorConstrainedImageSet> final_set=evolver.evolve(automaton,initial_set,evolution_time);
        ARIADNE_TEST_PRINT(final_set);
        HybridTaylorConstrainedImageSet expected_final_set(q2,Box(2, -0.0625,+0.0625, -0.0625,+0.0625));
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
    AtomicDiscreteLocation q1(1);
    AtomicDiscreteLocation q2(2);
    DiscreteEvent e12(12);
    //VectorAffineFunction dynamic(FMatrix(3,3, 0.,1.,0., 0.,0.,0., 0.,0.,0.), FVector(3, 0.0, -g, 1.0));
    //VectorAffineFunction reset(FMatrix(3,3, 1.0,0.0 ,0.0,-a,0.0, 0.0,0.0,1.0), FVector(3, 0.0,0.0,0.0));
    //VectorAffineFunction guard(FMatrix(1,3, -1.0,0.0,0.0), FVector(1, 0.0));
    VectorAffineFunction dynamic(FMatrix(2,2, 0.,1., 0.,0.), FVector(2, 0.0, -g));
    VectorAffineFunction reset(FMatrix(2,2, 1.0,0.0 ,0.0,-a), FVector(2, 0.0,0.0));
    VectorAffineFunction guard(FMatrix(1,2, -1.0,0.0), FVector(1, 0.0));

    /// Build the automaton
    MonolithicHybridAutomaton automaton;
    automaton.new_mode(q1,dynamic);
    automaton.new_mode(q2,dynamic);
    automaton.new_transition(e12,q1,q2,reset,guard,urgent);

    //TaylorConstrainedImageSet initial_enclosure(Box(3, x0-r0,x0+r0, -r0,+r0, 0.0,0.0));
    TaylorConstrainedImageSet initial_enclosure(Box(2, x0-r0,x0+r0, -r0,+r0));
    HybridTaylorConstrainedImageSet initial_set(q1,initial_enclosure);

    HybridEvolver evolver;
    evolver.verbosity=evolver_verbosity;
    evolver.parameters().maximum_step_size=0.125;

    ARIADNE_TEST_PRINT(automaton);
    ARIADNE_TEST_PRINT(initial_set);

    // Test continuous evolution without any jumps
    HybridTime evolution_time(1.5,2);
    ARIADNE_TEST_PRINT(evolution_time);
    Orbit<HybridTaylorConstrainedImageSet> orbit=evolver.orbit(automaton,initial_set,evolution_time);
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

    const AtomicDiscreteLocation location1(1);
    const AtomicDiscreteLocation location2(2);
    const DiscreteEvent event3(3);
    const DiscreteEvent event4(4);

    typedef TaylorConstrainedImageSet EnclosureType;
    typedef HybridBasicSet<TaylorConstrainedImageSet> HybridEnclosureType;

    // Set up the evolution parameters and grid
    Float step_size(0.5);
    Float enclosure_radius(0.25);

    EvolutionParameters parameters;
    parameters.maximum_enclosure_radius=enclosure_radius;
    parameters.maximum_step_size=step_size;

    // Set up the evaluators
    HybridEvolver evolver(parameters);
    evolver.verbosity = evolver_verbosity;


    // Make a hybrid automaton for the Van der Pol equation
    MonolithicHybridAutomaton automaton=system();
    ARIADNE_TEST_PRINT(automaton);

    // Define the initial box
    Box initial_box(2, -0.01,0.01, 0.49,0.51);
    cout << "initial_box=" << initial_box << endl;
    EnclosureType initial_set=TaylorConstrainedImageSet(initial_box);
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
    AtomicDiscreteLocation q1,q2;
    DiscreteEvent e;
    HybridEvolver evolver;
    RealScalarFunction z,o,x,y;
    RealScalarFunction x0,y0,t;
  public:
    TestHybridEvolver();
    MonolithicHybridAutomaton make_hybrid_automaton(const RealScalarFunction& guard);

    void test();
    void test_transverse_linear_crossing();
    void test_transverse_cubic_crossing();
    void test_transverse_cube_root_crossing();

};

TestHybridEvolver::TestHybridEvolver()
    : evolver()
{
    // Set up convenience variables
    q1=AtomicDiscreteLocation(1);
    q2=AtomicDiscreteLocation(2);
    e=DiscreteEvent(3);

    z=RealScalarFunction::constant(2,0.0);
    o=RealScalarFunction::constant(2,1.0);
    x=RealScalarFunction::variable(2,0);
    y=RealScalarFunction::variable(2,1);
    x0=RealScalarFunction::variable(3,0);
    y0=RealScalarFunction::variable(3,1);
    t=RealScalarFunction::variable(3,2);
}

MonolithicHybridAutomaton TestHybridEvolver::make_hybrid_automaton(const RealScalarFunction& guard)
{
    MonolithicHybridAutomaton system;
    system.new_mode(q1,RealVectorFunction(join(o,z)));
    system.new_mode(q2,RealVectorFunction(join(z,o)));
    system.new_transition(e,q1,q2,IdentityFunction(2),RealVectorFunction(1u,guard),true);
    return system;
}

void TestHybridEvolver::test_transverse_linear_crossing()
{
    Float r=1.0/8;
    Float tol=1e-5;
    RealScalarFunction guard=x+y/2-1;
    MonolithicHybridAutomaton system=make_hybrid_automaton(guard);
    Box initial_box(2, -r,+r, -r,+r);
    HybridTaylorConstrainedImageSet initial_set(q1,initial_box);
    HybridTime evolution_time(2.0,3);

    ListSet<HybridTaylorConstrainedImageSet> evolved_set=evolver.evolve(system,initial_set,evolution_time,UPPER_SEMANTICS);

    RealScalarFunction ct=-guard; // Crossing time
    RealVectorFunction f=join(x+ct,y+2-ct);
    Vector<Interval> tolerance(2,Interval(-tol,+tol));
    TaylorConstrainedImageSet expected_evolved_set(f,initial_box);
    ARIADNE_TEST_BINARY_PREDICATE(refines,expected_evolved_set.models(),evolved_set[q2][0].models());
    ARIADNE_TEST_BINARY_PREDICATE(refines,evolved_set[q2][0].models(),expected_evolved_set.models()+tolerance);

    plot("test_hybrid_evolution-transverse_linear_crossing",Box(2, -1.0,3.0, -1.0,3.0),
         Colour(0,0,1),evolved_set[q2][0]);
}

void TestHybridEvolver::test_transverse_cubic_crossing()
{
    Float r=1.0/8;
    Float tol=1e-5;
    RealScalarFunction guard=x-(1+y/2+y*y*y);
    MonolithicHybridAutomaton system=make_hybrid_automaton(guard);
    Box initial_box(2, -r,+r, -r,+r);
    HybridTaylorConstrainedImageSet initial_set(q1,initial_box);
    HybridTime evolution_time(2.0,3);

    ListSet<HybridTaylorConstrainedImageSet> evolved_set=evolver.evolve(system,initial_set,evolution_time,UPPER_SEMANTICS);

    RealScalarFunction ct=-guard; // Crossing time

    RealVectorFunction f=join(x+ct,y+2-ct);
    Vector<Interval> tolerance(2,Interval(-tol,+tol));
    TaylorConstrainedImageSet expected_evolved_set(f,initial_box);
    ARIADNE_TEST_BINARY_PREDICATE(refines,expected_evolved_set.models(),evolved_set[q2][0].models());
    ARIADNE_TEST_BINARY_PREDICATE(refines,evolved_set[q2][0].models(),expected_evolved_set.models()+tolerance);
    plot("test_hybrid_evolution-transverse_cubic_crossing",Box(2, -1.0,3.0, -1.0,3.0),
         Colour(0,0,1),evolved_set[q2][0]);
}

void TestHybridEvolver::test_transverse_cube_root_crossing()
{
    Float r=1.0/32;
    Float tol=1e-5;
    RealScalarFunction guard=((x-1)*(x-1)+1.0)*(x-1)-y-1./64;
    MonolithicHybridAutomaton system=make_hybrid_automaton(guard);
    Box initial_box(2, -r,+r, -r,+r);
    HybridTaylorConstrainedImageSet initial_set(q1,initial_box);
    HybridTime evolution_time(2.0,3);

    RealScalarFunction ct=y-pow(y,3)+3*pow(y,5)-12*pow(y,7)+55*pow(y,9)-273*pow(y,11)+1-x;
    RealVectorFunction f=join(x+ct,y+2-ct);
    Vector<Interval> tolerance(2,Interval(-tol,+tol));
    TaylorConstrainedImageSet expected_evolved_set(f,initial_box);

    ListSet<HybridTaylorConstrainedImageSet> evolved_set=evolver.evolve(system,initial_set,evolution_time,UPPER_SEMANTICS);

    //ARIADNE_TEST_BINARY_PREDICATE(refines,expected_evolved_set.models(),evolved_set[q2][0].models());
    ARIADNE_TEST_BINARY_PREDICATE(refines,evolved_set[q2][0].models(),expected_evolved_set.models()+tolerance);
}

void TestHybridEvolver::test() {
    ARIADNE_TEST_CALL(test_transverse_linear_crossing());
    ARIADNE_TEST_CALL(test_transverse_cubic_crossing());
    test_transverse_cube_root_crossing();
    ARIADNE_TEST_CALL(test_transverse_cube_root_crossing());
}
*/

int main(int argc, const char* argv[])
{
    int evolver_verbosity=get_verbosity(argc,argv);

    DRAWING_METHOD = AFFINE_DRAW; DRAWING_ACCURACY = 2u;

    GeneralHybridEvolver evolver;
    evolver.verbosity=evolver_verbosity;
    //TestSimpleHybridEvolver(evolver,"general_hybrid_evolver").test_unwind();
    //TestSimpleHybridEvolver(evolver,"general_hybrid_evolver").test_transverse_only();
    TestSimpleHybridEvolver(evolver,"general_hybrid_evolver").test_all();

    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

