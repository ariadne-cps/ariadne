/***************************************************************************
 *            test_hybrid_io_automaton.cc
 *
 *  Copyright  2010 Davide Bresolin
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

#include "numeric.h"
#include "space.h"
#include "discrete_state.h"
#include "expression.h"
#include "assignment.h"
#include "hybrid_io_automaton.h"
#include "hybrid_automaton.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

class TestHybridIOAutomaton
{
  public:
    HybridIOAutomaton tank,
                      valve,
                      controller,
                      system;
    RealVariable x,     // water level
                 y;     // valve aperture

    TestHybridIOAutomaton();
                      
    void test();
    void test_tank_definition();
    void test_valve_definition();
    void test_controller_definition();
    void test_composition();
    void test_conversion();
};

TestHybridIOAutomaton::TestHybridIOAutomaton()
    :   x("x"), y("y")
{
}


void TestHybridIOAutomaton::test_tank_definition() 
{
    tank = HybridIOAutomaton("watertank");
    
    // The water tank has one input var (the valve aperture) and one output var (the water level)
    tank.add_input_var(y);
    ARIADNE_TEST_ASSERT(tank.has_input_var(y));
    tank.add_output_var(x);
    ARIADNE_TEST_ASSERT(tank.has_output_var(x));
    // x cannot be both an output and internal var in tank
    ARIADNE_TEST_FAIL(tank.add_internal_var(x));
    
    // Only one state with no transitions and no invariants
    Float a = 0.02;
    Float b = 0.3;
    RealExpression dyn = - a * x + b * y;
    DiscreteState flow("flow");
    tank.new_mode(flow);
    ARIADNE_TEST_ASSERT(tank.has_mode(flow));
    tank.set_dynamics(flow, x, dyn);
    // setting the dynamic of an input var should cause an exception
    ARIADNE_TEST_FAIL(tank.set_dynamics(flow, y, dyn));
    
    ARIADNE_TEST_PRINT(tank);
    
}

void TestHybridIOAutomaton::test_valve_definition() 
{
    valve = HybridIOAutomaton("valve");
    
    // The valve has one output var (the valve aperture)
    valve.add_output_var(y);
    ARIADNE_TEST_ASSERT(valve.has_output_var(y));
    // Two input events (open and close) and one internal event
    DiscreteEvent e_open("open");
    valve.add_input_event(e_open);
    ARIADNE_TEST_ASSERT(valve.has_input_event(e_open));
    DiscreteEvent e_close("close");
    valve.add_input_event(e_close);
    ARIADNE_TEST_ASSERT(valve.has_input_event(e_close));
    DiscreteEvent e_idle("idle");
    valve.add_internal_event(e_idle);
    ARIADNE_TEST_ASSERT(valve.has_internal_event(e_idle));
    // setting open as an output event should cause an exception
    ARIADNE_TEST_FAIL(valve.add_output_event(e_open));
    
    // Three states:
    // Idle (valve either fully close or fully open)
    Float T = 4.0;
    RealExpression dynidle = 0.0;
    DiscreteState idle("idle");
    valve.new_mode(idle);
    ARIADNE_TEST_ASSERT(valve.has_mode(idle));
    valve.set_dynamics(idle, y, dynidle);
    // setting the dynamic of an unknown variable cause an exception
    ARIADNE_TEST_FAIL(valve.set_dynamics(idle, x, dynidle));
 
    // Opening (valve is opening)
    DiscreteState opening("opening");
    valve.new_mode(opening);
    ARIADNE_TEST_ASSERT(valve.has_mode(opening));
    RealExpression dynopening = 1.0/T;
    valve.set_dynamics(opening, y, dynopening);

    // Closing (valve is closing)
    DiscreteState closing("closing");
    valve.new_mode(closing);
    ARIADNE_TEST_ASSERT(valve.has_mode(closing));
    RealExpression dynclosing = -1.0/T;
    valve.set_dynamics(closing, y, dynclosing);
    
    // Transitions
    //
    // the reset function is always the identity y' = y.
    std::map< RealVariable, RealExpression> res;
    res[y] = y;
    // when open is received, go to opening
    valve.new_transition(e_open, idle, opening, false);
    ARIADNE_TEST_ASSERT(valve.has_transition(e_open,idle));
    ARIADNE_TEST_TRY(valve.set_reset(e_open, idle, y, y));
    valve.new_transition(e_open, opening, opening, res, false);
    ARIADNE_TEST_ASSERT(valve.has_transition(e_open,opening));
    valve.new_transition(e_open, closing, opening, res, false);
    ARIADNE_TEST_ASSERT(valve.has_transition(e_open,closing));
     // when closed is received, go to closing
    valve.new_unforced_transition(e_close, idle, closing, res);
    ARIADNE_TEST_ASSERT(valve.has_transition(e_close,idle));
    valve.new_unforced_transition(e_close, opening, closing, res);
    ARIADNE_TEST_ASSERT(valve.has_transition(e_close,opening));
    valve.new_unforced_transition(e_close, closing, closing, res);
    ARIADNE_TEST_ASSERT(valve.has_transition(e_close,closing));
    // when the valve is fully open go from opening to idle
    RealExpression y_geq_one = y - 1.0;
    valve.new_forced_transition(e_idle, opening, idle, res, y_geq_one);
    ARIADNE_TEST_ASSERT(valve.has_transition(e_idle,opening));
    // when the valve is fully close go from closing to idle
    RealExpression y_leq_zero = - y;
    valve.new_forced_transition(e_idle, closing, idle, res, y_leq_zero);
    ARIADNE_TEST_ASSERT(valve.has_transition(e_idle,closing));
    // creating a forced transition with an input event should cause an exception
    ARIADNE_TEST_FAIL(valve.new_forced_transition(e_open, idle, idle, y_geq_one));
    // creating a transition with an input event and a custom activation should cause an exception
    ARIADNE_TEST_FAIL(valve.new_unforced_transition(e_open, idle, idle, y_geq_one));
    
    ARIADNE_TEST_PRINT(valve);
   
}

void TestHybridIOAutomaton::test_controller_definition() 
{
    controller = HybridIOAutomaton("controller");
 
    // The valve has one input var (the water level)
    controller.add_input_var(x);
    ARIADNE_TEST_ASSERT(controller.has_input_var(x));
    // Two output events (open and close)
    DiscreteEvent e_open("open");
    controller.add_output_event(e_open);
    ARIADNE_TEST_ASSERT(controller.has_output_event(e_open));
    DiscreteEvent e_close("close");
    controller.add_output_event(e_close);
    ARIADNE_TEST_ASSERT(controller.has_output_event(e_close));
    // setting open as an input or internal event should cause an exception
    ARIADNE_TEST_FAIL(controller.add_input_event(e_open));
    ARIADNE_TEST_FAIL(controller.add_internal_event(e_open));
    
    // Two states:
    // Rising (water level is increasing)
    DiscreteState rising("rising");
    controller.new_mode(rising);
    ARIADNE_TEST_ASSERT(controller.has_mode(rising));
    // setting the dynamic of an input variable cause an exception
    ARIADNE_TEST_FAIL(controller.set_dynamics(rising, x, RealExpression(2.0*x)));
 
     // Falling (water level is decreasing)
    DiscreteState falling("falling");
    controller.new_mode(falling);
    ARIADNE_TEST_ASSERT(controller.has_mode(falling));

    // Transitions
    // when the water is greater than hmax, send a close command
    Float hmax = 8.0;
    Float delta = 0.05;
    RealExpression x_geq_hmax = x - hmax + delta;
    controller.new_unforced_transition(e_close, rising, falling, x_geq_hmax);
    ARIADNE_TEST_ASSERT(controller.has_transition(e_close,rising));
    // Add the invariant x < hmax + delta to rising
    RealExpression x_leq_hmax = x - hmax - delta;
    controller.new_invariant(rising, x_leq_hmax);
    
    // when the water is lower than hmin, send a open command
    Float hmin = 5.5;
    RealExpression x_leq_hmin = hmin + delta - x;
    controller.new_unforced_transition(e_open, falling, rising, x_leq_hmin);
    ARIADNE_TEST_ASSERT(controller.has_transition(e_close,rising));
    // Add the invariant x > hmin - delta to falling
    RealExpression x_geq_hmin = hmin - delta - x;
    controller.new_invariant(falling, x_geq_hmin);
 
    ARIADNE_TEST_PRINT(controller);
   
}

void TestHybridIOAutomaton::test_composition() 
{
    HybridIOAutomaton valvetank;
    ARIADNE_TEST_TRY(valvetank = compose("valvecontr", valve, tank, 
                                         DiscreteState("idle"), DiscreteState("flow")));
    ARIADNE_TEST_PRINT(valvetank);
    ARIADNE_TEST_TRY(system = compose("system", valvetank, controller,
                                      DiscreteState("idle,flow"), DiscreteState("rising")));
    ARIADNE_TEST_PRINT(system);
                     
}

void TestHybridIOAutomaton::test_conversion() 
{
    HybridAutomaton ha;
    RealSpace spc;
    
    ARIADNE_TEST_FAIL(make_monolithic_automaton(tank));
    ARIADNE_TEST_FAIL(make_monolithic_automaton(valve));
    ARIADNE_TEST_FAIL(make_monolithic_automaton(controller));
    
    ARIADNE_TEST_TRY(make_lpair(ha,spc) = make_monolithic_automaton(system));
    
    ARIADNE_TEST_PRINT(ha);
    ARIADNE_TEST_PRINT(spc);

}



void TestHybridIOAutomaton::test() 
{
    ARIADNE_TEST_CALL(test_tank_definition());
    ARIADNE_TEST_CALL(test_valve_definition());
    ARIADNE_TEST_CALL(test_controller_definition());
    ARIADNE_TEST_CALL(test_composition());
    ARIADNE_TEST_CALL(test_conversion());
}


int main()
{
    TestHybridIOAutomaton().test();
    return ARIADNE_TEST_FAILURES;
}

