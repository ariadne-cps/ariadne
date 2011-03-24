/***************************************************************************
 *            watertank-compositional-hysteresis.cc
 *
 *  Copyright  2010  Luca Geretti
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

#include "ariadne.h"

using namespace Ariadne;

int main(int argc,char *argv[])
{
	int analyserVerbosity = 1;
	if (argc > 1)
		analyserVerbosity = atoi(argv[1]);

    /// Set the system parameters
	RealConstant a("a",0.02);
	RealConstant b("b",Interval(0.3,0.32863));
	RealConstant T("T",4.0);
	RealConstant hmin("hmin",Interval(5.0,6.0)); // 5.5;
	RealConstant hmax("hmax",Interval(7.5,8.5)); // 8.0;
	RealConstant delta("Delta",0.1);

    // System variables
    RealVariable x("x");    // water level
    RealVariable y("y");    // valve aperture

    // Create the tank automaton

	    HybridIOAutomaton tank("tank");

		// States
		DiscreteState flow("flow");  

		// Add the input/output variables
    	tank.add_input_var(y);
    	tank.add_output_var(x);
		
		// Only one state with no transitions and no invariants
		RealExpression dyn = - a * x + b * y;
		tank.new_mode(flow);
		tank.set_dynamics(flow, x, dyn);

	// Create the valve automaton

		HybridIOAutomaton valve("valve");

		// States
		DiscreteState idle("idle");
		DiscreteState opening("opening");
		DiscreteState closing("closing");
		
		// The valve has one output var (the valve aperture)
		valve.add_output_var(y);

		// Two input events (open and close) and one internal event
		DiscreteEvent e_open("open");
		valve.add_input_event(e_open);
		DiscreteEvent e_close("close");
		valve.add_input_event(e_close);
		DiscreteEvent e_idle("idle");
		valve.add_internal_event(e_idle);

		// Three states:
		// Idle (valve either fully closed or fully opened)
		RealExpression dynidle = 0.0;
		valve.new_mode(idle);
		//valve.new_invariant(idle, -y);
		//valve.new_invariant(idle, y-1.0);
		valve.set_dynamics(idle, y, dynidle);
		// Opening (valve is opening)
		valve.new_mode(opening);
		//valve.new_invariant(opening, -y);
		//valve.new_invariant(opening, y-1.0);
		RealExpression dynopening = 1.0/T;
		valve.set_dynamics(opening, y, dynopening);
		// Closing (valve is closing)
		valve.new_mode(closing);
		//valve.new_invariant(closing, -y);
		//valve.new_invariant(closing, y-1.0);
		RealExpression dynclosing = -1.0/T;
		valve.set_dynamics(closing, y, dynclosing);

		// Transitions

		// the identity y' = y.
		std::map< RealVariable, RealExpression> reset_y_identity;
		reset_y_identity[y] = y;
		std::map< RealVariable, RealExpression> reset_y_one;
		reset_y_one[y] = 1.0;
		std::map< RealVariable, RealExpression> reset_y_zero;
		reset_y_zero[y] = 0.0;

		// when open is received, go to opening
		valve.new_unforced_transition(e_open, idle, opening, reset_y_identity);
		valve.new_unforced_transition(e_open, opening, opening, reset_y_identity);
		//valve.new_unforced_transition(e_open, closing, opening, res);
		 // when closed is received, go to closing
		valve.new_unforced_transition(e_close, idle, closing, reset_y_identity);
		//valve.new_unforced_transition(e_close, opening, closing, res);
		valve.new_unforced_transition(e_close, closing, closing, reset_y_identity);
		// when the valve is fully opened go from opening to idle
		RealExpression y_geq_one = y - 1.0;
		valve.new_forced_transition(e_idle, opening, idle, reset_y_identity, y_geq_one);
		// when the valve is fully closed go from closing to idle
		RealExpression y_leq_zero = - y;
		valve.new_forced_transition(e_idle, closing, idle, reset_y_identity, y_leq_zero);

	// Create the controller automaton

	    HybridIOAutomaton controller("controller");

		// States
		DiscreteState rising("rising");
		DiscreteState falling("falling");
 
		// The valve has one input var (the water level)
		controller.add_input_var(x);
		// Two output events (open and close)
		controller.add_output_event(e_open); 
		controller.add_output_event(e_close);
		
		// Accessible constants
		controller.register_accessible_constant(delta);
		controller.register_accessible_constant(hmin);
		controller.register_accessible_constant(hmax);

		// Two states:
		// Rising (water level is increasing)
		controller.new_mode(rising);
		 // Falling (water level is decreasing)
		controller.new_mode(falling);

		// Transitions
		// when the water is greater than hmax, send a close command
		RealExpression x_geq_hmax = x - hmax + delta;
		controller.new_unforced_transition(e_close, rising, falling, x_geq_hmax);
		// Add the invariant x < hmax + delta to rising
		RealExpression x_leq_hmax = x - hmax - delta;
		controller.new_invariant(rising, x_leq_hmax);
		
		// when the water is lower than hmin, send a open command
		RealExpression x_leq_hmin = hmin + delta - x;
		controller.new_unforced_transition(e_open, falling, rising, x_leq_hmin);
		// Add the invariant x > hmin - delta to falling
		RealExpression x_geq_hmin = hmin - delta - x;
		controller.new_invariant(falling, x_geq_hmin);

	/// Compose the automata
	HybridIOAutomaton tank_valve = compose("tank,valve",tank,valve,flow,idle);
	HybridIOAutomaton system_io = compose("watertank-io",tank_valve,controller,DiscreteState("flow,idle"),rising);

	cout << system_io << "\n\n";

	/// Create the monolithic automaton
	HybridAutomaton system("watertank-comp-hy");
	RealSpace space;
	make_lpair<HybridAutomaton,RealSpace>(system,space) = make_monolithic_automaton(system_io);

	cout << system << "\n\n";

	// Verification information

	// The initial values
	HybridImageSet initial_set;
	initial_set[DiscreteState("flow,idle,rising")] = Box(2, 6.0,6.0, 1.0,1.0);

	HybridBoxes domain;
	domain[DiscreteState("flow,opening,rising")] = Box(2,4.5,6.5,-0.1,1.1);
	domain[DiscreteState("flow,closing,falling")] = Box(2,7.0,9.0,-0.1,1.1);
	domain[DiscreteState("flow,idle,falling")] = Box(2,5.0,9.0,-0.1,0.1);
	domain[DiscreteState("flow,idle,rising")] = Box(2,5.0,9.0,0.9,1.1);

	// The safe region
	HybridBoxes safe_box = bounding_boxes(system.state_space(),Box(2, 5.25, 8.25, -std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));

	/// Verification

	// Create an evolver and analyser objects, then set their verbosity
	HybridEvolver evolver;
	evolver.verbosity = 0;
	HybridReachabilityAnalyser analyser(evolver);
	analyser.verbosity = analyserVerbosity;
	analyser.parameters().enable_lower_pruning = true;
	analyser.parameters().lowest_maximum_grid_depth = 0;
	analyser.parameters().highest_maximum_grid_depth = 6;
	Verifier verifier(analyser);

	RealConstantSet parameters;
	parameters.insert(RealConstant("hmin",Interval(5.0,6.0)));
	parameters.insert(RealConstant("hmax",Interval(7.5,8.5)));

	SystemVerificationInfo verInfo(system, initial_set, domain, safe_box);

	std::list<ParametricOutcome> results = verifier.parametric_safety(verInfo, parameters);
	draw(system.name(),results);
}
