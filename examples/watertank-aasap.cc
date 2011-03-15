/***************************************************************************
 *            watertank-aasap.cc
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
	RealConstant b("b",0.3);
	RealConstant Ta("Ta",4.0);
	RealConstant T("T",0.1);
	RealConstant hmin("hmin",5.5); 
	RealConstant hmax("hmax",8.0);

    // System variables
    RealVariable x("x");    // water level
    RealVariable y("y");    // valve aperture
    RealVariable t_out("t_out");    // valve aperture

	// The parameter to modify, its interval and the tolerance
    RealConstant Delta("Delta",Interval(0.0,0.0));
    RealConstant parameter = Delta;
	Float tolerance = 1e-2;

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
		DiscreteEvent e_open("OPEN");
		valve.add_input_event(e_open);
		DiscreteEvent e_close("CLOSE");
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
		RealExpression dynopening = 1.0/Ta;
		valve.set_dynamics(opening, y, dynopening);
		// Closing (valve is closing)
		valve.new_mode(closing);
		//valve.new_invariant(closing, -y);
		//valve.new_invariant(closing, y-1.0);
		RealExpression dynclosing = -1.0/Ta;
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
		valve.new_unforced_transition(e_open, closing, opening, reset_y_identity);
		 // when closed is received, go to closing
		valve.new_unforced_transition(e_close, idle, closing, reset_y_identity);
		valve.new_unforced_transition(e_close, opening, closing, reset_y_identity);
		valve.new_unforced_transition(e_close, closing, closing, reset_y_identity);
		// when the valve is fully opened go from opening to idle
		RealExpression y_geq_one = y - 1.0;
		valve.new_forced_transition(e_idle, opening, idle, reset_y_identity, y_geq_one);
		// when the valve is fully closed go from closing to idle
		RealExpression y_leq_zero = - y;
		valve.new_forced_transition(e_idle, closing, idle, reset_y_identity, y_leq_zero);

	// Create the evaluator automaton

		    HybridIOAutomaton evaluator("evaluator");

			// States
			DiscreteState deep("deep");
			DiscreteState shallow("shallow");

			// The evaluator has two output events
			DiscreteEvent e_high("HIGH");
			evaluator.add_output_event(e_high);
			DiscreteEvent e_low("LOW");
			evaluator.add_output_event(e_low);

			// An high event has been issued
			evaluator.new_mode(deep);
			// A low event has been issued;
			evaluator.new_mode(shallow);

			// Transitions
			// When the water is greater than hmax, send a high event
			RealExpression x_geq_hmax = x - hmax;
			evaluator.new_forced_transition(e_high, shallow, deep, x_geq_hmax);
			// When the water is lower than hmin, send a open command
			RealExpression x_leq_hmin = hmin - x;
			evaluator.new_forced_transition(e_low, deep, shallow, x_leq_hmin);


	// Create the controller automaton

	    HybridIOAutomaton controller("controller");

		// States
		DiscreteState increase("increase");
		DiscreteState decrease("decrease");
		DiscreteState nothing("nothing");

		// Involved variables
		controller.add_internal_var(t_out);
 
		// Two input events (high and low)
		controller.add_input_event(e_high);
		controller.add_input_event(e_low);
		// Two output events (open and close)
		controller.add_output_event(e_open); 
		controller.add_output_event(e_close);
		
		// Two states:
		// The controller is about to increase the level by opening the valve
		controller.new_mode(increase);
		controller.set_dynamics(increase, t_out, 1.0);
		// The controller is about to decrease the level by closing the valve
		controller.new_mode(decrease);
		controller.set_dynamics(decrease, t_out, 1.0);
		// The controller does nothing
		controller.new_mode(nothing);
		controller.set_dynamics(nothing, t_out, 0.0);

		// Resets

			// Reset t_out to zero
			std::map< RealVariable, RealExpression> reset_t_out_zero;
			reset_t_out_zero[t_out] = 0.0;

		// Transitions

		// When the water level is high, wait before issuing the close command
		controller.new_unforced_transition(e_high, nothing, decrease, reset_t_out_zero);
		// When the water level is low, wait before issuing the open command
		controller.new_unforced_transition(e_low, nothing, increase, reset_t_out_zero);
		// When the clock period has expired, issue the open command
		controller.new_forced_transition(e_open, increase, nothing, reset_t_out_zero, t_out-T);
		// When the clock period has expired, issue the close command
		controller.new_forced_transition(e_close, decrease, nothing, reset_t_out_zero, t_out-T);

	/// Compose the automata
	HybridIOAutomaton tank_valve = compose("tank,valve",tank,valve,flow,idle);
	HybridIOAutomaton tank_valve_evaluator = compose("tank,valve,evaluator",tank_valve,evaluator,DiscreteState("flow,idle"),shallow);
	HybridIOAutomaton system_io = compose("watertank-aasap",tank_valve_evaluator,controller,DiscreteState("flow,idle,shallow"),nothing);

	/// Create the monolithic automaton
	HybridAutomaton system;
	RealSpace space;
	make_lpair<HybridAutomaton,RealSpace>(system,space) = make_monolithic_automaton(system_io);

	/// Add access to the Delta constant
	//system.register_accessible_constant(Delta);

	// Verification information

	// The initial values
	HybridImageSet initial_set;
	initial_set[DiscreteState("flow,idle,shallow,nothing")] = Box(3, 0.0,0.0, 6.0,6.0, 1.0,1.0);

	// The safe region
	HybridBoxes safe_box = bounding_boxes(system.state_space(),Box(3, -std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
																      5.25, 8.25,
																	  -std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));


	HybridBoxes domain = bounding_boxes(system.state_space(),Box(3,-0.1,0.2,4.0,10.0,-0.5,1.5));
/*

	// The domain
	HybridBoxes domain;

	// The level is rising
	domain[DiscreteState("flow,idle,shallow,nothing")] = Box(3, -0.1,0.2, 4.5,8.5, 0.9,1.1);
	// A CLOSE command will be issued
	domain[DiscreteState("flow,idle,deep,decrease")] = Box(3, -0.1,0.2, 7.5,9.0, 0.1,1.1);
	// The valve is closing
	domain[DiscreteState("flow,closing,deep,nothing")] = Box(3, -0.1,0.2, 7.5,9.5, -0.1,1.1);
	// The level is falling
	domain[DiscreteState("flow,idle,deep,nothing")] = Box(3, -0.1,0.2, 5.0,8.5, -0.1,0.1);
	// An OPEN command will be issued
	domain[DiscreteState("flow,idle,shallow,increase")] = Box(3, -0.1,0.2, 4.5,6.0, -0.1,0.1);
	// The valve is opening
	domain[DiscreteState("flow,opening,shallow,nothing")] = Box(3, -0.1,0.2, 4.0,6.5, -0.1,1.1);

	// Spurious: the valve is closing but an OPEN command will be issued
	domain[DiscreteState("flow,closing,shallow,increase")] = Box(3, -0.1,0.2, 7.5,9.0, -0.1,1.1);
	// Spurious: the valve is opening but a CLOSE command will be issued
	domain[DiscreteState("flow,opening,deep,decrease")] = Box(3, -0.1,0.2, 4.5,6.0, -0.1,1.1);


*/
	/// Verification

	// Create an evolver and analyser objects, then set their verbosity
	HybridEvolver evolver;
	evolver.verbosity = 0;
	HybridReachabilityAnalyser analyser(evolver);
	analyser.verbosity = analyserVerbosity;
	evolver.parameters().enable_subdivisions = false;
	evolver.parameters().enable_set_model_reduction = false;
	analyser.parameters().enable_lower_pruning = true;
	analyser.parameters().lowest_maximum_grid_depth = 0;
	analyser.parameters().highest_maximum_grid_depth = 8;
	analyser.parameters().transient_time = 1e10;
	analyser.parameters().transient_steps = 1;
	analyser.parameters().lock_to_grid_time = 1e10;		
	analyser.parameters().lock_to_grid_steps = 1;
	analyser.plot_verify_results = false;
	analyser.free_cores = 0;
	Verifier verifier(analyser);

	// The resulting safe and unsafe intervals
	Interval safe_int, unsafe_int;

	SystemVerificationInfo verInfo(system, initial_set, domain, safe_box);
    verifier.verify_iterative(verInfo);

    
    HybridGridTreeSet reach = analyser.statistics().upper().reach;

	// Assigns local variables
	Figure fig;
	array<uint> xy(2,1,2);

	fig.set_projection_map(ProjectionFunction(xy,3));
	fig.set_bounding_box(Box(2,4.0,9.5,-0.05,1.05));

	// Appends the set, with the desired fill color
	fig.set_fill_colour(Colour(0.9,0.9,0.9));
	draw(fig,reach);
	fig.set_fill_colour(Colour(0.6,0.6,0.6));
	draw(fig,Box(3,0.0,0.0,3.5,5.25,-0.2,1.2));
	draw(fig,Box(3,0.0,0.0,8.25,10.0,-0.2,1.2));
	fig.write("watertank-aasap.png");


}
