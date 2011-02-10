/***************************************************************************
 *            coconut-final-demo.cc
 *
 *  Copyright  2010  Luca Geretti, Davide Bresolin
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

template<class SET> void plot(const char* filename, const int& xaxis, const int& yaxis, const int& numVariables, const Box& bbox, const SET& set) {

    Figure fig; 
    array<uint> xy(2,xaxis,yaxis);
    fig.set_projection_map(ProjectionFunction(xy,numVariables)); 
    fig.set_bounding_box(bbox); 
    fig.set_fill_colour(Colour(0.0,0.5,1.0));
    draw(fig,set); 
    fig.write(filename);

} 

int main(int argc, char** argv) 
{
    if(argc != 2) {
        std::cout << "USAGE: sast-run filename.cif" << std::endl;
        return 1;
    }

    std::cout << "Importing CIF description from " << argv[1] << "..." ;

    /// System constants
	RealConstant C1("C1",-0.4);			// The constant variation rate of the voltage difference
	RealConstant C2("C2",0.05);		// The step upwards/downwards that the voltage difference should perform
	RealConstant T1("T1",0.01);			// Time required for the step transition of the voltage difference
	RealConstant T("T",0.01);			// Delay for the emission of output events of the controller
	RealConstant H("H",0.001);			// Threshold for evaluating the positive-negative and negative-positive transitions

    // System variables
	RealVariable k1("k1");	// Clock for the voltage difference transitions
	RealVariable t_out("t_out"); // Clock for the controller output events
	RealVariable vd("vd");	// Voltage difference

	// System events
	DiscreteEvent up("UP");		// The core voltage must go up
	DiscreteEvent down("DOWN");	// The core voltage must go down
	DiscreteEvent n2p("N2P");	// The voltage difference goes from negative to positive
	DiscreteEvent p2n("P2N");	// The voltage difference goes from positive to negative

    // Create the voltage difference evaluator automaton

	    HybridIOAutomaton vdiffeval("vdiffeval");

		// States
		DiscreteState pos("Pos"); // Positive difference
		DiscreteState neg("Neg"); // Negative difference

		// Involved variables
    	vdiffeval.add_input_var(vd);

		// Involved events
		vdiffeval.add_output_event(p2n);
		vdiffeval.add_output_event(n2p);
		
		// Modes

			// Pos		
			vdiffeval.new_mode(pos);
			// Neg
			vdiffeval.new_mode(neg);

		// Transitions

			// From Pos to Neg if vd <= -H
			vdiffeval.new_forced_transition(p2n, pos, neg, -vd-H);
			// From Neg to Pos if vd >= H
			vdiffeval.new_forced_transition(n2p, neg, pos, vd-H);

	// Create the environment automaton

		HybridIOAutomaton environment("environment");

		// States
		DiscreteState idiff("Idiff"); // The voltage difference is idle
		DiscreteState rdiff("Rdiff"); // The voltage difference is rising
		DiscreteState fdiff("Fdiff"); // The voltage difference is falling
		
		// Internal events
		DiscreteEvent rdiff2idiff("Rdiff2Idiff");
		DiscreteEvent fdiff2idiff("Fdiff2Idiff");

		// Involved variables
		environment.add_output_var(vd);
		environment.add_internal_var(k1);

		// Involved events
		environment.add_internal_event(rdiff2idiff);
		environment.add_internal_event(fdiff2idiff);
		environment.add_input_event(up);
		environment.add_input_event(down);

		// Dynamics
		RealExpression d_vd_idle = C1; 			// The voltage difference is idle
		RealExpression d_vd_rising = C1+C2/T1;		// The voltage difference is rising
		RealExpression d_vd_falling = C1-C2/T1;	// The voltage difference is falling
		RealExpression d_k1_idle = 0.0;				// The clock is idle
		RealExpression d_k1_running = 1.0;			// The clock is running
		
		// Resets
		std::map< RealVariable, RealExpression> reset_k1_identity; // Keep both
		reset_k1_identity[vd] = vd;
		reset_k1_identity[k1] = k1;
		std::map< RealVariable, RealExpression> reset_k1_zero; // Reset k1, keep vc
		reset_k1_zero[vd] = vd;
		reset_k1_zero[k1] = 0.0;

		// Modes

			// Idiff
			environment.new_mode(idiff);
			environment.new_invariant(idiff, k1 - T1);
			environment.new_invariant(idiff, -k1);
			environment.set_dynamics(idiff, vd, d_vd_idle);
			environment.set_dynamics(idiff, k1, d_k1_idle);
			// Rdiff
			environment.new_mode(rdiff);
			environment.new_invariant(rdiff, k1 - T1);
			environment.new_invariant(rdiff, -k1);
			environment.set_dynamics(rdiff, vd, d_vd_rising);
			environment.set_dynamics(rdiff, k1, d_k1_running);
			// Fdiff
			environment.new_mode(fdiff);
			environment.new_invariant(fdiff, k1 - T1);
			environment.new_invariant(fdiff, -k1);
			environment.set_dynamics(fdiff, vd, d_vd_falling);
			environment.set_dynamics(fdiff, k1, d_k1_running);

		// Transitions

			// From Idiff to Rdiff as soon as the Up event is activated
			environment.new_unforced_transition(up, idiff, rdiff, reset_k1_zero);
			// From Idiff to Fdiff as soon as the Down event is activated
			environment.new_unforced_transition(down, idiff, fdiff, reset_k1_zero);

			// From Rdiff to Idiff as soon as k1 >= T1			
			environment.new_forced_transition(rdiff2idiff, rdiff, idiff, reset_k1_identity, k1-T1);
			// From Rdiff to Rdiff as soon as the Up event is activated
			environment.new_unforced_transition(up, rdiff, rdiff, reset_k1_zero);
			// From Rdiff to Fdiff as soon as the Down event is activated
			environment.new_unforced_transition(down, rdiff, fdiff, reset_k1_zero);

			// From Fdiff to Idiff as soon as k1 >= T1
			environment.new_forced_transition(fdiff2idiff, fdiff, idiff, reset_k1_identity, k1-T1);
			// From Fdiff to Fdiff as soon as the Down event is activated
			environment.new_unforced_transition(down, fdiff, fdiff, reset_k1_zero);
			// From Fdiff to Rdiff as soon as the Up event is activated
			environment.new_unforced_transition(up, fdiff, rdiff, reset_k1_zero);

	// Create the controller automaton

		HybridIOAutomaton controller("controller");

		// States
		DiscreteState decr("Decr");			// Decrease the vd level
		DiscreteState incr("Incr");			// Increase the vd level

		// Involved variables
		controller.add_internal_var(t_out);

		// Involved events
		controller.add_input_event(p2n);
		controller.add_input_event(n2p);
		controller.add_output_event(up);
		controller.add_output_event(down);

		// Resets

			// Resets t_out to zero
			std::map< RealVariable, RealExpression> reset_t_out_zero;
			reset_t_out_zero[t_out] = 0.0;
			// Keeps t_out
			std::map< RealVariable, RealExpression> reset_t_out_identity;
			reset_t_out_identity[t_out] = t_out;

		// Modes

			// Decr
			controller.new_mode(decr);
			controller.new_invariant(decr, t_out - T);
			controller.new_invariant(decr, -t_out);
			controller.set_dynamics(decr, t_out, 1.0);
			// Incr
			controller.new_mode(incr);
			controller.new_invariant(incr, t_out - T);
			controller.new_invariant(incr, -t_out);
			controller.set_dynamics(incr, t_out, 1.0);
	
		// Transitions

			// From Decr to Decr as soon as t_out >= T
			controller.new_forced_transition(down, decr, decr, reset_t_out_zero, t_out - T);
			// From Decr to Incr as soon as the p2n event is activated
			controller.new_unforced_transition(p2n, decr, incr, reset_t_out_identity);
			// From Decr to Decr as soon as the n2p event is activated
			controller.new_unforced_transition(n2p, decr, incr, reset_t_out_identity);
			
			// From Incr to Incr as soon as t_out >= T
			controller.new_forced_transition(up, incr, incr, reset_t_out_zero, t_out - T);
			// From Incr to Decr as soon as the n2p event is activated
			controller.new_unforced_transition(n2p, incr, decr, reset_t_out_identity);
			// From Incr to Incr as soon as the p2n event is activated
			controller.new_unforced_transition(p2n, incr, decr, reset_t_out_identity);

    std::cout << " done." << std::endl;
    
    std::cout << "Making the composition of the network of hybrid automata...";

	/// Compose the automata
	HybridIOAutomaton vdiffeval_environment = compose("vdiffeval,environment",vdiffeval,environment,pos,idiff);

	HybridIOAutomaton system_io = compose("pss-revised",vdiffeval_environment,controller,DiscreteState("Pos,Idiff"),decr);

	/// Create the monolithic automaton
	HybridAutomaton system;
	RealSpace space;
	make_lpair<HybridAutomaton,RealSpace>(system,space) = make_monolithic_automaton(system_io);
	
	std::cout << " done." << std::endl;

	std::cout << std::endl << "Safety property to verify: Always(-2.0*C2 <= vd <= 2.0*C2)." << std::endl;

    // std::cout << std::endl << "Starting verification loop." << std::endl;
    
	// Verification information

	// Variables order: k1, t_out vd

	// The initial values
	HybridImageSet initial_set;
	initial_set[DiscreteState("Neg,Idiff,Incr")] = Box(3,  0.0,0.0, // k1
							                               0.0,0.0, // t_out
	   						                             -0.01,-0.01 // vd
							                          );

	// The safe region
	HybridBoxes safe_box = bounding_boxes(system.state_space(),Box(3,
				   -std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
				   -std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
				   -0.1, 0.1
			   ));

	// The domain
	HybridBoxes domain = bounding_boxes(system.state_space(),Box(3,
								 -0.01,0.11,
								 -0.01,0.11,
								-0.2,0.2));



	// Create an evolver and analyser objects, then set their verbosity
	HybridEvolver evolver;
	evolver.verbosity = 0;		
	HybridReachabilityAnalyser analyser(evolver);
	analyser.verbosity = 3;
	evolver.parameters().enable_subdivisions = false;
	evolver.parameters().enable_set_model_reduction = false;
	analyser.parameters().enable_lower_pruning = true;
	analyser.parameters().lowest_maximum_grid_depth = 5;
	analyser.parameters().highest_maximum_grid_depth = 10;
	analyser.parameters().transient_time = 1e10;
	analyser.parameters().transient_steps = 8;
	analyser.parameters().lock_to_grid_time = 1e10;		
	analyser.parameters().lock_to_grid_steps = 1;
	analyser.plot_verify_results = false;
	analyser.free_cores = 0;

	// The resulting safe and unsafe intervals
	Interval safe_int, unsafe_int;

	// Perform the analysis
	SystemVerificationInfo verInfo(system,initial_set,domain,safe_box);
	tribool result = analyser.verify_iterative(verInfo);
	
	std::cout << "Done." << std::endl;
	
	std::cout << std::endl << "Plotting results...";

    HybridGridTreeSet reach = analyser.statistics().upper().reach;
    HybridGridTreeSet limits(reach.grid());
    limits.adjoin_outer_approximation(bounding_boxes(system.state_space(),Box(3,
								 -0.01,0.11,
								 -0.0,0.11,
								-0.2,0.2)), 7);
    reach.restrict(limits);
    
    Figure fig; 
    array<uint> xy(2,1,2);
    fig.set_projection_map(ProjectionFunction(xy,3)); 
    fig.set_x_axis_label("t_out");
    fig.set_y_axis_label("vd");
    fig.set_line_colour(Colour(0,0,0));
    fig.set_bounding_box(Box(2, -0.001,0.015, -0.15,0.15)); 
    fig.set_fill_colour(Colour(0.95,0.3,0.3));
    draw(fig, Box(3, 0.0,0.0, -0.1,0.1, -0.3, -0.1));
    draw(fig, Box(3, 0.0,0.0, -0.1,0.1, 0.1, 0.3));
    fig.set_fill_colour(Colour(0.0,0.5,1.0));
    draw(fig,reach); 
    fig.write("automaton_reach.png");	
	
	std::cout << " done." << std::endl;
	
	std::cout << std::endl << "Safety property Always(-2.0*C2 <= vd <= 2.0*C2) is " << result << "."  << std::endl;
	
	if(definitely(result))
	    return 0;

    return 1;
}
