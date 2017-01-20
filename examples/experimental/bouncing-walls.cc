/*****************************************************************************************************
 *            bouncing-walls.cc
 *
 *  Copyright  2016  Luca Geretti
 *
 * A simple automaton where a point bounces between two extremes.
 *
 *****************************************************************************************************/

#include "ariadne.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    /// Constants
    int VERBOSITY = 1;
	if (argc > 1)
		VERBOSITY = atoi(argv[1]);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
	HybridAutomaton automaton("bouncing-walls");

    // Parameters
    RealConstant velocity("velocity",1.0_dec); // Velocity in modulus
    RealConstant width("width",1.0_dec); // Width of the cut

    /// Modes

    DiscreteLocation scanning("scanning");

    // Variables

    RealVariable x("x"); // X position
    RealVariable vx("vx"); // X velocity

    DiscreteEvent switch_left("switch_left");
    DiscreteEvent switch_right("switch_right");

	automaton.new_mode(scanning, {dot(x)=vx,dot(vx)=0});

	automaton.new_guard(scanning,switch_right,x<=0,impact);
	automaton.new_update(scanning,switch_right,scanning,{next(x)=0,next(vx)=velocity});
	automaton.new_transition(scanning,switch_left,scanning,{next(x)=width,next(vx)=-velocity},x>=width);
	//automaton.new_transition(scanning,switch_right,scanning,{next(x)=0,next(vx)=velocity},x<=0);

    // Create a GeneralHybridEvolver object
    GeneralHybridEvolver evolver(automaton);
    evolver.verbosity = VERBOSITY;

    // Set the evolution parameters
    evolver.configuration().set_maximum_enclosure_radius(0.5);
    evolver.configuration().set_maximum_step_size(0.01);
    std::cout << evolver.configuration() << std::endl;

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolver::EnclosureType HybridEnclosureType;
    typedef GeneralHybridEvolver::OrbitType OrbitType;
    typedef GeneralHybridEvolver::EnclosureListType EnclosureListType;

    HybridSet initial_set(scanning,{vx==-velocity,x==Real(0.5)});

    HybridTime evolution_time(10.0,10);

    std::cout << "Computing orbit... " << std::flush;

    OrbitType orbit = evolver.orbit(initial_set,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;
}

