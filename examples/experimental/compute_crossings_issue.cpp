#include "ariadne.hpp"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    /// Constants
    int VERBOSITY = 1;
	if (argc > 1)
		VERBOSITY = atoi(argv[1]);

    /// Build the Hybrid System

	RealConstant velocity("velocity",0.092_dec);
	RealConstant L("L",0.00025_dec);
	RealConstant x0("x0",0.0023_dec);

    /// Create a HybridAutomaton object
	AtomicHybridAutomaton automaton("compute_crossings_issue");

    /// Modes

    AtomicDiscreteLocation far("far");
    AtomicDiscreteLocation close("close");

    // Variables

    RealVariable x("x"); // X position

	automaton.new_mode(far, {dot(x)=velocity});
	automaton.new_mode(close, {dot(x)=velocity});

    DiscreteEvent comes("comes");

	RealExpression sqr_distance = sqr(x-x0);

	automaton.new_transition(far,comes,close,{next(x)=x},sqr_distance<=L*L);

    // Create a GeneralHybridEvolver object
    GeneralHybridEvolver evolver(automaton);
    evolver.verbosity = VERBOSITY;

    // Set the evolution parameters
    evolver.configuration().set_maximum_enclosure_radius(0.5);
    evolver.configuration().set_maximum_step_size(0.0001);

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolver::EnclosureType HybridEnclosureType;
    typedef GeneralHybridEvolver::OrbitType OrbitType;
    typedef GeneralHybridEvolver::EnclosureListType EnclosureListType;

    HybridSet initial_set({automaton|far},{x==0});

    HybridTime evolution_time(0.5,2);

    std::cout << "Computing orbit... " << std::flush;

    OrbitType orbit = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    std::cout << "done." << std::endl;

}
