/*****************************************************************************************************
 *            circle-unstable.cc
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the behavior of a circle.
 *
 *****************************************************************************************************/

#include "ariadne.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    double A = 2.0;

    /// Constants
    double EVOL_TIME = 100.0/A;   /// Evolution time
    float MAX_ENCL_WIDTH = 0.1;   /// Maximum enclosure width
    float MAX_STEP_SIZE = 1e-2 / A;     /// Maximum step size
    int VERBOSITY = 1;              /// Verbosity of the HybridEvolver
	if (argc > 1)
		VERBOSITY = atoi(argv[1]);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    AtomicHybridAutomaton automaton("exponential");

    /// Create the discrete states
    AtomicDiscreteLocation work("work");

    RealVariable x("x");
    Real A_real(A);

	automaton.new_mode(work,{dot(x) = - A_real*(x + x*x/2)});

    GeneralHybridEvolver evolver(automaton);
    evolver.verbosity=VERBOSITY;

    /// Set the evolution parameters
    evolver.configuration().set_maximum_enclosure_radius(MAX_ENCL_WIDTH);
    evolver.configuration().set_maximum_step_size(MAX_STEP_SIZE);
    std::cout <<  evolver.configuration() << std::endl;

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolver::EnclosureType EnclosureType;
    typedef GeneralHybridEvolver::EnclosureListType EnclosureListType;
    typedef GeneralHybridEvolver::OrbitType OrbitType;

    HybridSet initial_set(automaton|work,{1<=x<=1});
    HybridTime evolution_time(EVOL_TIME,4);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(initial_set,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    //std::cout << "Orbit="<<orbit<<std::endl;
    std::cout << "Orbit.final()="<<orbit.final()<<std::endl;
    plot("exponential-tx",Axes2d(0.0,TimeVariable(),EVOL_TIME, -1.1,x,1.1), Colour(0.0,0.5,1.0), orbit);
}
