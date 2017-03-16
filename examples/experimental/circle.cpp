/*****************************************************************************************************
 *            circle-unstable.cpp
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the behavior of a circle.
 *
 *****************************************************************************************************/

#include "ariadne.hpp"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    /// Constants
    double EVOL_TIME = 32;   /// Evolution time
    float MAX_ENCL_WIDTH = 0.1;   /// Maximum enclosure width
    float MAX_STEP_SIZE = 4e-3;     /// Maximum step size
    int VERBOSITY = 1;              /// Verbosity of the HybridEvolver
	if (argc > 1)
		VERBOSITY = atoi(argv[1]);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    AtomicHybridAutomaton automaton("circle");

    /// Create the discrete states
    AtomicDiscreteLocation work("work");

    RealVariable x("x");
    RealVariable y("y");

	automaton.new_mode(work,{dot(x)=-y,dot(y)=x});

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

    Real e(0);//100.0/1024/1024);
    HybridSet initial_set(automaton|work,{1-e<=x<=1+e,y.in(0-e,0+e)});
    HybridTime evolution_time(EVOL_TIME,4);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(initial_set,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    //std::cout << "Orbit="<<orbit<<std::endl;
    std::cout << "Orbit.final()="<<orbit.final()<<std::endl;
    plot("circle-orbit",Axes2d(-1.1,x,1.1, -1.1,y,1.1), Colour(0.0,0.5,1.0), orbit);
    plot("circle-tx",Axes2d(0.0,TimeVariable(),EVOL_TIME, -1.1,x,1.1), Colour(0.0,0.5,1.0), orbit);
    plot("circle-ty",Axes2d(0.0,TimeVariable(),EVOL_TIME, -1.1,y,1.1), Colour(0.0,0.5,1.0), orbit);
}
