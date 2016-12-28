/*****************************************************************************************************
 *            laser.cc
 *
 *  Testing the laser cutting system.
 *
 *  Copyright  2016  Luca Geretti
 *
 *****************************************************************************************************/

#include "ariadne.h"

#include "laser/skin-exposure.h"
#include "laser/laser-trajectory.h"
#include "laser/cutting-depth.h"
#include "laser/skin-temperature.h"
#include "common/timer.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    /// Constants
    int VERBOSITY = 1;
	if (argc > 1)
		VERBOSITY = atoi(argv[1]);

    /// Build the Hybrid System

    /// Get the automata
	AtomicHybridAutomaton timer = getTimer();
    AtomicHybridAutomaton laser_trajectory = getLaserTrajectory();
    AtomicHybridAutomaton exposure = getSkinExposure();
    AtomicHybridAutomaton skin_temperature = getSkinTemperature();
    AtomicHybridAutomaton cutting_depth = getCuttingDepth();

    CompositeHybridAutomaton laser_system({laser_trajectory,timer,exposure,skin_temperature,cutting_depth});
    std::cout << "laser_system:\n" << laser_system << "\n";

    // Compute the system evolution

    // Create a GeneralHybridEvolver object
    GeneralHybridEvolver evolver(laser_system);
    evolver.verbosity = VERBOSITY;

    // Set the evolution parameters
    evolver.configuration().set_maximum_enclosure_radius(0.5);
    evolver.configuration().set_maximum_step_size(0.000002);
    std::cout << evolver.configuration() << std::endl;

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolver::EnclosureType HybridEnclosureType;
    typedef GeneralHybridEvolver::OrbitType OrbitType;
    typedef GeneralHybridEvolver::EnclosureListType EnclosureListType;

	AtomicDiscreteLocation work("work");
	AtomicDiscreteLocation scanning("scanning");
	AtomicDiscreteLocation far("far");
	AtomicDiscreteLocation varying("varying");
	AtomicDiscreteLocation idle("idle");
	DiscreteLocation initial_location={timer|work,laser_trajectory|scanning,exposure|far,skin_temperature|varying,cutting_depth|idle};
    RealVariable t("t"), p("p"), z("z"), zi("zi"), T("T"), vx("vx"), x("x");
    HybridSet initial_set(initial_location,{t==0,p==0,z==0,zi==0,T==37,vx==Real(-0.092),x==Real(0.0033)});

    HybridTime evolution_time(0.05,7);

    std::cout << "Computing orbit... " << std::flush;

    OrbitType orbit = evolver.orbit(initial_set,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << orbit << std::endl;

}
