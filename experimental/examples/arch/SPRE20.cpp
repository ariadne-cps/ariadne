/***************************************************************************
 *            SPRE20.cpp
 *
 *  Copyright  2020  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "SPRE20.hpp"
#include <cstring>

using namespace Ariadne;

Int main(Int argc, const char* argv[])
{
    Logger::set_verbosity(get_verbosity(argc,argv));

    auto problem = make_space_rendezvous_problem();

    auto initial_set = problem.initial_set;
    auto system = problem.system;
    HybridConstraintSet safe_set = problem.safe_set;

    ARIADNE_LOG_PRINTLN("Space Rendezvous system:");

    DiscreteLocation initial_location = initial_set.location();
    RealSpace initial_space = system.state_space()[initial_location];
    BoundedConstraintSet initial_constraint_set = initial_set.euclidean_set(initial_location,initial_space);

    MaximumError max_err=1e-3;
    TaylorSeriesIntegrator integrator(max_err,Order(3u));

    GeneralHybridEvolver evolver(system);
    evolver.set_integrator(integrator);
    evolver.configuration().set_maximum_enclosure_radius(12.0);
    evolver.configuration().set_maximum_step_size(1.0);
    evolver.configuration().set_maximum_spacial_error(1e-3);
    evolver.configuration().set_enable_subdivisions(true);

    StopWatch sw;

    ARIADNE_LOG_PRINTLN("Computing orbit...");
    HybridTime evolution_time(200.0,3);
    auto orbit=evolver.orbit(initial_set,evolution_time,Semantics::UPPER);

    StringVariable spacecraft("spacecraft");
    StringConstant approaching("approaching");
    StringConstant rendezvous("rendezvous");
    StringConstant aborting("aborting");

    ARIADNE_LOG_PRINTLN("Checking properties...");
    Nat num_ce = 0;
    for (auto reach : orbit.reach()) {
        if (reach.location() == DiscreteLocation(spacecraft|rendezvous) and not(definitely(safe_set.covers(reach.bounding_box())))) {
            ARIADNE_LOG_PRINTLN_AT(1,"Found counterexample in location " << reach.location() << " with bounding box " << reach.bounding_box() << ", unsafe");
            ++num_ce;
        }
        if (reach.location() == DiscreteLocation(spacecraft|aborting) and not(definitely(safe_set.separated(reach.bounding_box())))) {
            ARIADNE_LOG_PRINTLN_AT(1,"Found counterexample in location " << reach.location() << " with bounding box " << reach.bounding_box() << ", unsafe");
            ++num_ce;
        }
    }
    if (num_ce>0) ARIADNE_LOG_PRINTLN("Number of counterexamples: " << num_ce);

    sw.click();
    ARIADNE_LOG_PRINTLN("Done in " << sw.elapsed() << " seconds.");

    RealVariable t("t"), x("x"), y("y"), vx("vx"), vy("vy");

    ARIADNE_LOG_PRINTLN("Plotting...");
    plot("SPRE20",{-1000<=x<=200,-450<=y<=0},Colour(1.0,0.75,0.5),orbit.reach());
    ARIADNE_LOG_PRINTLN("File SPRE20.png written.");
}
