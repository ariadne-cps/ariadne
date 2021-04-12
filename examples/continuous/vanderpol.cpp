/***************************************************************************
 *            vanderpol.cpp
 *
 *  Copyright  2017-20  Luca Geretti
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

#include <cstdarg>
#include "ariadne.hpp"
#include "utility/stopwatch.hpp"

using namespace Ariadne;

int main(int argc, const char* argv[])
{
    ARIADNE_LOG_SET_VERBOSITY(get_verbosity(argc,argv));

    ARIADNE_LOG_PRINTLN("van der Pol oscillator");

    RealConstant mu("mu",1.0_dec);
    RealVariable x("x"), y("y");

    VectorField dynamics({dot(x)=y, dot(y)= mu*y*(1-sqr(x))-x});

    VectorFieldSimulator simulator(dynamics);
    simulator.configuration().set_step_size(0.02);

    MaximumError max_err=1e-6;

    TaylorPicardIntegrator integrator(max_err);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.02);
    evolver.configuration().set_maximum_spacial_error(1e-6);
    ARIADNE_LOG_PRINTLN(evolver.configuration());

    Real x0(1.40_dec);
    Real y0(2.40_dec);
    Real eps_x0 = 15/100_q;
    Real eps_y0 = 5/100_q;

    RealExpressionBoundedConstraintSet initial_set({x0-eps_x0<=x<=x0+eps_x0,y0-eps_y0<=y<=y0+eps_y0});

    ARIADNE_LOG_PRINTLN("Initial set: " << initial_set);
    Real evolution_time(7);

    StopWatch sw;
    ARIADNE_LOG_PRINTLN("Computing simulation...");
    ARIADNE_LOG_RUN_AT(1,auto simulation = simulator.orbit(initial_set,evolution_time));
    sw.click();
    ARIADNE_LOG_PRINTLN_AT(1,"Done in " << sw.elapsed() << " seconds.");

    ARIADNE_LOG_PRINTLN("Plotting...");;
    LabelledFigure sfig=LabelledFigure({-2.5<=x<=2.5,-3<=y<=3});
    sfig.draw(simulation.curve());
    sfig.write("vanderpol_simulation");

    sw.reset();
    ARIADNE_LOG_PRINTLN("Computing evolution... ");
    ARIADNE_LOG_RUN_AT(1,auto evolution = evolver.orbit(evolver.enclosure(initial_set),evolution_time));
    sw.click();
    ARIADNE_LOG_PRINTLN_AT(1,"Done in " << sw.elapsed() << " seconds.");

    ARIADNE_LOG_PRINTLN("Plotting...");;
    LabelledFigure efig=LabelledFigure({-2.5<=x<=2.5,-3<=y<=3});
    efig.draw(evolution.reach());
    efig.write("vanderpol_evolution");
}
