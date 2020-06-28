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
    Logger::configuration().set_verbosity(get_verbosity(argc,argv));

    ARIADNE_LOG_PRINTLN("van der Pol oscillator");

    RealConstant mu("mu",1.0_dec);
    RealVariable x("x"), y("y");

    VectorField dynamics({dot(x)=y, dot(y)= mu*y*(1-sqr(x))-x});

    MaximumError max_err=1e-6;

    TaylorPicardIntegrator integrator(max_err);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.02);
    evolver.configuration().set_maximum_spacial_error(1e-6);
    ARIADNE_LOG_PRINTLN(evolver.configuration());

    Real x0(1.40);
    Real y0(2.40);
    Real eps_x0 = 15/100_q;
    Real eps_y0 = 5/100_q;

    RealVariablesBox initial_set({x0-eps_x0<=x<=x0+eps_x0,y0-eps_y0<=y<=y0+eps_y0});

    ARIADNE_LOG_PRINTLN("Initial set: " << initial_set);
    Real evolution_time(7.0);

    StopWatch sw;
    ARIADNE_LOG_PRINTLN("Computing orbit... ");
    ARIADNE_LOG_RUN_AT(1,auto orbit = evolver.orbit(evolver.enclosure(initial_set),evolution_time,Semantics::UPPER));
    sw.click();
    ARIADNE_LOG_PRINTLN_AT(1,"Done in " << sw.elapsed() << " seconds.");

    sw.reset();
    ARIADNE_LOG_PRINTLN("Plotting...");
    Axes2d axes(-2.5<=x<=2.5,-3.0<=y<=3.0);
    LabelledFigure fig=LabelledFigure(axes);
    fig << line_colour(0.0,0.0,0.0);
    fig << line_style(false);
    fig << fill_colour(0.5,0.5,0.5);
    fig << fill_colour(1.0,0.75,0.5);
    fig.draw(orbit.reach());
    fig.write("vanderpol");

    sw.click();
    ARIADNE_LOG_PRINTLN_AT(1,"Done in " << sw.elapsed() << " seconds.");
}
