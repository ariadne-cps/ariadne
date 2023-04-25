/***************************************************************************
 *            constrained_vanderpol.cpp
 *
 *  Copyright  2023  Luca Geretti
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

#include "helper/stopwatch.hpp"
#include "ariadne_main.hpp"

void ariadne_main()
{
    RealConstant mu("mu",1);
    RealVariable x("x"), y("y");

    VectorField dynamics({dot(x)=y, dot(y)= mu*y*(1-sqr(x))-x});

    auto approximate_configuration = Configuration<VectorFieldEvolver>().
        set_maximum_step_size(0.1);

    AffineIntegrator approximate_integrator(2,1);
    VectorFieldEvolver approximate_evolver(dynamics,approximate_configuration,approximate_integrator);
    CONCLOG_PRINTLN_VAR(approximate_evolver.configuration())

    auto rigorous_configuration = Configuration<VectorFieldEvolver>().
            set_maximum_enclosure_radius(1.0).
            set_maximum_step_size(0.02).
            set_maximum_spacial_error(1e-6);

    StepMaximumError max_err=1e-6;
    TaylorPicardIntegrator rigorous_integrator(max_err);
    VectorFieldEvolver rigorous_evolver(dynamics,rigorous_configuration,rigorous_integrator);
    CONCLOG_PRINTLN_VAR(rigorous_evolver.configuration())

    Real x0 = 1.40_dec;
    Real y0 = 2.40_dec;
    Real eps_x0 = 0.15_dec;
    Real eps_y0 = 0.05_dec;

    RealExpressionBoundedConstraintSet initial_set({x0-eps_x0<=x<=x0+eps_x0,y0-eps_y0<=y<=y0+eps_y0});

    CONCLOG_PRINTLN("Initial set: " << initial_set)
    Real evolution_time = 7;

    Helper::Stopwatch<std::chrono::milliseconds> sw;
    CONCLOG_PRINTLN("Computing approximate evolution...")
    auto approximate_orbit = approximate_evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Done in " << sw.elapsed_seconds() << " seconds.")

    CONCLOG_PRINTLN("Plotting...")
    LabelledFigure fig=LabelledFigure({-2.5<=x<=2.5,-3<=y<=3});
    fig.draw(approximate_orbit);
    fig.write("vanderpol_approximate");

    sw.restart();
    CONCLOG_PRINTLN("Computing rigorous evolution... ")
    auto rigorous_evolution = rigorous_evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Done in " << sw.elapsed_seconds() << " seconds.")

    CONCLOG_PRINTLN("Plotting...")
    fig.clear();
    fig.draw(rigorous_evolution);
    fig.write("vanderpol_rigorous");
}
