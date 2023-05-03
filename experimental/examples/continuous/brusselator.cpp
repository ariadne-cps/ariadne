/***************************************************************************
 *            brusselator.cpp
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

#include "ariadne_main.hpp"
#include "helper/stopwatch.hpp"

using namespace Ariadne;

void ariadne_main()
{
    RealVariable x1("x1"),x2("x2"),x3("x3"),x4("x4"),x5("x5"),x6("x6"),x7("x7");
    RealConstant a("a",1),b("b",1.5_dec);
    VectorField dynamics({dot(x1)=a+sqr(x1)*x2-b*x1-x1, dot(x2)=b*x1-sqr(x1)*x2});

    StepMaximumError max_err=1e-6;
    //TaylorPicardIntegrator integrator(max_err);
    GradedTaylorSeriesIntegrator integrator(max_err,ThresholdSweeper<FloatDP>(DoublePrecision(),1e-8),LipschitzTolerance(0.5_x),MaximumTemporalOrder(5));

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.02);
    evolver.configuration().set_maximum_spacial_error(1e-6);
    CONCLOG_PRINTLN(evolver.configuration())

    RealExpressionBoundedConstraintSet initial_set({{0.8_dec<=x1<=1},{0<=x2<=0.2_dec}});

    CONCLOG_PRINTLN_VAR(initial_set);
    Real evolution_time(3);

    Stopwatch<Milliseconds> sw;
    CONCLOG_PRINTLN("Computing orbit...")
    auto orbit = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    sw.click();
    CONCLOG_PRINTLN("Done in " << sw.elapsed_seconds() << " seconds.");

    CONCLOG_PRINTLN("Plotting...")
    CONCLOG_PRINTLN_AT(1,"Plotting...")
    LabelledFigure fig(Axes2d{{0<=x1<=2},{0<=x2<=2}});
    fig << orbit;
    fig.write("brusselator");
    sw.click();
    CONCLOG_PRINTLN("Done in " << sw.elapsed_seconds() << " seconds.");
}
