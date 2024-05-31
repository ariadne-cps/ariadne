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
    RealVariable x("x"), y("y");
    VectorField dynamics({dot(x)=-y-1.5_dec*pow(x,2)-0.5_dec*pow(x,3)-0.5_dec,dot(y)=3*x-y});

    Real e1=5/100_q; Real e2=7/100_q;
    RealExpressionBoundedConstraintSet initial_set({1-e1<=x<=1+e1,1-e2<=y<=1+e2});

    Real evolution_time = 0.81_dec;

    TaylorPicardIntegrator integrator(Configuration<TaylorPicardIntegrator>().set_step_maximum_error(1e-6));

    auto configuration = Configuration<VectorFieldEvolver>().
            set_maximum_enclosure_radius(1.0).
            set_maximum_step_size(0.02).
            set_maximum_spacial_error(1e-6).
            set_integrator(integrator);

    VectorFieldEvolver evolver(dynamics,configuration);

    Stopwatch<Milliseconds> sw;
    CONCLOG_PRINTLN("Computing orbit...")
    auto orbit = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    sw.click();
    CONCLOG_PRINTLN("Done in " << sw.elapsed_seconds() << " seconds.");

    CONCLOG_PRINTLN("Plotting...")
    CONCLOG_PRINTLN_AT(1,"Plotting...")
    LabelledFigure fig(Axes2d{{-1<=x<=2},{-1<=y<=2}});
    fig << orbit;
    fig.write("brusselator");
    sw.click();
    CONCLOG_PRINTLN("Done in " << sw.elapsed_seconds() << " seconds.");
}
