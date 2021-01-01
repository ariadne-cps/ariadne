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
#include "concurrency/task_parameter.hpp"
#include "concurrency/task_parameter_space.hpp"

using namespace Ariadne;

int main(int argc, const char* argv[])
{
    Logger::configuration().set_verbosity(get_verbosity(argc,argv));
    Logger::configuration().set_theme(TT_THEME_DARK);
    Logger::configuration().set_thread_name_printing_policy(ThreadNamePrintingPolicy::BEFORE);

    ARIADNE_LOG_PRINTLN("van der Pol oscillator");

    RealConstant mu("mu",1);
    RealVariable x("x"), y("y");

    VectorField dynamics({dot(x)=y, dot(y)=mu*y*(1-sqr(x))-x});

    MaximumError max_err=1e-6;
    StartingStepSizeNumRefinements num_ref=2;
    TaylorPicardIntegrator integrator(max_err,
                                      ThresholdSweeper<FloatDP>(DoublePrecision(),max_err.value()/1024),
                                      LipschitzConstant(0.5),
                                      num_ref);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(1e20);
    evolver.configuration().set_maximum_spacial_error(1e-6);
    ARIADNE_LOG_PRINTLN(evolver.configuration());

    evolver.set_runner(SharedPointer<typename VectorFieldEvolver::RunnerType>(new VectorFieldEvolverFlowStepConcurrentRunner(evolver)));

    Real x0 = 1.4_dec;
    Real y0 = 2.4_dec;
    Real eps_x0 = 15/100_q;
    Real eps_y0 = 5/100_q;

    EnclosureConfiguration enclosure_config(evolver.function_factory());
    enclosure_config.set_reconditioning_num_blocks(4);
    auto initial_set = evolver.enclosure({x0-eps_x0<=x<=x0+eps_x0,y0-eps_y0<=y<=y0+eps_y0},enclosure_config);
    ARIADNE_LOG_PRINTLN("Initial set: " << initial_set);
    ARIADNE_LOG_PRINTLN(initial_set.configuration());

    Real evolution_time = 7;

    StopWatch sw;
    ARIADNE_LOG_PRINTLN("Computing orbit... ");
    ARIADNE_LOG_RUN_AT(1,auto orbit = evolver.orbit(initial_set,evolution_time,Semantics::UPPER));
    sw.click();
    ARIADNE_LOG_PRINTLN_AT(1,"Done in " << sw.elapsed() << " seconds.");

    sw.reset();
    ARIADNE_LOG_PRINTLN("Plotting...");
    LabelledFigure fig({-2.5<=x<=2.5,-3<=y<=3});
    fig << fill_colour(1.0,0.75,0.5);
    fig.draw(orbit.reach());
    fig.write("vanderpol");

    sw.click();
    ARIADNE_LOG_PRINTLN_AT(1,"Done in " << sw.elapsed() << " seconds.");
}
