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

#include "utility/stopwatch.hpp"
#include "function/taylor_function.hpp"
#include "function/constraint.hpp"
#include "dynamics/enclosure.hpp"
#include "ariadne_main.hpp"

void ariadne_main()
{
    CONCLOG_PRINTLN("van der Pol oscillator");

    RealConstant mu("mu",1);
    RealVariable x("x"), y("y");

    VectorField dynamics({dot(x)=y, dot(y)= mu*y*(1-sqr(x))-x});

    VectorFieldSimulator simulator(dynamics);
    simulator.configuration().set_step_size(0.02);

    StepMaximumError max_err=1e-6;

    GradedTaylorPicardIntegrator integrator(max_err,order=5);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.02);
    evolver.configuration().set_maximum_spacial_error(1e-6);
    // Temporary diagnostic: isolate the cost of graded propagation from
    // enclosure reconditioning.
    evolver.configuration().set_enable_reconditioning(false);
    CONCLOG_PRINTLN(evolver.configuration());

    Real x0 = 1.40_dec;
    Real y0 = 2.40_dec;
    Real eps_x0 = 0.15_dec;
    Real eps_y0 = 0.05_dec;

    RealExpressionBoundedConstraintSet initial_set({x0-eps_x0<=x<=x0+eps_x0,y0-eps_y0<=y<=y0+eps_y0});

    CONCLOG_PRINTLN("Initial set: " << initial_set);
    Real evolution_time = 7;

    Stopwatch<Milliseconds> sw;
    CONCLOG_PRINTLN("Computing simulation...");
    auto simulation = simulator.orbit(initial_set,evolution_time);
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Done in " << sw.elapsed_seconds() << " seconds.");

    CONCLOG_PRINTLN("Plotting...");;
    LabelledFigure fig=LabelledFigure({-2.5<=x<=2.5,-3<=y<=3});
    fig.draw(simulation);
    fig.write("vanderpol_simulation");

    sw.restart();
    CONCLOG_PRINTLN("Computing evolution... ");

    // Temporary diagnostic: compare the two Picard integrators on exactly
    // the same initial box and requested step before running the full evolver.
    auto diagnostic_initial_box = cast_exact_box(initial_set.euclidean_set(dynamics.state_space()).bounding_box());
    StepSizeType diagnostic_step=0.02_dy;
    TaylorPicardIntegrator diagnostic_taylor_picard(max_err);
    GradedTaylorPicardIntegrator diagnostic_graded_picard(max_err,order=5);

    Stopwatch<Microseconds> diagnostic_sw;
    std::cerr << "[vanderpol] TaylorPicard direct probe" << std::endl;
    auto diagnostic_taylor_flow = diagnostic_taylor_picard.flow_step(
        dynamics.function(),diagnostic_initial_box,suggest(diagnostic_step));
    diagnostic_sw.click();
    std::cerr << "[vanderpol] TaylorPicard time_us=" << diagnostic_sw.duration().count()
              << " error=" << diagnostic_taylor_flow.error() << std::endl;

    diagnostic_sw.restart();
    std::cerr << "[vanderpol] GradedTaylorPicard direct probe" << std::endl;
    auto diagnostic_graded_flow = diagnostic_graded_picard.flow_step(
        dynamics.function(),diagnostic_initial_box,suggest(diagnostic_step));
    diagnostic_sw.click();
    std::cerr << "[vanderpol] GradedTaylorPicard time_us=" << diagnostic_sw.duration().count()
              << " error=" << diagnostic_graded_flow.error() << std::endl;

    // Temporary synchronous propagation benchmark.  This deliberately avoids
    // VectorFieldEvolver/DynamicWorkload and measures flow construction and
    // enclosure propagation separately for a short prefix of the trajectory.
    LabelledEnclosure diagnostic_enclosure(
        initial_set.euclidean_set(dynamics.state_space()),dynamics.state_space(),
        EnclosureConfiguration(integrator.function_factory()));
    diagnostic_enclosure.set_auxiliary(dynamics.auxiliary_space(),dynamics.auxiliary_mapping());

    TimeStepType diagnostic_time(0u);
    for(Nat diagnostic_step_index=0; diagnostic_step_index!=20; ++diagnostic_step_index) {
        auto const& sf=diagnostic_enclosure.state_function();
        auto const& sf_taylor=dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(sf.reference());
        SizeType state_nnz=0;
        for(SizeType i=0; i!=sf_taylor.size(); ++i) { state_nnz+=sf_taylor[i].number_of_nonzeros(); }
        auto const state_error_before=sf.error();
        SizeType const state_params_before=diagnostic_enclosure.number_of_parameters();

        auto box=cast_exact_box(diagnostic_enclosure.euclidean_set().bounding_box());
        diagnostic_sw.restart();
        auto flow=integrator.flow_step(dynamics.function(),box,suggest(diagnostic_step));
        diagnostic_sw.click();
        auto const& flow_taylor=dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(flow.reference());
        SizeType flow_nnz=0;
        for(SizeType i=0; i!=flow_taylor.size(); ++i) { flow_nnz+=flow_taylor[i].number_of_nonzeros(); }
        auto const flow_us=diagnostic_sw.duration().count();

        StepSizeType actual_step=static_cast<StepSizeType>(flow.domain()[flow.argument_size()-1u].upper_bound());
        diagnostic_sw.restart();
        diagnostic_enclosure.apply_fixed_evolve_step(flow,actual_step);
        diagnostic_sw.click();

        std::cerr << "[SyncProfile] step=" << diagnostic_step_index
                  << " t=" << diagnostic_time
                  << " params_before=" << state_params_before
                  << " state_nnz_before=" << state_nnz
                  << " state_error_before=" << state_error_before
                  << " flow_us=" << flow_us
                  << " flow_nnz=" << flow_nnz
                  << " flow_error=" << flow.error()
                  << " evolve_us=" << diagnostic_sw.duration().count()
                  << " next_error=" << diagnostic_enclosure.state_function().error()
                  << std::endl;
        diagnostic_time+=TimeStepType(actual_step);
    }

    std::cerr << "[vanderpol] starting graded Taylor-Picard evolution" << std::endl;
    auto evolution = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Done in " << sw.elapsed_seconds() << " seconds.");

    CONCLOG_PRINTLN("Plotting...");
    fig.clear();
    fig.draw(evolution);
    fig.write("vanderpol_evolution");
}
