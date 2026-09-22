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

    GradedTaylorPicardIntegrator integrator(max_err,order=5,step_sweep_threshold=1e-12);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.02);
    evolver.configuration().set_maximum_spacial_error(1e-6);
    CONCLOG_PRINTLN(evolver.configuration());

    Real x0 = 1.40_dec;
    Real y0 = 2.30_dec;
    Real eps_x0 = 0.15_dec;
    Real eps_y0 = 0.05_dec;

    RealExpressionBoundedConstraintSet initial_set({x0-eps_x0<=x<=x0+eps_x0,y0-eps_y0<=y<=y0+eps_y0});

    CONCLOG_PRINTLN("Initial set: " << initial_set);

    // Local single-step comparison before exercising the evolver.  Both
    // integrators use fixed spatial/temporal order 5; the only algorithmic
    // difference is the diagonal state-domain preconditioning.
    {
        ExactBoxType const benchmark_box=cast_exact_box(
            initial_set.euclidean_set(dynamics.state_space()).bounding_box());
        ThresholdSweeper<FloatDP> benchmark_sweeper(DoublePrecision(),1e-12);

        GradedTaylorSeriesIntegrator graded_integrator(
            step_maximum_error=1e-3,benchmark_sweeper,lipschitz_tolerance=0.5_x,
            minimum_spacial_order=5,minimum_temporal_order=5,
            maximum_spacial_order=5,maximum_temporal_order=5);

        PreconditionedGradedTaylorSeriesIntegrator preconditioned_integrator(
            step_maximum_error=1e-3,benchmark_sweeper,lipschitz_tolerance=0.5_x,
            minimum_spacial_order=5,minimum_temporal_order=5,
            maximum_spacial_order=5,maximum_temporal_order=5);

        const StepSizeType benchmark_steps[] = {
            StepSizeType(0.02_dy), StepSizeType(0.01_dy),
            StepSizeType(0.005_dy), StepSizeType(0.0025_dy)
        };

        auto run_step = [&](const char* method,
                            IntegratorInterface const& candidate,
                            StepSizeType const& step)
        {
            Stopwatch<Microseconds> step_sw;
            bool completed=true;
            String failure;
            FlowStepModelType flow;

            try {
                step_sw.restart();
                flow=candidate.flow_step(dynamics.function(),benchmark_box,step);
                step_sw.click();
            } catch(const std::exception& e) {
                step_sw.click();
                completed=false;
                failure=e.what();
            }

            std::cerr << "[PreconditionedSeriesStep]"
                      << " method=" << method
                      << " requested_step=" << step
                      << " completed=" << completed
                      << " time_us=" << step_sw.duration().count();

            if(completed) {
                auto const& taylor=
                    dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(
                        flow.reference());
                SizeType nnz=0u;
                for(SizeType i=0u; i!=taylor.size(); ++i) {
                    nnz+=taylor[i].number_of_nonzeros();
                }

                auto at_zero=partial_evaluate(
                    flow,flow.argument_size()-1u,StepSizeType(0_dy));
                auto identity=factory(at_zero).create_identity();
                auto initial_defect=at_zero-identity;
                auto at_end=partial_evaluate(
                    flow,flow.argument_size()-1u,step);

                std::cerr << " nnz=" << nnz
                          << " error=" << flow.error()
                          << " component_errors=" << flow.errors()
                          << " initial_defect_range=" << initial_defect.range()
                          << " initial_defect_error=" << initial_defect.error()
                          << " endpoint_range=" << at_end.range()
                          << " flow_range=" << flow.range();
            } else {
                std::cerr << " failure=\"" << failure << "\"";
            }
            std::cerr << std::endl;
        };

        for(StepSizeType const& step : benchmark_steps) {
            run_step("GradedTaylorSeries",graded_integrator,step);
            run_step("PreconditionedGradedTaylorSeries",preconditioned_integrator,step);
        }
    }

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

    auto evolution = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Done in " << sw.elapsed_seconds() << " seconds.");

    CONCLOG_PRINTLN("Plotting...");
    fig.clear();
    fig.draw(evolution);
    fig.write("vanderpol_evolution");
}
