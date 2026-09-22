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


namespace {

ValidatedVectorMultivariateFunctionPatch affine_parameterisation(
    const ValidatedFunctionPatchFactory& factory,
    const ExactBoxType& unit_domain,
    const ExactBoxType& physical_box)
{
    SizeType const d=physical_box.size();
    ValidatedVectorMultivariateFunctionPatch result=
        factory.create_zeros(d,unit_domain);
    for(SizeType i=0u; i!=d; ++i) {
        FloatDP const c=physical_box[i].midpoint().raw();
        FloatDP const r=physical_box[i].radius().upper().raw();
        result[i]=factory.create_constant(unit_domain,c)
            + factory.create_coordinate(unit_domain,i)*FloatDPBounds(r);
    }
    return result;
}

ValidatedVectorMultivariateFunctionPatch flowstar_normalise_mapping(
    const ValidatedFunctionPatchFactory& factory,
    const ValidatedVectorMultivariateFunctionPatch& state,
    ExactBoxType& physical_box)
{
    // Match Flow*'s normalisation more closely: take the Taylor-model
    // constant as the new centre, remove it, then scale each centred
    // component by the magnitude of its validated range.  This differs from
    // centring the bounding box, which can move the centre using nonlinear
    // range information and is not what Flow* does.
    auto const& state_taylor =
        dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(
            state.reference());

    ExactBoxType const domain=state.domain();
    SizeType const d=state_taylor.size();
    physical_box=ExactBoxType(d);
    ValidatedVectorMultivariateFunctionPatch result=
        factory.create_zeros(d,domain);

    for(SizeType i=0u; i!=d; ++i) {
        FloatDP const c=state_taylor.model(i).value().raw();
        ValidatedScalarMultivariateFunctionPatch centred=
            state[i]-FloatDPBounds(c);
        FloatDP const r=cast_exact(mag(centred.range()));

        physical_box[i]=ExactIntervalType(sub(down,c,r),add(up,c,r));
        if(r==FloatDP(0,dp)) {
            result[i]=factory.create_zero(domain);
        } else {
            result[i]=centred/FloatDPBounds(r);
        }
    }
    return result;
}

}

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

    // Compare the two existing series-based local-flow constructors before
    // choosing the implementation base for the preconditioned integrator.
    //
    // This is deliberately a single-step benchmark: it measures the quality
    // and cost of constructing the local flow itself, without reconditioning
    // or long-horizon wrapping effects.  The TaylorSeries case uses one total
    // degree.  The GradedTaylorSeries case uses separate spatial/temporal
    // degree budgets, as intended by that implementation.
    {
        ExactBoxType const benchmark_box=cast_exact_box(
            initial_set.euclidean_set(dynamics.state_space()).bounding_box());
        const StepSizeType benchmark_steps[] = {
            StepSizeType(0.02_dy), StepSizeType(0.01_dy),
            StepSizeType(0.005_dy), StepSizeType(0.0025_dy)
        };

        ThresholdSweeper<FloatDP> series_sweeper(DoublePrecision(),1e-12);
        TaylorSeriesIntegrator series_integrator(
            series_sweeper,lipschitz_tolerance=0.5_x,order=5);
        GradedTaylorSeriesIntegrator graded_series_integrator(
            step_maximum_error=1e-3,series_sweeper,lipschitz_tolerance=0.5_x,
            minimum_spacial_order=1,minimum_temporal_order=5,
            maximum_spacial_order=5,maximum_temporal_order=5);

        for(SizeType step_case=0u; step_case!=4u; ++step_case) {
            StepSizeType const step=benchmark_steps[step_case];

            auto run_series_step = [&](const char* method, IntegratorInterface const& candidate) {
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

                std::cerr << "[SeriesStepBenchmark]"
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
                    std::cerr << " nnz=" << nnz
                              << " error=" << flow.error()
                              << " component_errors=" << flow.errors()
                              << " range=" << flow.range();
                } else {
                    std::cerr << " failure=\"" << failure << "\"";
                }
                std::cerr << std::endl;
            };

            run_series_step("TaylorSeries",series_integrator);
            run_series_step("GradedTaylorSeries",graded_series_integrator);
        }
    }

    // Reference Ariadne reconditioning at a fixed physical calendar.
    {
        const StepSizeType benchmark_steps[] = {
            StepSizeType(0.005_dy), StepSizeType(0.0025_dy)
        };
        const Nat benchmark_num_steps[] = {1400u,2800u};
        const Nat benchmark_period_steps[] = {8u,16u}; // 0.04 physical time

        for(SizeType step_case=0u; step_case!=2u; ++step_case) {
            GradedTaylorPicardIntegrator benchmark_integrator(
                step_maximum_error=1e-3,order=5,step_sweep_threshold=1e-12);
            benchmark_integrator.set_maximum_error_refinement_iterations(2u);
            benchmark_integrator.set_diagnostics(false);

            LabelledEnclosure benchmark_enclosure(
                initial_set.euclidean_set(dynamics.state_space()),dynamics.state_space(),
                EnclosureConfiguration(benchmark_integrator.function_factory()));
            benchmark_enclosure.set_auxiliary(
                dynamics.auxiliary_space(),dynamics.auxiliary_mapping());

            StepSizeType const benchmark_step=benchmark_steps[step_case];
            Nat const num_steps=benchmark_num_steps[step_case];
            Nat const period_steps=benchmark_period_steps[step_case];
            Nat reconditionings=0u;
            bool completed=true;
            Nat failed_step=num_steps;
            String failure;

            for(Nat step_index=0u; step_index!=num_steps; ++step_index) {
                auto box=cast_exact_box(
                    benchmark_enclosure.euclidean_set().bounding_box());
                try {
                    auto flow=benchmark_integrator.flow_step(
                        dynamics.function(),box,benchmark_step);
                    benchmark_enclosure.apply_fixed_evolve_step(flow,benchmark_step);
                    if(((step_index+1u)%period_steps)==0u) {
                        benchmark_enclosure.recondition();
                        ++reconditionings;
                    }
                } catch(const std::exception& e) {
                    completed=false;
                    failed_step=step_index;
                    failure=e.what();
                    break;
                }
            }

            auto const final_box=benchmark_enclosure.euclidean_set().bounding_box();
            double final_width=0.0;
            std::vector<double> component_widths(final_box.size());
            for(SizeType i=0u; i!=final_box.size(); ++i) {
                double const w=final_box[i].width().get_d();
                component_widths[i]=w;
                if(w>final_width) { final_width=w; }
            }
            std::cerr << "[StandardRecondition]"
                      << " requested_step=" << benchmark_step
                      << " interval=0.04"
                      << " completed=" << completed
                      << " completed_steps=" << (completed ? num_steps : failed_step)
                      << " reconditionings=" << reconditionings
                      << " final_width=" << final_width
                      << " component_widths=[" << component_widths[0]
                      << "," << component_widths[1] << "]"
                      << " final_box=" << final_box;
            if(not completed) { std::cerr << " failure=\"" << failure << "\""; }
            std::cerr << std::endl;
        }
    }

    // Flow*-style local-initial-set normalisation.
    //
    // Flow* keeps two Taylor-model layers: a local preconditioned flow and a
    // normalised map from the original parameters into the local coordinates.
    // At the end of every step it composes the two once, recentres the resulting
    // local initial set, scales each component to [-1,1], and carries the
    // normalised Taylor map to the next step.  No parameter dependence is
    // discarded here.
    {
        const StepSizeType steps[] = {
            StepSizeType(0.02_dy), StepSizeType(0.01_dy),
            StepSizeType(0.005_dy), StepSizeType(0.0025_dy)
        };
        const Nat num_steps_values[] = {350u,700u,1400u,2800u};

        for(SizeType step_case=0u; step_case!=4u; ++step_case) {
            GradedTaylorPicardIntegrator local_integrator(
                step_maximum_error=1e-3,order=5,step_sweep_threshold=1e-12);
            local_integrator.set_maximum_error_refinement_iterations(2u);
            local_integrator.set_diagnostics(false);

            auto const& factory=local_integrator.function_factory();
            StepSizeType const step=steps[step_case];
            Nat const num_steps=num_steps_values[step_case];

            ExactBoxType local_box=cast_exact_box(
                initial_set.euclidean_set(dynamics.state_space()).bounding_box());
            ExactBoxType const unit_domain(
                local_box.size(),ExactIntervalType(-1,+1));

            // r_0 is the identity on the normalised initial parameters.
            ValidatedVectorMultivariateFunctionPatch normalised_map=
                factory.create_identity(unit_domain);
            ValidatedVectorMultivariateFunctionPatch state=
                affine_parameterisation(factory,unit_domain,local_box);

            bool completed=true;
            Nat failed_step=num_steps;
            String failure;
            SizeType max_mapping_nnz=0u;
            Stopwatch<Microseconds> local_sw;
            long long flow_us=0;
            long long compose_us=0;
            long long normalise_us=0;

            for(Nat step_index=0u; step_index!=num_steps; ++step_index) {
                try {
                    local_sw.restart();
                    auto flow=local_integrator.flow_step(
                        dynamics.function(),local_box,step);
                    local_sw.click();
                    flow_us+=local_sw.duration().count();

                    ValidatedVectorMultivariateFunctionPatch flow_model=flow;
                    auto flow_end=partial_evaluate(
                        flow_model,flow_model.argument_size()-1u,step);

                    // P_l(y,h): compose the physical flow with x=c_l+S_l*y.
                    auto local_initial_map=
                        affine_parameterisation(factory,unit_domain,local_box);
                    auto preconditioned_end=compose(flow_end,local_initial_map);

                    // X_{l+1}(s)=P_l(r_l(s),h), preserving all parameter
                    // correlations before the normalisation for the next step.
                    local_sw.restart();
                    state=compose(preconditioned_end,normalised_map);
                    local_sw.click();
                    compose_us+=local_sw.duration().count();

                    auto const& state_taylor =
                        dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(
                            state.reference());
                    SizeType mapping_nnz=0u;
                    for(SizeType i=0u; i!=state_taylor.size(); ++i) {
                        mapping_nnz+=state_taylor[i].number_of_nonzeros();
                    }
                    if(mapping_nnz>max_mapping_nnz) { max_mapping_nnz=mapping_nnz; }

                    local_sw.restart();
                    normalised_map=
                        flowstar_normalise_mapping(factory,state,local_box);
                    local_sw.click();
                    normalise_us+=local_sw.duration().count();
                } catch(const std::exception& e) {
                    completed=false;
                    failed_step=step_index;
                    failure=e.what();
                    break;
                }
            }

            auto const final_box=state.codomain().bounding_box();
            double final_width=0.0;
            std::vector<double> component_widths(final_box.size());
            for(SizeType i=0u; i!=final_box.size(); ++i) {
                double const w=final_box[i].width().get_d();
                component_widths[i]=w;
                if(w>final_width) { final_width=w; }
            }

            std::cerr << "[FlowstarLikePrecondition]"
                      << " requested_step=" << step
                      << " completed=" << completed
                      << " completed_steps=" << (completed ? num_steps : failed_step)
                      << " max_mapping_nnz=" << max_mapping_nnz
                      << " flow_us=" << flow_us
                      << " compose_us=" << compose_us
                      << " normalise_us=" << normalise_us
                      << " final_error=" << state.error()
                      << " mapping_error=" << normalised_map.error()
                      << " final_width=" << final_width
                      << " component_widths=[" << component_widths[0]
                      << "," << component_widths[1] << "]"
                      << " final_box=" << final_box;
            if(not completed) { std::cerr << " failure=\"" << failure << "\""; }
            std::cerr << std::endl;
        }
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
