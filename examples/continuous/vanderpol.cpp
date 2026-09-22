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

    // Chen-aligned fixed-step benchmark.  This deliberately calls the exact
    // StepSizeType overload rather than suggest(...): Ariadne may try a smaller
    // step internally only to prove the flow bound, but then reports
    // IncompleteFlowException instead of silently accepting that smaller step.
    {
        GradedTaylorPicardIntegrator chen_integrator(
            step_maximum_error=1e-3,order=5,step_sweep_threshold=1e-12);
        chen_integrator.set_diagnostics(false);
        LabelledEnclosure chen_enclosure(
            initial_set.euclidean_set(dynamics.state_space()),dynamics.state_space(),
            EnclosureConfiguration(chen_integrator.function_factory()));
        chen_enclosure.set_auxiliary(dynamics.auxiliary_space(),dynamics.auxiliary_mapping());

        StepSizeType const chen_step=StepSizeType(0.02_dy);
        Nat const chen_num_steps=350u;
        Nat chen_reconditionings=0u;
        SizeType chen_max_state_nnz=0u;
        SizeType chen_max_reach_nnz=0u;
        long long chen_flow_us=0;
        long long chen_reach_us=0;
        long long chen_evolve_us=0;
        long long chen_recondition_us=0;
        Stopwatch<Microseconds> chen_sw;
        bool chen_completed=true;
        Nat chen_failed_step=chen_num_steps;
        String chen_failure;

        for(Nat step_index=0u; step_index!=chen_num_steps; ++step_index) {
            auto const errors_before_recondition=chen_enclosure.state_function().errors();
            bool did_recondition=false;
            if(possibly(chen_enclosure.state_function().error() > 1e-6_pr)) {
                chen_sw.restart();
                chen_enclosure.recondition();
                chen_sw.click();
                chen_recondition_us+=chen_sw.duration().count();
                ++chen_reconditionings;
                did_recondition=true;
            }
            auto const errors_after_recondition=chen_enclosure.state_function().errors();
            auto const state_errors_before_flow=chen_enclosure.state_function().errors();

            auto const& state_taylor =
                dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(
                    chen_enclosure.state_function().reference());
            SizeType state_nnz=0u;
            for(SizeType i=0; i!=state_taylor.size(); ++i) {
                state_nnz+=state_taylor[i].number_of_nonzeros();
            }
            if(state_nnz>chen_max_state_nnz) { chen_max_state_nnz=state_nnz; }

            auto box=cast_exact_box(chen_enclosure.euclidean_set().bounding_box());
            try {
                chen_sw.restart();
                auto flow=chen_integrator.flow_step(dynamics.function(),box,chen_step);
                chen_sw.click();
                chen_flow_us+=chen_sw.duration().count();
                auto const flow_errors=flow.errors();

                auto const actual_step=
                    static_cast<StepSizeType>(flow.domain()[flow.argument_size()-1u].upper_bound());
                if(actual_step!=chen_step) {
                    chen_completed=false;
                    chen_failed_step=step_index;
                    std::stringstream msg;
                    msg << "returned step " << actual_step << " instead of " << chen_step;
                    chen_failure=msg.str();
                    break;
                }

                LabelledEnclosure reach_enclosure=chen_enclosure;
                chen_sw.restart();
                reach_enclosure.apply_full_reach_step(flow);
                chen_sw.click();
                chen_reach_us+=chen_sw.duration().count();

                auto const& reach_taylor =
                    dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(
                        reach_enclosure.state_function().reference());
                SizeType reach_nnz=0u;
                for(SizeType i=0; i!=reach_taylor.size(); ++i) {
                    reach_nnz+=reach_taylor[i].number_of_nonzeros();
                }
                if(reach_nnz>chen_max_reach_nnz) { chen_max_reach_nnz=reach_nnz; }
                auto const reach_errors=reach_enclosure.state_function().errors();

                chen_sw.restart();
                chen_enclosure.apply_fixed_evolve_step(flow,chen_step);
                chen_sw.click();
                chen_evolve_us+=chen_sw.duration().count();

                auto const evolved_errors=chen_enclosure.state_function().errors();
                auto const evolved_box=chen_enclosure.euclidean_set().bounding_box();
                double evolved_width=0.0;
                for(SizeType i=0; i!=evolved_box.size(); ++i) {
                    double const component_width=evolved_box[i].width().get_d();
                    if(component_width>evolved_width) { evolved_width=component_width; }
                }
                std::cerr << "[ChenStep]"
                          << " step=" << (step_index+1u)
                          << " time=" << ((step_index+1u)*0.02)
                          << " reconditioned=" << did_recondition
                          << " errors_before_recondition=" << errors_before_recondition
                          << " errors_after_recondition=" << errors_after_recondition
                          << " state_errors_before_flow=" << state_errors_before_flow
                          << " flow_errors=" << flow_errors
                          << " reach_errors=" << reach_errors
                          << " evolved_errors=" << evolved_errors
                          << " evolved_width=" << evolved_width
                          << " evolved_box=" << evolved_box
                          << std::endl;
            } catch(const std::exception& e) {
                chen_completed=false;
                chen_failed_step=step_index;
                chen_failure=e.what();
                break;
            }
        }

        auto const chen_final_box=chen_enclosure.euclidean_set().bounding_box();
        double chen_final_width=0.0;
        for(SizeType i=0; i!=chen_final_box.size(); ++i) {
            double const component_width=chen_final_box[i].width().get_d();
            if(component_width>chen_final_width) { chen_final_width=component_width; }
        }

        std::cerr << "[ChenFixedStep]"
                  << " completed=" << chen_completed
                  << " requested_step=" << chen_step
                  << " completed_steps=" << (chen_completed ? chen_num_steps : chen_failed_step)
                  << " failed_step=" << (chen_completed ? -1 : static_cast<long long>(chen_failed_step))
                  << " reconditionings=" << chen_reconditionings
                  << " max_state_nnz=" << chen_max_state_nnz
                  << " max_reach_nnz=" << chen_max_reach_nnz
                  << " flow_us=" << chen_flow_us
                  << " reach_us=" << chen_reach_us
                  << " evolve_us=" << chen_evolve_us
                  << " recondition_us=" << chen_recondition_us
                  << " final_error=" << chen_enclosure.state_function().error()
                  << " final_width=" << chen_final_width
                  << " final_box=" << chen_final_box;
        if(not chen_completed) {
            std::cerr << " failure=\"" << chen_failure << "\"";
        }
        std::cerr << std::endl;
    }

    // Compare graded Taylor-Picard cutoff values on the full trajectory using
    // synchronous enclosure propagation.  Each case is capped in measured work
    // so that cutoff=0 cannot make the benchmark impractically long.
    // Chen-aligned state domain and integration parameters are used here,
    // but this sweep is intentionally adaptive and is therefore a performance
    // diagnostic, not a direct reproduction of Chen's fixed-step experiment.
    const double cutoff_values[] = {1e-14,1e-13,1e-12,1e-11,1e-10};
    const long long cutoff_work_limit_us = 15000000;
    for(double cutoff_value : cutoff_values) {
        GradedTaylorPicardIntegrator sweep_integrator(
            max_err,order=5,StepSweepThreshold(ApproximateDouble(cutoff_value)));
        LabelledEnclosure sweep_enclosure(
            initial_set.euclidean_set(dynamics.state_space()),dynamics.state_space(),
            EnclosureConfiguration(sweep_integrator.function_factory()));
        sweep_enclosure.set_auxiliary(dynamics.auxiliary_space(),dynamics.auxiliary_mapping());

        TimeStepType sweep_time(0u);
        Nat sweep_steps=0u;
        Nat sweep_reconditionings=0u;
        SizeType max_state_nnz=0u;
        SizeType max_reach_nnz=0u;
        long long flow_us_total=0;
        long long reach_us_total=0;
        long long evolve_us_total=0;
        long long recondition_us_total=0;
        Stopwatch<Microseconds> sweep_sw;

        while(possibly(sweep_time < TimeStepType(7u)) && sweep_steps<2000u
              && flow_us_total+reach_us_total+evolve_us_total+recondition_us_total < cutoff_work_limit_us) {
            if(possibly(sweep_enclosure.state_function().error() > 1e-6_pr)) {
                sweep_sw.restart();
                sweep_enclosure.recondition();
                sweep_sw.click();
                recondition_us_total+=sweep_sw.duration().count();
                ++sweep_reconditionings;
            }

            auto const& state_taylor =
                dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(
                    sweep_enclosure.state_function().reference());
            SizeType state_nnz=0u;
            for(SizeType i=0; i!=state_taylor.size(); ++i) {
                state_nnz+=state_taylor[i].number_of_nonzeros();
            }
            if(state_nnz>max_state_nnz) { max_state_nnz=state_nnz; }

            auto box=cast_exact_box(sweep_enclosure.euclidean_set().bounding_box());
            sweep_sw.restart();
            auto flow=sweep_integrator.flow_step(dynamics.function(),box,suggest(StepSizeType(0.02_dy)));
            sweep_sw.click();
            flow_us_total+=sweep_sw.duration().count();

            StepSizeType actual_step=
                static_cast<StepSizeType>(flow.domain()[flow.argument_size()-1u].upper_bound());

            LabelledEnclosure reach_enclosure=sweep_enclosure;
            sweep_sw.restart();
            reach_enclosure.apply_full_reach_step(flow);
            sweep_sw.click();
            reach_us_total+=sweep_sw.duration().count();
            auto const& reach_taylor =
                dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(
                    reach_enclosure.state_function().reference());
            SizeType reach_nnz=0u;
            for(SizeType i=0; i!=reach_taylor.size(); ++i) {
                reach_nnz+=reach_taylor[i].number_of_nonzeros();
            }
            if(reach_nnz>max_reach_nnz) { max_reach_nnz=reach_nnz; }

            sweep_sw.restart();
            sweep_enclosure.apply_fixed_evolve_step(flow,actual_step);
            sweep_sw.click();
            evolve_us_total+=sweep_sw.duration().count();

            sweep_time+=TimeStepType(actual_step);
            ++sweep_steps;
        }

        bool const completed=not possibly(sweep_time < TimeStepType(7u));
        auto const final_box=sweep_enclosure.euclidean_set().bounding_box();
        double final_width=0.0;
        for(SizeType i=0; i!=final_box.size(); ++i) {
            double const component_width=final_box[i].width().get_d();
            if(component_width>final_width) { final_width=component_width; }
        }
        std::cerr << "[CutoffSweep]"
                  << " cutoff=" << cutoff_value
                  << " completed=" << completed
                  << " t=" << sweep_time
                  << " steps=" << sweep_steps
                  << " reconditionings=" << sweep_reconditionings
                  << " max_state_nnz=" << max_state_nnz
                  << " max_reach_nnz=" << max_reach_nnz
                  << " flow_us=" << flow_us_total
                  << " reach_us=" << reach_us_total
                  << " evolve_us=" << evolve_us_total
                  << " recondition_us=" << recondition_us_total
                  << " final_error=" << sweep_enclosure.state_function().error()
                  << " final_width=" << final_width
                  << " final_box=" << final_box
                  << std::endl;
    }

    // The standalone first-step cutoff sweep is intentionally disabled here.
    // Cutoffs below Chen's 1e-12 can make a single Taylor-model composition
    // dominate the run without adding information to the fixed-step comparison.

    std::cerr << "[vanderpol] starting graded Taylor-Picard evolution" << std::endl;
    auto evolution = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Done in " << sw.elapsed_seconds() << " seconds.");

    CONCLOG_PRINTLN("Plotting...");
    fig.clear();
    fig.draw(evolution);
    fig.write("vanderpol_evolution");
}
