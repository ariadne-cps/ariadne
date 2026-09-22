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

    // Compare graded Taylor-Picard cutoff values on the full trajectory using
    // synchronous enclosure propagation.  Each case is capped in measured work
    // so that cutoff=0 cannot make the benchmark impractically long.
    const double cutoff_values[] = {0.0,1e-14,1e-13,1e-12,1e-11,1e-10};
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
            max_state_nnz=max(max_state_nnz,state_nnz);

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
            max_reach_nnz=max(max_reach_nnz,reach_nnz);

            sweep_sw.restart();
            sweep_enclosure.apply_fixed_evolve_step(flow,actual_step);
            sweep_sw.click();
            evolve_us_total+=sweep_sw.duration().count();

            sweep_time+=TimeStepType(actual_step);
            ++sweep_steps;
        }

        bool const completed=not possibly(sweep_time < TimeStepType(7u));
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
                  << " final_box=" << sweep_enclosure.euclidean_set().bounding_box()
                  << std::endl;
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
