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
#include "function/affine_model.hpp"
#include "function/constraint.hpp"
#include "geometry/zonotope.hpp"
#include "dynamics/enclosure.hpp"
#include "ariadne_main.hpp"


namespace {

LabelledEnclosure affine_precondition(const LabelledEnclosure& enclosure)
{
    auto const& state_taylor =
        dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(
            enclosure.state_function().reference());

    SizeType const state_dimension=state_taylor.size();
    SizeType const parameter_dimension=state_taylor.argument_size();

    Vector<FloatDP> centre(state_dimension,FloatDP(dp));
    Matrix<FloatDP> generators(state_dimension,parameter_dimension,FloatDP(dp));
    Vector<FloatDP> error(state_dimension,FloatDP(dp));

    for(SizeType i=0u; i!=state_dimension; ++i) {
        AffineModel<ValidatedTag,FloatDP> affine(state_taylor.model(i));
        centre[i]=affine.value();
        for(SizeType j=0u; j!=parameter_dimension; ++j) {
            generators[i][j]=affine.gradient(j);
        }
        error[i]=FloatDP(affine.error().raw());
    }

    // Convert nonlinear/uniform errors into explicit generators, then retain
    // one correlated generator block and one residual block.  For a
    // two-dimensional state this yields c + A*y + r with four normalized
    // parameters: two correlated affine directions and two independent
    // residual directions.
    Zonotope affine_set=error_free_over_approximation(
        Zonotope(centre,generators,error));
    Zonotope preconditioned=cascade_over_approximation(affine_set,2u);

    ExactBoxType domain(
        preconditioned.number_of_generators(),ExactIntervalType(-1,+1));
    auto const& factory=enclosure.configuration().function_factory();
    ValidatedVectorMultivariateFunctionPatch state=
        factory.create_zeros(state_dimension,domain);

    for(SizeType i=0u; i!=state_dimension; ++i) {
        state[i]=factory.create_constant(domain,preconditioned.centre()[i]);
        for(SizeType j=0u; j!=preconditioned.number_of_generators(); ++j) {
            state[i]=state[i]
                + factory.create_coordinate(domain,j)
                  * FloatDPBounds(preconditioned.generators()[i][j]);
        }
    }

    LabelledEnclosure result(
        Enclosure(domain,state,enclosure.configuration()),
        enclosure.state_space());
    result.set_auxiliary(enclosure.auxiliary_space(),enclosure.auxiliary_mapping());
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

    // Controlled comparison at identical physical reconditioning times.
    // The standard path uses Enclosure::recondition(); the affine path replaces
    // the current representation by c+A*y+r using the affine part of the current
    // Taylor model and explicit residual generators.
    {
        const StepSizeType benchmark_steps[] = {
            StepSizeType(0.005_dy), StepSizeType(0.0025_dy)
        };
        const Nat benchmark_num_steps[] = {1400u,2800u};
        const Nat benchmark_period_steps[] = {8u,16u}; // 0.04 physical time

        for(SizeType method=0u; method!=2u; ++method) {
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
                Nat preconditionings=0u;
                SizeType max_parameters=benchmark_enclosure.number_of_parameters();
                SizeType max_state_nnz=0u;
                bool completed=true;
                Nat failed_step=num_steps;
                String failure;

                Stopwatch<Microseconds> benchmark_sw;
                long long flow_us=0;
                long long evolve_us=0;
                long long precondition_us=0;

                for(Nat step_index=0u; step_index!=num_steps; ++step_index) {
                    auto box=cast_exact_box(
                        benchmark_enclosure.euclidean_set().bounding_box());
                    try {
                        benchmark_sw.restart();
                        auto flow=benchmark_integrator.flow_step(
                            dynamics.function(),box,benchmark_step);
                        benchmark_sw.click();
                        flow_us+=benchmark_sw.duration().count();

                        auto const actual_step=static_cast<StepSizeType>(
                            flow.domain()[flow.argument_size()-1u].upper_bound());
                        if(actual_step!=benchmark_step) {
                            completed=false;
                            failed_step=step_index;
                            std::stringstream msg;
                            msg << "returned step " << actual_step
                                << " instead of " << benchmark_step;
                            failure=msg.str();
                            break;
                        }

                        benchmark_sw.restart();
                        benchmark_enclosure.apply_fixed_evolve_step(flow,benchmark_step);
                        benchmark_sw.click();
                        evolve_us+=benchmark_sw.duration().count();

                        auto const& state_taylor =
                            dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(
                                benchmark_enclosure.state_function().reference());
                        SizeType state_nnz=0u;
                        for(SizeType i=0u; i!=state_taylor.size(); ++i) {
                            state_nnz+=state_taylor[i].number_of_nonzeros();
                        }
                        if(state_nnz>max_state_nnz) { max_state_nnz=state_nnz; }

                        if(((step_index+1u)%period_steps)==0u) {
                            benchmark_sw.restart();
                            if(method==0u) {
                                benchmark_enclosure.recondition();
                            } else {
                                benchmark_enclosure=affine_precondition(benchmark_enclosure);
                            }
                            benchmark_sw.click();
                            precondition_us+=benchmark_sw.duration().count();
                            ++preconditionings;
                            if(benchmark_enclosure.number_of_parameters()>max_parameters) {
                                max_parameters=benchmark_enclosure.number_of_parameters();
                            }
                        }
                    } catch(const std::exception& e) {
                        completed=false;
                        failed_step=step_index;
                        failure=e.what();
                        break;
                    }
                }

                auto const final_box=
                    benchmark_enclosure.euclidean_set().bounding_box();
                double final_width=0.0;
                std::vector<double> component_widths(final_box.size());
                for(SizeType i=0u; i!=final_box.size(); ++i) {
                    double const component_width=final_box[i].width().get_d();
                    component_widths[i]=component_width;
                    if(component_width>final_width) { final_width=component_width; }
                }

                std::cerr << (method==0u
                              ? "[StandardRecondition]"
                              : "[AffinePrecondition]")
                          << " requested_step=" << benchmark_step
                          << " interval=0.04"
                          << " completed=" << completed
                          << " completed_steps=" << (completed ? num_steps : failed_step)
                          << " preconditionings=" << preconditionings
                          << " max_parameters=" << max_parameters
                          << " max_state_nnz=" << max_state_nnz
                          << " flow_us=" << flow_us
                          << " evolve_us=" << evolve_us
                          << " precondition_us=" << precondition_us
                          << " final_error=" << benchmark_enclosure.state_function().error()
                          << " final_width=" << final_width
                          << " component_widths=["
                          << component_widths[0] << "," << component_widths[1] << "]"
                          << " final_box=" << final_box;
                if(not completed) {
                    std::cerr << " failure=\"" << failure << "\"";
                }
                std::cerr << std::endl;
            }
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
