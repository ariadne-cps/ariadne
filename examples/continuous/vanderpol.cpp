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
    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),1e-12);

    GradedTaylorSeriesIntegrator integrator(
        max_err,sweeper,lipschitz_tolerance=0.5_x,
        minimum_spacial_order=5,minimum_temporal_order=5,
        maximum_spacial_order=5,maximum_temporal_order=5);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.02);
    evolver.configuration().set_maximum_spacial_error(1e-6);
    evolver.configuration().set_enable_reconditioning(false);

    Real x0 = 1.40_dec;
    Real y0 = 2.30_dec;
    Real eps_x0 = 0.15_dec;
    Real eps_y0 = 0.05_dec;
    RealExpressionBoundedConstraintSet initial_set({
        x0-eps_x0<=x<=x0+eps_x0,
        y0-eps_y0<=y<=y0+eps_y0
    });

    auto orbit=evolver.orbit(initial_set,Real(0.02_dec),Semantics::UPPER);
    std::cerr << "[GradedRemainderDiagnosticSummary]"
              << " reach_sets=" << orbit.reach().size()
              << " intermediate_sets=" << orbit.intermediate().size();
    if(!orbit.final().empty()) {
        auto const& final_set=orbit.final()[0];
        std::cerr << " final_error=" << final_set.state_function().error()
                  << " final_radius=" << final_set.radius()
                  << " final_box=" << final_set.euclidean_set().bounding_box();
    }
    std::cerr << std::endl;

}
