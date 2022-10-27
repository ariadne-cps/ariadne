/***************************************************************************
 *            test_vector_field_evolver.cpp
 *
 *  Copyright  2006-21  Pieter Collins
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

#include <fstream>
#include <iostream>

#include "config.hpp"
#include "utility/tuple.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "function/constraint.hpp"
#include "dynamics/enclosure.hpp"
#include "geometry/box.hpp"
#include "geometry/list_set.hpp"
#include "solvers/integrator.hpp"
#include "symbolic/expression_set.hpp"
#include "dynamics/orbit.hpp"
#include "dynamics/vector_field_evolver.hpp"
#include "io/figure.hpp"
#include "io/command_line_interface.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

class TestVectorFieldEvolver {
public:
    Void test() const {
        ARIADNE_TEST_CALL(test_single_trajectory());
        ARIADNE_TEST_CALL(test_failure());
        ARIADNE_TEST_CALL(test_subdivide_initially());
        ARIADNE_TEST_CALL(test_subdivide_along_evolution());
    }

    Void test_single_trajectory() const {

        typedef VectorField::EnclosureType EnclosureType;

        // Set up the evolution parameters and grid
        Real time(2.0_dec);
        ExactDouble step_size(0.5_x);
        ExactDouble enclosure_radius(0.25_x);

        ThresholdSweeper<FloatDP> sweeper(DoublePrecision(), 1e-8_pr);

        // Set up the evaluators
        TaylorPicardIntegrator picard_integrator(step_maximum_error = 1e-7_pr, sweeper, lipschitz_tolerance = 0.5_x,
                                                 minimum_temporal_order = 0, maximum_temporal_order = 8);
        // Set up the evaluators
        GradedTaylorSeriesIntegrator series_integrator(step_maximum_error = 1e-7_pr, sweeper, lipschitz_tolerance = 0.5_x,
                                                       minimum_spacial_order = 1, minimum_temporal_order = 4,
                                                       maximum_spacial_order = 3, maximum_temporal_order = 8);

        IntegratorInterface &integrator = picard_integrator;

        ARIADNE_TEST_PRINT(integrator);

        // Set up the vector field
        Real mu = Dyadic(0.5_x);
        RealVariable x("x"), v("v");

        VectorField vanderpol({dot(x) = v, dot(v) = mu * (1 - x * x) * v - x});
        ARIADNE_TEST_PRINT(vanderpol);

        VectorFieldEvolver evolver(vanderpol, integrator);
        evolver.configuration().set_maximum_enclosure_radius(enclosure_radius);
        evolver.configuration().set_maximum_step_size(step_size);
        ARIADNE_TEST_PRINT(evolver.configuration());

        // Define the initial set
        RealExpressionBoundedConstraintSet initial_set({1.01_dec <= x <= 1.02_dec, 0.51_dec <= v <= 0.52_dec});
            initial_set=RealExpressionBoundedConstraintSet({1.01_dec <= x <= 1.02_dec, 0.5_dec <= v <= 0.5_dec});
            time=0.25_dec;

        Semantics semantics = Semantics::LOWER;

        // Compute the reachable sets
        Orbit<EnclosureType> orbit = evolver.orbit(initial_set, time, semantics);
        ARIADNE_TEST_PRINT(orbit);

        LabelledFigure fig(Axes2d(-1.0 <= x <= +21, -1.125 <= v <= +1.125));
        fig << line_style(true) << fill_colour(cyan) << orbit.reach();
        fig << fill_colour(magenta) << orbit.intermediate();
        fig << fill_colour(red) << orbit.final();

        // Define the initial set
        RealExpressionBoundedConstraintSet initial_point({1.0_dec <= x <= 1.0_dec, 0.5_dec <= v <= 0.5_dec});
            initial_point=RealExpressionBoundedConstraintSet({1.01_dec <= x <= 1.02_dec, 0.5_dec <= v <= 0.5_dec});
            time=0.25_dec;
        Orbit<EnclosureType> orbit = evolver.orbit(initial_set, time, semantics);

        fig.write("test_vector_field_evolver-vdp");
    }

    Void test_failure() const {
        // The system in this test is stiff and is expected to fail.

        RealVariable x("x"), y("y");

        // Set up the evolution parameters and grid
        Real time(0.5_dec);
        ExactDouble maximum_step_size(0.015625_pr);
        ExactDouble minimum_step_size(0.00097656_pr);
        ExactDouble enclosure_radius(0.25_x);


        ThresholdSweeper<FloatDP> sweeper(DoublePrecision(), 1e-8_pr);

        // Set up the evaluators
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-6_pr, sweeper, lipschitz_tolerance = 0.5_x,
                                          minimum_temporal_order = 0, maximum_temporal_order = 8);

        VectorField fail_vf({dot(x)=1,dot(y)=y*y*100});
        VectorFieldEvolver evolver(fail_vf, integrator);
        evolver.configuration().set_maximum_enclosure_radius(enclosure_radius);
        evolver.configuration().set_maximum_step_size(maximum_step_size);
//        evolver.configuration().set_minimum_step_size(minimum_step_size);

        RealVariablesBox initial_box({x==0,y==1});

        time = 1.5_dec;

        // Compute the reachable sets
        evolver.orbit(initial_box, time, Semantics::UPPER);
//        ARIADNE_TEST_FAIL(evolver.orbit(initial_box, time, Semantics::UPPER));
    }

    Void test_subdivide_initially() const {

        RealConstant mu("mu",1.0_dec);
        RealVariable x("x"), y("y");

        VectorField dynamics({dot(x)=y, dot(y)= mu*y*(1-sqr(x))-x});

        StepMaximumError max_err=1e-8;

        TaylorPicardIntegrator integrator(max_err);

        VectorFieldEvolver evolver(dynamics,integrator);
        evolver.configuration().set_maximum_enclosure_radius(0.1);
        evolver.configuration().set_maximum_step_size(0.1);
        evolver.configuration().set_maximum_spacial_error(1e-5);
        evolver.configuration().set_enable_subdivisions(true);

        Real x0(1.40_dec);
        Real y0(2.40_dec);
        Real eps_x0 = 15/100_q;
        Real eps_y0 = 5/100_q;

        RealExpressionBoundedConstraintSet initial_set({x0-eps_x0<=x<=x0+eps_x0,y0-eps_y0<=y<=y0+eps_y0});

        Real evolution_time(0.2_dec);

        auto evolution = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);

        ARIADNE_TEST_ASSERT(evolution.final().size()>1);
    }

    Void test_subdivide_along_evolution() const {

        RealVariable x("x"), y("y");

        VectorField dynamics({dot(x)=1, dot(y)= sqr(x)});

        StepMaximumError max_err=1e-8;

        TaylorPicardIntegrator integrator(max_err);

        VectorFieldEvolver evolver(dynamics,integrator);
        evolver.configuration().set_maximum_enclosure_radius(0.1);
        evolver.configuration().set_maximum_step_size(0.1);
        evolver.configuration().set_maximum_spacial_error(1e-5);
        evolver.configuration().set_enable_subdivisions(true);

        Real x0(1.40_dec);
        Real y0(2.40_dec);
        Real eps_x0 = 5/100_q;
        Real eps_y0 = 5/100_q;

        RealExpressionBoundedConstraintSet initial_set({x0-eps_x0<=x<=x0+eps_x0,y0-eps_y0<=y<=y0+eps_y0});

        Real evolution_time(3);

        auto orbit_upper = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);

        LabelledFigure fig(Axes2d(1 <= x <= 5, 2 <= y <= 30));
        fig << orbit_upper.reach();
        fig.write("test_vector_field_evolver-subdivisions");

        auto orbit_lower = evolver.orbit(initial_set,evolution_time,Semantics::LOWER);

        ARIADNE_TEST_ASSERT(orbit_upper.reach().size() > orbit_lower.reach().size());
    }
};

Int main(Int argc, const char* argv[])
{
    if (not CommandLineInterface::instance().acquire(argc,argv)) return -1;

    TestVectorFieldEvolver().test();
    return ARIADNE_TEST_FAILURES;
}
