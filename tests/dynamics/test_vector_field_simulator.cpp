/***************************************************************************
 *            test_vector_field_simulator.cpp
 *
 *  Copyright  2006-20  Luca Geretti
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
#include "function/constraint.hpp"
#include "dynamics/enclosure.hpp"
#include "geometry/box.hpp"
#include "geometry/list_set.hpp"
#include "solvers/integrator.hpp"
#include "symbolic/expression_set.hpp"
#include "dynamics/vector_field_simulator.hpp"
#include "dynamics/orbit.hpp"
#include "output/graphics.hpp"
#include "output/logging.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

class TestContinuousEvolution
{
  public:
    Void test() const {
        ARIADNE_TEST_CALL(test_one_vdp_cycle());
    }

    Void test_one_vdp_cycle() const {

        typedef VectorFieldSimulator::ApproximatePointType PointType;

        Real mu=Dyadic(0.5_x);
        RealVariable x("x"), v("v");

        VectorField vanderpol({dot(x)=v,dot(v)=mu*(1-x*x)*v-x});
        ARIADNE_TEST_PRINT(vanderpol);

        // Define the initial point
        RealPoint initial_point({-1.5_dec,1.0_dec});
        LabelledRealPoint pt2({x=-1.5_dec,v=1});

        Real time = 6.3_dec;

        VectorFieldSimulator simulator(vanderpol);
        simulator.configuration().set_step_size(0.05);

        ARIADNE_TEST_PRINT(simulator.configuration());

        Orbit<PointType> orbit = simulator.orbit(initial_point,time);

        ARIADNE_TEST_ASSERT(distance(orbit.curve().end()->second,PointType(initial_point,double_precision)).raw() <= 0.02);

        GraphicsBoundingBoxType bx({FloatDPApproximateInterval(-2.5,2.5),FloatDPApproximateInterval(-3,3)});
        Figure fig(bx,Projection2d(2,0,1));
        fig << orbit.curve();
        fig.write("test_vector_field_simulator");
    }
};

Int main()
{
    TestContinuousEvolution().test();
    return ARIADNE_TEST_FAILURES;
}


/*
Void TestContinuousEvolution::test() const
{
    typedef VectorField::EnclosureType EnclosureType;

    // Set up the evolution parameters and grid
    Real time(2.0_dec);
    ExactDouble step_size(0.5_x);
    ExactDouble enclosure_radius(0.25_x);

    // Set up the vector field
    Real mu=Dyadic(0.5_x);
    RealVariable x("x"), v("v");

    VectorField vanderpol({dot(x)=v,dot(v)=mu*(1-x*x)*v-x});
    ARIADNE_TEST_PRINT(vanderpol);

    // Define the initial box
    RealVariablesBox initial_box({x.in(1.01_dec,1.02_dec),v.in(0.51_dec,0.52_dec)});

    VectorFieldSimulator simulator();
    simulator.set_step_size(0.01);

    // Over-approximate the initial set by a grid cell
    TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp,1e-8));
    EnclosureType initial_set(initial_box,vanderpol.state_space(),EnclosureConfiguration(function_factory));
    ARIADNE_TEST_PRINT(initial_set);

    Semantics semantics=Semantics::UPPER;

    // Compute the reachable sets
    Orbit<EnclosureType> orbit = evolver.orbit(initial_set,time,semantics);
    ARIADNE_TEST_PRINT(orbit);

    LabelledFigure fig(Axes2d(-1.0<=x<=+21,-1.125<=v<=+1.125));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("test_vector_field_simulator");
}

*/