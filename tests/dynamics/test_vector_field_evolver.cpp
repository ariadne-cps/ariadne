/***************************************************************************
 *            test_continuous_evolution.cpp
 *
 *  Copyright  2006-20  Pieter Collins
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
#include "dynamics/vector_field_evolver.hpp"
#include "dynamics/orbit.hpp"
#include "output/graphics.hpp"
#include "output/logging.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

class TestContinuousEvolution
{
  public:
    Void test() const;
    Void failure_test() const;
};

Int main()
{
    ARIADNE_TEST_CALL(TestContinuousEvolution().test());
    //ARIADNE_TEST_CALL(TestContinuousEvolution().failure_test());
    return ARIADNE_TEST_FAILURES;
}


Void TestContinuousEvolution::test() const
{
    // cout << __PRETTY_FUNCTION__ << endl;

    typedef VectorField::EnclosureType EnclosureType;

    // Set up the evolution parameters and grid
    Real time(2.0_dec);
    ExactDouble step_size(0.5_x);
    ExactDouble enclosure_radius(0.25_x);

    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),1e-8_pr);

    // Set up the evaluators
    TaylorPicardIntegrator picard_integrator(maximum_error=1e-4_pr,sweeper,lipschitz_constant=0.5_x,
                                             step_maximum_error=1e-6_pr,minimum_temporal_order=0,maximum_temporal_order=8);
    // Set up the evaluators
    GradedTaylorSeriesIntegrator series_integrator(maximum_error=1e-4_pr,sweeper,lipschitz_constant=0.5_x,step_maximum_error=1e-6_pr,
                                             minimum_spacial_order=1,minimum_temporal_order=4,maximum_spacial_order=3,maximum_temporal_order=8);

    IntegratorInterface& integrator=picard_integrator;

    ARIADNE_TEST_PRINT(integrator);

    // cout << "initial_box=" << initial_box << endl;

    // Set up the vector field
    Real mu=Dyadic(0.5_x);
    RealVariable x("x"), v("v");

    VectorField vanderpol({dot(x)=v,dot(v)=mu*(1-x*x)*v-x});
    ARIADNE_TEST_PRINT(vanderpol);

    // Define the initial set
    RealExpressionBoundedConstraintSet initial_set({1.01_dec<=x<=1.02_dec,0.51_dec<=v<=0.52_dec});

    VectorFieldEvolver evolver(vanderpol,integrator);
    evolver.configuration().set_maximum_enclosure_radius(enclosure_radius);
    evolver.configuration().set_maximum_step_size(step_size);

    // Over-approximate the initial set by a grid cell
    TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp,1e-8));
    EnclosureType initial_encl(initial_set.euclidean_set(vanderpol.state_space()),vanderpol.state_space(),EnclosureConfiguration(function_factory));
    ARIADNE_TEST_PRINT(initial_encl);

    Semantics semantics=Semantics::UPPER;

    // Compute the reachable sets
    Orbit<EnclosureType> orbit = evolver.orbit(initial_encl,time,semantics);
    ARIADNE_TEST_PRINT(orbit);

    LabelledFigure fig(Axes2d(-1.0<=x<=+21,-1.125<=v<=+1.125));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_encl;
    fig.write("test_vector_field_evolver-vdp");
}

Void TestContinuousEvolution::failure_test() const
{
    // The systems in this test are stiff and are expected to fail.

    typedef VectorField::EnclosureType EnclosureType;

    // Set up the evolution parameters and grid
    Real time(0.5_dec);
    ExactDouble step_size(0.01_pr);
    ExactDouble enclosure_radius(0.25_x);

    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),1e-8_pr);

    // Set up the evaluators
    TaylorPicardIntegrator integrator(maximum_error=1e-6_pr,sweeper,lipschitz_constant=0.5_x,
                                      step_maximum_error=1e-8_pr,minimum_temporal_order=0,maximum_temporal_order=6);

    // Define the initial box
    ExactBoxType initial_box = ExactBoxType{{0.0_x,0.0_x},{0.9_pr,0.9_pr}};

    // cout << "initial_box=" << initial_box << endl;

    // Set up the vector field for the first test
    Real p = 200;
    RealVariable x("x"), y("y");

    VectorField failone_vf({dot(x)=1,dot(y)=-p*y+p});

    VectorFieldEvolver evolverone(failone_vf,integrator);
    evolverone.configuration().set_maximum_enclosure_radius(enclosure_radius);
    evolverone.configuration().set_maximum_step_size(step_size);

    TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp,1e-10));
    EnclosureType initial_set(initial_box,failone_vf.state_space(),EnclosureConfiguration(function_factory));
    // cout << "initial_set=" << initial_set << endl << endl;

    Semantics semantics=Semantics::UPPER;

    // Compute the reachable sets
    Orbit<EnclosureType> orbit = evolverone.orbit(initial_set,time,semantics);

    cout << "\norbit.final=\n" << orbit.final() << endl << endl;
    cout << "final set radius="<< orbit.final()[0].radius() << endl;

    ARIADNE_TEST_COMPARE(orbit.final()[0].radius(),>,0.5);

    // Print the intial, evolve and reach sets
    // cout << "Plotting sets" << endl;
    // cout << "evolve_set=" << hybrid_evolve_set << endl;
    // cout << "reach_set=" << hybrid_reach_set << endl;
    std::cout << "Plotting..."<< std::flush;
    LabelledFigure fig(Variables2d(x,y),widen(orbit.reach().bounding_box(),+0.25_x));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("test_vector_field_evolver-failone");
/*
    Figure fig; fig.set_bounding_box( widen(orbit.reach().bounding_box(),+0.25_x));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach().euclidean_set();
    fig << fill_colour(magenta) << orbit.intermediate().euclidean_set();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("test_continuous_evolution-failone");
*/
    // Set up the vector field for the second test
    p = Rational(1,10);
    EffectiveScalarMultivariateFunction z3=EffectiveScalarMultivariateFunction::zero(3);
    EffectiveScalarMultivariateFunction o3=EffectiveScalarMultivariateFunction::constant(3,1.0_q);
    EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(3,0);
    EffectiveScalarMultivariateFunction x1=EffectiveScalarMultivariateFunction::coordinate(3,1);
    EffectiveScalarMultivariateFunction x2=EffectiveScalarMultivariateFunction::coordinate(3,1);
    EffectiveVectorMultivariateFunction failtwo={o3,x1*x2/p,z3};
    VectorField failtwo_vf(failtwo);

    VectorFieldEvolver evolvertwo(failtwo_vf,integrator);
    evolvertwo.configuration().set_maximum_enclosure_radius(enclosure_radius);
    evolvertwo.configuration().set_maximum_step_size(step_size);

    ExactBoxType initial_box2 = ExactBoxType{{0.0_x,0.0_x},{1.0_x,1.0_x},{1.0_x,1.0_x}};
    initial_set = EnclosureType(initial_box2,failtwo_vf.state_space(),EnclosureConfiguration(function_factory));

    time = 1.5_dec;

    // Compute the reachable sets
    orbit = evolvertwo.orbit(initial_set,time,semantics);

    cout << "\norbit.final=\n" << orbit.final() << endl << endl;
    cout << "final set radius="<< orbit.final()[0].radius() << endl;

    ARIADNE_TEST_COMPARE(orbit.final()[0].radius(),>,1.0);

    // Print the intial, evolve and reach sets
    // cout << "Plotting sets" << endl;
    // cout << "evolve_set=" << hybrid_evolve_set << endl;
    // cout << "reach_set=" << hybrid_reach_set << endl;
    std::cout << "Plotting..." << std::flush;
    fig.clear(); fig.set_bounding_box(orbit.reach().bounding_box());
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("test_vector_field_evolver-failtwo");
    std::cout << std::endl;
}

