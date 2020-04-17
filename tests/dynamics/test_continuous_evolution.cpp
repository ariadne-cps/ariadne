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
#include "dynamics/vector_field_evolver.hpp"
#include "dynamics/orbit.hpp"
#include "output/graphics.hpp"
#include "output/logging.hpp"

#include "function/user_function.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

/// This function diverges heavily
struct FailOne : VectorMultivariateFunctionData<2,2,1> {
    template<class R, class A, class P> static Void
    compute(R& r, const A& x, const P& p) {
          r[0] = 1;
          r[1] = -p[0] * x[1] + p[0];
    }
};

/// This function diverges heavily
struct FailTwo : VectorMultivariateFunctionData<3,3,1> {
    template<class R, class A, class P> static Void
    compute(R& r, const A& x, const P& p) {
          r[0] = 1;
          r[1] = x[1] * x[2] / p[0];
          r[2] = 0;
    }
};


class TestContinuousEvolution
{
  public:
    Void test() const;
    Void failure_test() const;
};

Int main()
{
    //std::cerr<<"SKIPPED "; return 1;
    ARIADNE_TEST_CALL(TestContinuousEvolution().test());
    //ARIADNE_TEST_CALL(TestContinuousEvolution().failure_test());
    return ARIADNE_TEST_FAILURES;
}



Void TestContinuousEvolution::test() const
{
    // cout << __PRETTY_FUNCTION__ << endl;

    typedef Enclosure EnclosureType;

    // Set up the evolution parameters and grid
    Real time(2.0_dec);
    double step_size(0.5);
    double enclosure_radius(0.25);

    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),1e-8);

    // Set up the evaluators
    TaylorPicardIntegrator picard_integrator(maximum_error=1e-4,sweeper,lipschitz_constant=0.5,
                                             step_maximum_error=1e-6,minimum_temporal_order=0,maximum_temporal_order=8);
    // Set up the evaluators
    GradedTaylorSeriesIntegrator series_integrator(maximum_error=1e-4,sweeper,lipschitz_constant=0.5,step_maximum_error=1e-6,
                                             minimum_spacial_order=1,minimum_temporal_order=4,maximum_spacial_order=3,maximum_temporal_order=8);

    IntegratorInterface& integrator=picard_integrator;

    ARIADNE_TEST_PRINT(integrator);

    // Define the initial box
    ExactBoxType initial_box(2);
    initial_box[0]=ExactIntervalType(1.01,1.02);
    initial_box[1]=ExactIntervalType(0.51,0.52);

    // cout << "initial_box=" << initial_box << endl;

    // Set up the vector field
    Real mu=Dyadic(0.5);
    EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(2,0);
    EffectiveScalarMultivariateFunction xp=EffectiveScalarMultivariateFunction::coordinate(2,1);
    EffectiveVectorMultivariateFunction vdp={x,mu*(1-x*x)*xp-x};

    VectorField vanderpol(vdp);
    ARIADNE_TEST_PRINT(vanderpol);

    VectorFieldEvolver evolver(vanderpol,integrator);
    evolver.configuration().set_maximum_enclosure_radius(enclosure_radius);
    evolver.configuration().set_maximum_step_size(step_size);

    // Over-approximate the initial set by a grid cell
    TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp,1e-8));
    EnclosureType initial_set(initial_box,function_factory);
    ARIADNE_TEST_PRINT(initial_set);

    Semantics semantics=Semantics::UPPER;

    // Compute the reachable sets
    Orbit<EnclosureType> orbit = evolver.orbit(initial_set,time,semantics);
    ARIADNE_TEST_PRINT(orbit);

    // Print the intial, evolve and reach sets
    // cout << "Plotting sets" << endl;
    // cout << "evolve_set=" << hybrid_evolve_set << endl;
    // cout << "reach_set=" << hybrid_reach_set << endl;
    Figure fig;
    fig.set_bounding_box(ExactBoxType{{-1.0,+21.0},{-1.125,+1.125}});
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("test_continuous_evolution-vdp");

}

Void TestContinuousEvolution::failure_test() const
{
    // The systems in this test are stiff and are expected to fail.

    // cout << __PRETTY_FUNCTION__ << endl;

    typedef Enclosure EnclosureType;

    // Set up the evolution parameters and grid
    Real time(0.5_dec);
    double step_size(0.01);
    double enclosure_radius(0.25);

    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),1e-8);

    // Set up the evaluators
    TaylorPicardIntegrator integrator(maximum_error=1e-6,sweeper,lipschitz_constant=0.5,
                                      step_maximum_error=1e-8,minimum_temporal_order=0,maximum_temporal_order=6);

    // Define the initial box
    ExactBoxType initial_box = ExactBoxType{{0.0,0.0},{0.9,0.9}};

    // cout << "initial_box=" << initial_box << endl;

    // Set up the vector field for the first test
    Real p = 200;
    EffectiveScalarMultivariateFunction o=EffectiveScalarMultivariateFunction::constant(2,1.0_q);
    EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(2,0);
    EffectiveScalarMultivariateFunction y=EffectiveScalarMultivariateFunction::coordinate(2,1);
    EffectiveVectorMultivariateFunction failone={o,-p*y+p};
    VectorField failone_vf(failone);

    VectorFieldEvolver evolverone(failone_vf,integrator);
    evolverone.configuration().set_maximum_enclosure_radius(enclosure_radius);
    evolverone.configuration().set_maximum_step_size(step_size);

    TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp,1e-10));
    EnclosureType initial_set(initial_box,function_factory);
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
    Figure fig; fig.set_bounding_box( widen(orbit.reach().bounding_box(),+0.25_x));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("test_continuous_evolution-failone");

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

    ExactBoxType initial_box2 = ExactBoxType{{0.0,0.0},{1.0,1.0},{1.0,1.0}};
    initial_set = EnclosureType(initial_box2,function_factory);

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
    fig.set_projection(2,0,1);
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("test_continuous_evolution-failtwo");
    std::cout << std::endl;

}

