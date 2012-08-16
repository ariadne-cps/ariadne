/***************************************************************************
 *            test_continuous_evolution.cc
 *
 *  Copyright  2006-8  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <fstream>
#include <iostream>

#include "config.h"
#include "tuple.h"
#include "vector.h"
#include "matrix.h"
#include "function.h"
#include "taylor_function.h"
#include "enclosure.h"
#include "box.h"
#include "list_set.h"
#include "integrator.h"
#include "vector_field_evolver.h"
#include "orbit.h"
#include "graphics.h"
#include "textplot.h"
#include "logging.h"

#include "user_function.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

/// This function diverges heavily
struct FailOne : VectorFunctionData<2,2,1> {
    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p) {
          r[0] = 1.0;
          r[1] = -p[0] * x[1] + p[0];
    }
};

/// This function diverges heavily
struct FailTwo : VectorFunctionData<3,3,1> {
    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p) {
          r[0] = 1.0;
          r[1] = x[1] * x[2]/p[0];
          r[2] = 0.0;
    }
};


class TestContinuousEvolution
{
  public:
    void test() const;
    void failure_test() const;
};

int main()
{
    //std::cerr<<"SKIPPED "; return 1;
    ARIADNE_TEST_CALL(TestContinuousEvolution().test());
    //ARIADNE_TEST_CALL(TestContinuousEvolution().failure_test());
    return ARIADNE_TEST_FAILURES;
}



void TestContinuousEvolution::test() const
{
    // cout << __PRETTY_FUNCTION__ << endl;

    typedef Enclosure EnclosureType;

    // Set up the evolution parameters and grid
    Float time(3.0);
    double step_size(0.5);
    double enclosure_radius(0.25);

    // Set up the evaluators
    TaylorPicardIntegrator picard_integrator(maximum_error=1e-4,sweep_threshold=1e-8,lipschitz_constant=0.5,
                                             step_maximum_error=1e-6,step_sweep_threshold=1e-10,maximum_temporal_order=8);
    // Set up the evaluators
    TaylorSeriesIntegrator series_integrator(maximum_error=1e-4,sweep_threshold=1e-8,lipschitz_constant=0.5,
                                             step_maximum_error=1e-6,step_sweep_threshold=1e-10,
                                             minimum_spacial_order=1,minimum_temporal_order=4,maximum_spacial_order=3,maximum_temporal_order=8);

    IntegratorInterface& integrator=picard_integrator;

    ARIADNE_TEST_PRINT(integrator);

    // Define the initial box
    Box initial_box(2);
    initial_box[0]=Interval(1.01,1.02);
    initial_box[1]=Interval(0.51,0.52);

    // cout << "initial_box=" << initial_box << endl;

    // Set up the vector field
    Real mu=0.5;
    RealScalarFunction x=RealScalarFunction::coordinate(2,0);
    RealScalarFunction xp=RealScalarFunction::coordinate(2,1);
    RealVectorFunction vdp((x,mu*(1-x*x)*xp-x));

    VectorField vanderpol(vdp);
    ARIADNE_TEST_PRINT(vanderpol);

    VectorFieldEvolver evolver(vanderpol,integrator);
    evolver.configuration().maximum_enclosure_radius(enclosure_radius);
    evolver.configuration().maximum_step_size(step_size);

    // Over-approximate the initial set by a grid cell
    TaylorFunctionFactory function_factory(ThresholdSweeper(1e-10));
    EnclosureType initial_set(initial_box,function_factory);
    ARIADNE_TEST_PRINT(initial_set);

    Semantics semantics=UPPER_SEMANTICS;

    // Compute the reachable sets
    Orbit<EnclosureType> orbit = evolver.orbit(initial_set,time,semantics);
    ARIADNE_TEST_PRINT(orbit);

    // Print the intial, evolve and reach sets
    // cout << "Plotting sets" << endl;
    // cout << "evolve_set=" << hybrid_evolve_set << endl;
    // cout << "reach_set=" << hybrid_reach_set << endl;
    Figure fig;
    fig.set_bounding_box(Box{{-1.0,+21.0},{-1.125,+1.125}});
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("test_continuous_evolution-vdp");

}

void TestContinuousEvolution::failure_test() const
{
    // The systems in this test are stiff and are expected to fail.

    // cout << __PRETTY_FUNCTION__ << endl;

    typedef Enclosure EnclosureType;

    // Set up the evolution parameters and grid
    Float time(0.5);
    double step_size(0.01);
    double enclosure_radius(0.25);

    // Set up the evaluators
    TaylorPicardIntegrator integrator(maximum_error=1e-6,sweep_threshold=1e-8,lipschitz_constant=0.5,
                                      step_maximum_error=1e-8,step_sweep_threshold=1e-10,maximum_temporal_order=6);

    // Define the initial box
    Box initial_box = Box{{0.0,0.0},{0.9,0.9}};

    // cout << "initial_box=" << initial_box << endl;

    // Set up the vector field for the first test
    Vector<Real> p = {200.0};
    VectorUserFunction<FailOne> failone(p);
    VectorField failone_vf(failone);

    VectorFieldEvolver evolverone(failone_vf,integrator);
    evolverone.configuration().maximum_enclosure_radius(enclosure_radius);
    evolverone.configuration().maximum_step_size(step_size);

    TaylorFunctionFactory function_factory(ThresholdSweeper(1e-10));
    EnclosureType initial_set(initial_box,function_factory);
    // cout << "initial_set=" << initial_set << endl << endl;

    Semantics semantics=UPPER_SEMANTICS;

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
    Figure fig; fig.set_bounding_box( orbit.reach().bounding_box() + Box(2,Interval(-0.25,+0.25)) );
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(magenta) << orbit.intermediate();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("test_continuous_evolution-failone");

    // Set up the vector field for the second test
    p[0] = 0.1;
    VectorUserFunction<FailTwo> failtwo(p);
    VectorField failtwo_vf(failtwo);

    VectorFieldEvolver evolvertwo(failtwo_vf,integrator);
    evolvertwo.configuration().maximum_enclosure_radius(enclosure_radius);
    evolvertwo.configuration().maximum_step_size(step_size);

    Box initial_box2 = Box{{0.0,0.0},{1.0,1.0},{1.0,1.0}};
    initial_set = EnclosureType(initial_box2,function_factory);

    time = 1.5;

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

