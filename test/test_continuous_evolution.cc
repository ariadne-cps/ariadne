/***************************************************************************
 *      test_continuous_evolution.cc
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

#include "tuple.h"
#include "vector.h"
#include "matrix.h"
#include "function.h"
#include "taylor_function.h"
#include "taylor_set.h"
#include "box.h"
#include "list_set.h"
#include "evolution_parameters.h"
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

    typedef TaylorConstrainedImageSet EnclosureType;

    // Set up the evolution parameters and grid
    Float time(3.0);
    double step_size(0.5);
    double enclosure_radius(0.25);

    EvolutionParameters parameters;
    parameters.maximum_enclosure_radius=enclosure_radius;
    parameters.maximum_step_size=step_size;

    // Set up the evaluators
    TaylorIntegrator integrator(14,1e-4,1e-12);
    VectorFieldEvolver evolver(parameters,integrator);

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


    // Over-approximate the initial set by a grid cell
    EnclosureType initial_set(initial_box);
    // cout << "initial_set=" << initial_set << endl << endl;

    Semantics semantics=UPPER_SEMANTICS;

    // Compute the reachable sets
    ListSet<EnclosureType> evolve_set,reach_set;
    Orbit<EnclosureType> orbit = evolver.orbit(vanderpol,initial_set,time,semantics);

    ARIADNE_TEST_PRINT(orbit);

    // Print the intial, evolve and reach sets
    // cout << "Plotting sets" << endl;
    // cout << "evolve_set=" << hybrid_evolve_set << endl;
    // cout << "reach_set=" << hybrid_reach_set << endl;
    Figure fig;
    fig.set_bounding_box(Box(2, -1.0,+21.0, -1.125,+1.125));
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

    typedef TaylorConstrainedImageSet EnclosureType;

    // Set up the evolution parameters and grid
    Float time(0.5);
    double step_size(0.01);
    double enclosure_radius(0.25);

    EvolutionParameters parameters;
    parameters.maximum_enclosure_radius=enclosure_radius;
    parameters.maximum_step_size=step_size;

    // Set up the evaluators
    TaylorIntegrator integrator(6,1e-6,1e-10);
    VectorFieldEvolver evolver(parameters,integrator);

    // Define the initial box
    Box initial_box(2, 0.0,0.0, 0.9,0.9 );

    // cout << "initial_box=" << initial_box << endl;

    // Set up the vector field for the first test
    Vector<Float> p(1, 200.0);
    VectorUserFunction<FailOne> failone(p);
    VectorField failone_vf(failone);

    EnclosureType initial_set(initial_box);
    // cout << "initial_set=" << initial_set << endl << endl;

    Semantics semantics=UPPER_SEMANTICS;

    // Compute the reachable sets
    Orbit<EnclosureType> orbit = evolver.orbit(failone_vf,initial_set,time,semantics);

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

    Box initial_box2(3, 0.0,0.0, 1.0,1.0, 1.0,1.0);
    initial_set = EnclosureType(initial_box2);

    time = 1.5;

    // Compute the reachable sets
    orbit = evolver.orbit(failtwo_vf,initial_set,time,semantics);

    cout << "\norbit.final=\n" << orbit.final() << endl << endl;
    cout << "final set radius="<< orbit.final()[0].radius() << endl;

    ARIADNE_TEST_COMPARE(orbit.final()[0].radius(),>,1.0);

    // Print the intial, evolve and reach sets
    // cout << "Plotting sets" << endl;
    // cout << "evolve_set=" << hybrid_evolve_set << endl;
    // cout << "reach_set=" << hybrid_reach_set << endl;
    std::cout << "Plotting..." << std::flush;
    fig.clear(); fig.set_bounding_box(orbit.reach().bounding_box());
    fig.set_projection_map(ProjectionFunction(2,3,0));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(red) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("test_continuous_evolution-failtwo");
    std::cout << std::endl;

}

