/***************************************************************************
 *            test_discrete_evolution.cpp
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
#include "function/taylor_model.hpp"
#include "algebra/differential.hpp"
#include "function/constraint.hpp"
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "function/formula.hpp"
#include "dynamics/enclosure.hpp"
#include "geometry/box.hpp"
#include "geometry/list_set.hpp"
#include "dynamics/map.hpp"
#include "dynamics/map_evolver.hpp"
#include "output/graphics.hpp"
#include "output/logging.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

class TestMapEvolver
{
  public:
    Void test() const;
};

Int main()
{
    ARIADNE_TEST_CALL(TestMapEvolver().test());
    return ARIADNE_TEST_FAILURES;
}

Void TestMapEvolver::test() const
{
    typedef DoublePrecision PR;
    DoublePrecision pr;

    typedef Enclosure EnclosureType;

    // Define the initial box
    ExactBoxType initial_box(2);
    initial_box[0]=ExactIntervalType(1.01,1.03);
    initial_box[1]=ExactIntervalType(0.51,0.53);

    ARIADNE_TEST_PRINT(initial_box);

    // Set up the map field
    // The Henon map \f$(x,y)\mapsto(a-x^2+by,x)
    Real a=Dyadic(1.5); Real b=Dyadic(0.375);
    EffectiveVectorMultivariateFunction henon;
    {
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction y=EffectiveScalarMultivariateFunction::coordinate(2,1);
        henon = { a-x*x+b*y, x };
    }
    ARIADNE_TEST_PRINT(henon);

    // Function evaluation sanity check
    Vector<FloatApproximation<PR>> p={{a,b},pr};
    Vector<FloatApproximation<PR>> x={{0.5,0.25},pr};
    Vector<FloatApproximation<PR>> hx={p[0]-x[0]*x[0]+x[1]*p[1], x[0]};
    ARIADNE_TEST_EQUAL(henon.evaluate(x),hx);
    Matrix<FloatApproximation<PR>> dhx={{-2*x[0],p[1]},{1.0_approx,0.0_approx}};
    ARIADNE_TEST_EQUAL(henon.jacobian(x),dhx);


    // Function evaluation sanity check
    ARIADNE_TEST_PRINT(initial_box);
    ARIADNE_TEST_PRINT(image(initial_box,henon));
    ARIADNE_TEST_PRINT(jacobian_range(henon,cast_vector(initial_box)));



    // Over-approximate the initial set by a grid cell
    EnclosureType initial_set(initial_box,TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,1e-10)));
    ARIADNE_TEST_PRINT(initial_set);

    // Set up the evolution parameters and grid
    IteratedMap::TimeType time(3);
    double enclosure_radius(0.25);

    // Set up the evaluators
    MapEvolver evolver(henon);
    evolver.configuration().set_maximum_enclosure_radius(enclosure_radius);

    // Compute the reachable sets
    ListSet<EnclosureType> evolve_set,reach_set;
    //evolve_set = evolver.evolve(initial_set,time);
    reach_set = evolver.reach(initial_set,time);
    ARIADNE_TEST_PRINT(initial_set.bounding_box());
    //cout << "evolve_bounding_boxes=" << evolve_set.bounding_boxes() << endl;
    ARIADNE_TEST_PRINT(reach_set.bounding_boxes());

    // Print the intial, evolve and reach sets
    Figure fig;
    fig.set_bounding_box({{-4,2},{-3,3}});
    fig << line_style(true) << fill_colour(cyan) << reach_set;
    fig << fill_colour(yellow) << evolve_set;
    fig << fill_colour(blue) << initial_set;
    fig.write("test_discrete_evolution-henon");
}
