/***************************************************************************
 *            test_iterated_map_evolver.cpp
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
#include "dynamics/orbit.hpp"
#include "geometry/box.hpp"
#include "geometry/list_set.hpp"
#include "dynamics/iterated_map.hpp"
#include "dynamics/iterated_map_evolver.hpp"
#include "output/graphics.hpp"
#include "output/logging.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

class TestIteratedMapEvolver
{
  public:
    Void test() const;
};

Int main()
{
    ARIADNE_TEST_CALL(TestIteratedMapEvolver().test());
    return ARIADNE_TEST_FAILURES;
}

Void TestIteratedMapEvolver::test() const
{
    typedef DoublePrecision PR;
    DoublePrecision pr;

    // Set up the map field
    // The Henon map \f$(x,y)\mapsto(a-x^2+by,x)
    Real a=1.5_dyadic; Real b=0.375_dyadic;
    RealVariable x("x"), y("y");
    IteratedMap henon({ next(x)=a-x*x+b*y, next(y)=x });
    ARIADNE_TEST_PRINT(henon);

    // Define the initial set
    RealVariablesBox initial_box({1.01_dec<=x<=1.03_dec,0.51_dec<=y<=0.53_dec});

    ARIADNE_TEST_PRINT(initial_box);

    // Function evaluation sanity check
    Vector<FloatApproximation<PR>> p={{a,b},pr};
    Vector<FloatApproximation<PR>> xa={{0.5_x,0.25_x},pr};
    Vector<FloatApproximation<PR>> hxa={p[0]-xa[0]*xa[0]+xa[1]*p[1], xa[0]};
    ARIADNE_TEST_EQUAL(henon.update_function().evaluate(xa),hxa);
    Matrix<FloatApproximation<PR>> dhxa={{-2*xa[0],p[1]},{{1.0,pr},{0.0,pr}}};
    ARIADNE_TEST_EQUAL(henon.update_function().jacobian(xa),dhxa);

    // Set up the evolution parameters and grid
    IteratedMap::TimeType time(3);

    // Set up the evaluators
    IteratedMapEvolver evolver(henon);
    evolver.configuration().set_maximum_enclosure_radius(0.25);

    ARIADNE_TEST_PRINT(evolver.configuration());

    auto initial_set = evolver.enclosure(initial_box,EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,1e-10))));
    ARIADNE_TEST_PRINT(initial_set);

    // Compute the orbit
    auto orbit = evolver.orbit(initial_set,time);

    // Print the initial, evolve and reach sets
    LabelledFigure fig(Axes2d(-10,x,10, -10,y,10));
    fig << line_style(true) << fill_colour(cyan) << orbit.reach();
    fig << fill_colour(yellow) << orbit.final();
    fig << fill_colour(blue) << initial_set;
    fig.write("test_iterated_map_evolver");
}
