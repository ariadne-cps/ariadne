/***************************************************************************
 *            test_discrete_evolution.cpp
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
#include "geometry/enclosure.hpp"
#include "geometry/box.hpp"
#include "geometry/list_set.hpp"
#include "dynamics/map.hpp"
#include "dynamics/map_evolver.hpp"
#include "output/graphics.hpp"
#include "utility/logging.hpp"

#include "test.hpp"

using namespace Ariadne;
using namespace std;

class TestMapEvolver
{
  public:
    Void test() const;
};

Int main()
{
    TestMapEvolver().test();
    return ARIADNE_TEST_FAILURES;
}

Void TestMapEvolver::test() const
{
    cout << __PRETTY_FUNCTION__ << endl;

    DoublePrecision pr;

    typedef Enclosure EnclosureType;

    // Define the initial box
    ExactBoxType initial_box(2);
    initial_box[0]=ExactIntervalType(1.01,1.03);
    initial_box[1]=ExactIntervalType(0.51,0.53);

    cout << "initial_box=" << initial_box << endl;

    // Set up the map field
    // The Henon map \f$(x,y)\mapsto(a-x^2+by,x)
    Real a=Dyadic(1.5); Real b=Dyadic(0.375);
    EffectiveVectorFunction henon;
    {
        EffectiveScalarFunction x=EffectiveScalarFunction::coordinate(2,0);
        EffectiveScalarFunction y=EffectiveScalarFunction::coordinate(2,1);
        henon = { a-x*x+b*y, x };
    }
    cout << "henon_function=" << henon << endl;

    //VectorUserFunction evaluation sanity check
    Vector<ApproximateNumericType> p={{a,b},pr};
    Vector<ApproximateNumericType> x={{0.5,0.25},pr};
    Vector<ApproximateNumericType> hx={p[0]-x[0]*x[0]+x[1]*p[1], x[0]};
    ARIADNE_TEST_EQUAL(henon.evaluate(x),hx);
    Matrix<ApproximateNumericType> dhx={{-2*x[0],p[1]},{1.0_approx,0.0_approx}};
    ARIADNE_TEST_EQUAL(henon.jacobian(x),dhx);


    //VectorUserFunction evaluation sanity check
    cout << "apply(henon," << initial_box << ") " << flush; cout << " = " << apply(henon,initial_box) << endl;
    cout << "jacobian_range(henon," << initial_box << ") = " << jacobian_range(henon,initial_box) << endl;



    // Over-approximate the initial set by a grid cell
    EnclosureType initial_set(initial_box,TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,1e-10)));
    cout << "initial_set=" << initial_set << endl << endl << endl;

    // Set up the evolution parameters and grid
    IteratedMap::TimeType time(3);
    double enclosure_radius(0.25);

    // Set up the evaluators
    MapEvolver evolver(henon);
    evolver.configuration().maximum_enclosure_radius(enclosure_radius);

    // Compute the reachable sets
    ListSet<EnclosureType> evolve_set,reach_set;
    //evolve_set = evolver.evolve(initial_set,time);
    reach_set = evolver.reach(initial_set,time);
    cout << "initial_bounding_box=" << initial_set.bounding_box() << endl;
    //cout << "evolve_bounding_boxes=" << evolve_set.bounding_boxes() << endl;
    cout << "reach_bounding_boxes=" << reach_set.bounding_boxes() << endl;

    // Print the intial, evolve and reach sets
    Figure fig;
    fig << line_style(true) << fill_colour(cyan) << reach_set;
    fig << fill_colour(yellow) << evolve_set;
    fig << fill_colour(blue) << initial_set;
    fig.write("test_discrete_evolution-henon");
}
