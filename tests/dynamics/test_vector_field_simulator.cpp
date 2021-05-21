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
#include "io/figure.hpp"
#include "io/command_line_interface.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

class TestVectorFieldEvolver
{
  public:
    Void test() const {
        ARIADNE_TEST_CALL(test_one_vdp_cycle());
    }

    Void test_one_vdp_cycle() const {

        typedef VectorFieldSimulator::ApproximatePointType PointType;

        Real mu=Dyadic(0.5_x);
        RealVariable x("x"), y("y"), z("z");

        VectorField vanderpol({dot(x)=y,dot(y)=mu*(1-x*x)*y-x},{let(z)=sqrt(sqr(x)+sqr(y))});
        ARIADNE_TEST_PRINT(vanderpol);

        RealExpressionBoundedConstraintSet initial_set({x==-1.5_dec,y==1});

        Real time = 6.3_dec;

        VectorFieldSimulator simulator(vanderpol);
        simulator.configuration().set_step_size(0.05);

        ARIADNE_TEST_PRINT(simulator.configuration());

        Orbit<PointType> orbit = simulator.orbit(initial_set,time);

        auto final_pt_xy = project(orbit.curve().end()->second,Projection2d(3,0,1));
        auto initial_pt_xy = Point<FloatDPApproximation>(initial_set.euclidean_set(vanderpol.state_space()).bounding_box().midpoint(),dp);
        ARIADNE_TEST_ASSERT(distance(final_pt_xy,initial_pt_xy).raw() <= 0.02);

        auto orbit2 = simulator.orbit(RealVariablesBox({x==-1.5_dec,y==1}),time);
        auto final_pt_xy_2 = project(orbit2.curve().end()->second,Projection2d(3,0,1));
        ARIADNE_TEST_ASSERT(distance(final_pt_xy,final_pt_xy_2).raw() == 0);

        LabelledFigure fig({-2.5<=x<=2.5,-3<=y<=3});
        fig << orbit.curve();
        fig.write("test_vector_field_simulator_xy");

        LabelledFigure fig2({0<=TimeVariable()<=6.3_dec,0<=z<=2.5_dec});
        fig2 << orbit.curve();
        fig2.write("test_vector_field_simulator_tz");
    }
};

Int main(Int argc, const char* argv[])
{
    if (not CommandLineInterface::instance().acquire(argc,argv)) return -1;
    TestVectorFieldEvolver().test();
    return ARIADNE_TEST_FAILURES;
}