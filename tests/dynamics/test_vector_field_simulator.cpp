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
        RealVariable x("x"), y("y");

        VectorField vanderpol({dot(x)=y,dot(y)=mu*(1-x*x)*y-x});
        ARIADNE_TEST_PRINT(vanderpol);

        LabelledRealPoint initial_point({x=-1.5_dec,y=1});

        Real time = 6.3_dec;

        VectorFieldSimulator simulator(vanderpol);
        simulator.configuration().set_step_size(0.05);

        ARIADNE_TEST_PRINT(simulator.configuration());

        Orbit<PointType> orbit = simulator.orbit(initial_point,time);

        ARIADNE_TEST_ASSERT(distance(orbit.curve().end()->second,PointType(initial_point,double_precision)).raw() <= 0.02);

        LabelledFigure fig(Axes2d(-2.5<=x<=2.5,-3.0<=y<=3.0));
        fig << orbit.curve();
        fig.write("test_vector_field_simulator");
    }
};

Int main(Int argc, const char* argv[])
{
    ARIADNE_LOG_SET_VERBOSITY(get_verbosity(argc,argv));
    TestContinuousEvolution().test();
    return ARIADNE_TEST_FAILURES;
}