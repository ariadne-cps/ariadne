/***************************************************************************
 *            test_vector_field_simulator.cpp
 *
 *  Copyright  2006-20  Luca Geretti, Mirko Albanese
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
#include "io/figure.hpp"
#include "io/command_line_interface.hpp"
#include "dynamics/orbit.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

class TestVectorFieldSimulator
{
  public:
    Void test() const {
        ARIADNE_TEST_CALL(test_multiple_real_expr_vdp_cycle());
        ARIADNE_TEST_CALL(test_multiple_real_box_vdp_cycle());
    }

    Void test_multiple_real_expr_vdp_cycle() const {

        typedef VectorFieldSimulator::ApproximateListPointType ListPointType;

        Real mu=Dyadic(0.5_x);
        RealVariable x("x"), y("y"), z("z");

        VectorField vanderpol({dot(x)=y,dot(y)=mu*(1-x*x)*y-x},{let(z)=sqrt(sqr(x)+sqr(y))});
        ARIADNE_TEST_PRINT(vanderpol);

        RealExpressionBoundedConstraintSet initial_set({-2.0_dec<=x<=-1.5_dec,0.5_dec<=y<=1});

        Real time = 6.3_dec;

        VectorFieldSimulator simulator(vanderpol);
        simulator.configuration().set_step_size(0.05);
        simulator.configuration().set_num_subdivisions(1);

        ARIADNE_TEST_PRINT(simulator.configuration());

        Orbit<ListPointType> orbit = simulator.orbit(initial_set,time);

        /*
        auto final_pt_xy = project(orbit.curve()[orbit.num_curves()-1],Projection2d(3,0,1));
        auto initial_pt_xy = Point<FloatDPApproximation>(initial_set.euclidean_set(vanderpol.state_space()).bounding_box().midpoint(),dp);
        ARIADNE_TEST_ASSERT(distance(final_pt_xy,initial_pt_xy).raw() <= 0.02_dec);

        auto orbit2 = simulator.orbit(RealVariablesBox({x==-1.5_dec,y==1}),time);
        auto final_pt_xy_2 = project(orbit2.curve()[orbit.num_curves()-1],Projection2d(3,0,1));
        ARIADNE_TEST_ASSERT(distance(final_pt_xy,final_pt_xy_2).raw() == 0);*/

        LabelledFigure fig({-2.5_dec<=x<=2.5_dec,-2.5_dec<=y<=2.5_dec});
        fig << orbit;
        fig.write("test_real_expr_vector_field_simulator_xy");

        LabelledFigure fig2({0<=TimeVariable()<=6.3_dec,0<=z<=2.5_dec});
        fig2 << orbit;
        fig2.write("test_real_expr_vector_field_simulator_tz");

    }

    Void test_multiple_real_box_vdp_cycle() const{
        typedef VectorFieldSimulator::ApproximateListPointType ListPointType;
      
        Real mu=Dyadic(0.5_x);
        RealVariable x("x"), y("y"), z("z");

        VectorField vanderpol({dot(x)=y,dot(y)=mu*(1-x*x)*y-x},{let(z)=sqrt(sqr(x)+sqr(y))});
        ARIADNE_TEST_PRINT(vanderpol);

        RealVariablesBox initial_set({-2.0_dec<=x<=-1.5_dec,0.5_dec<=y<=1});
        
        Real time = 6.3_dec;

        VectorFieldSimulator simulator(vanderpol);
        simulator.configuration().set_step_size(0.05);
        simulator.configuration().set_num_subdivisions(0);

        ARIADNE_TEST_PRINT(simulator.configuration());

        Orbit<ListPointType> orbit = simulator.orbit(initial_set,time);

        LabelledFigure fig({-2.5<=x<=2.5,-3<=y<=3});
        fig << orbit;
        fig.write("test_real_box_list_vector_field_simulator_xy");

        LabelledFigure fig2({0<=TimeVariable()<=6.3_dec,0<=z<=2.5_dec});
        fig2 << orbit;
        fig2.write("test_real_box_list_vector_field_simulator_tz");
    }
};

Int main(Int argc, const char* argv[])
{
    if (not CommandLineInterface::instance().acquire(argc,argv)) return -1;
    TestVectorFieldSimulator().test();
    return ARIADNE_TEST_FAILURES;
}
