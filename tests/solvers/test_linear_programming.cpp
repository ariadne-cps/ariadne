/***************************************************************************
 *            test_linear_programming.cpp
 *
 *  Copyright  2009  Pieter Collins
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

#include <iostream>
#include <fstream>

#include "config.hpp"
#include "../test.hpp"

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "solvers/linear_programming.hpp"


using namespace std;
using namespace Ariadne;

class TestInteriorPointSolver
{
    InteriorPointSolver* optimiser;
  public:
    TestInteriorPointSolver(InteriorPointSolver& o) : optimiser(&o) { }

    Void test() {
        ARIADNE_TEST_CALL(test_validate_feasibility());
        //ARIADNE_TEST_CALL(test_optimization());
        //ARIADNE_TEST_CALL(test_feasibility());
        ARIADNE_TEST_CALL(test_constrained_feasibility());
    }

    Void test_validate_feasibility() {
        // A feasible instance
        RawFloatVector xl={0.0,0.0,0.0};
        RawFloatVector xu={inf,inf,inf};
        FloatMatrix A={{1.0,0.0,1.0}, {0.0,1.0,2.0}};
        RawFloatVector b={1.0,1.0};
        RawFloatVector x={0.81,0.51,0.21};
        RawFloatVector y={-1.0,-0.5};
        ARIADNE_ASSERT(definitely(optimiser->validate_feasibility(xl,xu,A,b,x,y)));

        b=RawFloatVector{1.0, -1.0};
        y=RawFloatVector{-0.5, -1.0};
        ARIADNE_ASSERT(definitely(!optimiser->validate_feasibility(xl,xu,A,b,x,y)));
    }


    Void test_feasibility() {
        // A feasible instance
        RawFloatVector xl={0.0,0.0,0.0};
        RawFloatVector xu={inf,inf,inf};
        FloatMatrix A={{1.0,0.0,1.0},{0.0,1.0,2.0}};
        RawFloatVector b={1.0,1.0};
        ARIADNE_ASSERT(definitely(optimiser->feasible(xl,xu,A,b)));

        b=RawFloatVector{1.0,-1.0};
        ARIADNE_ASSERT(definitely(!optimiser->feasible(xl,xu,A,b)));
    }


    Void test_constrained_feasibility() {
        FloatMatrix A={{1.0,0.0,1.0},{0.0,1.0,2.0}};
        RawFloatVector b={1.0,1.0};
        RawFloatVector xl={0.0,0.0,0.0};
        RawFloatVector xu={4.0,2.0,3.0};

        ARIADNE_ASSERT(definitely(optimiser->feasible(xl,xu,A,b)));

        b=RawFloatVector{1.0,-1.0};
        ARIADNE_ASSERT(definitely(!optimiser->feasible(xl,xu,A,b)));

        xu=RawFloatVector{+inf,+inf,3.0};
    }


    Void test_optimization() {
        FloatMatrix A={{1.0,0.0,1.0},{0.0,1.0,2.0}};
        RawFloatVector b={1.0,1.0};
        RawFloatVector c={1.0,0.5,-0.75};
        RawFloatVector xl={0.0,0.0,0.0};
        RawFloatVector xu={+inf,+inf,3.0};

        ARIADNE_TEST_PRINT(optimiser->minimise(c,xl,xu,A,b));
    }

};


Int main(Int argc, const char* argv[])
{
    auto verbosity = get_verbosity(argc,argv);

    InteriorPointSolver interior_point_optimiser;
    interior_point_optimiser.verbosity=verbosity;
    TestInteriorPointSolver(interior_point_optimiser).test();

    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

