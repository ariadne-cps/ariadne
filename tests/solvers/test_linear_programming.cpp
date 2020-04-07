/***************************************************************************
 *            test_linear_programming.cpp
 *
 *  Copyright  2009-20  Pieter Collins
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
        RawFloatDPVector xl={0.0,0.0,0.0};
        RawFloatDPVector xu={inf,inf,inf};
        RawFloatDPMatrix A={{1.0,0.0,1.0}, {0.0,1.0,2.0}};
        RawFloatDPVector b={1.0,1.0};
        RawFloatDPVector x={0.81,0.51,0.21};
        RawFloatDPVector y={-1.0,-0.5};
        ARIADNE_ASSERT(definitely(optimiser->validate_feasibility(xl,xu,A,b,x,y)));

        b=RawFloatDPVector{1.0, -1.0};
        y=RawFloatDPVector{-0.5, -1.0};
        ARIADNE_ASSERT(definitely(!optimiser->validate_feasibility(xl,xu,A,b,x,y)));
    }


    Void test_feasibility() {
        // A feasible instance
        RawFloatDPVector xl={0.0,0.0,0.0};
        RawFloatDPVector xu={inf,inf,inf};
        RawFloatDPMatrix A={{1.0,0.0,1.0},{0.0,1.0,2.0}};
        RawFloatDPVector b={1.0,1.0};
        ARIADNE_ASSERT(definitely(optimiser->feasible(xl,xu,A,b)));

        b=RawFloatDPVector{1.0,-1.0};
        ARIADNE_ASSERT(definitely(!optimiser->feasible(xl,xu,A,b)));
    }


    Void test_constrained_feasibility() {
        RawFloatDPMatrix A={{1.0,0.0,1.0},{0.0,1.0,2.0}};
        RawFloatDPVector b={1.0,1.0};
        RawFloatDPVector xl={0.0,0.0,0.0};
        RawFloatDPVector xu={4.0,2.0,3.0};

        ARIADNE_ASSERT(definitely(optimiser->feasible(xl,xu,A,b)));

        b=RawFloatDPVector{1.0,-1.0};
        ARIADNE_ASSERT(definitely(!optimiser->feasible(xl,xu,A,b)));

        xu=RawFloatDPVector{+inf,+inf,3.0};
    }


    Void test_optimization() {
        RawFloatDPMatrix A={{1.0,0.0,1.0},{0.0,1.0,2.0}};
        RawFloatDPVector b={1.0,1.0};
        RawFloatDPVector c={1.0,0.5,-0.75};
        RawFloatDPVector xl={0.0,0.0,0.0};
        RawFloatDPVector xu={+inf,+inf,3.0};

        ARIADNE_TEST_PRINT(optimiser->minimise(c,xl,xu,A,b));
    }

};


Int main(Int argc, const char* argv[])
{
    auto verb = get_verbosity(argc,argv);

    InteriorPointSolver interior_point_optimiser;
    interior_point_optimiser.verbosity=verb;
    TestInteriorPointSolver(interior_point_optimiser).test();

    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

