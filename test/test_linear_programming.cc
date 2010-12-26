/***************************************************************************
 *            test_linear_programming.cc
 *
 *  Copyright  2009  Pieter Collins
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

#include <iostream>
#include <fstream>

#include "test.h"

#include "config.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "linear_programming.h"


using namespace std;
using namespace Ariadne;

class TestInteriorPointSolver
{
    InteriorPointSolver* optimiser;
  public:
    TestInteriorPointSolver(InteriorPointSolver& o) : optimiser(&o) { }

    void test() {
        ARIADNE_TEST_CALL(test_validate_feasibility());
        //ARIADNE_TEST_CALL(test_optimization());
        //ARIADNE_TEST_CALL(test_feasibility());
        ARIADNE_TEST_CALL(test_constrained_feasibility());
    }

    void test_validate_feasibility() {
        // A feasible instance
        FloatVector xl(3, 0.0);
        FloatVector xu(3, infty);
        FloatMatrix A(2,3, 1.0,0.0,1.0, 0.0,1.0,2.0);
        FloatVector b(2, 1.0,1.0);
        FloatVector x(3, 0.81,0.51,0.21);
        FloatVector y(2, -1.0,-0.5);
        ARIADNE_ASSERT(definitely(optimiser->validate_feasibility(xl,xu,A,b,x,y)));

        b=FloatVector(2, 1.0, -1.0);
        y=FloatVector(2, -0.5, -1.0);
        ARIADNE_ASSERT(definitely(!optimiser->validate_feasibility(xl,xu,A,b,x,y)));
    }


    void test_feasibility() {
        // A feasible instance
        FloatVector xl(3, 0.0);
        FloatVector xu(3, infty);
        FloatMatrix A(2,3, 1.0,0.0,1.0, 0.0,1.0,2.0);
        FloatVector b(2, 1.0,1.0);
        ARIADNE_ASSERT(definitely(optimiser->feasible(xl,xu,A,b)));

        b=FloatVector(2, 1.0, -1.0);
        ARIADNE_ASSERT(definitely(!optimiser->feasible(xl,xu,A,b)));
    }


    void test_constrained_feasibility() {
        FloatMatrix A(2,3, 1.0,0.0,1.0, 0.0,1.0,2.0);
        FloatVector b(2, 1.0,1.0);
        FloatVector xl(3, 0.0,0.0,0.0);
        FloatVector xu(3, 4.0,2.0,3.0);

        ARIADNE_ASSERT(definitely(optimiser->feasible(xl,xu,A,b)));

        b=FloatVector(2, 1.0, -1.0);
        ARIADNE_ASSERT(definitely(!optimiser->feasible(xl,xu,A,b)));

        xu=FloatVector(3, +infty,+infty,3.0);
    }


    void test_optimization() {
        FloatMatrix A(2,3, 1.0,0.0,1.0, 0.0,1.0,2.0);
        FloatVector b(2, 1.0,1.0);
        FloatVector c(3, 1.0,0.5,-0.75);
        FloatVector xl(3, 0.0,0.0,0.0);
        FloatVector xu(3, +infty,+infty,3.0);

        ARIADNE_TEST_PRINT(optimiser->minimise(c,xl,xu,A,b));
    }

};


int main(int argc, const char* argv[])
{
    uint optimiser_verbosity = 0;
    if(argc>1) { optimiser_verbosity=atoi(argv[1]); }

    InteriorPointSolver interior_point_optimiser;
    interior_point_optimiser.verbosity=optimiser_verbosity;
    TestInteriorPointSolver(interior_point_optimiser).test();

    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}

