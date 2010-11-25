/***************************************************************************
 *            test_constraint_solver.cc
 *
 *  Copyright  2010  Pieter Collins
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

#include "numeric.h"
#include "vector.h"
#include "function.h"
#include "constraint.h"
#include "constraint_solver.h"
#include "box.h"


using namespace std;
using namespace Ariadne;

class TestConstraintSolver
{
    uint verbosity;
  public:
    TestConstraintSolver(uint v) : verbosity(v) { }

    void test() {
        ARIADNE_TEST_CALL(test_empty_reduce_inequality());
        ARIADNE_TEST_CALL(test_empty_reduce_equality());
        ARIADNE_TEST_CALL(test_empty_reduce_mixed());
        ARIADNE_TEST_CALL(test_empty_hull_reduce());
        ARIADNE_TEST_CALL(test_empty_box_reduce());
        ARIADNE_TEST_CALL(test_hull_reduce());
        ARIADNE_TEST_CALL(test_box_reduce());
        ARIADNE_TEST_CALL(test_monotone_reduce());
        ARIADNE_TEST_CALL(test_split());
    }

    void test_empty_reduce_inequality() {
        List<RealScalarFunction> x=RealScalarFunction::coordinates(2);
        Box D(2, 0.0,1.0, 0.0,1.0);
        List<NonlinearConstraint> c;
        c.append(4.0<=2*x[0]+x[1]);

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.reduce(D,c));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.empty());
    }

    void test_empty_reduce_equality() {
        List<RealScalarFunction> x=RealScalarFunction::coordinates(2);
        Box D(2, 0.0,1.0, 0.0,1.0);
        List<NonlinearConstraint> c;
        c.append(2*x[0]+x[1]==4);

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.reduce(D,c.function(),c.bounds()));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.empty());
    }

    void test_empty_reduce_mixed() {
        List<RealScalarFunction> x=RealScalarFunction::coordinates(2);
        Box D(2, 0.0,0.25, 0.0, 2.0);
        List<NonlinearConstraint> c;
        c.append(x[1]<=1);
        c.append(x[0]+x[1]==2);

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.reduce(D,c));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.empty());
    }

    void test_empty_hull_reduce() {
        List<RealScalarFunction> x=RealScalarFunction::coordinates(2);
        Box D(2, 0.0,0.25, 0.0, 2.0);
        List<NonlinearConstraint> c;
        c.append(x[1]<=1);
        c.append(x[0]+x[1]==2);

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[0]));
        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[1]));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.empty());
    }

    void test_empty_box_reduce() {
        List<RealScalarFunction> x=RealScalarFunction::coordinates(2);
        Box D(2, 0.0,0.25, 0.0, 2.0);
        List<NonlinearConstraint> c;
        c.append(x[1]<=1);
        c.append(x[0]+x[1]==2);

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c[0],0));
        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c[1],0));
        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c[0],1));
        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c[1],1));
        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[0]));
        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[1]));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.empty());
    }

    void test_hull_reduce() {
        List<RealScalarFunction> x=RealScalarFunction::coordinates(2);
        Box D(2, 0.0,2.0, 0.0,2.0);
        List<NonlinearConstraint> c;
        c.append(-2.0<=2*x[0]+x[1]<=1.0);

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[0]));
        ARIADNE_TEST_EQUAL(D,Box(2,0.0,0.5, 0.0,1.0));
    }

    void test_box_reduce() {
        List<RealScalarFunction> x=RealScalarFunction::coordinates(2);
        Box D(2, 0.0,2.0, 0.0,2.0);
        NonlinearConstraint c=(-2.0<=2*x[0]+x[1]<=1.0);

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c,0));
        ARIADNE_TEST_EQUAL(D,Box(2,0.0,0.75, 0.0,2.0));
        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c,1));
        ARIADNE_TEST_EQUAL(D,Box(2,0.0,0.75, 0.0,1.25));
    }


    void test_monotone_reduce() {
        List<RealScalarFunction> x=RealScalarFunction::coordinates(2);
        Box D(2, 0.0,2.0, 0.0,2.0);
        NonlinearConstraint c=(-2.0<=2*x[0]+x[1]<=1.0);

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c,0));
        ARIADNE_TEST_EQUAL(D,Box(2,0.0,0.75, 0.0,2.0));
        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c,1));
        ARIADNE_TEST_EQUAL(D,Box(2,0.0,0.75, 0.0,1.25));
    }

    void test_split() {
        ARIADNE_TEST_WARN("test_split: Not implemented");
    }
};

int main(int argc, const char* argv[]) {
    TestConstraintSolver(get_verbosity(argc,argv)).test();
    return ARIADNE_TEST_FAILURES;
}

