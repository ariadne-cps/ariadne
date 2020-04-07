/***************************************************************************
 *            test_solvers.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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
#include <sstream>
#include <string>

#include "config.hpp"

#include "solvers/solver.hpp"
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "algebra/vector.hpp"
#include "algebra/algebra.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/space.hpp"
#include "function/formula.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

typedef Vector<DyadicInterval> DyadicIntervalVector;

class TestSolver
{
  private:
    std::unique_ptr<SolverInterface> solver;
    StringType solver_class_name;
  public:
    TestSolver(const SolverInterface& s,const char* n)
        : solver(s.clone()), solver_class_name(n) { }

    Int test() {
        ARIADNE_TEST_PRINT(*solver);
        ARIADNE_TEST_CALL(test_solve());
        ARIADNE_TEST_CALL(test_implicit());
        ARIADNE_TEST_CALL(test_scalar_implicit());
        return 0;
    }

    Void test_solve() {
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(1,0);
        ExactBoxType d({ExactIntervalType(0.0,1.0)});
        EffectiveVectorMultivariateFunction f({(x*x+1)*x-1});
        FloatDPBoundsVector p=solver->solve(f,d);
        ARIADNE_TEST_BINARY_PREDICATE(contains,ExactIntervalType(0.6823,0.6824),p[0]);
    }

    Void test_implicit() {
        //TaylorModelAccuracy::set_default_sweep_threshold(1e-12);

        // FIXME: Should be able to use numbers yielding exact results for p,r
        EffectiveScalarMultivariateFunction aa=EffectiveScalarMultivariateFunction::coordinate(1,0);
        EffectiveScalarMultivariateFunction a=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(2,1);
        EffectiveScalarMultivariateFunction bb;
        ExactBoxType p,r;
        EffectiveVectorMultivariateFunction f;
        ValidatedVectorMultivariateFunctionModelDP h;
        EffectiveVectorMultivariateFunction e;
        FloatDPValue tol;

        // Test solution of x-a=0. This should be very easy to solve.
        p=ExactBoxType({ExactIntervalType(-0.25,0.25)});
        r=ExactBoxType({ExactIntervalType(-2.0,2.0)});
        f=EffectiveVectorMultivariateFunction({x-a});
        ARIADNE_TEST_PRINT(f);
        h=solver->implicit(f,p,r);
        ARIADNE_TEST_PRINT(h);
        e=EffectiveVectorMultivariateFunction(1u,aa);
        ARIADNE_TEST_PRINT(e);
        ARIADNE_TEST_COMPARE(norm(h-e),<,1e-8);

        // Test solution of 4x^2+x-4-a=0 on [0.875,1.125]. There is a unique solution with positive derivative.
        p=ExactBoxType({ExactIntervalType(0.875,1.125)});
        r=ExactBoxType({ExactIntervalType(0.25,1.25)});
        f=EffectiveVectorMultivariateFunction({(x*x+1)*x-a});
        ARIADNE_TEST_PRINT(f);
        h=solver->implicit(f,p,r);
        ARIADNE_TEST_PRINT(h);
        bb=EffectiveScalarMultivariateFunction(aa-DyadicInterval(p[0]).midpoint())/DyadicInterval(p[0]).radius();
        Decimal a0(0.682328), a1(0.0521547), a2(-0.0023232), a3(0.000147778);
        e=EffectiveVectorMultivariateFunction( { a0+bb*(a1+bb*(a2+bb*a3)) } );
        ARIADNE_TEST_PRINT(e);
        ARIADNE_TEST_COMPARE(norm(h-e),<,1e-4);

        // Test solution of 4x^2+x-4-a=0 on [-0.25,0.25]. There is a unique solution with positive derivative.
        p=ExactBoxType({ExactIntervalType(-0.25,0.25)});
        r=ExactBoxType({ExactIntervalType(0.25,2.0)});
        f=EffectiveVectorMultivariateFunction({4*x+x*x-a-4});
        ARIADNE_TEST_PRINT(f);
        h=solver->implicit(f,p,r);
        ARIADNE_TEST_PRINT(h);
        bb=EffectiveScalarMultivariateFunction(aa-DyadicInterval(p[0]).midpoint())/DyadicInterval(p[0]).radius();
        Decimal c0(0.828427), c1(0.0441942), c2(-0.000345267), c3(0.00000539468);
        e=EffectiveVectorMultivariateFunction( { c0+bb*(c1+bb*(c2+bb*c3)) } );
        ARIADNE_TEST_PRINT(e);
        ARIADNE_TEST_COMPARE(norm(h-e),<,1e-4);

        // Test solution of x-2*a=0 on [-1,+1], taking values in [-1,+1].
        // There is at most one solution, but this lies partially outside the range.
        // Should obtain PartialSolutionException
        p=ExactBoxType({ExactIntervalType(-1,1)});
        r=ExactBoxType({ExactIntervalType(-1,1)});
        f=EffectiveVectorMultivariateFunction({x-2*a});
        ARIADNE_TEST_PRINT(f);
        try {
            h=solver->implicit(f,p,r);
            ARIADNE_TEST_NOTIFY(solver_class_name<<" silently returns partially defined vector implicit function.");
        }
        catch(const SolverException& exc) {
            ARIADNE_TEST_THROWS(solver->implicit(f,p,r),SolverException);
            ARIADNE_TEST_NOTIFY(solver_class_name<<" throws error on partially defined vector implicit function.");
        }

    }

    Void test_scalar_implicit() {
        //TaylorModelAccuracy::set_default_sweep_threshold(1e-12);

        EffectiveScalarMultivariateFunction aa=EffectiveScalarMultivariateFunction::coordinate(1,0);
        EffectiveScalarMultivariateFunction a=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(2,1);
        // Test solution of x-2*a=0 on [-1,+1], taking values in [-1,+1]. There is at most one solution.
        // Uses scalar implicit
        ExactBoxType p; ExactIntervalType r;
        EffectiveScalarMultivariateFunction e,f,s; // s is unscaling functions
        ValidatedScalarMultivariateFunctionModelDP h;

        ARIADNE_TEST_PRINT(*solver);

        // Test solution of 4x^2+x-4-a=0 on [0.875,1.125]. There is a unique solution with positive derivative.
        p=ExactBoxType({ExactIntervalType(0.875,1.125)});
        r=ExactIntervalType(0.25,1.25);
        f=EffectiveScalarMultivariateFunction((x*x+1)*x-a);
        ARIADNE_TEST_PRINT(f);
        h=solver->implicit(f,p,r);
        ARIADNE_TEST_PRINT(h);
        s=EffectiveScalarMultivariateFunction(aa-DyadicInterval(p[0]).midpoint())/DyadicInterval(p[0]).radius();
        Decimal a0(0.682328), a1(0.0521547), a2(-0.0023232), a3(0.000147778);
        e=EffectiveScalarMultivariateFunction( a0+s*(a1+s*(a2+s*a3)) );
        ARIADNE_TEST_PRINT(e);
        ARIADNE_TEST_COMPARE(mag((h-e).range()),<,1e-4);

        // Test solution of x-2*a=0 on [-1,+1], taking values in [-1,+1].
        // There is at most one solution, but this lies partially outside the range.
        // Should obtain PartialSolutionException
        p=ExactBoxType({ExactIntervalType(-1,1)});
        r=ExactIntervalType(-1,1);
        f=EffectiveScalarMultivariateFunction(x-2*a);
        ARIADNE_TEST_PRINT(f);
        ValidatedScalarMultivariateFunctionModelDP g=ValidatedScalarMultivariateTaylorFunctionModelDP(product(p,r),f,ThresholdSweeper<FloatDP>(dp,1e-12));
        ARIADNE_TEST_PRINT(g);
        try {
            h=solver->implicit(g,p,r);
            ARIADNE_TEST_NOTIFY(solver_class_name<<" silently returns partially defined scalar implicit function.");
        }
        catch(const SolverException& exc) {
            ARIADNE_TEST_THROWS(solver->implicit(g,p,r),SolverException);
            ARIADNE_TEST_NOTIFY(solver_class_name<<" throws error on partially defined scalar implicit function.");
        }

    }

};

#include "algebra/differential.hpp"

Int main(Int argc, const char **argv) {

    unsigned int verb=get_verbosity(argc,argv);

    IntervalNewtonSolver interval_newton_solver(maximum_error=1e-5,maximum_number_of_steps=12);
    interval_newton_solver.verbosity=verb;
    TestSolver(interval_newton_solver,"IntervalNewtonSolver").test();

    KrawczykSolver krawczyk_solver(maximum_error=1e-5,maximum_number_of_steps=12);
    ARIADNE_TEST_PRINT(krawczyk_solver.function_factory());
    krawczyk_solver.verbosity=verb;
    TestSolver(krawczyk_solver,"KrawczykSolver").test();

    FactoredKrawczykSolver factored_krawczyk_solver(maximum_error=1e-5,maximum_number_of_steps=12);
    factored_krawczyk_solver.verbosity=verb;
    TestSolver(factored_krawczyk_solver,"FactoredKrawczykSolver").test();

    std::cerr<<"INCOMPLETE "<<std::flush;
    return ARIADNE_TEST_FAILURES;
}
