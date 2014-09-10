/***************************************************************************
 *            test_nonlinear_programming.cc
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

#include "config.h"

#include "test.h"

#include "numeric.h"
#include "vector.h"
#include "function.h"
#include "nonlinear_programming.h"
#include "box.h"


using namespace std;
using namespace Ariadne;

class TestOptimiser
{
  private:
    std::unique_ptr<OptimiserInterface> optimiser;
  public:
    TestOptimiser(const OptimiserInterface& opt)
        : optimiser(opt.clone()) { }

    void test() {
        ARIADNE_TEST_CALL(test_feasibility_check());
        ARIADNE_TEST_CALL(test_unconstrained_optimisation());
        ARIADNE_TEST_CALL(test_constrained_optimisation());
        ARIADNE_TEST_CALL(test_equality_constrained_optimisation());
        ARIADNE_TEST_CALL(test_linear_feasibility());
        ARIADNE_TEST_CALL(test_nonlinear_feasibility());
        ARIADNE_TEST_CALL(test_nonlinear_equality_feasibility());
    }

    void test_unconstrained_optimisation() {
        // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(2);
        EffectiveScalarFunction x0s = sqr(x[0]);
        EffectiveScalarFunction f(x0s*(12+x0s*(Decimal(6.3)+x0s))+6*x[1]*(x[1]-x[0]));
        ARIADNE_TEST_PRINT(f);
        EffectiveVectorFunction g(0u,2u);
        ARIADNE_TEST_PRINT(g);
        Box D=Box{{-1.0,2.0},{-3.0,5.0}};
        Box C=Box{};
        ARIADNE_TEST_PRINT(Ariadne::make_tuple(f,D,g,C));

        IntervalVector x_optimal=optimiser->minimise(f,D,g,C);
        ARIADNE_TEST_BINARY_PREDICATE(subset,x_optimal,D);
        ARIADNE_TEST_BINARY_PREDICATE(subset,g(x_optimal),C);
        //ARIADNE_TEST_LESS(norm(x_optimal),1e-8);
    }

    void test_equality_constrained_optimisation() {
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(2);
        EffectiveScalarFunction f=(sqr(x[0])+sqr(x[1]));
        ARIADNE_TEST_PRINT(f);
        ExactFloat a(1.5); ExactFloat b(0.25);
        EffectiveVectorFunction g={a+x[0]+2*x[1]+b*x[0]*x[1]};
        ARIADNE_TEST_PRINT(g);
        IntervalVector C={{0.0,0.0}};
        Box D=Box{{-1.0,2.0},{-3.0,5.0}};
        ARIADNE_TEST_PRINT(Ariadne::make_tuple(f,D,g,C));

        IntervalVector x_optimal=optimiser->minimise(f,D,g,C);
        ARIADNE_TEST_BINARY_PREDICATE(subset,x_optimal,D);
        ARIADNE_TEST_LESS(norm(g(x_optimal)),1e-7);
    }

    void test_constrained_optimisation() {
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(3);
        EffectiveScalarFunction x0s = sqr(x[0]);
        EffectiveScalarFunction f = x0s*(12+x0s*(Decimal(6.3)+x0s))+6*x[1]*(x[1]-x[0])+x[2];
        ARIADNE_TEST_PRINT(f);
        //EffectiveVectorFunction g( (x[0]-1, x[0]+x[1]*x[1], x[1]*x[1]) );
        Box D = Box{{-1.0,2.0},{-3.0,5.0},{-3.0,5.0}};
        ARIADNE_TEST_PRINT(D);
        EffectiveVectorFunction g = {2*x[1]+x[0], x[0]+x[1]*x[1]-ExactFloat(0.875)};
        ARIADNE_TEST_PRINT(g);
        Box C = Box{{0.0,inf},{0.0,inf}};
        ARIADNE_TEST_PRINT(C);

        IntervalVector x_optimal=optimiser->minimise(f,D,g,C);
        ARIADNE_TEST_BINARY_PREDICATE(subset,x_optimal,D);
        ARIADNE_TEST_BINARY_PREDICATE(subset,g(x_optimal),C);
        //ARIADNE_TEST_LESS(norm(x_optimal),1e-6);
    }

    void test_mixed_constrained_optimisation() {
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(3);
        EffectiveScalarFunction f(+(sqr(x[0])+sqr(x[1])+x[1]*x[2]));
        ARIADNE_TEST_PRINT(f);
        Box D = Box{{-1.0,2.0},{-3.0,5.0},{1.25,2.25}};
        ARIADNE_TEST_PRINT(D);
        EffectiveScalarFunction g = x[0]*x[1]-x[0]*ExactFloat(1.25);
        EffectiveVectorFunction h = {ExactFloat(1.5)+x[0]+2*x[1]+ExactFloat(0.25)*x[0]*x[1]};
        EffectiveVectorFunction gh=join(g,h);
        ARIADNE_TEST_PRINT(gh);
        Box C = Box {{-1.0,-0.5},{0.0,0.0}};
        ARIADNE_TEST_PRINT(C);

        IntervalVector x_optimal=optimiser->minimise(f,D,gh,C);
        ARIADNE_TEST_LESS(norm(h(x_optimal)),1e-8);
    }

    void test_linear_feasibility() {
        // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(2);
        EffectiveVectorFunction g=EffectiveVectorFunction(1u, 2*x[0]+x[1]);
        ARIADNE_TEST_PRINT(g);
        Box D = Box{{0.0,2.0},{0.0,2.0}};
        Box C = Box{{-2.0,1.0}};

        ARIADNE_TEST_ASSERT(optimiser->feasible(D,g,C));
        C=Box{{1.0,1.5}};
        ARIADNE_TEST_ASSERT(optimiser->feasible(D,g,C));
        D=Box{{1.0,1.5},{0.5,1.0}};
        ARIADNE_TEST_ASSERT(!optimiser->feasible(D,g,C));
    }

    void test_nonlinear_feasibility() {
        // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(2);
        EffectiveVectorFunction g = {2*x[0]+x[1]+x[0]*x[1]/8};
        ARIADNE_TEST_PRINT(g);
        Box D = Box{{0.0,2.0},{0.0,2.0}};
        Box C = Box{{-2.0,1.0}};

        ARIADNE_TEST_ASSERT(optimiser->feasible(D,g,C));
        C=Box{{1.0,1.5}};
        ARIADNE_TEST_ASSERT(optimiser->feasible(D,g,C));
        D=Box{{1.0,1.5},{0.5,1.0}};
        ARIADNE_TEST_ASSERT(!optimiser->feasible(D,g,C));
    }

    void test_nonlinear_equality_feasibility() {
        // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(2);
        EffectiveVectorFunction h = { 2*x[0]-x[1]+x[0]*x[1]/8 };
        ARIADNE_TEST_PRINT(h);
        Box D = Box{{0.0,2.0},{0.0,2.0}};
        Box C = Box{{0.0,0.0}};

        ARIADNE_TEST_ASSERT(optimiser->feasible(D,h,C));
    }

    void test_feasibility_check() {
        EffectiveVectorFunction x=EffectiveVectorFunction::identity(2);
        ARIADNE_TEST_CONSTRUCT( EffectiveVectorFunction, g, ({sqr(x[0])+2*sqr(x[1])-1}) );
        ARIADNE_TEST_CONSTRUCT( IntervalVector, D, ({{-1.0, 1.0},{-1.0,1.0}}) );
        ARIADNE_TEST_CONSTRUCT( IntervalVector, C, ({{0.0,0.0}}) );

        ARIADNE_TEST_CONSTRUCT( IntervalVector, X1, ({{0.30,0.40},{0.60,0.70}}) );
        ARIADNE_TEST_ASSERT( definitely(optimiser->contains_feasible_point(D,g,C,X1)) );

        // The following test fails since it is difficult to find the feasible
        // point in the box.
        ARIADNE_TEST_CONSTRUCT( IntervalVector, X2, ({{0.30,0.40},{0.65,0.65}}) );
        ARIADNE_TEST_ASSERT( optimiser->contains_feasible_point(D,g,C,X2) );

        ARIADNE_TEST_CONSTRUCT( IntervalVector, X3, ({{0.30,0.40},{0.65,0.68}}) );
        ARIADNE_TEST_ASSERT( definitely(optimiser->contains_feasible_point(D,g,C,X3)) );

        ARIADNE_TEST_CONSTRUCT(FloatVector, x2, ({0.35,0.655}) );
        ARIADNE_TEST_ASSERT( optimiser->validate_feasibility(D,g,C,x2) );
    }

};

int main(int argc, const char* argv[]) {
    uint optimiser_verbosity = get_verbosity(argc,argv);

    NonlinearInfeasibleInteriorPointOptimiser nlio;
    nlio.verbosity=optimiser_verbosity;
    TestOptimiser(nlio).test();
    return ARIADNE_TEST_FAILURES;
    NonlinearInteriorPointOptimiser nlo;
    nlo.verbosity=optimiser_verbosity;
    TestOptimiser(nlo).test();

    ApproximateOptimiser appo;
    appo.verbosity=optimiser_verbosity;
    TestOptimiser(appo).test_nonlinear_equality_feasibility();

    IntervalOptimiser ivlo;
    ivlo.verbosity=optimiser_verbosity;
    //TestOptimiser(ivlo).test_nonlinear_equality_feasibility();
    return ARIADNE_TEST_FAILURES;
}

