/***************************************************************************
 *            test_nonlinear_programming.cpp
 *
 *  Copyright  2010-20  Pieter Collins
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
#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/taylor_model.hpp"
#include "function/formula.hpp"
#include "solvers/nonlinear_programming.hpp"
#include "geometry/box.hpp"
#include "io/command_line_interface.hpp"

using namespace std;
using namespace Ariadne;


class TestFeasibilityChecker
{
  private:
    std::unique_ptr<FeasibilityCheckerInterface> feasibility_checker;
    DoublePrecision pr;
  public:
    TestFeasibilityChecker(const FeasibilityCheckerInterface& fc)
        : feasibility_checker(fc.clone()) { }

    Void test() {
        ARIADNE_TEST_CALL(test_feasibility_check());
    }

    Void test_feasibility_check() {
        EffectiveVectorMultivariateFunction x=EffectiveVectorMultivariateFunction::identity(2);
        ARIADNE_TEST_CONSTRUCT( EffectiveVectorMultivariateFunction, g, ({sqr(x[0])+2*sqr(x[1])-1}) );
        ARIADNE_TEST_CONSTRUCT( ExactBoxType, D, ({{-1.0_x, 1.0_x},{-1.0_x,1.0_x}}) );
        ARIADNE_TEST_CONSTRUCT( ExactBoxType, C, ({{0.0_x,0.0_x}}) );

        ARIADNE_TEST_CONSTRUCT( ValidatedFeasibilityProblem, p, (D,g,C) );

        ARIADNE_TEST_CONSTRUCT( FloatDPVector, x1, ({1,0},pr) );
        ARIADNE_TEST_ASSERT( definitely(feasibility_checker->is_feasible_point(p,x1)) );
        ARIADNE_TEST_CONSTRUCT( FloatDPVector, x2, ({0.46875_x,0.62500_x},pr) );
        ARIADNE_TEST_ASSERT( not possibly(feasibility_checker->is_feasible_point(p,x2)) );

        ARIADNE_TEST_CONSTRUCT( FloatDPApproximation, e2l, (0.000732421875,pr) );
        ARIADNE_TEST_ASSERT( not probably(feasibility_checker->almost_feasible_point(p,x2,e2l)) );
        ARIADNE_TEST_CONSTRUCT( FloatDPApproximation, e2u, (0.001220703125,pr) );
        ARIADNE_TEST_ASSERT( probably(feasibility_checker->almost_feasible_point(p,x2,e2u)) );

        ARIADNE_TEST_CONSTRUCT( FloatDPBoundsVector, x3, ({{0.296875_x,0.406250_x},{0.593750_x,0.703125_x}},pr) );
        ARIADNE_TEST_CONSTRUCT( FloatDPVector, y3, (1, {1,pr}) );
        ARIADNE_TEST_CONSTRUCT( UpperBoxType, X3, ({{0.296875_x,0.406250_x},{0.593750_x,0.703125_x}}) );
        ARIADNE_TEST_ASSERT( definitely(feasibility_checker->contains_feasible_point(p,X3)) );
        ARIADNE_TEST_ASSERT( definitely(feasibility_checker->check_feasibility(p,x3,y3)) );
        ARIADNE_TEST_ASSERT( feasibility_checker->validate_feasibility(p,x3) );
        ARIADNE_TEST_ASSERT( feasibility_checker->validate_feasibility(p.g,x3) );
        ARIADNE_TEST_ASSERT( not feasibility_checker->validate_infeasibility(p,y3) );
        ARIADNE_TEST_ASSERT( not feasibility_checker->validate_infeasibility(p,x3,y3) );

        ARIADNE_TEST_CONSTRUCT( FloatDPBoundsVector, x4, ({{0.328125_x,0.406250_x},{0.671875_x,0.703125_x}},pr) );
        ARIADNE_TEST_CONSTRUCT( FloatDPVector, y4, (1, {1,pr}) );
        ARIADNE_TEST_CONSTRUCT( UpperBoxType, X4, ({{0.328125_x,0.406250_x},{0.671875_x,0.703125_x}}) );
        ARIADNE_TEST_PRINT( feasibility_checker->contains_feasible_point(p,X4) );
        ARIADNE_TEST_PRINT( feasibility_checker->check_feasibility(p,x4,y4) );
        ARIADNE_TEST_ASSERT( not definitely(feasibility_checker->contains_feasible_point(p,X4)) );
        ARIADNE_TEST_ASSERT( not definitely(feasibility_checker->contains_feasible_point(p,X4)) );
        ARIADNE_TEST_ASSERT( not feasibility_checker->validate_feasibility(p,x4) );
        ARIADNE_TEST_ASSERT( not feasibility_checker->validate_feasibility(p.g,x4) );
        ARIADNE_TEST_ASSERT( not definitely(feasibility_checker->check_feasibility(p,x4,y4)) );

        ARIADNE_TEST_CONSTRUCT( FloatDPBoundsVector, X5, ({{0.296875_x,0.406250_x},{0.656250_x,0.656250_x}},pr) );
        ARIADNE_TEST_ASSERT( feasibility_checker->contains_feasible_point(p,UpperBoxType(X5.array())) );

        ARIADNE_TEST_CONSTRUCT( FloatDPBoundsVector, X6, ({{0.296875_x,0.406250_x},{0.656250_x,0.687500_x}},pr) );
        ARIADNE_TEST_ASSERT( definitely(feasibility_checker->contains_feasible_point(p,UpperBoxType(X6.array()))) );

        ARIADNE_TEST_CONSTRUCT(FloatDPVector, x7, ({0.343750_x,0.656250_x},pr) );
        ARIADNE_TEST_ASSERT( not feasibility_checker->validate_feasibility(p,x7) );
    }

};

class TestOptimiser
{
  private:
    std::unique_ptr<OptimiserInterface> optimiser;
    DoublePrecision pr;
  public:
    TestOptimiser(const OptimiserInterface& opt)
        : optimiser(opt.clone()) { }

    Void test() {
        ARIADNE_TEST_CALL(test_unconstrained_optimisation());
        ARIADNE_TEST_CALL(test_constrained_optimisation());
        ARIADNE_TEST_CALL(test_equality_constrained_optimisation());
        ARIADNE_TEST_CALL(test_approximate_inequality_constrained_optimisation());
        ARIADNE_TEST_CALL(test_linear_feasibility());
        ARIADNE_TEST_CALL(test_nonlinear_feasibility());
        ARIADNE_TEST_CALL(test_nonlinear_equality_feasibility());
    }

    Void test_unconstrained_optimisation() {
        // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        EffectiveScalarMultivariateFunction x0s = sqr(x[0]);
        EffectiveScalarMultivariateFunction f(x0s*(12+x0s*(Decimal(6.3)+x0s))+6*x[1]*(x[1]-x[0]));
        ARIADNE_TEST_PRINT(f);
        EffectiveVectorMultivariateFunction g(0u,2u);
        ARIADNE_TEST_PRINT(g);
        ExactBoxType D=ExactBoxType{{-1,2},{-3,5}};
        ExactBoxType C=ExactBoxType{};
        ARIADNE_TEST_PRINT(Ariadne::make_tuple(f,D,g,C));

        FloatDPBoundsVector x_optimal=optimiser->minimise(f,D,g,C);
        ARIADNE_TEST_BINARY_PREDICATE(element,x_optimal,D);
        ARIADNE_TEST_BINARY_PREDICATE(element,g(x_optimal),C);
        //ExactDouble required_accuracy=1e-8_pr;
        //ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
    }

    Void test_equality_constrained_optimisation() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        EffectiveScalarMultivariateFunction f=(sqr(x[0])+sqr(x[1]));
        ARIADNE_TEST_PRINT(f);
        Real a(1.5_x); Real b(0.25_x);
        EffectiveVectorMultivariateFunction g={a+x[0]+2*x[1]+b*x[0]*x[1]};
        ARIADNE_TEST_PRINT(g);
        ExactBoxType C={{0.0_x,0.0_x}};
        ExactBoxType D=ExactBoxType{{-1.0_x,2.0_x},{-3.0_x,5.0_x}};
        ARIADNE_TEST_PRINT(Ariadne::make_tuple(f,D,g,C));

        ExactDouble required_accuracy=1e-7_pr;
        FloatDPBoundsVector x_optimal=optimiser->minimise(f,D,g,C);
        ARIADNE_TEST_BINARY_PREDICATE(element,x_optimal,D);
        ARIADNE_TEST_LESS(norm(g(x_optimal)),required_accuracy);
    }

    Void test_approximate_inequality_constrained_optimisation() {
        List<ApproximateScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(3);
        ApproximateScalarMultivariateFunction x0s = sqr(x[0]);
        ApproximateScalarMultivariateFunction f = x0s*(12+x0s*(Decimal(6.3)+x0s))+6*x[1]*(x[1]-x[0])+x[2];
        ARIADNE_TEST_PRINT(f);
        //EffectiveVectorMultivariateFunction g( (x[0]-1, x[0]+x[1]*x[1], x[1]*x[1]) );
        ApproximateBoxType D = ApproximateBoxType{{-1.0_x,2.0_x},{-3.0_x,5.0_x},{-3.0_x,5.0_x}};
        ARIADNE_TEST_PRINT(D);
        ApproximateVectorMultivariateFunction g = {2*x[1]+x[0], x[0]+x[1]*x[1]-Real(0.875_x)};
        ARIADNE_TEST_PRINT(g);
        ApproximateBoxType C = ApproximateBoxType{{0.0_x,4.0_x},{0.0_x,4.0_x}};
        ARIADNE_TEST_PRINT(C);

        ARIADNE_TEST_PRINT(optimiser);
        FloatDPApproximationVector x_optimal=optimiser->minimise(f,D,g,C);
        ARIADNE_TEST_PRINT(x_optimal);
        ARIADNE_TEST_BINARY_PREDICATE(element,x_optimal,D);
        ARIADNE_TEST_BINARY_PREDICATE(element,g(x_optimal),C);
        //ExactDouble required_accuracy=1e-6_pr;
        //ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
    }

    Void test_constrained_optimisation() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(3);
        EffectiveScalarMultivariateFunction x0s = sqr(x[0]);
        EffectiveScalarMultivariateFunction f = x0s*(12+x0s*(Decimal(6.3)+x0s))+6*x[1]*(x[1]-x[0])+x[2];
        ARIADNE_TEST_PRINT(f);
        //EffectiveVectorMultivariateFunction g( (x[0]-1, x[0]+x[1]*x[1], x[1]*x[1]) );
        ExactBoxType D = ExactBoxType{{-1.0_x,2.0_x},{-3.0_x,5.0_x},{-3.0_x,5.0_x}};
        ARIADNE_TEST_PRINT(D);
        EffectiveVectorMultivariateFunction g = {2*x[1]+x[0], x[0]+x[1]*x[1]-Real(0.875_x)};
        ARIADNE_TEST_PRINT(g);
        ExactBoxType C = ExactBoxType{{0.0_x,inf},{0.0_x,inf}};
        ARIADNE_TEST_PRINT(C);

        FloatDPBoundsVector x_optimal=optimiser->minimise(f,D,g,C);
        ARIADNE_TEST_BINARY_PREDICATE(element,x_optimal,D);
        ARIADNE_TEST_BINARY_PREDICATE(element,g(x_optimal),C);
        //ExactDouble required_accuracy=1e-6_pr;
        //ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
    }

    Void test_mixed_constrained_optimisation() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(3);
        EffectiveScalarMultivariateFunction f(+(sqr(x[0])+sqr(x[1])+x[1]*x[2]));
        ARIADNE_TEST_PRINT(f);
        ExactBoxType D = ExactBoxType{{-1.0_x,2.0_x},{-3.0_x,5.0_x},{1.25_x,2.25_x}};
        ARIADNE_TEST_PRINT(D);
        EffectiveScalarMultivariateFunction g = x[0]*x[1]-x[0]*Real(1.25_x);
        EffectiveVectorMultivariateFunction h = {Real(1.5_x)+x[0]+2*x[1]+Real(0.25_x)*x[0]*x[1]};
        EffectiveVectorMultivariateFunction gh=join(g,h);
        ARIADNE_TEST_PRINT(gh);
        ExactBoxType C = ExactBoxType {{-1.0_x,-0.5_x},{0.0_x,0.0_x}};
        ARIADNE_TEST_PRINT(C);

        FloatDPBoundsVector x_optimal=optimiser->minimise(f,D,gh,C);
        ExactDouble required_accuracy=1e-8_pr;
        ARIADNE_TEST_LESS(norm(h(x_optimal)),required_accuracy);
    }

    Void test_linear_feasibility() {
        // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        EffectiveVectorMultivariateFunction g=EffectiveVectorMultivariateFunction(1u, 2*x[0]+x[1]);
        ARIADNE_TEST_PRINT(g);
        ExactBoxType D = ExactBoxType{{0.0_x,2.0_x},{0.0_x,2.0_x}};
        ExactBoxType C = ExactBoxType{{-2.0_x,1.0_x}};

        ARIADNE_TEST_ASSERT(optimiser->feasible(D,g,C));
        C=ExactBoxType{{1.0_x,1.5_x}};
        ARIADNE_TEST_ASSERT(optimiser->feasible(D,g,C));
        D=ExactBoxType{{1.0_x,1.5_x},{0.5_x,1.0_x}};
        ARIADNE_TEST_ASSERT(!optimiser->feasible(D,g,C));
    }

    Void test_nonlinear_feasibility() {
        // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        EffectiveVectorMultivariateFunction g = {2*x[0]+x[1]+x[0]*x[1]/8};
        ARIADNE_TEST_PRINT(g);
        ExactBoxType D = ExactBoxType{{0.0_x,2.0_x},{0.0_x,2.0_x}};
        ExactBoxType C = ExactBoxType{{-2.0_x,1.0_x}};

        ARIADNE_TEST_ASSERT(optimiser->feasible(D,g,C));
        C=ExactBoxType{{1.0_x,1.5_x}};
        ARIADNE_TEST_ASSERT(optimiser->feasible(D,g,C));
        D=ExactBoxType{{1.0_x,1.5_x},{0.5_x,1.0_x}};
        ARIADNE_TEST_ASSERT(!optimiser->feasible(D,g,C));
    }

    Void test_nonlinear_equality_feasibility() {
        // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        EffectiveVectorMultivariateFunction h = { 2*x[0]-x[1]+x[0]*x[1]/8 };
        ARIADNE_TEST_PRINT(h);
        ExactBoxType D = ExactBoxType{{0.0_x,2.0_x},{0.0_x,2.0_x}};
        ExactBoxType C = ExactBoxType{{0.0_x,0.0_x}};

        ARIADNE_TEST_ASSERT(optimiser->feasible(D,h,C));
    }

};

Int main(Int argc, const char* argv[]) {
    if (not CommandLineInterface::instance().acquire(argc,argv)) return -1;

    FeasibilityChecker fc;
    TestFeasibilityChecker(fc).test();

    InteriorPointOptimiser nlo;
    ARIADNE_TEST_CLASS("InteriorPointOptimiser",TestOptimiser(nlo));

    InfeasibleInteriorPointOptimiser nlio;
    ARIADNE_TEST_CLASS("InfeasibleInteriorPointOptimiser",TestOptimiser(nlio));
/*
    ApproximateOptimiser appo;
    ARIADNE_TEST_CLASS("ApproximateOptimiser",TestOptimiser(appo));

    IntervalOptimiser ivlo;
    ARIADNE_TEST_CLASS("IntervalOptimiser",TestOptimiser(ivlo));
*/
    return ARIADNE_TEST_FAILURES;
}

