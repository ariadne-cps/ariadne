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

decltype(auto) operator<(FloatDPApproximation x, ApproximateDouble y) {
    return x < FloatDPApproximation(y,dp);
}

bool consistent(FloatDPBounds x, Rational y) {
    return x.lower_raw()<=y && y <= x.upper_raw();
}
bool consistent(FloatDPBoundsVector x, RationalVector y) {
    ARIADNE_PRECONDITION(x.size()==y.size());
    for (SizeType i=0; i!=x.size(); ++i) {
        if (not consistent(x[i],y[i])) {
            return false;
        }
    }
    return true;
}
template<class X, class Y> decltype(auto) inconsistent(X const& x, Y const& y) {
    return not consistent(x,y);
}

template<class A, class X, class E> decltype(auto) within(A a, X x, E e) {
    return decide(abs(a-x)<e);
}

class TestOptimiser
{
  private:
    std::unique_ptr<OptimiserInterface> optimiser;
    DoublePrecision pr;
  public:
    TestOptimiser(const OptimiserInterface& opt)
        : optimiser(opt.clone()) { }

    Void test() {
        ARIADNE_TEST_CALL(test_nonconvex_one_dimensional_optimisation());
        ARIADNE_TEST_CALL(test_nonconvex_unconstrained_optimisation());
#warning
        return;
        ARIADNE_TEST_CALL(test_feasibility_check());
        ARIADNE_TEST_CALL(test_unconstrained_optimisation());
        ARIADNE_TEST_CALL(test_nonconvex_one_dimensional_optimisation());
        ARIADNE_TEST_CALL(test_constrained_optimisation());
        ARIADNE_TEST_CALL(test_equality_constrained_optimisation());
        ARIADNE_TEST_CALL(test_nonconvex_unconstrained_optimisation());
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

    Void test_nonconvex_one_dimensional_optimisation() {
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

    Void test_nonconvex_unconstrained_optimisation() {
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction y=EffectiveScalarMultivariateFunction::coordinate(2,1);
        EffectiveScalarMultivariateFunction f=x*(1-x)+2*y*y+x*y+y/4;
//        EffectiveScalarMultivariateFunction f=x*x*(x-1)+2*y*y+x*y+y/4;
        ARIADNE_TEST_PRINT(f);
        ExactBoxType D=ExactBoxType{{-1.0_x,3.0_x},{-1.0_x,1.0_x}};
        ARIADNE_TEST_PRINT(D);
        EffectiveVectorMultivariateFunction g(0,2);
        ExactBoxType C;
        ARIADNE_TEST_PRINT(Ariadne::make_tuple(f,D,g,C));
        //FloatDPBoundsVector expected_x_critical=FloatDPBoundsVector({5,-2},dp)/12;
        RationalVector expected_x_critical=RationalVector({5,-2})/12;
        FloatDPBoundsVector expected_x_optimal({3,-0.8125_x},dp);

        ExactDouble required_accuracy=1e-2_pr;
        ARIADNE_TEST_CONSTRUCT(FloatDPBoundsVector,x_optimal,(optimiser->minimise(f,D,g,C)));
        ARIADNE_TEST_BINARY_PREDICATE(element,x_optimal,D);
        ARIADNE_TEST_BINARY_PREDICATE(consistent,x_optimal,expected_x_optimal);
        //ARIADNE_TEST_LESS(error(x_optimal),required_accuracy);

        ARIADNE_TEST_PRINT(FloatDPBoundsVector(expected_x_critical,dp));
        ARIADNE_TEST_BINARY_PREDICATE(inconsistent,x_optimal,expected_x_critical);
        ARIADNE_TEST_COMPARE(error(x_optimal-expected_x_critical),>=,required_accuracy);

        ARIADNE_TEST_PRINT(f(expected_x_optimal));
        ARIADNE_TEST_PRINT(f(x_optimal));
        ARIADNE_TEST_PRINT(gradient(f,x_optimal));
        ARIADNE_TEST_PRINT(hessian(f,x_optimal));


#warning
/*
        EuclideanDomain dom(2);
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(dom);
        EffectiveScalarMultivariateFunction f(sqr(x[0])*(x[0]+1)-5*x[1]*x[2]+sqr(x[1]));
        ARIADNE_TEST_PRINT(f);
        ApproximateBoxType D = ApproximateBoxType{{-1.0_x,2.0_x},{-3.0_x,5.0_x}};
        ARIADNE_TEST_PRINT(D);
        ApproximateBoxType C = {C0,C1};
        ARIADNE_TEST_PRINT(C);

//        g = {g1};
//        C = {C1};

        FloatDPApproximationVector  x_optimal=optimiser->minimise(f,D,g,C).x;
        ApproximateDouble required_accuracy=1e-3_pr;
        ARIADNE_TEST_BINARY_PREDICATE(element,x_optimal,D);
        ARIADNE_TEST_BINARY_PREDICATE(element,g0(x_optimal),C0);
//        ARIADNE_TEST_BINARY_PREDICATE(element,g1(x_optimal),C1);
        ARIADNE_TEST_TERNARY_PREDICATE(within,g1(x_optimal),c1,required_accuracy);
*/
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

    Void test_feasibility_check() {
        EffectiveVectorMultivariateFunction x=EffectiveVectorMultivariateFunction::identity(2);
        ARIADNE_TEST_CONSTRUCT( EffectiveVectorMultivariateFunction, g, ({sqr(x[0])+2*sqr(x[1])-1}) );
        ARIADNE_TEST_CONSTRUCT( ExactBoxType, D, ({{-1.0_x, 1.0_x},{-1.0_x,1.0_x}}) );
        ARIADNE_TEST_CONSTRUCT( ExactBoxType, C, ({{0.0_x,0.0_x}}) );

        ARIADNE_TEST_CONSTRUCT( FloatDPBoundsVector, X1, ({{0.296875_x,0.406250_x},{0.593750_x,0.703125_x}},pr) );
        ARIADNE_TEST_ASSERT( definitely(optimiser->contains_feasible_point(D,g,C,X1)) );

        // The following test fails since it is difficult to find the feasible
        // point in the box.
        ARIADNE_TEST_CONSTRUCT( FloatDPBoundsVector, X2, ({{0.296875_x,0.406250_x},{0.656250_x,0.656250_x}},pr) );
        ARIADNE_TEST_ASSERT( optimiser->contains_feasible_point(D,g,C,X2) );

        ARIADNE_TEST_CONSTRUCT( FloatDPBoundsVector, X3, ({{0.296875_x,0.406250_x},{0.656250_x,0.687500_x}},pr) );
        ARIADNE_TEST_ASSERT( definitely(optimiser->contains_feasible_point(D,g,C,X3)) );

        ARIADNE_TEST_CONSTRUCT(FloatDPVector, x2, ({0.343750_x,0.656250_x},pr) );
        ARIADNE_TEST_ASSERT( optimiser->validate_feasibility(D,g,C,x2) );
    }


};

class TestApproximateOptimiser
{
  private:
    std::unique_ptr<ApproximateOptimiserInterface> optimiser;
    DoublePrecision pr;
  public:
    TestApproximateOptimiser(const ApproximateOptimiserInterface& opt)
        : optimiser(opt.clone()) { }

    Void test_mixed_constrained_optimisation() {
        EuclideanDomain dom(3);
        List<ApproximateScalarMultivariateFunction> x=ApproximateScalarMultivariateFunction::coordinates(dom);
        ApproximateScalarMultivariateFunction f(+(sqr(x[0])+sqr(x[1])+x[1]*x[2]));
        ARIADNE_TEST_PRINT(f);
        ApproximateBoxType D = ApproximateBoxType{{-1.0_x,2.0_x},{-3.0_x,5.0_x},{1.25_x,2.25_x}};
        ARIADNE_TEST_PRINT(D);
        ApproximateScalarMultivariateFunction g0 = x[0]*x[1]-x[0]*Real(1.25_x);
        ApproximateScalarMultivariateFunction g1 = Real(1.5_x)+x[0]+2*x[1]+Real(0.25_x)*x[0]*x[1];
        ApproximateVectorMultivariateFunction g = {g0,g1};
        ARIADNE_TEST_PRINT(g);
        ApproximateIntervalType C0 = {-1.0_x,-0.5_x};
        ApproximateIntervalType C1 = {0.0_x,0.0_x};
        auto c1 = C1.lower_bound();
        ApproximateBoxType C = {C0,C1};
        ARIADNE_TEST_PRINT(C);

//        g = {g1};
//        C = {C1};

        FloatDPApproximationVector  x_optimal=optimiser->minimise(f,D,g,C).x;
        ApproximateDouble required_accuracy=1e-3_pr;
        ARIADNE_TEST_BINARY_PREDICATE(element,x_optimal,D);
        ARIADNE_TEST_BINARY_PREDICATE(element,g0(x_optimal),C0);
//        ARIADNE_TEST_BINARY_PREDICATE(element,g1(x_optimal),C1);
        ARIADNE_TEST_TERNARY_PREDICATE(within,g1(x_optimal),c1,required_accuracy);
    }
};


Int main(Int argc, const char* argv[]) {
    if (not CommandLineInterface::instance().acquire(argc,argv)) return -1;

    // C0BranchAndBoundOptimiser bbo(1./1024);
    // TestOptimiser(bbo).test();

    PrimalSlackPenaltyBarrierFunctionApproximateOptimiser pbappo;
    TestApproximateOptimiser(pbappo).test_mixed_constrained_optimisation();

    NonlinearInfeasibleInteriorPointOptimiser nlio;
    TestOptimiser(nlio).test();

    return ARIADNE_TEST_FAILURES;

    NonlinearInteriorPointOptimiser nlo;
    TestOptimiser(nlo).test();

    ApproximateOptimiser appo;
    TestOptimiser(appo).test_nonlinear_equality_feasibility();

    IntervalOptimiser ivlo;
    TestOptimiser(ivlo).test_nonlinear_equality_feasibility();
    return ARIADNE_TEST_FAILURES;
}

