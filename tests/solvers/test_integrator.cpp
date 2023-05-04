/***************************************************************************
 *            test_integrator.cpp
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

#include "solvers/integrator.hpp"
#include "algebra/vector.hpp"
#include "algebra/algebra.hpp"
#include "function/polynomial.hpp"
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "function/formula.hpp"
#include "io/command_line_interface.hpp"

#include "pronest/configurable.tpl.hpp"
#include "pronest/configuration_property.tpl.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

struct UnsafeReal : Real { UnsafeReal(double d) : Real(operator""_x(d)) { } };

class TestIntegrator
{
    typedef Vector<ExactIntervalType> ExactIntervalVectorType;
  private:
    std::unique_ptr<IntegratorInterface> integrator_ptr;
    std::unique_ptr<BounderInterface> bounder_ptr;
  public:
    TestIntegrator(IntegratorInterface const& i)
        : integrator_ptr(i.clone()), bounder_ptr(new EulerBounder(Configuration<EulerBounder>()))
    { }

    Int test() {
        ARIADNE_TEST_PRINT(integrator_ptr.operator->());
        ARIADNE_TEST_PRINT(*integrator_ptr);
        ARIADNE_TEST_CALL(test_constant_derivative());
        ARIADNE_TEST_CALL(test_quadratic_flow());
        ARIADNE_TEST_CALL(test_linear());
        ARIADNE_TEST_CALL(test_spiral());
        ARIADNE_TEST_CALL(test_logistic());
        ARIADNE_TEST_CALL(test_time_variant());
        ARIADNE_TEST_CALL(test_time_variant_with_parameters());
        return 0;
    }

    Int test_affine() {
        ARIADNE_TEST_CALL(test_constant_derivative());
        ARIADNE_TEST_CALL(test_quadratic_flow());
        ARIADNE_TEST_CALL(test_linear());
        ARIADNE_TEST_CALL(test_spiral());
        return 0;
    }

    Void test_constant_derivative() {

        EffectiveScalarMultivariateFunction o=EffectiveScalarMultivariateFunction::constant(2,1);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction y=EffectiveScalarMultivariateFunction::coordinate(2,1);
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(3,0);
        EffectiveScalarMultivariateFunction y0=EffectiveScalarMultivariateFunction::coordinate(3,1);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(3,2);

        EffectiveVectorMultivariateFunction f={o*2,o*3};
        ARIADNE_TEST_PRINT(f);
        ExactBoxType d={{0.0_x,1.0_x},{-0.5_x,1.5_x}};
        StepSizeType h=0.25_x;
        ValidatedVectorMultivariateFunctionPatch flow=integrator_ptr->flow_step(f,d,h);
        EffectiveVectorMultivariateFunction expected_flow={x0+2*t,y0+3*t};
        ARIADNE_TEST_PRINT(flow);
        ARIADNE_TEST_PRINT(expected_flow);
        ARIADNE_TEST_PRINT(flow.errors());
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow-expected_flow),1e-8_pr);
    }

    Void test_quadratic_flow() {

        EffectiveScalarMultivariateFunction o=EffectiveScalarMultivariateFunction::constant(2,1);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction y=EffectiveScalarMultivariateFunction::coordinate(2,1);
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(3,0);
        EffectiveScalarMultivariateFunction y0=EffectiveScalarMultivariateFunction::coordinate(3,1);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(3,2);

        EffectiveVectorMultivariateFunction f={o,x};
        ExactBoxType d={ExactIntervalType(0.0_x,1.0_x),ExactIntervalType(-0.5_x,1.5_x)};
        StepSizeType h=0.25_x;
        ValidatedVectorMultivariateFunctionPatch flow=integrator_ptr->flow_step(f,d,h);
        EffectiveVectorMultivariateFunction expected_flow={x0+t,y0+x0*t+t*t/2};
        ARIADNE_TEST_PRINT(f);
        ARIADNE_TEST_PRINT(flow);
        ARIADNE_TEST_PRINT(expected_flow);
        ARIADNE_TEST_PRINT(flow.errors());
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow-expected_flow),1e-8_pr);
    }

    Void test_linear() {

        EffectiveScalarMultivariateFunction o=EffectiveScalarMultivariateFunction::constant(2,1);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction y=EffectiveScalarMultivariateFunction::coordinate(2,1);
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(3,0);
        EffectiveScalarMultivariateFunction y0=EffectiveScalarMultivariateFunction::coordinate(3,1);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(3,2);

        EffectiveVectorMultivariateFunction f={x,-y};
        ExactBoxType d={ExactIntervalType(-0.25_x,0.25_x),ExactIntervalType(-0.25_x,0.25_x)};
        StepSizeType h=0.25_x;
        ValidatedVectorMultivariateFunctionPatch flow=integrator_ptr->flow_step(f,d,h);
        EffectiveVectorMultivariateFunction expected_flow={x0*(1+t+t*t/2+t*t*t/6+t*t*t*t/24),y0*(1-t+t*t/2-t*t*t/6+t*t*t*t/24)};
        ARIADNE_TEST_PRINT(f);
        ARIADNE_TEST_PRINT(flow);
        ARIADNE_TEST_PRINT(expected_flow);
        ARIADNE_TEST_PRINT(flow.errors());
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow-expected_flow),1e-4_pr);
    };

    Void test_spiral() {

        EffectiveScalarMultivariateFunction o=EffectiveScalarMultivariateFunction::constant(2,1);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction y=EffectiveScalarMultivariateFunction::coordinate(2,1);
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(3,0);
        EffectiveScalarMultivariateFunction y0=EffectiveScalarMultivariateFunction::coordinate(3,1);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(3,2);

        Real half(0.5_x);
        EffectiveVectorMultivariateFunction f={-half*x-y,x-half*y};
        ExactBoxType d={ExactIntervalType(0.75_x,1.25_x),ExactIntervalType(-0.25_x,0.25_x)};
        StepSizeType h=0.25_x;
        ValidatedVectorMultivariateFunctionPatch flow=integrator_ptr->flow_step(f,d,h);
        EffectiveVectorMultivariateFunction expected_flow={exp(-half*t)*(x0*cos(t)-y0*sin(t)),exp(-half*t)*(x0*sin(t)+y0*cos(t))};
        ARIADNE_TEST_PRINT(f);
        ARIADNE_TEST_PRINT(flow);
        ARIADNE_TEST_PRINT(expected_flow);
        ARIADNE_TEST_PRINT(flow.errors());
        ARIADNE_TEST_PRINT(flow-expected_flow);
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow-expected_flow),1e-3_pr);

    };

    // Equation: dx/dt = x*(1-x)
    // Solution:  x = 1/(1+(1/x0-1)*exp(-t))
    Void test_logistic() {

        EffectiveScalarMultivariateFunction o=EffectiveScalarMultivariateFunction::constant(1,1);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(1,0);
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(2,1);

        EffectiveVectorMultivariateFunction f={x*(o-x)};
        ExactBoxType d={ExactIntervalType(0.25_x,0.5_x)};
        StepSizeType h=0.5_x;
        ValidatedVectorMultivariateFunctionPatch flow=integrator_ptr->flow_step(f,d,h);
        ValidatedVectorMultivariateTaylorFunctionModelDP taylor_flow=dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP&>(flow.reference());
        EffectiveVectorMultivariateFunction expected_flow={1/(1+(1/x0-1)*exp(-t))};
        ARIADNE_TEST_PRINT(*integrator_ptr);
        ARIADNE_TEST_PRINT(f);
        ARIADNE_TEST_PRINT(flow);
        ARIADNE_TEST_PRINT(expected_flow);
        ARIADNE_TEST_PRINT(taylor_flow.errors());
        ARIADNE_TEST_PRINT(taylor_flow-expected_flow);
        ARIADNE_TEST_BINARY_PREDICATE(operator<,taylor_flow.error(),0.01_pr);
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(taylor_flow-expected_flow),0.014_pr);
    };

    Void test_time_variant() {

        // dx/dt=a*x+b*t+c; x=(x0+b/a*t0+(b/a+c)/a)*exp(a*(t-t0))-b/a*t-(b/a+c)/a
        // dx/dt=b*t+c; x=x0+c*(t-t0)+b/2*(t^2-t0^2)
        // dx/dt=c; x=x0+c*(t-t0)

        EffectiveScalarMultivariateFunction o=EffectiveScalarMultivariateFunction::constant(4,1);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(4,0);
        EffectiveScalarMultivariateFunction y=EffectiveScalarMultivariateFunction::coordinate(4,1);
        EffectiveScalarMultivariateFunction z=EffectiveScalarMultivariateFunction::coordinate(4,2);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(4,3);

        Rational a=2; Rational b=3; Rational c=5;
        EffectiveVectorMultivariateFunction f={a*x+b*t+c,b*t+c,c*o};
        //EffectiveVectorMultivariateFunction f={a*x+b*t,b*t+c};
        ExactBoxType d={ExactIntervalType(-0.5_x,1.5_x),ExactIntervalType(-0.5_x,2.5_x),ExactIntervalType(0.0_x,1.0_x)};
        StepSizeType t0=3.0_x;
        StepSizeType hsug=0.0625_x;
        Pair<StepSizeType,UpperBoxType> step_bounds = bounder_ptr->compute(f,d,t0,ExactBoxType(0u),hsug);
        StepSizeType h = step_bounds.first;
        UpperBoxType B = step_bounds.second;
        ValidatedVectorMultivariateFunctionPatch flow=integrator_ptr->flow_step(f,d,Interval<StepSizeType>(t0,t0+h),ExactBoxType(0u),B);
        EffectiveVectorMultivariateFunction expected_flow={(x+b/a*t0+(b/a+c)/a)*exp(a*(t-t0))-b/a*t-(b/a+c)/a,y+c*(t-t0)+b/2*(t*t-t0*t0),z+c*(t-t0)};
        ARIADNE_TEST_PRINT(f);
        ARIADNE_TEST_PRINT(flow);
        ARIADNE_TEST_PRINT(expected_flow);
        ARIADNE_TEST_PRINT(flow.errors());
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow-expected_flow),1e-4_pr);
    };

    Void test_time_variant_with_parameters() {
        // Equation: dx/dt=a*x+b*t+c;
        // Solution: x=c*exp(a*t)-b/a*t-b/a^2
        //     where c=(x0+b/a*t0+b/a^2)/exp(a*t0)
        //         since dx/dt=a*c*exp(a*t)-b/a
        //             a*x+b*t=a*c*exp(a*t)-b*t-b/a + b*t
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(5,0);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(5,1);
        EffectiveScalarMultivariateFunction a=EffectiveScalarMultivariateFunction::coordinate(5,2);
        EffectiveScalarMultivariateFunction b=EffectiveScalarMultivariateFunction::coordinate(5,3);
        EffectiveScalarMultivariateFunction c=EffectiveScalarMultivariateFunction::coordinate(5,4);
        EffectiveVectorMultivariateFunction f={a*x+b*t+c};
        ExactBoxType domx={ExactIntervalType(0.5_x,1.0_x)};
        StepSizeType t0=0.0_x;
        StepSizeType hsug=0.0625_x;
        ExactBoxType domp={ExactIntervalType{2.0_x,2.5_x},ExactIntervalType{-0.5_x,1.0_x},ExactIntervalType{0.5_x,1.0_x}};
        Pair<StepSizeType,UpperBoxType> step_bounds = bounder_ptr->compute(f,domx,t0,domp,hsug);
        StepSizeType h = step_bounds.first;
        UpperBoxType B = step_bounds.second;
        Interval<StepSizeType> domt(t0,t0+h);
        ValidatedVectorMultivariateFunctionPatch flow=integrator_ptr->flow_step(f,domx,domt,domp,B);
        EffectiveVectorMultivariateFunction expected_flow={(x+b/a*(t0)+(b/a+c)/a)*exp(a*(t-t0))-b/a*t-(b/a+c)/a};

        ARIADNE_TEST_PRINT(f);
        ARIADNE_TEST_PRINT(flow);
        ARIADNE_TEST_PRINT(expected_flow);
        ARIADNE_TEST_PRINT(flow-expected_flow);
        ARIADNE_TEST_PRINT(flow.errors());
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow-expected_flow),1e-4_pr);
    }
};

Int main(Int argc, const char* argv[]) {
    if (not CommandLineInterface::instance().acquire(argc,argv)) return -1;

    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),1e-10);

    TaylorPicardIntegrator taylor_picard_integrator(Configuration<TaylorPicardIntegrator>().set_sweeper(sweeper).set_step_maximum_error(1e-8));
    ARIADNE_TEST_CLASS("TaylorPicardIntegrator",TestIntegrator(taylor_picard_integrator))


    GradedTaylorPicardIntegrator unbounded_taylor_picard_integrator(Configuration<GradedTaylorPicardIntegrator>().set_sweeper(sweeper).set_step_maximum_error(1e-8).set_temporal_order(8));
    ARIADNE_TEST_CLASS("GradedTaylorPicardIntegrator",TestIntegrator(unbounded_taylor_picard_integrator))

    TaylorSeriesIntegrator taylor_series_integrator(Configuration<TaylorSeriesIntegrator>().set_sweeper(sweeper).set_order(6));
    ARIADNE_TEST_CLASS("TaylorSeriesIntegrator",TestIntegrator(taylor_series_integrator))

    //GradedTaylorSeriesIntegrator graded_taylor_series_integrator(Configuration<GradedTaylorSeriesIntegrator>().set_sweeper(sweeper).set_step_maximum_error(1e-7).set_maximum_temporal_order(8));
    //ARIADNE_TEST_CLASS("GradedTaylorSeriesIntegrator",TestIntegrator(graded_taylor_series_integrator))

    ARIADNE_PRINT_TEST_CASE_TITLE("AffineIntegrator");
    //AffineIntegrator affine_integrator(Configuration<AffineIntegrator>().set_sweeper(sweeper).set_spacial_order(2).set_temporal_order(4));
    //TestIntegrator(affine_integrator).test_affine();

    return ARIADNE_TEST_FAILURES;
}
