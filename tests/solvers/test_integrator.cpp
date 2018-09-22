/***************************************************************************
 *            test_integrator.cpp
 *
 *  Copyright  2008-10  Pieter Collins
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

#include "function/polynomial.hpp"
#include "solvers/integrator.hpp"
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "algebra/vector.hpp"
#include "algebra/algebra.hpp"
#include "function/formula.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

inline EffectiveScalarFunction operator^(EffectiveScalarFunction f, Int m) { return pow(f,m); }
inline EffectiveScalarFunction operator*(double c, EffectiveScalarFunction f) { return Real(c)*f; }

struct UnsafeReal : Real { UnsafeReal(double d) : Real(d) { } };

class TestIntegrator
{
    typedef Vector<ExactIntervalType> ExactIntervalVectorType;
  private:
    std::unique_ptr<IntegratorInterface> integrator_ptr;
    EffectiveScalarFunction o,x,y,x0,y0,t;
  public:
    TestIntegrator(const IntegratorInterface& i)
        : integrator_ptr(i.clone())
    {
        o=EffectiveScalarFunction::constant(2,1);
        x=EffectiveScalarFunction::coordinate(2,0);
        y=EffectiveScalarFunction::coordinate(2,1);
        x0=EffectiveScalarFunction::coordinate(3,0);
        y0=EffectiveScalarFunction::coordinate(3,1);
        t=EffectiveScalarFunction::coordinate(3,2);
    }

    Int test() {
        ARIADNE_TEST_PRINT(integrator_ptr.operator->());
        ARIADNE_TEST_PRINT(*integrator_ptr);
        ARIADNE_TEST_CALL(test_constant_derivative());
        ARIADNE_TEST_CALL(test_quadratic_flow());
        ARIADNE_TEST_CALL(test_linear());
        ARIADNE_TEST_CALL(test_spiral());
        ARIADNE_TEST_CALL(test_logistic());
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
        EffectiveVectorFunction f={o*2,o*3};
        ARIADNE_TEST_PRINT(f);
        ExactBoxType d={ExactIntervalType(0.0,1.0),ExactIntervalType(-0.5,1.5)};
        FloatDP h=0.25;
        ValidatedVectorFunctionModelDP flow=integrator_ptr->flow_step(f,d,h);
        EffectiveVectorFunction expected_flow={x0+2*t,y0+3*t};
        ARIADNE_TEST_PRINT(flow);
        ARIADNE_TEST_PRINT(expected_flow);
        ARIADNE_TEST_PRINT(flow.errors());
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow-expected_flow),1e-8);
    }

    Void test_quadratic_flow() {
        EffectiveVectorFunction f={o,x};
        ExactBoxType d={ExactIntervalType(0.0,1.0),ExactIntervalType(-0.5,1.5)};
        FloatDP h=0.25;
        ValidatedVectorFunctionModelDP flow=integrator_ptr->flow_step(f,d,h);
        EffectiveVectorFunction expected_flow={x0+t,y0+x0*t+t*t/2};
        ARIADNE_TEST_PRINT(f);
        ARIADNE_TEST_PRINT(flow);
        ARIADNE_TEST_PRINT(expected_flow);
        ARIADNE_TEST_PRINT(flow.errors());
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow-expected_flow),1e-8);
    }

    Void test_linear() {
        EffectiveVectorFunction f={x,-y};
        ExactBoxType d={ExactIntervalType(-0.25,0.25),ExactIntervalType(-0.25,0.25)};
        FloatDP h=0.25;
        ValidatedVectorFunctionModelDP flow=integrator_ptr->flow_step(f,d,h);
        EffectiveVectorFunction expected_flow={x0*(1+t+t*t/2+t*t*t/6+t*t*t*t/24),y0*(1-t+t*t/2-t*t*t/6+t*t*t*t/24)};
        ARIADNE_TEST_PRINT(f);
        ARIADNE_TEST_PRINT(flow);
        ARIADNE_TEST_PRINT(expected_flow);
        ARIADNE_TEST_PRINT(flow.errors());
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow-expected_flow),1e-4);
    };

    Void test_spiral() {
        Real half(0.5);
        EffectiveVectorFunction f={-half*x-y,x-half*y};
        ExactBoxType d={ExactIntervalType(0.75,1.25),ExactIntervalType(-0.25,0.25)};
        FloatDP h=0.25;
        ValidatedVectorFunctionModelDP flow=integrator_ptr->flow_step(f,d,h);
        EffectiveVectorFunction expected_flow={exp(-half*t)*(x0*cos(t)-y0*sin(t)),exp(-half*t)*(x0*sin(t)+y0*cos(t))};
        ARIADNE_TEST_PRINT(f);
        ARIADNE_TEST_PRINT(flow);
        ARIADNE_TEST_PRINT(expected_flow);
        ARIADNE_TEST_PRINT(flow.errors());
        ARIADNE_TEST_PRINT(flow-expected_flow);
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow-expected_flow),1e-3);

    };

    Void test_logistic() {
        o=EffectiveScalarFunction::constant(1,1);
        x=EffectiveScalarFunction::coordinate(1,0);
        x0=EffectiveScalarFunction::coordinate(2,0);
        t=EffectiveScalarFunction::coordinate(2,1);

        EffectiveScalarFunction flowf=
            -1.9364*(x0^6)*(t^12)+0.30657*(x0^5)*(t^12)-18.4777*(x0^8)*(t^11)+10.4865*(x0^7)*(t^11)-3.65151*(x0^6)*(t^11)
            +0.74626*(x0^5)*(t^11)-0.0792375*(x0^4)*(t^11)+0.00334345*(x0^3)*(t^11)
            -15.9276*(x0^8)*(t^10)+12.2428*(x0^7)*(t^10)-5.69044*(x0^6)*(t^10)+1.54717*(x0^5)*(t^10)-0.223181*(x0^4)*(t^10)+0.0136828*(x0^3)*(t^10)
            -0.000190697*(x0^2)*(t^10)+4.93651*(x0^9)*(t^9)-10.28*(x0^8)*(t^9)+11.5009*((x0^7))*(t^9)-7.4244*(x0^6)*(t^9)+2.75346*(x0^5)*(t^9)
            -0.546434*(x0^4)*(t^9)+0.0484458*(x0^3)*(t^9)-0.00117394*(x0^2)*(t^9)+0.996825*(x0^9)*(t^8)-4.48571*(x0^8)*(t^8)+8.22222*(x0^7)*(t^8)
            -7.84444*(x0^6)*(t^8)+4.11667*(x0^5)*(t^8)-1.14722*(x0^4)*(t^8)+0.147619*(x0^3)*(t^8)-0.00595238*(x0^2)*(t^8)
            -(x0^8)*(t^7)+4*(x0^7)*(t^7)-6.33333*(x0^6)*(t^7)+5*(x0^5)*(t^7)-2.025*(x0^4)*(t^7)+0.383333*(x0^3)*(t^7)-0.0251984*(x0^2)*(t^7)
            +0.000198413*x0*(t^7)+(x0^7)*(t^6)-3.5*(x0^6)*(t^6)+4.66667*(x0^5)*(t^6)-2.91667*(x0^4)*(t^6)+0.836111*(x0^3)*(t^6)-0.0875*(x0^2)*(t^6)
            +0.00138889*x0*(t^6)-(x0^6)*(t^5)+3*(x0^5)*(t^5)-3.25*(x0^4)*(t^5)+1.5*(x0^3)*(t^5)-0.258333*(x0^2)*(t^5)+0.00833333*x0*(t^5)+(x0^5)*(t^4)
            -2.5*(x0^4)*(t^4)+2.08333*(x0^3)*(t^4)-0.625*(x0^2)*(t^4)+0.0416667*x0*(t^4)-(x0^4)*(t^3)
            +2*(x0^3)*(t^3)-1.16667*(x0^2)*(t^3)+0.166667*x0*(t^3)+(x0^3)*(t^2)-1.5*(x0^2)*(t^2)+0.5*x0*(t^2)-(x0^2)*t+x0*t+x0;

        //EffectiveVectorFunction f=(x*(o-x),o);
        EffectiveVectorFunction f={x*(o-x)};
        //ExactIntervalVectorType d=(ExactIntervalType(0.25,0.5),ExactIntervalType(-0.25,0.25));
        ExactBoxType d={ExactIntervalType(0.25,0.5)};
        FloatDP h=0.5;
        //ExactIntervalVectorType d(1u,ExactIntervalType(-0.125,+0.125));
        //FloatDP h=0.125;
        ValidatedVectorFunctionModelDP flow=integrator_ptr->flow_step(f,d,h);
        ValidatedVectorTaylorFunctionModelDP taylor_flow=dynamic_cast<ValidatedVectorTaylorFunctionModelDP&>(flow.reference());
        //EffectiveVectorFunction expected_flow( (x0+x0*(1-x0)*t+x0*(1-x0)*(1-2*x0)/2*t*t, y0+t) );
        //EffectiveVectorFunction expected_flow(1u, (x0+x0*(1-x0)*t+x0*(1-x0)*(1-2*x0)/2*t*t) );
        EffectiveVectorFunction expected_flow(1u, flowf );
        ARIADNE_TEST_PRINT(*integrator_ptr);
        ARIADNE_TEST_PRINT(f);
        ARIADNE_TEST_PRINT(flow);
        ARIADNE_TEST_PRINT(expected_flow);
        ARIADNE_TEST_PRINT(taylor_flow.errors());
        ARIADNE_TEST_PRINT(taylor_flow-expected_flow);
        ARIADNE_TEST_BINARY_PREDICATE(operator<,taylor_flow.error(),0.01);
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(taylor_flow-expected_flow),0.01+0.004);
    };
};

Int main(Int argc, const char* argv[]) {
    auto verbosity = get_verbosity(argc,argv);

    ARIADNE_TEST_PRINT("Testing TaylorSeriesIntegrator");
    TaylorSeriesIntegrator taylor_series_integrator(
            maximum_error=1e-6,sweep_threshold=1e-10,lipschitz_constant=0.5,
            step_maximum_error=1e-8,step_sweep_threshold=1e-12,
            minimum_spacial_order=1,minimum_temporal_order=4,
            maximum_spacial_order=3,maximum_temporal_order=8);
    taylor_series_integrator.verbosity=verbosity;
    ARIADNE_TEST_PRINT(taylor_series_integrator);
    TestIntegrator(taylor_series_integrator).test();

    ARIADNE_TEST_PRINT("Testing TaylorPicardIntegrator");
    TaylorPicardIntegrator taylor_picard_integrator(
            maximum_error=1e-6,sweep_threshold=1e-10,lipschitz_constant=0.5,
            step_maximum_error=1e-8,step_sweep_threshold=1e-12, maximum_temporal_order=16);
    ARIADNE_TEST_PRINT(taylor_picard_integrator);
    taylor_picard_integrator.verbosity=verbosity;
    TestIntegrator(taylor_picard_integrator).test();

    ARIADNE_TEST_PRINT("Testing AffineIntegrator");
    AffineIntegrator affine_integrator(1e-6, 6);
    affine_integrator.verbosity=verbosity;
    //TestIntegrator(affine_integrator).test_affine();
    ARIADNE_TEST_WARN("AffineIntegrator does not work correctly.");

    std::cerr<<"INCOMPLETE "<<std::flush;
    return ARIADNE_TEST_FAILURES;
}
