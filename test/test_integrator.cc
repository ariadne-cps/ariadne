/***************************************************************************
 *            test_integrator.cc
 *
 *  Copyright  2008-10  Pieter Collins
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
#include <sstream>
#include <string>

#include "polynomial.h"
#include "integrator.h"
#include "function.h"
#include "taylor_function.h"
#include "vector.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

inline Vector<Interval> operator,(const Interval& ivl1, const Interval& ivl2) {
    Vector<Interval> r(2); r[0]=ivl1; r[1]=ivl2; return r; }

class TestIntegrator
{
    typedef Vector<Interval> IntervalVector;
  private:
    scoped_ptr<IntegratorInterface> integrator;
    RealScalarFunction o,x,y,x0,y0,t;
  public:
    TestIntegrator(const IntegratorInterface& i)
        : integrator(i.clone())
    {
        o=RealScalarFunction::constant(2,1.0);
        x=RealScalarFunction::coordinate(2,0);
        y=RealScalarFunction::coordinate(2,1);
        x0=RealScalarFunction::coordinate(3,0);
        y0=RealScalarFunction::coordinate(3,1);
        t=RealScalarFunction::coordinate(3,2);
    }

    int test() {
        ARIADNE_TEST_CALL(test_constant_derivative());
        ARIADNE_TEST_CALL(test_quadratic_flow());
        ARIADNE_TEST_CALL(test_linear());
        ARIADNE_TEST_CALL(test_spiral());
        ARIADNE_TEST_CALL(test_logistic());
        return 0;
    }

    void test_constant_derivative() {
        RealVectorFunction f=(o*2,o*3);
        IntervalVector d=(Interval(0.0,1.0),Interval(-0.5,1.5));
        Float h=0.25;
        VectorTaylorFunction flow=integrator->flow_step(f,d,h);
        RealVectorFunction expected_flow( (x0+2*t,y0+3*t) );
        ARIADNE_TEST_PRINT(f);
        ARIADNE_TEST_PRINT(flow);
        ARIADNE_TEST_PRINT(expected_flow);
        ARIADNE_TEST_PRINT(flow.errors());
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow-expected_flow),1e-8);
    }

    void test_quadratic_flow() {
        RealVectorFunction f=(o,x);
        IntervalVector d=(Interval(0.0,1.0),Interval(-0.5,1.5));
        Float h=0.25;
        VectorTaylorFunction flow=integrator->flow_step(f,d,h);
        RealVectorFunction expected_flow( (x0+t,y0+x0*t+t*t/2) );
        ARIADNE_TEST_PRINT(f);
        ARIADNE_TEST_PRINT(flow);
        ARIADNE_TEST_PRINT(expected_flow);
        ARIADNE_TEST_PRINT(flow.errors());
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow-expected_flow),1e-8);
    }

    void test_linear() {
        RealVectorFunction f=(x,-y);
        IntervalVector d=(Interval(-0.25,0.25),Interval(-0.25,0.25));
        Float h=0.25;
        VectorTaylorFunction flow=integrator->flow_step(f,d,h);
        RealVectorFunction expected_flow( (x0*(1+t+t*t/2+t*t*t/6+t*t*t*t/24),y0*(1-t+t*t/2-t*t*t/6+t*t*t*t/24)) );
        ARIADNE_TEST_PRINT(f);
        ARIADNE_TEST_PRINT(flow);
        ARIADNE_TEST_PRINT(expected_flow);
        ARIADNE_TEST_PRINT(flow.errors());
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow-expected_flow),1e-4);
    };

    void test_spiral() {
        RealVectorFunction f=(-0.5*x-y,x-0.5*y);
        IntervalVector d=(Interval(0.75,1.25),Interval(-0.25,0.25));
        Float h=0.25;
        VectorTaylorFunction flow=integrator->flow_step(f,d,h);
        RealVectorFunction expected_flow( (exp(-0.5*t)*(x0*cos(t)-y0*sin(t)),exp(-0.5*t)*(x0*sin(t)+y0*cos(t))) );
        ARIADNE_TEST_PRINT(f);
        ARIADNE_TEST_PRINT(flow);
        ARIADNE_TEST_PRINT(expected_flow);
        ARIADNE_TEST_PRINT(flow.errors());
        ARIADNE_TEST_PRINT((flow-expected_flow).sweep(GradedSweeper(3)));
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow-expected_flow),1e-3);

    };

    void test_logistic() {
        RealVectorFunction f=(x*(o-x),o);
        IntervalVector d=(Interval(0.25,0.5),Interval(-0.25,0.25));
        Float h=0.5;
        VectorTaylorFunction flow=integrator->flow_step(f,d,h);
        RealVectorFunction expected_flow( (x0+x0*(1-x0)*t+x0*(1-x0)*(1-2*x0)/2*t*t, y0+t) );
        ARIADNE_TEST_PRINT(f);
        ARIADNE_TEST_PRINT(flow);
        ARIADNE_TEST_PRINT(expected_flow);
        ARIADNE_TEST_PRINT(flow.errors());
        ARIADNE_TEST_PRINT(flow-expected_flow);
        ARIADNE_TEST_PRINT((flow-expected_flow).sweep(GradedSweeper(3)));
        ARIADNE_TEST_BINARY_PREDICATE(operator<,flow.error(),0.01);
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow-expected_flow),0.01+0.004);
    };
};

int main(int argc, const char **argv) {
    int verbosity=get_verbosity(argc,argv);

    ARIADNE_TEST_PRINT("Testing TaylorSeriesIntegrator");
    TaylorSeriesIntegrator taylor_series_integrator(2,6,1e-8,1e-16);
    taylor_series_integrator.verbosity=verbosity;
    TestIntegrator(taylor_series_integrator).test();

    ARIADNE_TEST_PRINT("Testing TaylorIntegrator");
    TaylorIntegrator taylor_integrator(16,1e-6,1e-16);
    taylor_integrator.verbosity=verbosity;
    TestIntegrator(taylor_integrator).test();

    ARIADNE_TEST_PRINT("Testing AffineIntegrator");
    AffineIntegrator affine_integrator(1e-6, 6);
    affine_integrator.verbosity=verbosity;
    TestIntegrator(affine_integrator).test();

    std::cerr<<"INCOMPLETE "<<std::flush;
    return ARIADNE_TEST_FAILURES;
}
