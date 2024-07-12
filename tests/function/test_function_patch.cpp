/***************************************************************************
 *            test_function_patch.cpp
 *
 *  Copyright  2024  Pieter Collins
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

#include "config.hpp"
#include "algebra/algebra.hpp"
#include "function/taylor_model.hpp"
#include "function/taylor_function.hpp"
#include "function/function_patch.hpp"
#include "function/domain.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/expression_patch.hpp"
#include "symbolic/expression_set.hpp"
#include "symbolic/space.hpp"

#include "../test.hpp"

using namespace std;
using namespace Ariadne;


class TestFunctionPatch
{
  public:
    Void test();
  private:
    Void test_expression_function();
};

Void TestFunctionPatch::test()
{
    ARIADNE_TEST_CALL(test_expression_function());
}


Void TestFunctionPatch::test_expression_function()
{
    RealVectorVariable x("x",2);
    RealVariable t("t");
    RealVariable u("u"), u0("u0"), u1("u1");

    // dx0/dt=u; dx1/dt=-x1+x0 = -x1+(x0i+u*t)
    // x0=x0i+u*t; x1=(x1i-x0i+u)*exp(-t)+ u*t +x0i-u
    Vector<RealExpression> phie={ x[0]+u*t, (x[1]-x[0]+u)*exp(-t)+u*t+x[0]-u};
    RealSpace spc={x[0],x[1],t,u};

    EffectiveVectorMultivariateFunction phif({x[0],x[1],t,u},{ x[0]+u*t, (x[1]-x[0]+u)*exp(-t)+u*t+x[0]-u});
    ARIADNE_TEST_PRINT(phif);

    Dyadic h(0.125_x);
    BoxDomainType xdom={{0.0_x,0.5_x},{0.25_x,0.75_x}};
    IntervalDomainType tdom={0,h};
    IntervalDomainType udom={-1.0_x,+1.0_x};

    BoxDomainType xtudom=product(xdom,tdom,udom);
    ThresholdSweeper<FloatDP> swp(dp,1e-8);

    EffectiveVectorMultivariateFunction psif({x[0],x[1],u}, call(phif,{x[0],x[1],h,u}));
    ARIADNE_TEST_PRINT(psif);
    EffectiveVectorMultivariateFunction psi2f({x[0],x[1],u0,u1}, call(phif,join( call(phif,{x[0],x[1],h,u0}),h,u1 ) ));
    ARIADNE_TEST_PRINT(psi2f);


    using ValidatedVectorExpressionPatch = ExpressionPatch<ValidatedTag,RealVector>;

    VectorVariableBoxDomainType xd(x,xdom);
    xd=x|xdom;
    VariableIntervalDomainType xd0(x[0],xdom[0]);
    VariableIntervalDomainType xd1(x[1],xdom[1]);
    VariableIntervalDomainType td(t,tdom);
    td=t|tdom;
    VariableIntervalDomainType ud(u,udom);
    VariableIntervalDomainType ud0(u0,udom);
    VariableIntervalDomainType ud1(u1,udom);

    ValidatedVectorMultivariateFunctionPatch phifp = ValidatedVectorMultivariateTaylorFunctionModelDP(xtudom,phif,swp);
    ARIADNE_TEST_PRINT(phifp);

    auto xtud=ValidatedVectorExpressionPatch({x|xdom,t|tdom,u|udom});
    ARIADNE_TEST_PRINT(xtud);
    auto xhud=ValidatedVectorExpressionPatch({xd[0],xd[1],h,ud});
    ARIADNE_TEST_PRINT(xhud);

    ValidatedVectorExpressionPatch phi2ep = phifp({phifp({x,t,u0}),t,u1});
    ARIADNE_TEST_PRINT(phi2ep);
    auto phi2fp = ValidatedVectorMultivariateFunctionPatch({xd,td,ud0,ud1},phifp({phifp({x,t,u0}),t,u1}));
    phi2fp = ValidatedVectorMultivariateFunctionPatch({x|xdom,t|tdom,u0|udom,u1|udom},phifp({phifp({x,t,u0}),t,u1}));
    ARIADNE_TEST_PRINT(phi2fp);

    ValidatedVectorExpressionPatch psiep = phifp({x,h,u});
    ARIADNE_TEST_PRINT(psiep);
    psiep = phifp({x[0],x[1],h,u});
    ARIADNE_TEST_PRINT(psiep);

    auto psifp = ValidatedVectorMultivariateFunctionPatch({xd,ud},psiep);
    psifp = ValidatedVectorMultivariateFunctionPatch({x|xdom,u|udom},psiep);
    ARIADNE_TEST_PRINT(psifp);

    auto psi2fp = ValidatedVectorMultivariateFunctionPatch({x|xdom,u0|udom,u1|udom},psifp({psifp({x,u0}),u1}));
    ARIADNE_TEST_PRINT(psi2fp);
}


Int main() {
    TestFunctionPatch().test();

    return ARIADNE_TEST_FAILURES;
}
