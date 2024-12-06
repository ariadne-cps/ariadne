/***************************************************************************
 *            test_expression_patch.cpp
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
#include "function/function_patch.hpp"
#include "function/domain.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/expression_patch.hpp"
#include "symbolic/expression_set.hpp"
#include "symbolic/space.hpp"

#include "../test.hpp"

using namespace std;
using namespace Ariadne;


class TestExpressionPatch {
    RealVectorVariable x;
    RealVariable t;
    RealVariable u,u0,u1;
    BoxDomainType xdom;
    IntervalDomainType tdom;
    IntervalDomainType udom;
  public:
    TestExpressionPatch();
    Void test();
  private:
    Void test_space_patch();
    Void test_syntax();
    Void test_expression_function_patch();
    Void test_construct_expression_patch();
    Void test_evaluate_function_patch();
    Void test_construct_function_patch();
    Void test_operations();
};

TestExpressionPatch::TestExpressionPatch()
    : x("x",2), t("t"), u("u"), u0("u0"), u1("u1")
    , xdom{{0.0_x,0.5_x},{0.25_x,0.75_x}}, tdom{0,0.125_x}, udom{-1.0_x,+1.0_x}
{
}

Void TestExpressionPatch::test()
{
    ARIADNE_TEST_CALL(test_space_patch());
    ARIADNE_TEST_CALL(test_syntax());
    ARIADNE_TEST_CALL(test_construct_expression_patch());
    ARIADNE_TEST_CALL(test_evaluate_function_patch());
    ARIADNE_TEST_CALL(test_construct_function_patch());
    ARIADNE_TEST_CALL(test_operations());
}



Void TestExpressionPatch::test_space_patch()
{
    // Constructors of SpacePatch
    //   SpacePatch(InitializerList<SpacePatch<T>> lst);
    //   SpacePatch(VariableIntervalDomainType vivl);
    //   SpacePatch(VectorVariableBoxDomainType vvbx);
    //   SpacePatch(List<RealVariable> vars, BoxDomainType dom);
    //   SpacePatch(List<RealVariable> vars, Map<RealVariable,IntervalDomainType> vivls);

    RealVariable t("t");
    IntervalDomainType tdom={0,0.5_dy};
    RealVectorVariable x("x",2);
    BoxDomainType xdom={{0.25_x,0.75_x},{0.375_x,1.125_x}};

    ARIADNE_TEST_CONSTRUCT(RealSpacePatch,spcptch1,(t|tdom));
    ARIADNE_TEST_CONSTRUCT(RealSpacePatch,spcptch2,(x|xdom));
    ARIADNE_TEST_CONSTRUCT(RealSpacePatch,spcptch3,({x[0],x[1],t},product(xdom,tdom)));
    ARIADNE_TEST_CONSTRUCT(RealSpacePatch,spcptch4,({x[0],x[1],t},{{x[0],xdom[0]},{x[1],xdom[1]},{t,tdom}}));
    ARIADNE_TEST_CONSTRUCT(RealSpacePatch,spcptch5,({t|tdom,spcptch2}));

    List<RealVariable> vars={x[0],x[1],t};
    Map<RealVariable,IntervalDomainType> vivls={{x[0],xdom[0]},{x[1],xdom[1]},{t,tdom}};
    ARIADNE_TEST_CONSTRUCT(RealSpacePatch,spcptch6,(vars,vivls));
}


Void TestExpressionPatch::test_expression_function_patch()
{
}

Void TestExpressionPatch::test_construct_expression_patch()
{
    // Constructors of ExpressionPatch
    //   ExpressionPatch(Map<RealVariable,IntervalDomainType> dom, RealExpression expr);
    //   ExpressionPatch(RealSpacePatch spcptch, RealExpression expr);
    //   ExpressionPatch(Vector<RealVariable> vars, BoxDomainType dom, Function<P,SIG> f);
    //   ExpressionPatch(RealSpacePatch spcptch, Function<P,SIG> f);
    //   ExpressionPatch(Vector<RealVariable> vars, FunctionPatch<P,SIG> fp);
    //   ExpressionPatch(RealSpace spc, FunctionPatch<P,SIG> fp);

    Dyadic h=0.125_dy;
    VectorVariableBoxDomainType xd(x,xdom);
    xd=x|xdom;
    VariableIntervalDomainType xd0(x[0],xdom[0]);
    VariableIntervalDomainType xd1(x[1],xdom[1]);
    VariableIntervalDomainType td(t,tdom);
    td=t|tdom;
    VariableIntervalDomainType ud(u,udom);
    VariableIntervalDomainType ud0(u0,udom);
    VariableIntervalDomainType ud1(u1,udom);
    ARIADNE_TEST_PRINT(xd);
    ARIADNE_TEST_PRINT(td);

    auto xtud=ValidatedVectorExpressionPatch({x|xdom,t|tdom,u|udom},{x[0],x[1],t,u});
    ARIADNE_TEST_PRINT(xtud);
    auto xhud=ValidatedVectorExpressionPatch({xd[0],xd[1],ud},{x[0],x[1],h,u});
    ARIADNE_TEST_PRINT(xhud);
}

Void TestExpressionPatch::test_evaluate_function_patch() {
    Dyadic h=0.125_dy;
    BoxDomainType xrdom={{0.125_x,0.375_x},{0.375_x,0.625_x}};
    BoxDomainType xtudom=product(xdom,tdom,udom);

    VectorVariableBoxDomainType xd(x,xdom);
    VariableIntervalDomainType td(t,tdom);
    VariableIntervalDomainType ud0(u0,udom);
    VariableIntervalDomainType ud1(u1,udom);

    EffectiveVectorMultivariateFunction phif({x[0],x[1],t,u},{ x[0]+u*t, (x[1]-x[0]+u)*exp(-t)+u*t+x[0]-u});
    ValidatedVectorMultivariateFunctionPatch phifp = ValidatedVectorMultivariateRestrictedFunction(phif,xtudom);
    ARIADNE_TEST_PRINT(phifp);
    ValidatedVectorMultivariateFunctionPatch phi1fp = restriction(phifp,product(xrdom,tdom,udom));

    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedVectorExpressionPatch,phiep,phifp({x[0],x[1],t,u0}));
    ARIADNE_TEST_ASSIGN(phiep,(phifp({x[0],x[1],t,u0})));
    ValidatedVectorExpressionPatch phi1ep = phi1fp({x[0],x[1],t,u});
    ARIADNE_TEST_PRINT(phi1ep);
    ValidatedVectorExpressionPatch phieptu({phi1ep,t|tdom,u|udom});
    ARIADNE_TEST_PRINT(phieptu);
    ValidatedVectorExpressionPatch phi2ep = phifp(phieptu);
    ARIADNE_TEST_PRINT(phi2ep);

    phi2ep = phifp({phi1fp({x,t,u0}),t|tdom,u|udom});
    ARIADNE_TEST_PRINT(phi2ep);


    ValidatedVectorExpressionPatch psi1ep = phi1fp({x,h,u0});
    ValidatedVectorExpressionPatch psi2ep = phifp({phi1fp({x,h,u0}),h,ud1});
    ARIADNE_TEST_PRINT(psi1ep);
    ARIADNE_TEST_PRINT(psi2ep);
}

Void TestExpressionPatch::test_construct_function_patch() {
    Dyadic h=0.125_dy;
    BoxDomainType xrdom={{0.125_x,0.375_x},{0.375_x,0.625_x}};
    BoxDomainType xtudom=product(xdom,tdom,udom);

    VectorVariableBoxDomainType xd(x,xdom);
    VariableIntervalDomainType td(t,tdom);
    VariableIntervalDomainType ud(u,udom);
    VariableIntervalDomainType ud0(u0,udom);
    VariableIntervalDomainType ud1(u1,udom);

    EffectiveVectorMultivariateFunction phif({x[0],x[1],t,u},{ x[0]+u*t, (x[1]-x[0]+u)*exp(-t)+u*t+x[0]-u});
    ValidatedVectorMultivariateFunctionPatch phifp = ValidatedVectorMultivariateRestrictedFunction(phif,xtudom);
    ValidatedVectorMultivariateFunctionPatch phi1fp = restriction(phifp,product(xrdom,tdom,udom));
    auto phi2fp = ValidatedVectorMultivariateFunctionPatch({xd,td,ud0,ud1},phifp({phi1fp({x,t,u0}),td,ud1}));
    phi2fp = ValidatedVectorMultivariateFunctionPatch({x|xdom,t|tdom,u0|udom,u1|udom},phifp({phi1fp({x,t,u0}),td,ud1}));
    ARIADNE_TEST_PRINT(phi2fp);

    ValidatedVectorExpressionPatch psiep = phifp({x,h,u});
    ARIADNE_TEST_PRINT(psiep);
    psiep = phifp({x[0],x[1],h,u});
    ARIADNE_TEST_PRINT(psiep);


    auto psifp = ValidatedVectorMultivariateFunctionPatch({xd,ud},psiep);
    psifp = ValidatedVectorMultivariateFunctionPatch({x|xdom,u|udom},psiep);
    ARIADNE_TEST_PRINT(psifp);

    auto psi1fp=restriction(psifp,product(xrdom,udom));
    auto psi2fp = ValidatedVectorMultivariateFunctionPatch({x|xdom,u0|udom,u1|udom},psifp({psi1fp({x,u0}),ud1}));
    //psi2fp = ValidatedVectorMultivariateFunctionPatch({x|xdom,u0|udom,u1|udom},psifp({psifp({xd,ud0}),ud1}));
    ARIADNE_TEST_PRINT(psi2fp);
}

Void TestExpressionPatch::test_syntax()
{
}

Void TestExpressionPatch::test_operations()
{
#warning  Add test constructing expression patch from constant vector expression.
    {
        RealVector v({5,23});
        std::cerr<<"v="<<v<<"\n";
        Expression<RealVector> ve(v);
        std::cerr<<"ve="<<ve<<"\n";
        Vector<RealExpression> ev(ve);
        Map<RealVariable,IntervalDomainType> doms;
        std::cerr<<"ev="<<ev<<"\n";
        ValidatedVectorExpressionPatch ep(doms,ev);
        std::cerr<<"ep="<<ep<<"\n";

    }

    Int n=5;
#warning RealConstant and Real are ambiguous (convert to RealExpression or ValidatedNumber).
    //RealConstant sc(5);
    //Real sc(5);
    ValidatedNumber sc(5);
    //RealVectorConstant vc({5,23});
    //RealVector vc({5,23});
    Vector<ValidatedNumber> vc({5,23});
    RealVariable x("x");
    RealVariable t("t");

    IntervalDomainType xdom={0.25_x,0.75_x};
    IntervalDomainType tdom={0,0.125_x};

    ValidatedScalarExpressionPatch sep({x|xdom},x);
    ValidatedVectorExpressionPatch vep({x|xdom,t|tdom},{x,t});

#warning Intermediate conversion from Vector<RealExpression> to FunctionPatch doesn't work
    add(vep,vc);

    add(sep,sep); sub(sep,sep); mul(sep,sep); div(sep,sep);
    add(sep,sc); sub(sep,sc); mul(sep,sc); div(sep,sc);
    add(sc,sep); sub(sc,sep); mul(sc,sep); div(sc,sep);
    //max(sep,sep); min(sep,sep); abs(sep);
    nul(sep); pos(sep); neg(sep); sqr(sep); hlf(sep); rec(sep); pow(sep,n);
    sqrt(sep); exp(sep); log(sep);
    sin(sep); cos(sep); tan(sep);

    sep+sep; sep-sep; sep*sep; sep/sep;
    sep+sc; sep-sc; sep*sc; sep/sc;
    sc+sep; sc-sep; sc*sep; sc/sep;

    nul(vep); pos(vep); neg(vep);
    add(vep,vep); sub(vep,vep); mul(sep,vep); mul(vep,sep); div(vep,sep);
    add(vep,vc); sub(vep,vc); mul(sep,vc); mul(vep,sc); div(vep,sc);
    add(vc,vep); sub(vc,vep); mul(sc,vep); mul(vc,sep); div(vc,sep);

    +vep; -vep;
    vep+vep; vep-vep; sep*vep; vep*sep; vep/sep;
    vep+vc; vep-vc; sep*vc; vep*sc; vep/sc;
    vc+vep; vc-vep; sc*vep; vc*sep; vc/sep;

#warning
//    join(sep,sep); join(sep,vep); join(vep,sep);
    join(vep,vep);
}

Int main() {
    TestExpressionPatch().test();

    return ARIADNE_TEST_FAILURES;
}
