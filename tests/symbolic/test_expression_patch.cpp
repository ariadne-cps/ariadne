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
#include "symbolic/expression_patch.hpp"

#include "algebra/algebra.hpp"
#include "function/function.decl.hpp"
#include "function/restricted_function.hpp"
#include "function/function_patch.hpp"
#include "function/taylor_function.hpp"
#include "function/domain.hpp"
#include "symbolic/expression.hpp"
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
    Void test_usage();
};

TestExpressionPatch::TestExpressionPatch()
    : x("x",2), t("t"), u("u"), u0("u0"), u1("u1")
    , xdom{{0.0_x,0.5_x},{0.25_x,0.75_x}}, tdom{0,0.125_x}, udom{-1.0_x,+1.0_x}
{
}

Void TestExpressionPatch::test()
{
#warning
    ARIADNE_TEST_CALL(test_usage()); return;
    ARIADNE_TEST_CALL(test_space_patch());
    ARIADNE_TEST_CALL(test_syntax());
    ARIADNE_TEST_CALL(test_construct_expression_patch());
    ARIADNE_TEST_CALL(test_evaluate_function_patch());
    ARIADNE_TEST_CALL(test_construct_function_patch());
    ARIADNE_TEST_CALL(test_operations());
    ARIADNE_TEST_CALL(test_usage());
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

    RealVariableDomainMap doms = {xd0,xd1,td,ud};

    auto x0x1tud=ValidatedVectorRestrictedExpression({xd0,xd1,td,ud},{x[0],x[1],t,u});

    auto xtud=ValidatedVectorRestrictedExpression({x|xdom,t|tdom,u|udom},{x[0],x[1],t,u});
    ARIADNE_TEST_PRINT(xtud);
    auto xhud=ValidatedVectorRestrictedExpression({xd[0],xd[1],ud},{x[0],x[1],h,u});
    ARIADNE_TEST_PRINT(xhud);
}

Void TestExpressionPatch::test_evaluate_function_patch() {
    Dyadic h=0.125_dy;
    RealConstant hc("h",h);
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
    ARIADNE_TEST_PRINT(phi1fp);


//    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedVectorRestrictedExpression,phire,phifp({x,t,u}));
    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedVectorRestrictedExpression,phire,phifp({x|xdom,t|tdom,u|udom}));
//    ARIADNE_TEST_ASSIGN(phire,phifp({x[0],x[1],t,u}));
//    ARIADNE_TEST_ASSIGN(phire,(phifp({x[0],x[1],t,u})));
//    ValidatedVectorRestrictedExpression phiretu = {phire,t,u};
//    ARIADNE_TEST_PRINT(phiretu);
    ValidatedVectorRestrictedExpression phiretur = {phire,t|tdom,u|udom};
    ARIADNE_TEST_PRINT(phiretur);

//    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedVectorRestrictedExpression,psi1ep,phifp({x,h,u0}));
//    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedVectorRestrictedExpression,psi2ep,phifp({psi1ep,hc,u1}));
//    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedVectorRestrictedExpression,psi1ep,phifp({x,h,u0}));
//    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedVectorRestrictedExpression,psi2ep,phifp({psi1ep,hc,u1}));
}

Void TestExpressionPatch::test_construct_function_patch() {
/*
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

    ValidatedVectorRestrictedExpression psiep = phifp({x,h,u});
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
*/
}


Void TestExpressionPatch::test_usage()
{
    RealVectorVariable x("x",2);
    RealVariable t("t");
    RealVariable u("u");

    BoxDomainType xd={{1,3},{2,5}};
    IntervalDomainType td={0,2};
    IntervalDomainType ud={-1,+1};
    Real h=2;
    RealConstant hc("h",2);

    ValidatedVectorMultivariateFunction phif = Function({x,t,u},{x[0]*exp(t)+u, (x[1]-1)*exp(t)});

    ValidatedVectorMultivariateFunctionPatch phifp = FunctionPatch({x|xd,t|td,u|ud},Expression<RealVector>({x[0]*exp(t)+u, (x[1]-1)*exp(t)}));

#warning
//    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedVectorRestrictedExpression, phiep, phifp({x,t,u}));
    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedVectorRestrictedExpression, phiep, phifp({x[0],x[1],t,u}));
//    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedVectorRestrictedExpression, psiep, phifp({x,h,u}));
    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedVectorRestrictedExpression, psiep, phifp({x[0],x[1],h,u}));
    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedVectorMultivariateFunctionPatch, psifp, FunctionPatch({x|xd,u|ud},psiep));

    RealVectorVariable x0("x0",2), x1("x1",2), x2("x2",2);
    RealVariable u0("u0"), u1("u1");



//    ValidatedVectorRestrictedExpression phi0ep = phifp({x0|xd,hc,u0|ud});
//    ValidatedVectorRestrictedExpression phi01ep = phifp({phifp({x0|xd,h,u0|ud}),h,u1|ud});
#warning Removing sub
//    ValidatedVectorMultivariateFunctionPatch g0fp = FunctionPatch({x0|xd,x1|xd,x2|xd,u0|ud,u1|ud}, phifp({x0,hc,u0})-x1);
//    ValidatedVectorMultivariateFunctionPatch g1fp = FunctionPatch({x0|xd,x1|xd,x2|xd,u0|ud,u1|ud}, phifp({x1,hc,u1})-x2);


    {
        BoxDomainType xd({{-1,+1},{-1,+1}});
        IntervalDomainType ud(-0.125_x,+0.125_x);
        IntervalDomainType td(0.0_x,+0.0625_x);

        ValidatedVectorMultivariateFunctionPatch phi0=ValidatedVectorMultivariateTaylorFunctionModelDP(
            BoxDomainType({{-1,+1},{-1,+1},{-0.125_dy,+0.125_dy},{0.0_dy,0.0625_dy}}),
            Function({x[0],x[1],u,t},{-x[1]*t+x[0],-u*t+x[1]}),
            ThresholdSweeper<FloatDP>(dp,1e-8) );

        ValidatedScalarRestrictedExpression zero({},u-u);
        ARIADNE_TEST_PRINT(zero);

        // Then phi0 is the function phi with domain x in x0Dom, u in uDom, t in t0Dom=[0,h]
        ARIADNE_TEST_PRINT(phi0({x[0],x[1],u,Real(h)}));
        ValidatedVectorMultivariateFunctionPatch psi0 = ValidatedVectorMultivariateFunctionPatch({x[0]|xd[0],x[1]|xd[1],u|ud}, phi0({x[0],x[1],u,Real(h)}));
        ARIADNE_TEST_PRINT(psi0);

        //( result_size=2, dom=[{-1.0000000000000000:1.0000000000000000},{-1.0000000000000000:1.0000000000000000},{-0.12500000000000000:0.12500000000000000},{0.:0.062500000000000000}], rng=[{-1.0627442:1.0627442},{-1.0078125:1.0078125}], f=[ { ~1.000~*x1*x3 +~1.000~*x0+/-0.000245}, { ~1.000~*x2*x3 +~1.000~*x1} ] )
    }

}

Void TestExpressionPatch::test_syntax()
{
}

Void TestExpressionPatch::test_operations()
{
    // Test constructing expression patch from constant vector expression.
    {
        RealVector v({5,23});
        std::cerr<<"v="<<v<<"\n";
        Expression<RealVector> ve(v);
        std::cerr<<"ve="<<ve<<"\n";
        Vector<RealExpression> ev(ve);
        std::cerr<<"ev="<<ev<<"\n";

        Variable<RealVector> x("x",1);
        ValidatedVectorMultivariateFunction f=Function(x,ev);
        std::cerr<<"f="<<f<<"\n";

        x=Variable<RealVector>("x",0);
        ARIADNE_TEST_FAIL(f=Function(x,ev));

        Map<RealVariable,IntervalDomainType> doms;
        x=Variable<RealVector>("x",1);
        doms[x[0]]=IntervalDomainType(-1,0);
        auto re=ValidatedVectorRestrictedExpression(doms,ev);
        std::cerr<<"re="<<re<<"\n";
        doms.clear();
        ARIADNE_TEST_FAIL(re=ValidatedVectorRestrictedExpression(doms,ev));
    }

    {
        RealVector v({5,23});
        std::cerr<<"v="<<v<<"\n";
        Variable<RealVector> x("x",2);
        Expression<RealVector> ve(v+x);
        std::cerr<<"ve="<<ve<<"\n";
        Vector<RealExpression> ev(ve);
        std::cerr<<"ev="<<ev<<"\n";

        ValidatedVectorMultivariateFunction f=Function(x,ev);
        std::cerr<<"f="<<f<<"\n";

        Map<RealVariable,IntervalDomainType> doms;
        doms[x[0]]=IntervalDomainType(-1,0);
        doms[x[1]]=IntervalDomainType(0,+1);
        auto re=ValidatedVectorRestrictedExpression(doms,ev);
        std::cerr<<"re="<<re<<"\n";
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

    ValidatedScalarRestrictedExpression sep({x|xdom},x);
    ValidatedVectorRestrictedExpression vep({x|xdom,t|tdom},{x,t});

#warning Intermediate conversion from Vector<RealExpression> to FunctionPatch does not work
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
