/***************************************************************************
 *            test_multifunction.cpp
 *
 *  Copyright  2020  Pieter Collins
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

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <fenv.h>

#include "config.hpp"
#include "algebra/sweeper.hpp"
#include "function/multifunction.hpp"
#include "function/taylor_function.hpp"
#warning
//    #include "function/scaled_function_patch.tpl.hpp"
#include "function/taylor_multifunction.hpp"
#include "numeric/numeric.hpp"

#include "symbolic/variable.hpp"
#include "symbolic/space.hpp"
#include "symbolic/expression.hpp"

#include "../test.hpp"

using namespace std;
using namespace Ariadne;

class TestMultifunction
{
    using ScalarIntervalFunctionModelType = ValidatedScalarMultivariateIntervalTaylorFunctionModel<FloatMP>;
    using VectorIntervalFunctionModelType = ValidatedVectorMultivariateIntervalTaylorFunctionModel<FloatMP>;
    using UpperInterval = Interval<FloatMPUpperBound>;
  public:
    Void test();
  private:
    Void test_evaluate();
    Void test_image();
    Void test_taylor_evaluate();
    Void test_interval_taylor_concept();
    Void test_interval_taylor_evaluate();
    Void test_function_set();
    Void test_inclusion_solutions();
};

Void TestMultifunction::test()
{
    Dyadic::set_default_writer(DecimalWriter());
    ARIADNE_TEST_CALL(test_evaluate());
    ARIADNE_TEST_CALL(test_image());
    ARIADNE_TEST_CALL(test_taylor_evaluate());
    ARIADNE_TEST_CALL(test_interval_taylor_evaluate());
    ARIADNE_TEST_CALL(test_function_set());
    ARIADNE_TEST_CALL(test_inclusion_solutions());
}

Void TestMultifunction::test_evaluate()
{
    RealVariable x0("x0"), x1("x1"), p("p");

    BoxDomainType dom({{-1,3},{-1,3}});
    BoxDomainType prms({{-1,1}});;

    EffectiveVectorMultivariateFunction f=make_function({x0,x1,p},{1.5_dy-x0*x0-0.375_dy*x1+p,x0});

    Sweeper<FloatDP> swp;

    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateParametrisedPatchMultifunction,mvfp,(f,prms));
    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateMultifunction,mvf,(mvfp));

    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateParametrisedMultifunctionModel<DP>,mvfmdp,(dom,f,prms,swp));
    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateParametrisedMultifunctionModel<DP>,mvfdp,(mvfmdp));

    Vector<ValidatedNumber> v(2,1);
    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedLocatedSet<RealVector>,res,mvfp(v));
    ARIADNE_TEST_ASSIGN(res,mvf(v));
    ARIADNE_TEST_ASSIGN(res,mvfmdp(v));
    ARIADNE_TEST_ASSIGN(res,mvfdp(v));
}

Void TestMultifunction::test_image() {

}

Void TestMultifunction::test_taylor_evaluate() {
    MP pr(128);
    ThresholdSweeper<FloatMP> swp(pr,1e-10);
    RealVariable x0("x0"), x1("x1"), p("p");
    Decimal c01=4_dec, c=9_dec;
    BoxDomainType dom({{-1,+3},{-1,+1}});
    BoxDomainType pdom({{-1,+1}});
    ValidatedVectorMultivariateFunction f=make_function({x0,x1,p},{c01*x0*x1+x0*p+c,x0+x1});
    Vector<ValidatedNumber> v({0,0.5_dec});
    //ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateFunctionPatch,fp,(product(dom,pdom),f));
    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateTaylorFunctionModel<FloatMP>,fm,(product(dom,pdom),f,swp));
    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateFunctionPatch,fp,(fm));
    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateParametrisedPatchMultifunction,mfp,(fp,pdom));
    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateParametrisedMultifunctionModel<MP>,mfm,(2u,fm));
    ARIADNE_TEST_ASSIGN_CONSTRUCT(auto,mfpv,mfp(v));
    ARIADNE_TEST_ASSIGN_CONSTRUCT(auto,mfmv,mfm(v));

    ValidatedConstrainedImageSet vcis(dom);
#warning
//    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedConstrainedImageSet,imvcis,image(vcis,mfp));
}

inline Interval<FloatMPUpperBound> create_zero(Interval<FloatMPUpperBound> const& ivl) {
    std::cerr<<"create_zero(ivl)\n";
    std::cerr<<"ivl="<<ivl<<", nul(ivl.lower_bound())="<<nul(ivl.lower_bound())<<"\n";
    return Interval<FloatMPUpperBound>(nul(ivl.lower_bound()),nul(ivl.upper_bound()));
}

Void TestMultifunction::test_interval_taylor_concept() {

    MP pr(128);
    BoxDomainType dom({{-1,+1},{-1,+1}});
    ThresholdSweeper<FloatMP> swp(pr,1e-10);
    UpperInterval c(pr);
    ScalarIntervalFunctionModelType f(dom,swp);
    Bounds<FloatMP> x(pr);
    ValidatedNumber y;

    +f; -f;
    f+f; f-f; f*f; f/f;
    f+c; f-c; f*c; f/c;
    c+f; c-f; c*f; c/f;
    f+x; f-x; f*x; f/x;
    x+f; x-f; x*f; x/f;
    f+y; f-y; f*y; f/y;
    y+f; y-f; y*f; y/f;

    nul(f); pos(f); neg(f); sqr(f); rec(f); pow(f,0u); pow(f,0);
    sqrt(f); exp(f); log(f); sin(f); cos(f); tan(f); atan(f);
}

Void TestMultifunction::test_interval_taylor_evaluate()
{
    MP pr(128);
    BoxDomainType dom({{0,1},{1,3}});
    ThresholdSweeper<FloatMP> swp(pr,1e-10);

    ScalarIntervalFunctionModelType x0=ScalarIntervalFunctionModelType::coordinate(dom,0,swp);
    ScalarIntervalFunctionModelType x1=ScalarIntervalFunctionModelType::coordinate(dom,1,swp);
    UpperInterval one(1,1,pr);
    UpperInterval c01(3,5,pr);
    UpperInterval c(7,11,pr);
    ScalarIntervalFunctionModelType xc=ScalarIntervalFunctionModelType::constant(dom,c,swp);

    Vector<ValidatedNumber> u({0.75_x,2.5_x});
    //Vector<ValidatedNumber> u({0.75_x,-0.5_x});
    Vector<FloatMPBounds> v(u,pr);
    Vector<FloatDPBounds> w(u,dp);

    ARIADNE_TEST_SAME(xc(v),c);
    ARIADNE_TEST_SAME(x0(v),one*v[0]);
    ARIADNE_TEST_SAME(x1(v),one*v[1]);
    ARIADNE_TEST_CONSTRUCT(ScalarIntervalFunctionModelType,sf,(c01*x0*x1+c));
    ARIADNE_TEST_SAME(sf(v),c01*v[0]*v[1]+c);
    //ARIADNE_TEST_SAME(sf(w),c01*w[0]*w[1]+c);
    ARIADNE_TEST_PRINT(sf(w));

    ARIADNE_TEST_CONSTRUCT(VectorIntervalFunctionModelType,vf,({c01*x0*x1+c,x0+x1}));
    ARIADNE_TEST_SAME(vf(v),Box<UpperInterval>({c01*v[0]*v[1]+c,one*v[0]+v[1]}));

    using IntervalFunctionSetType = ValidatedVectorMultivariateIntervalTaylorFunctionModelSet<FloatMP>;
    ARIADNE_TEST_CONSTRUCT(IntervalFunctionSetType,set,(dom,vf));

#warning
    ARIADNE_TEST_PRINT(vf.domain());
    ARIADNE_TEST_PRINT(vf.models());
    using VectorFunctionModelType = ValidatedVectorMultivariateTaylorFunctionModel<FloatMP>;
    ARIADNE_TEST_CONSTRUCT(VectorFunctionModelType,vfm,(explicitly_parametrise(vf,0.0625_x)));
    ARIADNE_TEST_PRINT(vfm.models());
    ARIADNE_TEST_CONSTRUCT(VectorFunctionModelType,vfm2,(explicitly_parametrise(vf,0.75_x)));
    ARIADNE_TEST_PRINT(vfm2.models());
}

Void TestMultifunction::test_function_set()
{
    using ScalarIntervalFunctionModelType = ValidatedScalarMultivariateIntervalTaylorFunctionModel<FloatMP>;
    using Ivl = Interval<FloatMPUpperBound>;

    MP pr(128);
    BoxDomainType dom({{-1,+1},{-1,+1}});
//    BoxDomainType dom({{0,1},{1,3}});
    ThresholdSweeper<FloatMP> swp(pr,1e-10);
    ScalarIntervalFunctionModelType x0=ScalarIntervalFunctionModelType::coordinate(dom,0,swp);
    ScalarIntervalFunctionModelType x1=ScalarIntervalFunctionModelType::coordinate(dom,1,swp);

    ScalarIntervalFunctionModelType sf=Ivl(3,5,pr)*x0*x1+Ivl(7,11,pr);

    ValidatedScalarMultivariateIntervalTaylorFunctionModelSet<FloatMP> fset=sf;
    FunctionSet<ValidatedTag,Real(RealVector),CompactSet> gset(fset);
    CompactSet<ValidatedTag,Real(RealVector)> gfset(fset);
}


Void TestMultifunction::test_inclusion_solutions()
{
    using ScalarIntervalFunctionModelType = ValidatedScalarMultivariateIntervalTaylorFunctionModel<FloatMP>;
    using VectorIntervalFunctionModelType = ValidatedVectorMultivariateIntervalTaylorFunctionModel<FloatMP>;
    using Ivl = Interval<FloatMPUpperBound>;
    using P=ValidatedTag;

    MP pr(128);
    BoxDomainType dom({{-1,+1},{-1,+1},{0,+1}});
//    BoxDomainType dom({{0,1},{1,3}});
    ThresholdSweeper<FloatMP> swp(pr,1e-10);
    ScalarIntervalFunctionModelType x0=ScalarIntervalFunctionModelType::coordinate(dom,0,swp);
    ScalarIntervalFunctionModelType x1=ScalarIntervalFunctionModelType::coordinate(dom,1,swp);
    ScalarIntervalFunctionModelType t=ScalarIntervalFunctionModelType::coordinate(dom,1,swp);

    VectorIntervalFunctionModelType vf={Ivl(3,5,pr)*x0*x1*+Ivl(7,11,pr),x1*cos(t)};

    Vector<ValidatedNumber> w0({.375_x,0.625_x});
    ValidatedNumber s(0.25_x);

    ValidatedIntervalTaylorMultiflow<FloatMP> phi(vf);
    ValidatedIntervalTaylorCurriedTrajectorySet<FloatMP> phi_w0=phi(w0);
    CompactSet<P,RealVector> phi_w0_s=phi_w0(s);

    //    Function<P,CompactSet<P,RealVector(Real)>(RealVector)> const& phif=phi;
    ValidatedCompactMultiflow::Interface const& phif=phi;
    CompactSet<P,RealVector(Real)> phif_w0=phif._call(w0);
    CompactSet<P,RealVector> phif_w0_s=phif_w0(s);

    ARIADNE_TEST_PRINT(phi_w0);
    ARIADNE_TEST_PRINT(phif_w0);
    ARIADNE_TEST_PRINT(phi_w0_s);
    ARIADNE_TEST_PRINT(phif_w0_s);

    FloatDPLowerBound(FloatMP(0,precision(128)),dp);

    BoxDomainType xdom({dom[0],{0.625_x,0.75_x}});
    BoxDomainType pdom({dom[2]});
    ARIADNE_TEST_PRINT(pdom);
    ARIADNE_TEST_PRINT(vf.argument_size());
    auto mf=ParametrisedMultifunction(vf,pdom);
    ARIADNE_TEST_PRINT(mf);
    ARIADNE_TEST_PRINT(mf.argument_size());
    auto mfimset=mf(w0);
    ARIADNE_TEST_PRINT(mfimset);
    auto imset=explicitly_parametrise(mfimset,0.0009765625_x);
    ARIADNE_TEST_PRINT(imset);
    ARIADNE_TEST_PRINT(imset.separated(xdom));

    mfimset=image(mf,xdom);
    ARIADNE_TEST_PRINT(mfimset);
    ARIADNE_TEST_PRINT(mfimset.separated(xdom));


}

Int main() {
    TestMultifunction().test();

    return ARIADNE_TEST_FAILURES;
}

