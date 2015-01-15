/***************************************************************************
 *            test_taylor_model.cc
 *
 *  Copyright 2008  Pieter Collins
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
#include <iomanip>
#include "config.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/multi_index.h"
#include "function/taylor_model.h"
#include "algebra/differential.h"
#include "function/function.h"
#include "function/polynomial.h"

#include "test.h"
using std::cout; using std::cerr; using std::endl;
using namespace Ariadne;

Vector<ExactFloat> v(Nat n, Nat i) { return Vector<ExactFloat>::unit(n,i); }
ValidatedTaylorModel ctm(Nat m, double c, Sweeper swp) { return ValidatedTaylorModel::constant(m,ExactFloat(c),swp); }
ValidatedTaylorModel ctm(Nat m, Sweeper swp) { return ValidatedTaylorModel::constant(m,ExactFloat(1.0),swp); }
//ValidatedTaylorModel tm(Nat m, Nat i, Sweeper swp) { return ValidatedTaylorModel::coordinate(m,i,swp); }


class TestTaylorModel
{
    typedef MultiIndex MI;
    typedef Expansion<Float> E;
    typedef Polynomial<Float> P;
    typedef ValidatedTaylorModel T;
  public:
    Sweeper swp;
  public:
    TestTaylorModel(Sweeper swp);
    Void test();
  private:
    Void test_concept();
    Void test_constructors();
    Void test_predicates();
    Void test_approximation();
    Void test_unscale();
    Void test_evaluate();
    Void test_arithmetic();
    Void test_range();
    Void test_functions();
    Void test_rescale();
    Void test_restrict();
    Void test_refinement();
    Void test_split();
    Void test_antiderivative();
    Void test_compose();
};


TestTaylorModel::TestTaylorModel(Sweeper sweeper)
    : swp(sweeper)
{
}

Void TestTaylorModel::test()
{
    std::cerr<<std::setprecision(17);
    std::cout<<std::setprecision(17);
    std::clog<<std::setprecision(17);

    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_predicates());
    ARIADNE_TEST_CALL(test_approximation());
    ARIADNE_TEST_CALL(test_unscale());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_range());
    ARIADNE_TEST_CALL(test_functions());
    ARIADNE_TEST_CALL(test_rescale());
    ARIADNE_TEST_CALL(test_restrict());
    ARIADNE_TEST_CALL(test_refinement());
    ARIADNE_TEST_CALL(test_split());
    ARIADNE_TEST_CALL(test_antiderivative());
    ARIADNE_TEST_CALL(test_compose());
}


Void TestTaylorModel::test_concept()
{
    const ExactFloat f=0;
    const ValidatedFloat i;
    const Vector<ExactFloat> vf;
    const Vector<ValidatedFloat> vi;
    const ValidatedTaylorModel  t(0,swp);
    ValidatedTaylorModel tr(0,swp);

    tr=t+f; tr=t-f; tr=t*f; tr=t/f;
    tr=f+t; tr=f-t; tr=f*t; tr=f/t;
    tr=t+i; tr=t-i; tr=t*i; tr=t/i;
    tr=i+t; tr=i-t; tr=i*t; tr=i/t;
    tr=t+t; tr=t-t; tr=t*t; tr=t/t;

    tr+=f; tr-=f; tr*=f; tr/=f;
    tr+=i; tr-=i; tr*=i; tr/=i;
    tr+=t; tr-=t;

    tr=exp(t); tr=log(t); tr=sqrt(t);
    tr=sin(t); tr=cos(t); tr=tan(t);
    //tr=asin(t); tr=acos(t); tr=atan(t);

    tr.sweep(); tr.clobber();

    evaluate(t,vi);
    t.domain(); t.range(); t.expansion(); t.error();

}

Void TestTaylorModel::test_constructors()
{
    ARIADNE_TEST_CONSTRUCT(ValidatedTaylorModel,tv1,(E(2,3, {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0}), 0.25, swp));
    ARIADNE_TEST_CONSTRUCT(ValidatedTaylorModel,tv2,(E({ {{0,0},1.0}, {{1,0},2.0}, {{0,1},3.0}, {{2,0},4.0}, {{1,1},5.0}, {{0,2},6.0}, {{3,0},7.0}, {{2,1},8.0}, {{1,2},9.0}, {{0,3},10.0} }), 0.25, swp));

    ARIADNE_ASSERT_EQUAL(tv1.value(),1.0);
    ARIADNE_ASSERT_EQUAL(tv1.error(),0.25);
    ARIADNE_ASSERT_EQUAL(tv1.norm(),55.25);

    ARIADNE_ASSERT_EQUAL(tv2,tv1);
}

Void TestTaylorModel::test_predicates()
{
    ValidatedTaylorModel tv1({{{0},1.00},{{1},2.00},{{2},3.00}}, 0.75, swp);
    ValidatedTaylorModel tv2({{{0},1.00},{{1},1.75},{{2},3.25}}, 0.25, swp);
    ValidatedTaylorModel tv3({{{0},1.125},{{1},1.75},{{2},3.25}}, 0.25, swp);
    ValidatedTaylorModel tv4({{{0},1.00},{{1},2.25},{{2},3.00},{{4},-0.25}}, 0.25, swp);

    ARIADNE_TEST_BINARY_PREDICATE(refines,tv1,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(refines,tv2,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(!refines,tv3,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(refines,tv4,tv1);
}

Void TestTaylorModel::test_approximation()
{
    ARIADNE_TEST_CONSTRUCT(ValidatedTaylorModel,tv2,(E(1,2, {1.0,2.0,3.0}),0.25, swp));
}


Void TestTaylorModel::test_unscale()
{

    if(unscale(ValidatedTaylorModel({{{0},3.0}},0.0,swp),ExactInterval(1.0)).codomain()!=ExactInterval(1.0)) {
        ARIADNE_TEST_WARN("Unscaling over singleton domain does not yield constant");
    }
}

Void TestTaylorModel::test_evaluate()
{
    Vector<ValidatedFloat> iv={{0.25,0.5},{-0.75,-0.5}};
    ValidatedTaylorModel tv({{{0,0},1.0},{{1,0},2.0},{{0,1},3.0},{{2,0},4.0},{{1,1},5.0},{{0,2},6.0}},0.25,swp);
    ARIADNE_TEST_EQUAL(evaluate(tv,iv),ValidatedFloat(-1,1));
}

Void TestTaylorModel::test_arithmetic()
{
    //Operations which can be performed exactly with floating-point arithmetic.
    ARIADNE_TEST_EQUAL(ValidatedTaylorModel(E(1,2, {1.0,-2.0,3.0}), 0.75,swp)+(-3), ValidatedTaylorModel(E(1,2, {-2.0,-2.0,3.0}), 0.75,swp));
    ARIADNE_TEST_EQUAL(ValidatedTaylorModel(E(1,2, {1.0,-2.0,3.0}), 0.75,swp)-(-3), ValidatedTaylorModel(E(1,2, {4.0,-2.0,3.0}), 0.75,swp));
    ARIADNE_TEST_EQUAL(ValidatedTaylorModel(E(1,2, {1.0,-2.0,3.0}), 0.75,swp)*(-3), ValidatedTaylorModel(E(1,2, {-3.0,6.0,-9.0}), 2.25,swp));
    ARIADNE_TEST_EQUAL(ValidatedTaylorModel(E(1,2, {1.0,-2.0,3.0}), 0.75,swp)/(-4), ValidatedTaylorModel(E(1,2, {-0.25,0.5,-0.75}), 0.1875,swp));
    ARIADNE_TEST_EQUAL(ValidatedTaylorModel(E(1,2, {1.0,-2.0,3.0}), 0.75,swp)+ValidatedFloat(-1,2), ValidatedTaylorModel(E(1,2, {1.5,-2.0,3.0}), 2.25,swp));
    ARIADNE_TEST_EQUAL(ValidatedTaylorModel(E(1,2, {1.0,-2.0,3.0}), 0.75,swp)-ValidatedFloat(-1,2), ValidatedTaylorModel(E(1,2, {0.5,-2.0,3.0}), 2.25,swp));
    ARIADNE_TEST_EQUAL(ValidatedTaylorModel(E(1,2, {1.0,-2.0,3.0}), 0.75,swp)*ValidatedFloat(-1,2), ValidatedTaylorModel(E(1,2, {0.5,-1.0,1.5}), 10.5,swp));
    ARIADNE_TEST_EQUAL(ValidatedTaylorModel(E(1,2, {1.0,-2.0,3.0}), 0.75,swp)/ValidatedFloat(0.25,2.0), ValidatedTaylorModel(E(1,2, {2.25,-4.5,6.75}), 13.5,swp));
    ARIADNE_TEST_EQUAL(+ValidatedTaylorModel(E(1,2, {1.0,-2.0,3.0}), 0.75,swp), ValidatedTaylorModel(E(1,2, {1.0,-2.0,3.0}), 0.75,swp));
    ARIADNE_TEST_EQUAL(-ValidatedTaylorModel(E(1,2, {1.0,-2.0,3.0}), 0.75,swp), ValidatedTaylorModel(E(1,2, {-1.0,2.0,-3.0}), 0.75,swp));
    ARIADNE_TEST_EQUAL(ValidatedTaylorModel(E(1,2, {1.0,-2.0,3.0}), 0.75,swp)+ValidatedTaylorModel(E(1,2, {3.0,2.0,-4.0}),0.5,swp), ValidatedTaylorModel(E(1,2, {4.0,0.0,-1.0}), 1.25,swp));
    ARIADNE_TEST_EQUAL(ValidatedTaylorModel(E(1,2, {1.0,-2.0,3.0}), 0.75,swp)-ValidatedTaylorModel(E(1,2, {3.0,2.0,-4.0}),0.5,swp), ValidatedTaylorModel(E(1,2, {-2.0,-4.0,7.0}), 1.25,swp));
    ARIADNE_TEST_EQUAL(ValidatedTaylorModel(E(1,2, {0.0,0.0,3.0}), 0.75,swp)*ValidatedTaylorModel(E(1,2, {3.0,2.0,-4.0}),0.5,swp), ValidatedTaylorModel(E(1,4, {0.0,0.0,9.0,6.0,-12.0}), 8.625,swp));
    ARIADNE_TEST_EQUAL(ValidatedTaylorModel(E(1,2, {1.0,-2.0,3.0}), 0.75,swp)*ValidatedTaylorModel(E(1,2, {3.0,2.0,-4.0}),0.5,swp), ValidatedTaylorModel(E(1,4, {3.0,-4.0,1.0,14.0,-12.0}), 10.125,swp));

    ValidatedTaylorModel tm_inf(Expansion<Float>(2),+inf,swp);
    if(isnan(numeric_cast<double>(numeric_cast<Float>((tm_inf * 0.0).error())))) {
        ARIADNE_TEST_WARN("Multiplying 0+/-inf by 0 yields 0+/-NaN");
    } else if((tm_inf * 0).error()==+infty) {
        ARIADNE_TEST_WARN("Multiplying 0+/-inf by 0 yields 0+/-inf");
    } else if((tm_inf * 0).error()==0) {
        ARIADNE_TEST_PRINT("Multiplying 0+/-inf by 0 yields 0+/-0");
    }
}

Void TestTaylorModel::test_range()
{
    ValidatedTaylorModel x0 = ValidatedTaylorModel::coordinate(2,0,swp);
    ValidatedTaylorModel x1 = ValidatedTaylorModel::coordinate(2,1,swp);

    // Test range of cubic, which should be exact
    ValidatedTaylorModel t1 = x0*x0*x0+x0;
    ARIADNE_TEST_ASSERT(t1.range()==ExactInterval(-2,+2));

    // Test range of quadratic, which could be exact, but need not be
    ValidatedTaylorModel t2 = x0*x0+x0;
    ARIADNE_TEST_BINARY_PREDICATE(refines,t2.range(),ExactInterval(-2,+2));
    ARIADNE_TEST_BINARY_PREDICATE(refines,ExactInterval(-0.25,+2),t2.range());
    if(make_exact_interval(t2.range())!=ExactInterval(-0.2578125,2.0)) {
        ARIADNE_TEST_WARN("ValidatedTaylorModel::range() not exact for quadratic functions."); }
}

Void TestTaylorModel::test_functions()
{
    ValidatedTaylorModel x(E(1,1, {0.0, 1.0}), 0.0, swp);
    ValidatedTaylorModel xz(E(1,1, {0.0, 0.5}), 0.0, swp);
    ValidatedTaylorModel xo(E(1,1, {1.0, 0.5}), 0.0, swp);

    ARIADNE_TEST_PRINT(exp(x));
    ARIADNE_TEST_PRINT(sin(x));
    ARIADNE_TEST_PRINT(cos(x));

    //Functions based on their natural defining points with variable dependence 0.5
    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(T(E(1,1,{0.0,1.0}),0.0,swp)),T(E(1,6, {1.0,1.00,0.500,0.1667,0.0417,0.0083,0.0014}),0.0004,swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,sin(x),T(E(1,6, {0.0,1.0000,0.0,-0.1667,0.0,0.0083,0.0}),0.0003,swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,cos(x),T(E(1,6, {1.0000,0.0,-0.5000,0.0,0.0417,0.0,-0.0014}),0.0003,swp));

    //Functions based on their natural defining points with variable dependence 0.5
    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(xz),ValidatedTaylorModel(E(1,6, {1.00000,0.50000,0.12500,0.02083,0.00260,0.00026,0.00002}), 0.00003, swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,sin(xz),ValidatedTaylorModel(E(1,6, {0.00000,0.50000,0.0000,-0.02083,0.00000,0.00026,0.00000}), 0.00003, swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,cos(xz),ValidatedTaylorModel(E(1,6, {1.00000,0.0000,-0.12500,0.00000,0.00260,0.0000,-0.00002}), 0.00003, swp));

    ARIADNE_TEST_BINARY_PREDICATE(refines,rec(xo),ValidatedTaylorModel(E(1,6,  {1.000000,-0.500000, 0.250000,-0.125000, 0.062500,-0.031250, 0.015625}), 0.018, swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,sqrt(xo),ValidatedTaylorModel(E(1,6, {1.000000, 0.250000,-0.031250, 0.007813,-0.002441, 0.000854,-0.000320}), 0.0003, swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,log(xo),ValidatedTaylorModel(E(1,6,  {0.000000, 0.500000,-0.125000, 0.041667,-0.015625, 0.006250,-0.002604}), 0.003, swp));

    // Test exponential based at log2
    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(T(E(1,1,{0.693147,0.5}),0.0,swp)),
                                  T(E(1,6, {2.00000,1.00000,0.25000,0.04166,0.00520,0.00052,0.00004}), 0.00006, swp));

}


Void TestTaylorModel::test_rescale()
{
}

Void TestTaylorModel::test_restrict()
{
}

Void TestTaylorModel::test_refinement()
{
    ValidatedTaylorModel x=ValidatedTaylorModel::coordinate(2,0,swp);
    ValidatedTaylorModel y=ValidatedTaylorModel::coordinate(2,1,swp);
    ValidatedTaylorModel e=ValidatedTaylorModel::error(2,1u,swp);

    // Test refinement with no roundoff errors
    ARIADNE_TEST_EQUAL(refinement(T(E(1,4, {1.0,-0.75,0.0,3.0,3.25}),0.5,swp),T(E(1,4, {1.0,0.0,0.25,2.0,3.0}),1.0,swp)),
        T(E(1,4, {1.0,-0.625,0.0,2.75,3.25}),0.50,swp));

    // Test refinement with roundoff errors
    ARIADNE_TEST_EQUAL(refinement(T(E(1,0, {2./3}),0.5,swp),T(E(1,0, {6./5}),0.25,swp)),
        T(E(1,0, {1.0583333333333331}),0.10833333333333339,swp));
}

Void TestTaylorModel::test_split()
{
    ValidatedTaylorModel x=ValidatedTaylorModel::coordinate(2,0,swp);
    ValidatedTaylorModel y=ValidatedTaylorModel::coordinate(2,1,swp);
    ValidatedTaylorModel z=ValidatedTaylorModel::zero(2,swp);
    ValidatedTaylorModel t=1+3*x+2*y-5*x*x-7*x*y+11*y*y;
    ValidatedTaylorModel es1=-1.75+4*x+5.5*y-1.25*x*x-3.5*x*y+11*y*y;
    ValidatedTaylorModel es2=1+1.5*x+2*y-1.25*x*x-3.5*x*y+11*y*y;
    ValidatedTaylorModel es3=1.25-1*x-1.5*y-1.25*x*x-3.5*x*y+11*y*y;

    ARIADNE_TEST_PRINT(t);
    ARIADNE_TEST_EQUAL(split(t,0,SplitPart::LOWER),es1);
    ARIADNE_TEST_EQUAL(split(t,0,SplitPart::MIDDLE),es2);
    ARIADNE_TEST_EQUAL(split(t,0,SplitPart::UPPER),es3);
}


Void TestTaylorModel::test_antiderivative()
{
    ValidatedFloat unit_interval(-1,+1);
    ValidatedTaylorModel tm=ValidatedTaylorModel::constant(2,1,swp);
    ValidatedTaylorModel atm=antiderivative(tm,1);

    ARIADNE_TEST_EQUAL(antiderivative(T(E({ {{0,0},2.0} }),0.,swp),0),T(E({ {{1,0},2.0} }),0.,swp));
    ARIADNE_TEST_EQUAL(antiderivative(T(E({ {{0,0},2.0} }),0.,swp),1),T(E({ {{0,1},2.0} }),0.,swp));
    ARIADNE_TEST_EQUAL(antiderivative(T(E({ {{1,0},3.0} }),0.,swp),0),T(E({ {{2,0},1.5} }),0.,swp));
    ARIADNE_TEST_EQUAL(antiderivative(T(E({ {{2,0},7.5} }),0.,swp),0),T(E({ {{3,0},2.5} }),0.,swp));
    ARIADNE_TEST_EQUAL(antiderivative(T(E({ {{2,4},7.5} }),0.,swp),0),T(E({ {{3,4},2.5} }),0.,swp));
    ARIADNE_TEST_EQUAL(antiderivative(T(E({ {{2,4},7.5} }),0.,swp),1),T(E({ {{2,5},1.5} }),0.,swp));

    // Test error control
    ValidatedTaylorModel x=ValidatedTaylorModel::coordinate(1,0,swp);
    ValidatedTaylorModel e=ValidatedTaylorModel::zero(1,swp)+ValidatedFloat(-1,+1);
    ARIADNE_TEST_EQUAL(antiderivative(2.0*x*x,0),0.66666666666666663*x*x*x+5.5511151231257827021e-17*e);
    ARIADNE_TEST_EQUAL(antiderivative(2.0*x*x+e,0),0.66666666666666663*x*x*x+1.0000000000000002*e);
    ARIADNE_TEST_EQUAL(antiderivative(T({{{2},2.0}},0.,swp),0),T({{{3},0.66666666666666663}},5.5511151231257827021e-17,swp));
    ARIADNE_TEST_EQUAL(antiderivative(T({{{2},2.0}},1.,swp),0),T({{{3},0.66666666666666663}},1.0000000000000002,swp));

    // Regression test
    T t1=T({ {{0,0},1.}, {{1,0},2.}, {{0,1},3.}, {{2,0},4.}, {{1,1},5.}, {{0,2},6.} }, 0., swp);
    T at1=T({ {{1,0},1.}, {{2,0},1.}, {{1,1},3.}, {{3,0},1.33333333333333333}, {{2,1},2.5}, {{1,2},6.} }, 1.1102230246251565404e-16, swp);
    ARIADNE_TEST_EQUAL(antiderivative(t1,0),at1);
}


Void TestTaylorModel::test_compose()
{

}



namespace Ariadne {
Vector<UpperInterval> range(const Vector<ValidatedTaylorModel>& tm) {
    Vector<UpperInterval> r(tm.size()); for(Nat i=0; i!=tm.size(); ++i) { r[i]=tm[i].range(); } return r; }
}

Int main() {
    ThresholdSweeper sweeper(1e-8);
    TestTaylorModel(sweeper).test();

    return ARIADNE_TEST_FAILURES;
}
