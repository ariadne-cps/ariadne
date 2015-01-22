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

ValidatedTaylorModel operator^(ValidatedTaylorModel tm, Nat m) { return pow(tm,m); }
ValidatedTaylorModel clobber(ValidatedTaylorModel tm) { tm.clobber(); return std::move(tm); }

template<class A> void display_difference(A const& res, A const& exp) {
    auto tol=exp.error();
    const_cast<A&>(exp).clobber();
    auto dif=res-exp;
    auto err=norm(dif);
    std::cerr<<"\nres="<<res<<"\nexp="<<exp<<"\ndif="<<dif<<"\nerr="<<err<<"\ntol="<<tol<<"\n";
}



class TestTaylorModel
{
    typedef MultiIndex MI;
    typedef Expansion<Float> E;
    typedef Polynomial<Float> P;
    typedef ValidatedTaylorModel T;
  public:
    Sweeper swp;
    ValidatedTaylorModel z,o,x,y,e;
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
    Void test_recondition();
};


TestTaylorModel::TestTaylorModel(Sweeper sweeper)
    : swp(sweeper)
    , z(ValidatedTaylorModel::zero(2,swp))
    , o(ValidatedTaylorModel::constant(2,1.0_exact,swp))
    , x(ValidatedTaylorModel::coordinate(2,0,swp))
    , y(ValidatedTaylorModel::coordinate(2,1,swp))
    , e(ValidatedTaylorModel::ball(2,1u,swp))
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
    ARIADNE_TEST_CALL(test_recondition());
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
    ARIADNE_TEST_BINARY_PREDICATE(refines,1+x+(x^2)/2+(x^3)/4,1+x+(x^2)/2+e/4);
    ARIADNE_TEST_BINARY_PREDICATE(not refines,1+x+(x^2)/2+(x^3)/6+e/pow(2.0_exact,31),1+x+(x^2)/2+e/6);
    ValidatedTaylorModel tv1=1+2*x+3*(x^2)+e*3/4;
    ValidatedTaylorModel tv2=1+x*7/4+(x^2)*13/4+e/4;
    ValidatedTaylorModel tv3=o*9/8+x*7/4+(x^2)*13/4+e/4;
    ValidatedTaylorModel tv4=1+x*9/4+(x^2)*3-(x^4)/4+e/4;

    ARIADNE_TEST_BINARY_PREDICATE(refines,tv1,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(refines,tv2,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(not refines,tv3,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(refines,tv4,tv1);

    ExactFloat h(0.5);
    ARIADNE_TEST_BINARY_PREDICATE(consistent,1+2*x+3*y+e*3/4,1+2*x+3*y+e*3/4);
    ARIADNE_TEST_BINARY_PREDICATE(consistent,1+e*3/4,2+e/4);
    ARIADNE_TEST_BINARY_PREDICATE(not consistent,x-(x^3)*2/3+e/(3.0_exact+pow(h,20)),z);

    ARIADNE_TEST_BINARY_PREDICATE(inconsistent,1-2*(x^2),e*7/8);
    ARIADNE_TEST_BINARY_PREDICATE(not inconsistent,1+2*x+3*y+e*3/4,1+2*x+3*y+e*3/4);
    ARIADNE_TEST_BINARY_PREDICATE(not inconsistent,1+e*3/4,2+e/4);

    if(not refines(1-(x^2)/2,e/2)) {
        ARIADNE_TEST_WARN("refines(ValidatedTaylorModel,ValidatedTaylorModel) may return false even if the first model refines the second.");
    }

    if(not consistent(x-(x^3)*2/3+e/3,z)) {
        ARIADNE_TEST_WARN("consistent(ValidatedTaylorModel,ValidatedTaylorModel) may return false even if models are consistent with representing the same function.");
    }

    if(not inconsistent(x-(x^3)*2/3+e/4,z)) {
        ARIADNE_TEST_WARN("inconsistent(ValidatedTaylorModel,ValidatedTaylorModel) may return false even if models are consistent with representing the same function.");
    }

}


Void TestTaylorModel::test_approximation()
{
    ARIADNE_TEST_CONSTRUCT(ValidatedTaylorModel,tv2,(E(1,2, {1.0,2.0,3.0}),0.25, swp));
}


Void TestTaylorModel::test_unscale()
{

    ExactInterval singleton(1.0,1.0);
    if(unscale(3*o,singleton).codomain()!=singleton) {
        ARIADNE_TEST_WARN("Unscaling over singleton domain does not yield constant");
    }
}

Void TestTaylorModel::test_evaluate()
{
    Vector<ValidatedFloat> iv={{0.25,0.5},{-0.75,-0.5}};
    ValidatedTaylorModel tv=1+2*x+3*y+4*(x^2)+5*(x*y)+6*(y^2)+e/2;
    ARIADNE_TEST_EQUAL(evaluate(tv,iv),ValidatedFloat(-1,1));
}


Void TestTaylorModel::test_arithmetic()
{
    {
        ValidatedTaylorModel x=ValidatedTaylorModel::coordinate(1,0,swp);
        ValidatedTaylorModel e=ValidatedTaylorModel::error(1,1.0_exact,swp);
        ARIADNE_TEST_EQUALS(1-2*x+3*(x^2)+e*3/4, ValidatedTaylorModel({{{0},1.0},{{1},-2.0},{{2},3.0}}, 0.75,swp));
    }

    //Operations which can be performed exactly with floating-point arithmetic.
    ARIADNE_TEST_EQUALS(ValidatedTaylorModel(E(1,2, {1.0,-2.0,3.0}),0.75,swp)+(-3), ValidatedTaylorModel(E(1,2, {-2.0,-2.0,3.0}), 0.75,swp));
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

    ValidatedTaylorModel o=ValidatedTaylorModel::constant(2,1.0_exact,swp);
    ValidatedTaylorModel x=ValidatedTaylorModel::coordinate(2,0,swp);
    ValidatedTaylorModel y=ValidatedTaylorModel::coordinate(2,1,swp);
    ValidatedTaylorModel e=ValidatedTaylorModel::ball(2,1.0_exact,swp);

    ValidatedTaylorModel t=1-2*x+3*y;

    ARIADNE_TEST_EQUAL((1-2*x+3*y)*(3+2*x-4*y),3-4*x+5*y-4*(x^2)+14*(x*y)-12*(y^2));
    ARIADNE_TEST_EQUAL(sqr(t),t*t);
    ARIADNE_TEST_EQUAL(pow(t,0),o);
    ARIADNE_TEST_EQUAL(pow(t,1),t);
    ARIADNE_TEST_EQUAL(pow(t,2),t*t);
    ARIADNE_TEST_EQUAL(pow(t,3),t*t*t);

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
    // Test range of cubic, which should be exact
    ValidatedTaylorModel t1 = x*x*x+x;
    ARIADNE_TEST_EQUALS(t1.range(),ExactInterval(-2,+2));

    //x^2+x = (x+1/2)^2-1/4; Range [-0.25,+2.0]
    // Test range of quadratic, which could be exact, but need not be
    ValidatedTaylorModel t2 = x*x+x;
    ARIADNE_TEST_BINARY_PREDICATE(refines,t2.range(),ExactInterval(-2,+2));
    ARIADNE_TEST_BINARY_PREDICATE(refines,ExactInterval(-0.25,+2),t2.range());
    if(make_exact_interval(t2.range())!=ExactInterval(-0.25,2.0)) {
        ARIADNE_TEST_WARN("ValidatedTaylorModel::range() not exact for quadratic functions."); }
}

Void TestTaylorModel::test_functions()
{
    ValidatedTaylorModel x=ValidatedTaylorModel::coordinate(1,0,swp);
    ValidatedTaylorModel hx=x/2;
    ValidatedTaylorModel ophx=1+x/2;
    ValidatedFloat e(-1,+1);

    ARIADNE_TEST_PRINT(exp(x));
    ARIADNE_TEST_PRINT(sin(x));
    ARIADNE_TEST_PRINT(cos(x));
    ARIADNE_TEST_PRINT(x.tolerance());

    // Expected tolerance based on sweeper characteristics
    static const ExactFloat xtol=ExactFloat(x.tolerance());
    static const ValidatedFloat tol=xtol*ValidatedFloat(-1,+1);
    ExactFloat LAXITY=1;

    // exp, sin and cos have error bound e^N/N!*(1+1/N), where e is bound for |x| N is the first term omitted
    // Error bound for rec is e^(N-1); log is e^(N-1)/N; sqrt is ???, where e is bound for |x-1|

    //Functions based on their natural defining points with variable dependence 1.0
    ValidatedTaylorModel expected_exp_x = 1+x+(x^2)/2+(x^3)/6+(x^4)/24+(x^5)/120+(x^6)/720+e/4410+tol;
    ValidatedTaylorModel expected_sin_x = x-(x^3)/6+(x^5)/120+e/4410+tol;
    ValidatedTaylorModel expected_cos_x = 1-(x^2)/2+(x^4)/24-(x^6)/720+e/35840+tol;
    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(x),expected_exp_x);
    ARIADNE_TEST_BINARY_PREDICATE(refines,sin(x),expected_sin_x);
    ARIADNE_TEST_BINARY_PREDICATE(refines,cos(x),expected_cos_x);

    //Functions based on their natural defining points with variable dependence 0.5
    ValidatedTaylorModel expected_exp_hx = 1+hx+(hx^2)/2+(hx^3)/6+(hx^4)/24+(hx^5)/120+(hx^6)/720+e/128/4410+tol;
    ValidatedTaylorModel expected_sin_hx = hx-(hx^3)/6+(hx^5)/120+e/128/4410+tol;
    ValidatedTaylorModel expected_cos_hx = 1-(hx^2)/2+(hx^4)/24-(hx^6)/720+e/128/4410+tol;
    // Uncomment the last line so that expected_cos_hx contains error term base on x^8 term rather than x^7 term.
    //ValidatedTaylorModel expected_cos_hx = 1-(hx^2)/2+(hx^4)/24-(hx^6)/720+e/256/35840+tol;

    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(hx),expected_exp_hx);
    ARIADNE_TEST_BINARY_PREDICATE(refines,sin(hx),expected_sin_hx);
    ARIADNE_TEST_BINARY_PREDICATE(refines,cos(hx),expected_cos_hx);

    ValidatedTaylorModel expected_rec_ophx = 1-hx+(hx^2)-(hx^3)+(hx^4)-(hx^5)+(hx^6)+e/64+tol;
    ValidatedTaylorModel expected_sqrt_ophx = 1+hx/2-(hx^2)/8+(hx^3)/16-(hx^4)*5/128+(hx^5)*7/256-(hx^6)*21/1024+e/64+tol;
    ValidatedTaylorModel expected_log_ophx = hx-(hx^2)/2+(hx^3)/3-(hx^4)/4+(hx^5)/5-(hx^6)/6+e/64/7+tol;
    ARIADNE_TEST_BINARY_PREDICATE(refines,rec(ophx),expected_rec_ophx);
    ARIADNE_TEST_BINARY_PREDICATE(refines,sqrt(ophx),expected_sqrt_ophx);
    ARIADNE_TEST_BINARY_PREDICATE(refines,log(ophx),expected_log_ophx);

    // Doubling formulae
    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(2*x),expected_exp_x^2);
    ARIADNE_TEST_BINARY_PREDICATE(refines,sin(2*x),2*expected_sin_x*expected_cos_x);
    ARIADNE_TEST_BINARY_PREDICATE(refines,cos(2*x),2*(expected_cos_x^2)-1);

    ARIADNE_TEST_BINARY_PREDICATE(refines,rec(3*ophx),expected_rec_ophx/3);
    ARIADNE_TEST_BINARY_PREDICATE(refines,sqrt(9*ophx),expected_sqrt_ophx*3);
    ARIADNE_TEST_BINARY_PREDICATE(refines,log(3*ophx),expected_log_ophx+log(3));

    Nat rn=3; ExactFloat c(2); ExactFloat r(rn);
    ARIADNE_TEST_PRINT(c);
    ARIADNE_TEST_PRINT(r);
    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(c+r*hx),exp(c)*pow(expected_exp_hx,rn)+tol);
    ARIADNE_TEST_BINARY_PREDICATE(refines,sin(c+hx),sin(c)*expected_cos_hx+cos(c)*expected_sin_hx+tol);
    ARIADNE_TEST_BINARY_PREDICATE(refines,cos(c+hx),cos(c)*expected_cos_hx-sin(c)*expected_sin_hx+tol);

    c=ExactFloat(10);
    ARIADNE_TEST_BINARY_PREDICATE(refines,sin(c+hx),sin(c)*expected_cos_hx+cos(c)*expected_sin_hx+LAXITY*tol);
    ARIADNE_TEST_BINARY_PREDICATE(refines,cos(c+hx),cos(c)*expected_cos_hx-sin(c)*expected_sin_hx+LAXITY*tol);

    // Test exponential based at log2; exp(log(2)+x/2)=2*exp(x/2)
    ExactFloat log2_apprx(0.693147);
    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(log2_apprx+x/2),
                                  2*(1+hx+(hx^2)/2+(hx^3)/6+(hx^4)/24+(hx^5)/120+(hx^6)/720)+e*6/100000);
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

    ARIADNE_TEST_BINARY_PREDICATE(refines,1+x+(x^2)/2+(x^3)/4,1+x+(x^2)/2+e/4);
    ARIADNE_TEST_BINARY_PREDICATE(not refines,1+x+(x^2)/2+(x^3)/6+e/pow(2.0_exact,31),1+x+(x^2)/2+e/6);

    // Test refinement with no roundoff errors
    ARIADNE_TEST_EQUAL(refinement(1-x*3/4+(x^3)*3+(x^4)*13/4+e/2,1+(x^2)/4+(x^3)*2+(x^4)*3+e),1-x*5/8+(x^3)*11/4+(x^4)*13/4+e/2);

    // Test refinement with roundoff errors
    ARIADNE_TEST_EQUAL(refinement(ExactFloat(2.0/3)*x+e/2,ExactFloat(6.0/5)*x+e/4),
        ExactFloat(1.05833333333333335)*x+ExactFloat(0.108333333333333393)*e);

    // Code below computes expected values for second test
    // Float xv=2./3; Float xe=1./2; Float yv=6./5; Float ye=1./4;
    // Float rl=sub_down(yv,ye); Float ru=add_up(xv,xe); Float rv=add_near(rl,ru)/2; Float re=sub_up(ru,rl)/2;
    // std::cerr << std::setprecision(18) << "xv="<<xv<<" yv="<<yv<<" rl="<<rl<<" ru="<<ru<<" rv="<<rv<<" re="<<re<<"\n";


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

    ARIADNE_TEST_EQUAL(antiderivative(2*o,0),2*x);
    ARIADNE_TEST_EQUAL(antiderivative(2*o,1),2*y);
    ARIADNE_TEST_EQUAL(antiderivative(3*x,0),(x^2)*3/2);
    ARIADNE_TEST_EQUAL(antiderivative(2*o+3*x,0),x*2+(x^2)*3/2);
    ARIADNE_TEST_EQUAL(antiderivative((x^2)*15/2,0),(x^3)*5/2);
    ARIADNE_TEST_EQUAL(antiderivative((x^2)*(y^4)*15/2,0),(x^3)*(y^4)*5/2);
    ARIADNE_TEST_EQUAL(antiderivative((x^2)*(y^4)*15/2,1),(x^2)*(y^5)*3/2);

    // Test error control
    ValidatedTaylorModel x=ValidatedTaylorModel::coordinate(1,0,swp);
    ValidatedTaylorModel e=ValidatedTaylorModel::zero(1,swp)+ValidatedFloat(-1,+1);
    ARIADNE_TEST_EQUAL(antiderivative(2.0*x*x,0),0.66666666666666663*x*x*x+5.5511151231257827021e-17*e);
    ARIADNE_TEST_EQUAL(antiderivative(2.0*x*x+e,0),0.66666666666666663*x*x*x+1.0000000000000002*e);
    ARIADNE_TEST_EQUAL(antiderivative(2*(x^2),0),ExactFloat(0.66666666666666663)*(x^3)+ExactFloat(5.5511151231257827021e-17)*e);
    ARIADNE_TEST_EQUAL(antiderivative(2*(x^2)+e,0),ExactFloat(0.66666666666666663)*(x^3)+ExactFloat(1.0000000000000002)*e);

    // Regression test
    ValidatedTaylorModel t1({ {{0,0},1.}, {{1,0},2.}, {{0,1},3.}, {{2,0},4.}, {{1,1},5.}, {{0,2},6.} }, 0., swp);
    ValidatedTaylorModel at1({ {{1,0},1.}, {{2,0},1.}, {{1,1},3.}, {{3,0},1.33333333333333333}, {{2,1},2.5}, {{1,2},6.} }, 1.1102230246251565404e-16, swp);
    ARIADNE_TEST_EQUAL(antiderivative(t1,0),at1);
}


Void TestTaylorModel::test_compose()
{
    ARIADNE_TEST_EQUAL(compose(2-x*x-y/4,{2-x*x-y/4,x}),-2-(x^4)-(y^2)/16+4*(x^2)+y-(x^2)*y/2-x/4);

}

Void TestTaylorModel::test_recondition()
{
    Sweeper swp;
    ARIADNE_TEST_CONSTRUCT( ValidatedTaylorModel, tm1, ({ {{0,0},2.0}, {{1,0},3.0}, {{0,1},5.0} }, 0.5, swp) );
    ARIADNE_TEST_CONSTRUCT( ValidatedTaylorModel, tm2, ({ {{0,0,1},0.5}, {{0,0,0},2.0}, {{1,0,0},3.0}, {{0,1,0},5.0} }, 0.0, swp) );
    ARIADNE_TEST_EQUAL(embed_error(tm1),tm2);
    ARIADNE_TEST_EQUAL(discard_variables(tm2,{2}),tm1);
}




Int main() {
    ThresholdSweeper sweeper(1e-8);
    TestTaylorModel(sweeper).test();

    return ARIADNE_TEST_FAILURES;
}
