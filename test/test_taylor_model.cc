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
#include "numeric.h"
#include "real.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "taylor_model.h"
#include "differential.h"
#include "function.h"
#include "models.h"

#include "test.h"
using namespace Ariadne;

Vector<Float> v(uint n, uint i) { return Vector<Float>::unit(n,i); }
TaylorModel ctm(uint m, double c) { return TaylorModel::constant(m,c); }
TaylorModel ctm(uint m) { return TaylorModel::constant(m,1.0); }
TaylorModel tm(uint m, uint i) { return TaylorModel::variable(m,i); }


class TestTaylorModel
{
    typedef MultiIndex MI;
    typedef Expansion<Float> E;
    typedef Polynomial<Float> P;
    typedef TaylorModel T;
  public:
    void test();
  private:
    void test_concept();
    void test_constructors();
    void test_predicates();
    void test_approximation();
    void test_evaluate();
    void test_arithmetic();
    void test_functions();
    void test_rescale();
    void test_restrict();
    void test_intersection();
    void test_split();
    void test_antiderivative();
    void test_compose();
    void test_flow();
    void test_solve();
    void test_implicit();
};


void TestTaylorModel::test()
{
    std::cerr<<std::setprecision(17);
    std::cout<<std::setprecision(17);
    std::clog<<std::setprecision(17);

    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_predicates());
    ARIADNE_TEST_CALL(test_approximation());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_functions());
    ARIADNE_TEST_CALL(test_rescale());
    ARIADNE_TEST_CALL(test_restrict());
    ARIADNE_TEST_CALL(test_intersection());
    ARIADNE_TEST_CALL(test_split());
    ARIADNE_TEST_CALL(test_antiderivative());
    ARIADNE_TEST_CALL(test_compose());
    ARIADNE_TEST_CALL(test_flow());
    ARIADNE_TEST_CALL(test_solve());
    ARIADNE_TEST_CALL(test_implicit());
}


void TestTaylorModel::test_concept()
{
    const Float f=0.0;
    const Interval i;
    const Vector<Float> vf;
    const Vector<Interval> vi;
    const TaylorModel  t;
    TaylorModel tr;

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

    tr.sweep(); tr.truncate(); tr.clean();

    t.evaluate(vi); evaluate(t,vi);
    t.domain(); t.range(); t.expansion(); t.error();

}

void TestTaylorModel::test_constructors()
{
    ARIADNE_TEST_CONSTRUCT(TaylorModel,tv1,(E(2,3, 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0), 0.25));

    ARIADNE_ASSERT_EQUAL(tv1.value(),1.0);
    ARIADNE_ASSERT_EQUAL(tv1.error(),0.25);
}

void TestTaylorModel::test_predicates()
{
    TaylorModel tv1(E(1,2, 1.00,2.00,3.00), 0.75);
    TaylorModel tv2(E(1,2, 1.00,1.75,3.25), 0.25);
    TaylorModel tv3(E(1,2, 1.125,1.75,3.25), 0.25);
    TaylorModel tv4(E(1,3, 1.00,2.25,3.00,-0.25), 0.25);

    ARIADNE_TEST_BINARY_PREDICATE(refines,tv1,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(refines,tv2,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(!refines,tv3,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(refines,tv4,tv1);
}

void TestTaylorModel::test_approximation()
{
    ARIADNE_TEST_CONSTRUCT(TaylorModel,tv2,(E(1,2,1.0,2.0,3.0),0.25));
}

void TestTaylorModel::test_evaluate()
{
    Vector<Interval> iv(2, 0.25,0.5, -0.75,-0.5);
    TaylorModel tv(E(2,2,1.0,2.0,3.0,4.0,5.0,6.0),0.25);
    ARIADNE_TEST_EQUAL(evaluate(tv,iv),Interval(-1,1));
}

void TestTaylorModel::test_arithmetic()
{
    //Operations which can be performed exactly with floating-point arithmetic.
    ARIADNE_TEST_EQUAL(TaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75)+(-3), TaylorModel(E(1,2, -2.0,-2.0,3.0), 0.75));
    ARIADNE_TEST_EQUAL(TaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75)-(-3), TaylorModel(E(1,2, 4.0,-2.0,3.0), 0.75));
    ARIADNE_TEST_EQUAL(TaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75)*(-3), TaylorModel(E(1,2, -3.0,6.0,-9.0), 2.25));
    ARIADNE_TEST_EQUAL(TaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75)/(-4), TaylorModel(E(1,2, -0.25,0.5,-0.75), 0.1875));
    ARIADNE_TEST_EQUAL(TaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75)+Interval(-1,2), TaylorModel(E(1,2, 1.5,-2.0,3.0), 2.25));
    ARIADNE_TEST_EQUAL(TaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75)-Interval(-1,2), TaylorModel(E(1,2, 0.5,-2.0,3.0), 2.25));
    ARIADNE_TEST_EQUAL(TaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75)*Interval(-1,2), TaylorModel(E(1,2, 0.5,-1.0,1.5), 10.5));
    ARIADNE_TEST_EQUAL(TaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75)/Interval(0.25,2.0), TaylorModel(E(1,2, 2.25,-4.5,6.75), 13.5));
    ARIADNE_TEST_EQUAL(+TaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75), TaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75));
    ARIADNE_TEST_EQUAL(-TaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75), TaylorModel(E(1,2, -1.0,2.0,-3.0), 0.75));
    ARIADNE_TEST_EQUAL(TaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75)+TaylorModel(E(1,2, 3.0,2.0,-4.0), 0.5), TaylorModel(E(1,2, 4.0,0.0,-1.0), 1.25));
    ARIADNE_TEST_EQUAL(TaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75)-TaylorModel(E(1,2, 3.0,2.0,-4.0), 0.5), TaylorModel(E(1,2, -2.0,-4.0,7.0), 1.25));
    ARIADNE_TEST_EQUAL(TaylorModel(E(1,2, 0.0,0.0,3.0), 0.75)*TaylorModel(E(1,2, 3.0,2.0,-4.0), 0.5), TaylorModel(E(1,4, 0.0,0.0,9.0,6.0,-12.0), 8.625));
    ARIADNE_TEST_EQUAL(TaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75)*TaylorModel(E(1,2, 3.0,2.0,-4.0), 0.5), TaylorModel(E(1,4, 3.0,-4.0,1.0,14.0,-12.0), 10.125));
}

void TestTaylorModel::test_functions()
{
    TaylorModel x(E(1,1, 0.0, 1.0), 0.0);
    TaylorModel xz(E(1,1, 0.0, 0.5), 0.0);
    TaylorModel xo(E(1,1, 1.0, 0.5), 0.0);

    ARIADNE_TEST_PRINT(exp(x));
    ARIADNE_TEST_PRINT(sin(x));
    ARIADNE_TEST_PRINT(cos(x));

    //Functions based on their natural defining points with variable dependence 0.5
    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(T(E(1,1,0.0,1.0),0.0)),T(E(1,6, 1.0,1.00,0.500,0.1667,0.0417,0.0083,0.0014), 0.0004));
    ARIADNE_TEST_BINARY_PREDICATE(refines,sin(x),T(E(1,6, 0.0,1.0000,0.0,-0.1667,0.0,0.0083,0.0), 0.0003));
    ARIADNE_TEST_BINARY_PREDICATE(refines,cos(x),T(E(1,6, 1.0000,0.0,-0.5000,0.0,0.0417,0.0,-0.0014), 0.0003));

    //Functions based on their natural defining points with variable dependence 0.5
    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(xz),TaylorModel(E(1,6, 1.00000,0.50000,0.12500,0.02083,0.00260,0.00026,0.00002), 0.00003));
    ARIADNE_TEST_BINARY_PREDICATE(refines,sin(xz),TaylorModel(E(1,6, 0.00000,0.50000,0.0000,-0.02083,0.00000,0.00026,0.00000), 0.00003));
    ARIADNE_TEST_BINARY_PREDICATE(refines,cos(xz),TaylorModel(E(1,6, 1.00000,0.0000,-0.12500,0.00000,0.00260,0.0000,-0.00002), 0.00003));

    ARIADNE_TEST_BINARY_PREDICATE(refines,rec(xo),TaylorModel(E(1,6,  1.000000,-0.500000, 0.250000,-0.125000, 0.062500,-0.031250, 0.015625), 0.018));
    ARIADNE_TEST_BINARY_PREDICATE(refines,sqrt(xo),TaylorModel(E(1,6, 1.000000, 0.250000,-0.031250, 0.007813,-0.002441, 0.000854,-0.000320), 0.0003));
    ARIADNE_TEST_BINARY_PREDICATE(refines,log(xo),TaylorModel(E(1,6,  0.000000, 0.500000,-0.125000, 0.041667,-0.015625, 0.006250,-0.002604), 0.003));

    // Test exponential based at log2
    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(T(E(1,1,0.693147,0.5))),
                                  T(E(1,6, 2.00000,1.00000,0.25000,0.04166,0.00520,0.00052,0.00004), 0.00006));

}


void TestTaylorModel::test_rescale()
{
}

void TestTaylorModel::test_restrict()
{
}

void TestTaylorModel::test_intersection()
{
    TaylorModel x=tm(2,0); TaylorModel y=tm(2,1); TaylorModel e=tm(2,0)*0+Interval(-1,1);

    // Test intersection with no roundoff errors
    ARIADNE_TEST_EQUAL(intersection(T(E(1,4, 1.0,-0.75,0.0,3.0,3.25),0.5),T(E(1,4, 1.0,0.0,0.25,2.0,3.0),1.0)),
        T(E(1,4, 1.0,-0.625,0.0,2.75,3.25),0.50));

    // Test intersection with roundoff errors
    ARIADNE_TEST_EQUAL(intersection(T(E(1,0, 2./3),0.5),T(E(1,0, 6./5),0.25)),
        T(E(1,0, 1.0583333333333331),0.10833333333333339));
}

void TestTaylorModel::test_split()
{
    TaylorModel x=tm(2,0); TaylorModel y=tm(2,1);
    TaylorModel t=1+3*x+2*y-5*x*x-7*x*y+11*y*y;
    TaylorModel es1=-1.75+4*x+5.5*y-1.25*x*x-3.5*x*y+11*y*y;
    TaylorModel es2=1.25-1*x-1.5*y-1.25*x*x-3.5*x*y+11*y*y;

    ARIADNE_TEST_PRINT(t);
    ARIADNE_TEST_EQUAL(split(t,0).first,es1);
    ARIADNE_TEST_EQUAL(split(t,0).second,es2);
    ARIADNE_TEST_EQUAL(split(t,0,false),es1);
    ARIADNE_TEST_EQUAL(split(t,0,true),es2);
}


void TestTaylorModel::test_antiderivative()
{
    Interval unit_interval(-1,+1);
    TaylorModel tm=TaylorModel::constant(2,1.0);
    TaylorModel atm=antiderivative(tm,1);

    ARIADNE_TEST_EQUAL(antiderivative(T(E(2,1, 0,0,2.0),0.),0),T(E(2,1, 1,0,2.0)));
    ARIADNE_TEST_EQUAL(antiderivative(T(E(2,1, 0,0,2.0),0.),1),T(E(2,1, 0,1,2.0)));
    ARIADNE_TEST_EQUAL(antiderivative(T(E(2,1, 1,0,3.0),0.),0),T(E(2,1, 2,0,1.5)));
    ARIADNE_TEST_EQUAL(antiderivative(T(E(2,1, 2,0,7.5),0.),0),T(E(2,1, 3,0,2.5)));
    ARIADNE_TEST_EQUAL(antiderivative(T(E(2,1, 2,4,7.5),0.),0),T(E(2,1, 3,4,2.5)));
    ARIADNE_TEST_EQUAL(antiderivative(T(E(2,1, 2,4,7.5),0.),1),T(E(2,1, 2,5,1.5)));

    // Test error control
    ARIADNE_TEST_EQUAL(antiderivative(T(E(1,1, 2,2.0),0.),0),T(E(1,1, 3,0.66666666666666663),5.5511151231257827021e-17));
    ARIADNE_TEST_EQUAL(antiderivative(T(E(1,1, 2,2.0),1.),0),T(E(1,1, 3,0.66666666666666663),1.0000000000000002));

    // Regression test for
    T t1=T(E(2,6, 0,0,1., 1,0,2., 0,1,3., 2,0,4., 1,1,5., 0,2,6.), 0.);
    T at1=T(E(2,6, 1,0,1., 2,0,1., 1,1,3., 3,0,1.33333333333333333, 2,1,2.5, 1,2,6.), 1.1102230246251565404e-16);
    ARIADNE_TEST_EQUAL(antiderivative(t1,0),at1);
}


void TestTaylorModel::test_compose()
{
}


namespace Ariadne {
Vector<TaylorModel> _implicit1(const Vector<TaylorModel>& f, uint n=4);
Vector<TaylorModel> _implicit2(const Vector<TaylorModel>& f, uint n=4);
Vector<TaylorModel> _implicit3(const Vector<TaylorModel>& f, uint n=4);
Vector<TaylorModel> _implicit4(const Vector<TaylorModel>& f, uint n=4);
Vector<TaylorModel> _implicit5(const Vector<TaylorModel>& f, uint n=4);
}

void TestTaylorModel::test_solve()
{
    TaylorModel f(E(1,3, 0,1.0, 1,4.0, 2,1.0),0.125);
    Interval x=solve(Vector<TaylorModel>(1u,f))[0];
    ARIADNE_TEST_PRINT(x);
    ARIADNE_TEST_PRINT(f.evaluate(Vector<Interval>(1u,x)));
}


void TestTaylorModel::test_implicit()
{
    TaylorModel f(E(2,5, 0,0,0.0000002, 1,0,1.000000000000003, 2,0,0.000000000000003, 0,1,4.000000000000001, 0,2,1.000000000000001),0.0);
    //TaylorModel ha=implicit_approx(Vector<TaylorModel>(1u,f),8)[0];
    //TaylorModel h1=_implicit1(Vector<TaylorModel>(1u,f),4)[0];
    //TaylorModel h2=_implicit2(Vector<TaylorModel>(1u,f),4)[0];
    //TaylorModel h3=_implicit3(Vector<TaylorModel>(1u,f),4)[0];
    //std::cerr<<"\n\nh1="<<h1<<"\nh2="<<h2<<"\nh3="<<h3<<std::endl;
    TaylorModel h2=_implicit2(Vector<TaylorModel>(1u,f))[0];
    TaylorModel h5=_implicit5(Vector<TaylorModel>(1u,f))[0];
    ARIADNE_TEST_PRINT(h2);
    ARIADNE_TEST_PRINT(h5);

    {
        //Test computation of sqrt(4+x)-2 on [-1,+1] by solving 4+x-(y+2)*(y+2)=0
        TaylorModel f2=4+tm(2,0)-(sqr(tm(2,1)+2));
        TaylorModel i2=tm(1,0);
        TaylorModel h2=implicit(f2);
        ARIADNE_TEST_BINARY_PREDICATE(operator<,mag(compose(f2,join(i2,h2)).range()),1e-6);
    }

    return;
    TaylorModel h=implicit(f);
    TaylorModel id(E(1,1, 1,1.0),0.0);
    TaylorModel z(1);
    TaylorModel c=compose(f,join(id,h));
    TaylorModel hh=implicit_step(f,h);
    // Compute the power series for square root
    TaylorModel s(1); double cc=2; MultiIndex a(1);
    for(int i=1; i!=24; ++i) { cc*=(((2*i-3)*0.25)/(2*i)); ++a[0]; s[a]=cc; }
    ARIADNE_TEST_PRINT(s);


    ARIADNE_TEST_PRINT(f);
    ARIADNE_TEST_PRINT(h);
    ARIADNE_TEST_PRINT(s);
    ARIADNE_TEST_PRINT(id);
    ARIADNE_TEST_PRINT(join(id,h));
    ARIADNE_TEST_PRINT(compose(f,join(id,h)));
    ARIADNE_TEST_PRINT(h-s);
    ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(c),1e-2);
    ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(h-s),1e-4);
    ARIADNE_TEST_BINARY_PREDICATE(refines,hh,h);
    ARIADNE_TEST_BINARY_PREDICATE(refines,z,c);
    ARIADNE_TEST_BINARY_PREDICATE(refines,s,h);
    Float he=h.error(); h.set_error(0);
    std::cerr<<"\n\n";
    std::cerr<<"hh="<<hh<<"\nh="<<h<<"\n";
    TaylorModel d=h-hh;
    std::cerr<<"h-hh="<<d<<"\n";
    std::cerr<<"|h-hh|="<<norm(h-hh)<<" he="<<he<<"\n\n\n";

}

namespace Ariadne {
Vector<Interval> range(const Vector<TaylorModel>& tm) {
    Vector<Interval> r(tm.size()); for(uint i=0; i!=tm.size(); ++i) { r[i]=tm[i].range(); } return r; }
}

void TestTaylorModel::test_flow()
{
    {
        // Test flow dx/dt=3/8 on an assymetric domain
        // phi=x0+3t/8 x0=0.25s+0.25
        Vector<TaylorModel> vf=ctm(1,3.0/8)*v(1,0);
        Vector<TaylorModel> d1=(0.25+0.25*tm(1,0))*v(1,0);
        Vector<TaylorModel> phi1=flow(vf,d1,2);
        Vector<TaylorModel> expected_phi1=(0.25*ctm(2)+0.25*tm(2,0)+0.375*tm(2,1))*v(1,0);
        ARIADNE_TEST_EQUAL(phi1,expected_phi1);
   }

    // Test flow dx/dt=1/4; dy/dt=y/4 on domain [0.5,0.5]x[0.5,0.5]
    //
    Vector<TaylorModel> vf=(ctm(2,2.0)*v(2,0)+tm(2,1)*v(2,1))*0.25;
    Vector<TaylorModel> dom=0.5*tm(2,0)*v(2,0)+0.5*tm(2,1)*v(2,1);
    uint o=6;

    Vector<TaylorModel> phi=flow(vf,dom,o);
    Vector<TaylorModel> next_phi=antiderivative(compose(vf,phi),2)+embed(dom,1u);
    ARIADNE_TEST_PRINT(vf);
    ARIADNE_TEST_PRINT(phi);
    ARIADNE_TEST_PRINT(next_phi);
    ARIADNE_TEST_BINARY_PREDICATE(operator<,(norm(Ariadne::range(phi-next_phi))),0.1);


}


int main() {
    TestTaylorModel().test();

    return ARIADNE_TEST_FAILURES;
}
