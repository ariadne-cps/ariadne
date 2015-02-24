/***************************************************************************
 *            test_validated_float.cc
 *
 *  Copyright  2006-14  Alberto Casagrande, Pieter Collins
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

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "config.h"
#include "numeric/rational.h"
#include "numeric/float.decl.h"
#include "numeric/float-user.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

namespace Ariadne {
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator==(Float64 const& x, Q const& q) { return Rational(x)==q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator!=(Float64 const& x, Q const& q) { return Rational(x)!=q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator<=(Float64 const& x, Q const& q) { return Rational(x)<=q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator>=(Float64 const& x, Q const& q) { return Rational(x)>=q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator< (Float64 const& x, Q const& q) { return Rational(x)< q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator> (Float64 const& x, Q const& q) { return Rational(x)> q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator==(FloatMP const& x, Q const& q) { return Rational(x)==q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator!=(FloatMP const& x, Q const& q) { return Rational(x)!=q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator<=(FloatMP const& x, Q const& q) { return Rational(x)<=q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator>=(FloatMP const& x, Q const& q) { return Rational(x)>=q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator< (FloatMP const& x, Q const& q) { return Rational(x)< q; }
template<class Q, EnableIf<IsSame<Q,Rational>> =dummy> Bool operator> (FloatMP const& x, Q const& q) { return Rational(x)> q; }

template<class PR> Bool models(Float<Lower,PR> x, Rational q) { return x.raw() <= q; }
template<class PR> Bool models(Float<Upper,PR> x, Rational q) { return x.raw() >= q; }
template<class PR> Bool models(Float<Bounded,PR> x, Rational q) { return x.lower_raw() <= q and x.upper_raw() >= q; }
template<class PR> Bool models(Float<Metric,PR> x, Rational q) { return x.error_raw() >= abs(Rational(x.value_raw())-q); }

template<> String class_name<Precision64>() { return "Precision64"; }
template<> String class_name<PrecisionMP>() { return "PrecisionMP"; }
} // namespace Ariadne

template<class PR>
class TestFloats
{
  public:
    static Rational to_rational(Float<Approximate,PR> x) { return Rational(x.raw()); }
    static Rational to_rational(Float<Lower,PR> x) { return Rational(x.raw()); }
    static Rational to_rational(Float<Upper,PR> x) { return Rational(x.raw()); }
    static Rational to_rational(Float<Exact,PR> x) { return Rational(x.raw()); }
};


template<class PR>
class TestDirectedFloats
    : public TestFloats<PR>
{
    typedef Float<Approximate,PR> ApproximateFloatType;
    typedef Float<Lower,PR> LowerFloatType;
    typedef Float<Upper,PR> UpperFloatType;
    typedef Float<Bounded,PR> BoundedFloatType;
    typedef Float<Metric,PR> MetricFloatType;
    typedef Float<Exact,PR> ExactFloatType;

  private:
    PR precision;
  public:
    TestDirectedFloats(PR pr) : precision(pr) { };
    Void test();
  private:
    Void test_concept();
    Void test_conversions();
    Void test_validation();
    Void test_rounded_arithmetic();
};

template<class PR> Void
TestDirectedFloats<PR>::test()
{
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_validation());
    ARIADNE_TEST_CALL(test_rounded_arithmetic());
}

template<class PR> Void
TestDirectedFloats<PR>::test_concept()
{
    Nat m;
    ApproximateFloatType ax(1);
    LowerFloatType lx(1);
    UpperFloatType ux(1);
    ExactFloatType ex(1);

    lx=+lx; lx=-ux; lx=lx+lx; lx=lx-ux; lx=lx*m; lx=lx/m;
    ex=nul(lx); lx=pos(lx); lx=neg(ux); lx=half(lx);
    lx=add(lx,lx); lx=sub(lx,ux);
    lx=sqrt(lx); lx=exp(lx); lx=log(lx); lx=atan(lx);
    lx=max(lx,lx); lx=min(lx,lx);

    ux=+ux; ux=-lx; ux=ux+ux; ux=ux-lx; ux=ux*m; ux=ux/m;
    ex=nul(ux); ux=pos(ux); ux=neg(lx); ux=half(ux);
    ux=add(ux,ux); ux=sub(ux,lx);
    ux=sqrt(ux); ux=exp(ux); ux=log(ux);
    ux=max(ux,ux); ux=min(ux,ux);
}

template<class PR> Void
TestDirectedFloats<PR>::test_conversions()
{
    Rational one=1;
    Rational five_thirds=5*one/3;
    Rational neg_five_thirds=-5*one/3;

    ARIADNE_TEST_COMPARE(LowerFloatType(five_thirds).raw(),<=,five_thirds);
    ARIADNE_TEST_COMPARE(LowerFloatType(neg_five_thirds).raw(),<=,neg_five_thirds);
    ARIADNE_TEST_COMPARE(UpperFloatType(five_thirds).raw(),>=,five_thirds);
    ARIADNE_TEST_COMPARE(UpperFloatType(neg_five_thirds).raw(),>=,neg_five_thirds);
}

template<class PR> Void
TestDirectedFloats<PR>::test_validation() {
    PR pr=precision;
    Rational one=1;
    Rational two=2;
    ARIADNE_TEST_ASSERT(refines(LowerFloatType(one,pr),LowerFloatType(-one,pr)));
    ARIADNE_TEST_ASSERT(refines(UpperFloatType(-one,pr),UpperFloatType(+one,pr)));
    ARIADNE_TEST_ASSERT(refines(UpperFloatType(-two,pr),UpperFloatType(-one,pr)));
    ARIADNE_TEST_ASSERT(refines(rec(UpperFloatType(-two,pr)),rec(UpperFloatType(-one,pr))));
}

template<class PR> Void
TestDirectedFloats<PR>::test_rounded_arithmetic() {
    PR pr=precision;
    Rational one=1;
    Rational third=one/3;
    Rational fifth=one/3;
    ARIADNE_TEST_COMPARE((LowerFloatType(third,pr)+LowerFloatType(fifth,pr)).raw(),<=,third+fifth);
    ARIADNE_TEST_COMPARE((LowerFloatType(third,pr)-UpperFloatType(fifth,pr)).raw(),<=,third-fifth);
    ARIADNE_TEST_ASSERT(refines(UpperFloatType(third+fifth,pr),UpperFloatType(third,pr)+UpperFloatType(fifth,pr)));
    ARIADNE_TEST_ASSERT(refines(LowerFloatType(third+fifth,pr),LowerFloatType(third,pr)+LowerFloatType(fifth,pr)));
}


template<class PR>
class TestMetricFloat
    : public TestFloats<PR>
{
    typedef RawFloat<PR> RawFloatType;
    typedef Float<Approximate,PR> ApproximateFloatType;
    typedef Float<Lower,PR> LowerFloatType;
    typedef Float<Upper,PR> UpperFloatType;
    typedef Float<Bounded,PR> BoundedFloatType;
    typedef Float<Metric,PR> MetricFloatType;
    typedef Float<Exact,PR> ExactFloatType;
  private:
    PR precision;
  public:
    TestMetricFloat(PR pr) : precision(pr) { };
    Void test();
  private:
    using TestFloats<PR>::to_rational;
    Void test_concept();
    Void test_conversions();
    Void test_validation();
    Void test_rounded_arithmetic();
};

template<class PR> Void
TestMetricFloat<PR>::test()
{
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_rounded_arithmetic());
}

template<class PR> Void
TestMetricFloat<PR>::test_conversions()
{
    Rational one=1;
    Rational third=one/3;
    PR pr=precision;
}

template<class PR> Void
TestMetricFloat<PR>::test_rounded_arithmetic()
{
    Rational one=1;
    Rational three=3;
    Rational five=5;
    Rational six=6;
    Rational third=one/three;
    Rational fifth=one/five;
    PR pr=precision;

    ARIADNE_TEST_BINARY_PREDICATE(models,MetricFloatType(third,pr)+MetricFloatType(fifth,pr),third+fifth);
    ARIADNE_TEST_BINARY_PREDICATE(models,MetricFloatType(third,pr)-MetricFloatType(fifth,pr),third-fifth);
    ARIADNE_TEST_BINARY_PREDICATE(models,MetricFloatType(third,pr)*MetricFloatType(fifth,pr),third*fifth);
    ARIADNE_TEST_BINARY_PREDICATE(models,MetricFloatType(third,pr)/MetricFloatType(fifth,pr),third/fifth);
//    ARIADNE_TEST_BINARY_PREDICATE(models,1/MetricFloatType(three,pr)/MetricFloatType(fifth,pr),1/three);

    ARIADNE_TEST_BINARY_PREDICATE(models,pow(MetricFloatType(three/five,pr),4u),pow(three/five,4u));
    ARIADNE_TEST_BINARY_PREDICATE(models,pow(MetricFloatType(three/five,pr),4),pow(three/five,4));
    ARIADNE_TEST_BINARY_PREDICATE(models,pow(MetricFloatType(three/five,pr),-4),pow(three/five,-4));
    ARIADNE_TEST_BINARY_PREDICATE(models,pow(MetricFloatType(three/five,pr),-7),pow(three/five,-7));
}



template<class PR>
class TestBoundedFloat
{
    typedef RawFloat<PR> RawFloatType;
    typedef Float<Bounded,PR> BoundedFloatType;
    typedef Float<Exact,PR> ExactFloatType;
  private:
    PR pr;
  public:
    TestBoundedFloat(PR prec) : pr(prec) { }
    Void test();
  private:
    Void test_concept();
    Void test_constructors();
    Void test_input();
    Void test_class();
    Void test_comparison();
    Void test_correct_rounded_arithmetic();
    Void test_accurate_rounded_arithmetic();
    Void test_exact_rounded_arithmetic();
    Void test_aliasing();
    Void test_monotone_functions();
    Void test_trigonometric_functions();
    Void regression_tests();
};


template<class PR> Void
TestBoundedFloat<PR>::test()
{
    ARIADNE_TEST_CALL(test_concept());
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_input());
    ARIADNE_TEST_CALL(test_class());
    ARIADNE_TEST_CALL(test_comparison());
    ARIADNE_TEST_CALL(test_correct_rounded_arithmetic());
    ARIADNE_TEST_CALL(test_accurate_rounded_arithmetic());
    ARIADNE_TEST_CALL(test_exact_rounded_arithmetic());
    ARIADNE_TEST_CALL(test_monotone_functions());
    ARIADNE_TEST_CALL(test_trigonometric_functions());
    ARIADNE_TEST_CALL(regression_tests());
}

template<class PR> Void
TestBoundedFloat<PR>::test_concept()
{
    BoundedFloatType::set_output_precision(17);

    Nat m=1;
    Int n=1;
    double d=1;
    ExactFloatType x(1);
    RawFloatType a,b;
    BoundedFloatType vx(1);
    BoundedFloatType rx(1);

    // Constructors
    rx=BoundedFloatType(); rx=BoundedFloatType(n); rx=BoundedFloatType(m); rx=BoundedFloatType(d); rx=BoundedFloatType(x); rx=BoundedFloatType(vx);
    rx=BoundedFloatType(n,n); rx=BoundedFloatType(m,m); rx=BoundedFloatType(d,d); rx=BoundedFloatType(a,b);
    rx=BoundedFloatType(n,m); rx=BoundedFloatType(m,d); rx=BoundedFloatType(d,n);

    // Assignment
    rx=n; rx=m; rx=x; rx=vx;

    // Exact operations
    rx=nul(vx); rx=pos(vx); rx=neg(vx); rx=half(vx); rx=sqr(vx); rx=rec(vx);

    rx=operator+(x,x); rx=operator+(x,vx); rx=operator+(vx,x); rx=operator+(vx,vx);
    rx=operator-(x,x); rx=operator-(x,vx); rx=operator-(vx,x); rx=operator-(vx,vx);
    rx=operator*(x,x); rx=operator*(x,vx); rx=operator*(vx,x); rx=operator*(vx,vx);
    rx=operator/(x,x); rx=operator/(x,vx); rx=operator/(vx,x); rx=operator/(vx,vx);

    // Arithmetic
    rx=add(x,x); rx=add(vx,vx);
    rx=sub(x,x); rx=sub(vx,vx);
    rx=mul(x,x); rx=mul(vx,vx);
    rx=div(x,x); rx=div(vx,vx);
    rx=pow(x,m); rx=pow(x,m);
    rx=pow(x,n); rx=pow(x,n);

    // Order
    rx=max(vx,vx); rx==min(vx,vx); rx=abs(vx);

    // Transcendental functions
    rx=sqrt(vx);
    rx=exp(vx);
    rx=log(vx);
    rx=sin(vx);
    rx=cos(vx);
    rx=tan(vx);
    rx=atan(vx);
}



// Test that interval arithmetic is rounded correctly,
// without paying attention to accuracy issues.
template<class PR> Void
TestBoundedFloat<PR>::test_correct_rounded_arithmetic()
{
    BoundedFloatType onethird=BoundedFloatType(1)/BoundedFloatType(3);
    ARIADNE_TEST_COMPARE( onethird.lower_raw(), < , onethird.upper_raw() );
    BoundedFloatType one_approx=onethird*BoundedFloatType(3);
    ARIADNE_TEST_COMPARE( one_approx.lower_raw(), < , 1.0 );
    ARIADNE_TEST_COMPARE( one_approx.upper_raw(), > , 1.0 );
}


// Test that interval arithmetic gives the most accurate rounded values
template<class PR> Void
TestBoundedFloat<PR>::test_accurate_rounded_arithmetic()
{
    const double min=std::numeric_limits<double>::min();
    const double eps=std::numeric_limits<double>::epsilon();

    ARIADNE_TEST_SAME(BoundedFloatType(1.5)+BoundedFloatType(min),BoundedFloatType(1.5,1.5+eps));
    ARIADNE_TEST_SAME(BoundedFloatType(1.5)-BoundedFloatType(min),BoundedFloatType(1.5-eps,1.5));
    ARIADNE_TEST_SAME(BoundedFloatType(1+eps,1+2*eps)*BoundedFloatType(1+eps,1+3*eps),BoundedFloatType(1+2*eps,1+6*eps));
    ARIADNE_TEST_SAME(BoundedFloatType(1)/BoundedFloatType(3),BoundedFloatType(0.33333333333333331,0.33333333333333337));
    ARIADNE_TEST_SAME(BoundedFloatType(2)/BoundedFloatType(5),BoundedFloatType(0.39999999999999997,0.40000000000000002));

    ARIADNE_TEST_SAME(BoundedFloatType(1.5)+ExactFloatType(min),BoundedFloatType(1.5,1.5+eps));
    ARIADNE_TEST_SAME(BoundedFloatType(1.5)-ExactFloatType(min),BoundedFloatType(1.5-eps,1.5));
    ARIADNE_TEST_SAME(BoundedFloatType(1+eps,1+2*eps)*ExactFloatType(1+eps),BoundedFloatType(1+2*eps,1+4*eps));
    ARIADNE_TEST_SAME(BoundedFloatType(1+3*eps,1+5*eps)/ExactFloatType(1+eps),BoundedFloatType(1+eps,1+4*eps));

    ARIADNE_TEST_SAME(ExactFloatType(min)-BoundedFloatType(1.5),BoundedFloatType(-1.5,eps-1.5));
    ARIADNE_TEST_SAME(ExactFloatType(1+5*eps)/BoundedFloatType(1+2*eps,1+3*eps),BoundedFloatType(1+eps,1+3*eps));

    ARIADNE_TEST_SAME(sqr(BoundedFloatType(1-eps,1+eps)),BoundedFloatType(1-4*eps/2,1+3*eps));

    ARIADNE_TEST_SAME(pow(BoundedFloatType(3,5),-1),BoundedFloatType(0.19999999999999998,0.33333333333333337));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(3,5),-2),BoundedFloatType(0.039999999999999986955,0.11111111111111114658));

    ARIADNE_TEST_SAME(rec(BoundedFloatType(1+2*eps,1+5*eps)),BoundedFloatType(1-10*eps/2,1-3*eps/2));

}


// Test that interval arithmetic gives exact values if possible
template<class PR> Void
TestBoundedFloat<PR>::test_exact_rounded_arithmetic()
{
    ARIADNE_TEST_SAME(BoundedFloatType(5,7)+BoundedFloatType(2,4),BoundedFloatType(7,11));
    ARIADNE_TEST_SAME(BoundedFloatType(5,7)-BoundedFloatType(2,6),BoundedFloatType(-1,5));

    ARIADNE_TEST_SAME(BoundedFloatType(5,7)*BoundedFloatType(2,4),BoundedFloatType(10,28));
    ARIADNE_TEST_SAME(BoundedFloatType(5,7)*BoundedFloatType(-2,4),BoundedFloatType(-14,28));
    ARIADNE_TEST_SAME(BoundedFloatType(5,7)*BoundedFloatType(-4,-2),BoundedFloatType(-28,-10));
    ARIADNE_TEST_SAME(BoundedFloatType(-7,5)*BoundedFloatType(2,4),BoundedFloatType(-28,20));
    ARIADNE_TEST_SAME(BoundedFloatType(-7,5)*BoundedFloatType(-2,4),BoundedFloatType(-28,20));
    ARIADNE_TEST_SAME(BoundedFloatType(-7,5)*BoundedFloatType(-4,-2),BoundedFloatType(-20,28));
    ARIADNE_TEST_SAME(BoundedFloatType(-7,-5)*BoundedFloatType(2,4),BoundedFloatType(-28,-10));
    ARIADNE_TEST_SAME(BoundedFloatType(-7,-5)*BoundedFloatType(-2,4),BoundedFloatType(-28,14));
    ARIADNE_TEST_SAME(BoundedFloatType(-7,-5)*BoundedFloatType(-4,-2),BoundedFloatType(10,28));

    ARIADNE_TEST_SAME(BoundedFloatType(5,7)/BoundedFloatType(2,4),BoundedFloatType(1.25,3.50));
    ARIADNE_TEST_SAME(BoundedFloatType(5,7)/BoundedFloatType(-4,-2),BoundedFloatType(-3.50,-1.25));
    ARIADNE_TEST_SAME(BoundedFloatType(-7,5)/BoundedFloatType(2,4),BoundedFloatType(-3.50,2.50));
    ARIADNE_TEST_SAME(BoundedFloatType(-7,5)/BoundedFloatType(-4,-2),BoundedFloatType(-2.50,3.5));
    ARIADNE_TEST_SAME(BoundedFloatType(-7,-5)/BoundedFloatType(2,4),BoundedFloatType(-3.50,-1.25));
    ARIADNE_TEST_SAME(BoundedFloatType(-7,-5)/BoundedFloatType(-4,-2),BoundedFloatType(1.25,3.50));

    ARIADNE_TEST_SAME(pow(BoundedFloatType(5,7),0u),BoundedFloatType(1,1));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-5,7),0u),BoundedFloatType(1,1));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,5),0u),BoundedFloatType(1,1));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,-5),0u),BoundedFloatType(1,1));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(5,7),1u),BoundedFloatType(5,7));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-5,7),1u),BoundedFloatType(-5,7));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,5),1u),BoundedFloatType(-7,5));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,-5),1u),BoundedFloatType(-7,-5));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(5,7),2u),BoundedFloatType(25,49));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-5,7),2u),BoundedFloatType(0,49));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,5),2u),BoundedFloatType(0,49));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,-5),2u),BoundedFloatType(25,49));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(5,7),3u),BoundedFloatType(125,343));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-5,7),3u),BoundedFloatType(-125,343));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,5),3u),BoundedFloatType(-343,125));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,-5),3u),BoundedFloatType(-343,-125));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(5,7),4u),BoundedFloatType(625,2401));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-5,7),4u),BoundedFloatType(0,2401));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,5),4u),BoundedFloatType(0,2401));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,-5),4u),BoundedFloatType(625,2401));

    ARIADNE_TEST_SAME(pow(BoundedFloatType(5,7),0),BoundedFloatType(1,1));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,5),0),BoundedFloatType(1,1));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,-5),0),BoundedFloatType(1,1));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(5,7),1),BoundedFloatType(5,7));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-5,7),1),BoundedFloatType(-5,7));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,5),1),BoundedFloatType(-7,5));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,-5),1),BoundedFloatType(-7,-5));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(5,7),2),BoundedFloatType(25,49));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-5,7),2),BoundedFloatType(0,49));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,5),2),BoundedFloatType(0,49));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,-5),2),BoundedFloatType(25,49));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(5,7),3),BoundedFloatType(125,343));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-5,7),3),BoundedFloatType(-125,343));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,5),3),BoundedFloatType(-343,125));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,-5),3),BoundedFloatType(-343,-125));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(5,7),4),BoundedFloatType(625,2401));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-5,7),4),BoundedFloatType(0,2401));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,5),4),BoundedFloatType(0,2401));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,-5),4),BoundedFloatType(625,2401));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(5,7),5),BoundedFloatType(3125,16807));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-5,7),5),BoundedFloatType(-3125,16807));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,5),5),BoundedFloatType(-16807,3125));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,-5),5),BoundedFloatType(-16807,-3125));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(5,7),7),BoundedFloatType(78125,823543));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-5,7),7),BoundedFloatType(-78125,823543));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,5),7),BoundedFloatType(-823543,78125));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-7,-5),7),BoundedFloatType(-823543,-78125));

    ARIADNE_TEST_SAME(pow(BoundedFloatType(2,4),-1),BoundedFloatType(0.25,0.5));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-4,-2),-1),BoundedFloatType(-0.5,-0.25));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(2,4),-2),BoundedFloatType(0.0625,0.25));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-4,-2),-2),BoundedFloatType(0.0625,0.25));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(2,4),-3),BoundedFloatType(0.015625,0.125));
    ARIADNE_TEST_SAME(pow(BoundedFloatType(-4,-2),-3),BoundedFloatType(-0.125,-0.015625));

    ARIADNE_TEST_SAME(rec(BoundedFloatType(2,4)),BoundedFloatType(0.25,0.50));
    ARIADNE_TEST_SAME(rec(BoundedFloatType(-4,-2)),BoundedFloatType(-0.50,-0.25));
}



template<> Void
TestBoundedFloat<Precision64>::test_constructors()
{
    Float64 zero=0;

    // Construct from pair
    BoundedFloatType xd1(Float64(1.125),Float64(2.25));
    ARIADNE_TEST_ASSERT(xd1.lower_raw()==1.125); ARIADNE_TEST_ASSERT(xd1.upper_raw()==2.25);

    // Default constructor
    BoundedFloatType xd2;
    if(xd2.lower_raw()>xd2.upper_raw()) {
        ARIADNE_TEST_WARN("BoundedFloatType default constructor returns an empty set.");
    } else {
        ARIADNE_TEST_BINARY_PREDICATE(same,xd2,BoundedFloatType(zero,zero));
    }

    // Constructor with approximations
    BoundedFloatType xd3(Rational(21,10));
    ARIADNE_TEST_COMPARE(Rational(xd3.lower_raw()),<,Rational(21,10));
    ARIADNE_TEST_COMPARE(Rational(xd3.upper_raw()),>,Rational(21,10));

    // Constructor from approximate values
    BoundedFloatType xd4(2.1,3.2);
    ARIADNE_TEST_COMPARE(xd4.lower_raw(),<=,2.1);
    ARIADNE_TEST_COMPARE(xd4.upper_raw(),>=,3.2);

    // Approximate constructor from a single value
    BoundedFloatType xd5(Rational(1,3));
    ARIADNE_TEST_COMPARE(Rational(xd5.lower_raw()),<,Rational(1,3));
    ARIADNE_TEST_COMPARE(Rational(xd5.upper_raw()),>,Rational(1,3));

    // Exact constructor from a single value
    BoundedFloatType xd6(Float64(1.25));
    ARIADNE_TEST_EQUAL(xd6.lower_raw(),Float64(1.25));
    ARIADNE_TEST_EQUAL(xd6.upper_raw(),Float64(1.25));
}

template<class PR> Void
TestBoundedFloat<PR>::test_constructors()
{
    RawFloatType zero=0;

    // Construct from pair
    BoundedFloatType xd1(RawFloatType(1.125),RawFloatType(2.25));
    ARIADNE_TEST_ASSERT(xd1.lower_raw()==1.125); ARIADNE_TEST_ASSERT(xd1.upper_raw()==2.25);

    // Default constructor
    BoundedFloatType xd2;
    if(xd2.lower_raw()>xd2.upper_raw()) {
        ARIADNE_TEST_WARN("BoundedFloatType default constructor returns an empty set.");
    } else {
        ARIADNE_TEST_BINARY_PREDICATE(same,xd2,BoundedFloatType(zero,zero));
    }

    // Constructor with approximations
    BoundedFloatType xd3(Rational(21,10));
    ARIADNE_TEST_COMPARE(Rational(xd3.lower_raw()),<,Rational(21,10));
    ARIADNE_TEST_COMPARE(Rational(xd3.upper_raw()),>,Rational(21,10));

    // Constructor from approximate values
    BoundedFloatType xd4(2.1,3.2);
    ARIADNE_TEST_COMPARE(xd4.lower_raw(),<=,2.1);
    ARIADNE_TEST_COMPARE(xd4.upper_raw(),>=,3.2);

    // Approximate constructor from a single value
    BoundedFloatType xd5(Rational(1,3));
    ARIADNE_TEST_COMPARE(Rational(xd5.lower_raw()),<,Rational(1,3));
    ARIADNE_TEST_COMPARE(Rational(xd5.upper_raw()),>,Rational(1,3));

    // Exact constructor from a single value
    BoundedFloatType xd6(RawFloatType(1.25));
    ARIADNE_TEST_EQUAL(xd6.lower_raw(),RawFloatType(1.25));
    ARIADNE_TEST_EQUAL(xd6.upper_raw(),RawFloatType(1.25));
}

template<class PR> Void TestBoundedFloat<PR>::test_class()
{
    // Test lower, upper, midpoint, radius, width
    RawFloatType one=1.0;
    RawFloatType two=2.0;

    // Tests for exact operations
    ARIADNE_TEST_EQUAL(BoundedFloatType(-0.25,0.50).lower().raw(),-0.25);
    ARIADNE_TEST_EQUAL(BoundedFloatType(-0.25,0.50).upper().raw(),0.5);
    ARIADNE_TEST_EQUAL(BoundedFloatType(-0.25,0.50).value().raw(),0.125);
    ARIADNE_TEST_EQUAL(BoundedFloatType(-0.25,0.50).error().raw(),0.375)

    // Tests for inexact operations
    ARIADNE_TEST_EQUAL(BoundedFloatType(-1./3,2./3).lower().raw(),-0.33333333333333331483);
    ARIADNE_TEST_EQUAL(BoundedFloatType(-1./3,2./3).upper().raw(),0.66666666666666662966);
    ARIADNE_TEST_EQUAL(BoundedFloatType(-1./3,2./3).value().raw(),0.16666666666666665741);
    ARIADNE_TEST_EQUAL(BoundedFloatType(-1./3,2./3).error().raw(),0.5)

    // Tests for inexact operations
    ARIADNE_TEST_EQUAL(BoundedFloatType(div_down(-one,3),div_up(two,3)).lower().raw(),-0.33333333333333337034);
    ARIADNE_TEST_EQUAL(BoundedFloatType(div_down(-one,3),div_up(two,3)).upper().raw(),0.66666666666666674068);
    ARIADNE_TEST_EQUAL(BoundedFloatType(div_down(-one,3),div_up(two,3)).value().raw(),0.16666666666666668517);
    ARIADNE_TEST_EQUAL(BoundedFloatType(div_down(-one,3),div_up(two,3)).error().raw(),0.50000000000000011102);
}

template<class PR> Void TestBoundedFloat<PR>::test_input()
{
    BoundedFloatType x;

    string input("[1.125,2.75] [0.4,0.6]");
    stringstream iss(input);

    iss >> x;
    ARIADNE_TEST_EQUALS(x.lower_raw(),Rational(9,8))
    ARIADNE_TEST_EQUALS(x.upper_raw(),Rational(11,4))

    iss >> x;
    // TODO: Make stream input produce rounded values
    //ARIADNE_TEST_COMPARE(x.lower_raw(),<,Rational(2,5))
    // ARIADNE_TEST_COMPARE(x.upper_raw(),>,Rational(3,5))
    if(not(x.lower_raw()<=Rational(2,5) and x.upper_raw()>=Rational(3,5))) {
        ARIADNE_TEST_WARN("Float<Bounded,"<<class_name<PR>()<<"> string constructor returns an approximate interval, not an outwardly rounded interval.");
    }
}

template<class PR> Void TestBoundedFloat<PR>::test_comparison() {
    // FIXME: If using Boost style interval tests, uncomment the line below
    // and comment out the line after
    //ARIADNE_TEST_ASSERT(indeterminate(ivld1==ivld2));
    BoundedFloatType ivl1(1.125,2.25);
    BoundedFloatType ivl2=ivl1;

    ARIADNE_TEST_ASSERT(!definitely(ivl1==ivl2));
    ARIADNE_TEST_ASSERT(possibly(ivl1==ivl2));
    BoundedFloatType& ivl1ref=ivl1;
    ivl1ref=BoundedFloatType(5.25,7.375);
    cout << "ivl1ref=" << ivl1ref << endl;
    ARIADNE_TEST_ASSERT(ivl1ref.lower_raw()==RawFloatType(5.25));
}

template<class PR> Void TestBoundedFloat<PR>::test_aliasing() {
    ExactFloatType ex2(1.5);
    ExactFloatType ex3(2.25);

    BoundedFloatType vx1;
    BoundedFloatType vx2(1.5,2.25);
    BoundedFloatType vx3(3.125,4.0625);

    // Check to make sure aliases are handled correctly
    vx1=vx3; vx1=vx2-vx1; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,BoundedFloatType(vx2-vx3));
    vx1=vx3; vx1=vx2*vx1; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,BoundedFloatType(vx2*vx3));
    vx1=vx2; vx1=vx1*vx3; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,BoundedFloatType(vx2*vx3));
    vx1=vx2; vx1=vx1*ex3; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,BoundedFloatType(vx2*ex3));
    vx1=vx3; vx1=ex2*vx1; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,BoundedFloatType(ex2*vx3));
    vx1=vx2; vx1=vx1/vx3; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,BoundedFloatType(vx2/vx3));
    vx1=vx2; vx1=vx1/ex3; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,BoundedFloatType(vx2/ex3));
    vx1=vx3; vx1=ex2/vx1; ARIADNE_TEST_BINARY_PREDICATE(same,vx1,BoundedFloatType(ex2/vx3));
}

template<class PR> Void TestBoundedFloat<PR>::test_monotone_functions()
{

    BoundedFloatType two(2.0);
    BoundedFloatType sqrttwo=sqrt(two);
    ARIADNE_TEST_PRINT(sqrttwo);
    ARIADNE_TEST_COMPARE(sqrttwo.lower_raw(),<=,1.4142135623730949);
    ARIADNE_TEST_COMPARE(sqrttwo.lower_raw(),> ,1.4142135623730947);
    ARIADNE_TEST_COMPARE(sqrttwo.upper_raw(),>=,1.4142135623730951);
    ARIADNE_TEST_COMPARE(sqrttwo.upper_raw(),< ,1.4142135623730954);

    BoundedFloatType one(1.0);
    BoundedFloatType expone=exp(one);
    ARIADNE_TEST_PRINT(expone);
    ARIADNE_TEST_COMPARE(expone.lower_raw(),<,2.71828182845905);
    ARIADNE_TEST_COMPARE(expone.lower_raw(),>,2.71828182845903);
    ARIADNE_TEST_COMPARE(expone.upper_raw(),>,2.71828182845904);
    ARIADNE_TEST_COMPARE(expone.upper_raw(),<,2.71828182845906);
    ARIADNE_TEST_ASSERT(expone.lower_raw()<expone.upper_raw());

    BoundedFloatType e(2.7182818284590451,2.7182818284590455);
    BoundedFloatType loge=log(e);
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_COMPARE(loge.lower_raw(),<,1.0);
    ARIADNE_TEST_COMPARE(loge.lower_raw(),>,0.9999999999998);
    ARIADNE_TEST_COMPARE(loge.upper_raw(),>,1.0);
    ARIADNE_TEST_COMPARE(loge.upper_raw(),<,1.000000000002);
}

template<class PR> Void TestBoundedFloat<PR>::test_trigonometric_functions()
{
    try {
        BoundedFloatType x(6.283185307179586,6.283185307179587);
        BoundedFloatType sinx=sin(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(sinx.lower_raw(),<,0.0);
        ARIADNE_TEST_COMPARE(sinx.lower_raw(),>,-1e-14);
        ARIADNE_TEST_COMPARE(sinx.upper_raw(),>,0.0);
        ARIADNE_TEST_COMPARE(sinx.upper_raw(),<,+1e-14);
        ARIADNE_TEST_ASSERT(sinx.lower_raw()<sinx.upper_raw());
    }
    catch(...) { }

    try {
        BoundedFloatType x(7.0685834705770345);
        BoundedFloatType sinx=sin(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(sinx.lower_raw(),<,0.7071067811866);
        ARIADNE_TEST_COMPARE(sinx.upper_raw(),>,0.7071067811865);
        ARIADNE_TEST_ASSERT(sinx.lower_raw()<sinx.upper_raw());
    }
    catch(...) { }

}

template<class PR> Void TestBoundedFloat<PR>::regression_tests() {
    RawFloatType inf=RawFloatType::inf();

    // Regression test; fails dramatically on certain types of rounding
    {
        BoundedFloatType x(1.5707963267948966,1.5707963267948968);
        BoundedFloatType cosx=cos(x);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_COMPARE(cosx.lower_raw(),<,0.0);
        ARIADNE_TEST_COMPARE(cosx.lower_raw(),>,-1e-14);
        ARIADNE_TEST_COMPARE(cosx.upper_raw(),>,0.0);
        ARIADNE_TEST_COMPARE(cosx.upper_raw(),<,+1e-14);
        ARIADNE_TEST_ASSERT(cosx.lower_raw()<cosx.upper_raw());
    }

    // Regression test for dividing by interval with lower endpoint -0.0 or upper endpoint +0.0

    ARIADNE_TEST_EQUAL((BoundedFloatType(1.0,2.0)/BoundedFloatType(-0.0,1.0)).upper_raw(),+inf);
    ARIADNE_TEST_EQUAL((BoundedFloatType(1.0,2.0)/BoundedFloatType(-1.0,+0.0)).lower_raw(),-inf);

    ARIADNE_TEST_EQUAL(rec(BoundedFloatType(-0.0,+1.0)).upper_raw(),+inf);
    ARIADNE_TEST_EQUAL(rec(BoundedFloatType(-1.0,+0.0)).lower_raw(),-inf);
}

Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestDirectedFloats<Precision64>(Precision64()).test();
    TestDirectedFloats<PrecisionMP>(PrecisionMP(64)).test();

    TestMetricFloat<Precision64>(Precision64()).test();
    TestMetricFloat<PrecisionMP>(PrecisionMP(64)).test();

    TestBoundedFloat<Precision64>(Precision64()).test();
    TestBoundedFloat<PrecisionMP>(PrecisionMP(64)).test();

    return ARIADNE_TEST_FAILURES;
}

