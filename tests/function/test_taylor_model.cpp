/***************************************************************************
 *            test_taylor_model.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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
#include <iomanip>
#include "config.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/multi_index.hpp"
#include "function/taylor_model.hpp"
#include "algebra/differential.hpp"
#include "function/function.hpp"
#include "function/polynomial.hpp"

#include "../test.hpp"

using std::cout; using std::cerr; using std::endl;
using namespace Ariadne;

extern template Ariadne::Nat Ariadne::Error<Ariadne::FloatDP>::output_places;
extern template Ariadne::Nat Ariadne::Error<Ariadne::FloatMP>::output_places;

inline Dyadic operator"" _exd (long double x) { return Dyadic(x); }

template<class T, class = decltype(declval<T>().clobber())> True has_clobber(int);
template<class T> False has_clobber(...);
template<class T> using HasClobber = decltype(has_clobber<T>(1));

template<class T, EnableIf<HasClobber<T>> =dummy> void do_clobber(T& t) { t=t.clobber(); }
template<class T, DisableIf<HasClobber<T>> =dummy> void do_clobber(T& t) { }

template<class T, class = decltype(declval<T>().norm())> decltype(declval<T>().norm()) norm(T& t, int=0) { return t.norm(); }
template<class T> T norm(T const& t) { return t; }

#define ARIADNE_TEST_REFINES(expression,expected)                         \
    {                                                                   \
        std::cout << "refines(" << #expression << "," << #expected << "): " << std::flush; \
        auto result = (expression); \
        Bool ok = decide(refines(result,(expected)));                       \
        if(ok) {                                                        \
            std::cout << "true\n" << std::endl;                         \
        } else {                                                        \
            ++ARIADNE_TEST_FAILURES;                                    \
            std::cout << "false\nERROR: refines(" << #expression << "," << #expected << "):\n"; \
            std::cout << "    " << #expression << "=" << (result) << "\n    #expected="<< (expected) << std::endl; \
            std::cerr << "ERROR: " << __FILE__ << ":" << __LINE__ << ": " << __PRETTY_FUNCTION__ << ": "; \
            std::cerr << "`refines(" << #expression << "," << #expected << ")' failed;\n"; \
            std::cerr << "result=" << (result) << "\nexpected="<< (expected) << std::endl; \
            auto exact_result=result; auto exact_expected=expected; clobber(exact_expected); clobber(exact_result); \
            auto difference=(exact_result-exact_expected); \
            std::cerr << "difference=" << difference << std::endl; \
            std::cerr<<"  norm(difference)="<<norm(difference)<<" result.error="<<(result).error()<<" expected.error="<<(expected).error()<<"\n"; \
        }                                                               \
    }                                                                   \

inline Vector<Dyadic> v(Nat n, Nat i) { return Vector<Dyadic>::unit(n,i); }
template<class F> ValidatedTaylorModel<F> ctm(Nat m, double c, Sweeper<F> swp) {
    typedef typename F::PrecisionType PR; return ValidatedTaylorModel<F>::constant(m,FloatValue<PR>(c),swp); }
template<class F> ValidatedTaylorModel<F> ctm(Nat m, Sweeper<F> swp) {
    typedef typename F::PrecisionType PR; return ValidatedTaylorModel<F>::constant(m,FloatValue<PR>(1.0),swp); }
//template<class F> ValidatedTaylorModel<F> tm(Nat m, Nat i, Sweeper swp) { return ValidatedTaylorModel<F>::coordinate(m,i,swp); }

// Allow addition and multiplication by double
template<class F, class D, EnableIf<IsSame<D,double>> =dummy> ValidatedTaylorModel<F> operator+(D d, ValidatedTaylorModel<F> tm) {
    typedef typename F::PrecisionType PR; return FloatValue<PR>(d)+tm; }
template<class F, class D, EnableIf<IsSame<D,double>> =dummy> ValidatedTaylorModel<F> operator-(D d, ValidatedTaylorModel<F> tm) {
    typedef typename F::PrecisionType PR; return FloatValue<PR>(d)-tm; }
template<class F, class D, EnableIf<IsSame<D,double>> =dummy> ValidatedTaylorModel<F> operator*(D d, ValidatedTaylorModel<F> tm) {
    typedef typename F::PrecisionType PR; return FloatValue<PR>(d)*tm; }

// Allow carat for power
template<class F> ValidatedTaylorModel<F> operator^(ValidatedTaylorModel<F> tm, Nat m) { return pow(tm,m); }
template<class F> ValidatedTaylorModel<F> clobber(ValidatedTaylorModel<F> tm) { tm.clobber(); return tm; }

template<class A> void display_difference(A const& res, A const& exp) {
    auto tol=exp.error();
    const_cast<A&>(exp).clobber();
    auto dif=res-exp;
    auto err=norm(dif);
    std::cerr<<"\nres="<<res<<"\nexp="<<exp<<"\ndif="<<dif<<"\nerr="<<err<<"\ntol="<<tol<<"\n";
}

template<class F> class TestTaylorModel
{
    typedef PrecisionType<F> PR;
    typedef MultiIndex MI;
    typedef Expansion<MI,F> E;
    typedef Polynomial<F> P;
    typedef ValidatedTaylorModel<F> T;
    typedef Expansion<MI,F> ExpansionType;
  public:
    Sweeper<F> swp; PR pr;
    ValidatedTaylorModel<F> x0,x1,x2,r,o;
  public:
    typedef F FloatType;
  public:
    TestTaylorModel(Sweeper<F> swp);
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


template<class F> TestTaylorModel<F>::TestTaylorModel(Sweeper<F> sweeper)
    : swp(sweeper), pr(sweeper.precision())
    , x0(ValidatedTaylorModel<F>::coordinate(2,0,swp))
    , x1(ValidatedTaylorModel<F>::coordinate(2,1,swp))
    , x2(ValidatedTaylorModel<F>::zero(2,swp))
    , r(ValidatedTaylorModel<F>::unit_ball(2,swp))
    , o(ValidatedTaylorModel<F>::constant(2,Dyadic(1),swp))
{
}

template<class F> Void TestTaylorModel<F>::test()
{
    std::cerr<<std::setprecision(18);
    std::cout<<std::setprecision(18);
    std::clog<<std::setprecision(18);
    FloatValue<PR>::set_output_places(18);
    FloatBounds<PR>::set_output_places(18);

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


template<class F> Void TestTaylorModel<F>::test_concept()
{
    const FloatValue<PR> f={0,pr};
    const FloatBounds<PR> i;
    const Vector<FloatValue<PR>> vf;
    const Vector<FloatBounds<PR>> vi;
    const ValidatedTaylorModel<F>  t(0,swp);
    ValidatedTaylorModel<F> tr(0,swp);

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

template<class F> Void TestTaylorModel<F>::test_constructors()
{


    ARIADNE_TEST_CONSTRUCT(ExpansionType,raw_expansion,({ {{0,0},1.0}, {{1,0},2.0}, {{0,1},3.0}, {{2,0},4.0}, {{1,1},5.0}, {{0,2},6.0}, {{3,0},7.0}, {{2,1},8.0}, {{1,2},9.0}, {{0,3},10.0} },pr));
    ARIADNE_TEST_CONSTRUCT(F,raw_error,(0.25,pr));
    ARIADNE_TEST_CONSTRUCT(ExpansionType,expansion,(raw_expansion));
    ARIADNE_TEST_CONSTRUCT(FloatError<PR>,error,(raw_error));

    ARIADNE_TEST_CONSTRUCT(ValidatedTaylorModel<F>,tv,(raw_expansion, raw_error, swp));
    ARIADNE_TEST_PRINT(tv.expansion());
    ARIADNE_TEST_PRINT(tv.error());
    ARIADNE_TEST_EQUAL(tv.value().precision(),swp.precision());

    ARIADNE_TEST_EQUALS(tv.value(),1.0);
    ARIADNE_TEST_EQUALS(tv.error(),0.25);
    ARIADNE_TEST_EQUALS(tv.norm(),55.25);

    ARIADNE_TEST_EQUALS((tv[{0,0}]),1.0);
    ARIADNE_TEST_EQUALS((tv[{1,0}]),2.0);
    ARIADNE_TEST_EQUALS((tv[{0,1}]),3.0);
    ARIADNE_TEST_EQUALS((tv[{2,1}]),8.0);

}

template<class F> Void TestTaylorModel<F>::test_predicates()
{
    ARIADNE_TEST_REFINES(1+x0+(x0^2)/2+(x0^3)/4,1+x0+(x0^2)/2+r/4);
    ARIADNE_TEST_BINARY_PREDICATE(not refines,1+x0+(x0^2)/2+(x0^3)/6+r/1048576,1+x0+(x0^2)/2+r/6);

    ARIADNE_TEST_CONSTRUCT(ValidatedTaylorModel<F>,tv1,(1+2*x0+3*(x0^2)+r*3/4));
    ARIADNE_TEST_CONSTRUCT(ValidatedTaylorModel<F>,tv2,(1+x0*7/4+(x0^2)*13/4+r/4));
    ARIADNE_TEST_CONSTRUCT(ValidatedTaylorModel<F>,tv3,(o*9/8+x0*7/4+(x0^2)*13/4+r/4));
    ARIADNE_TEST_CONSTRUCT(ValidatedTaylorModel<F>,tv4,(1+x0*9/4+(x0^2)*3-(x0^4)/4+r/4));

    ARIADNE_TEST_REFINES(tv1,tv1);
    ARIADNE_TEST_REFINES(tv2,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(not refines,tv3,tv1);
    ARIADNE_TEST_REFINES(tv4,tv1);

    FloatValue<PR> h(Dyadic(0.5),pr);
    ARIADNE_TEST_BINARY_PREDICATE(consistent,1+2*x0+3*x1+r*3/4,1+2*x0+3*x1+r*3/4);
    ARIADNE_TEST_BINARY_PREDICATE(consistent,1+r*3/4,2+r/4);
    ARIADNE_TEST_BINARY_PREDICATE(not consistent,x0-(x0^3)*2/3+r/(3.0_exd+pow(h,20)),x2);
    ARIADNE_TEST_BINARY_PREDICATE(inconsistent,1-2*(x0^2),r*7/8);
    ARIADNE_TEST_BINARY_PREDICATE(not inconsistent,1+2*x0+3*x1+r*3/4,1+2*x0+3*x1+r*3/4);
    ARIADNE_TEST_BINARY_PREDICATE(not inconsistent,1+r*3/4,2+r/4);

    if(not refines(1-(x0^2)/2,r/2)) {
        ARIADNE_TEST_WARN("refines(ValidatedTaylorModel<F>,ValidatedTaylorModel<F>) may return false even if the first model refines the second.");
    }

    if(not consistent(x0-(x0^3)*2/3+r/3,x2)) {
        ARIADNE_TEST_WARN("consistent(ValidatedTaylorModel<F>,ValidatedTaylorModel<F>) may return false even if models are consistent with representing the same function.");
    }

    if(not inconsistent(x0-(x0^3)*2/3+r/4,x2)) {
        ARIADNE_TEST_WARN("inconsistent(ValidatedTaylorModel<F>,ValidatedTaylorModel<F>) may return false even if models are consistent with representing the same function.");
    }

}


template<class F> Void TestTaylorModel<F>::test_approximation()
{
    ARIADNE_TEST_CONSTRUCT(ValidatedTaylorModel<F>,tv2,(E({{{0,0},1.0},{{1,0},2.0},{{0,1},3.0}}),0.25, swp));
}


template<class F> Void TestTaylorModel<F>::test_unscale()
{

    ExactIntervalType is_singleton(1.0,1.0);
    if(unscale(3*o,is_singleton).codomain()!=is_singleton) {
        ARIADNE_TEST_WARN("Unscaling over is_singleton domain does not yield constant");
    }
}

template<class F> Void TestTaylorModel<F>::test_evaluate()
{
    Vector<FloatBounds<PR>> iv={{0.25_exact,0.5_exact},{-0.75_exact,-0.5_exact}};
    ValidatedTaylorModel<F> tv=1+2*x0+3*x1+4*(x0^2)+5*(x0*x1)+6*(x1^2)+r/2;
    ARIADNE_TEST_EQUAL(evaluate(tv,iv),FloatBounds<PR>(-1,1));
}


template<class F> Void TestTaylorModel<F>::test_arithmetic()
{
    //Operations which can be performed exactly with floating-point arithmetic.
    ARIADNE_TEST_SAME(ValidatedTaylorModel<F>({{{0},1.0},{{1},-2.0},{{2},3.0}},0.75,swp)+(-3), ValidatedTaylorModel<F>({{{0},-2.0},{{1},-2.0},{{2},3.0}}, 0.75,swp));
    ARIADNE_TEST_SAME(ValidatedTaylorModel<F>({{{0},1.0},{{1},-2.0},{{2},3.0}},0.75,swp)-(-3), ValidatedTaylorModel<F>({{{0},4.0},{{1},-2.0},{{2},3.0}}, 0.75,swp));
    ARIADNE_TEST_SAME(ValidatedTaylorModel<F>({{{0},1.0},{{1},-2.0},{{2},3.0}},0.75,swp)*(-3), ValidatedTaylorModel<F>({{{0},-3.0},{{1},6.0},{{2},-9.0}}, 2.25,swp));
    ARIADNE_TEST_SAME(ValidatedTaylorModel<F>({{{0},1.0},{{1},-2.0},{{2},3.0}},0.75,swp)/(-4), ValidatedTaylorModel<F>({{{0},-0.25},{{1},0.5},{{2},-0.75}}, 0.1875,swp));
    ARIADNE_TEST_SAME(ValidatedTaylorModel<F>({{{0},1.0},{{1},-2.0},{{2},3.0}},0.75,swp)+FloatBounds<PR>(-1,2), ValidatedTaylorModel<F>({{{0},1.5},{{1},-2.0},{{2},3.0}}, 2.25,swp));
    ARIADNE_TEST_SAME(ValidatedTaylorModel<F>({{{0},1.0},{{1},-2.0},{{2},3.0}},0.75,swp)-FloatBounds<PR>(-1,2), ValidatedTaylorModel<F>({{{0},0.5},{{1},-2.0},{{2},3.0}}, 2.25,swp));
    ARIADNE_TEST_SAME(ValidatedTaylorModel<F>({{{0},1.0},{{1},-2.0},{{2},3.0}},0.75,swp)*FloatBounds<PR>(-1,2), ValidatedTaylorModel<F>({{{0},0.5},{{1},-1.0},{{2},1.5}}, 10.5,swp));
    ARIADNE_TEST_SAME(ValidatedTaylorModel<F>({{{0},1.0},{{1},-2.0},{{2},3.0}},0.75,swp)/FloatBounds<PR>(0.25,2.0), ValidatedTaylorModel<F>({{{0},2.25},{{1},-4.5},{{2},6.75}}, 13.5,swp));
    ARIADNE_TEST_SAME(+ValidatedTaylorModel<F>({{{0},1.0},{{1},-2.0},{{2},3.0}},0.75,swp), ValidatedTaylorModel<F>({{{0},1.0},{{1},-2.0},{{2},3.0}},0.75,swp));
    ARIADNE_TEST_SAME(-ValidatedTaylorModel<F>({{{0},1.0},{{1},-2.0},{{2},3.0}},0.75,swp), ValidatedTaylorModel<F>({{{0},-1.0},{{1},2.0},{{2},-3.0}}, 0.75,swp));
    ARIADNE_TEST_SAME(ValidatedTaylorModel<F>({{{0},1.0},{{1},-2.0},{{2},3.0}},0.75,swp)+ValidatedTaylorModel<F>({{{0},3.0},{{1},2.0},{{2},-4.0}},0.5,swp),
                      ValidatedTaylorModel<F>({{{0},4.0},{{1},0.0},{{2},-1.0}}, 1.25,swp));
    ARIADNE_TEST_SAME(ValidatedTaylorModel<F>({{{0},1.0},{{1},-2.0},{{2},3.0}},0.75,swp)-ValidatedTaylorModel<F>({{{0},3.0},{{1},2.0},{{2},-4.0}},0.5,swp),
                      ValidatedTaylorModel<F>({{{0},-2.0},{{1},-4.0},{{2},7.0}}, 1.25,swp));
    ARIADNE_TEST_SAME(ValidatedTaylorModel<F>({{{0},0.0},{{1},0.0},{{2},3.0}}, 0.75,swp)*ValidatedTaylorModel<F>({{{0},3.0},{{1},2.0},{{2},-4.0}},0.5,swp),
                      ValidatedTaylorModel<F>({{{0},0.0},{{1},0.0},{{2},9.0},{{3},6.0},{{4},-12.0}}, 8.625,swp));
    ARIADNE_TEST_SAME(ValidatedTaylorModel<F>({{{0},1.0},{{1},-2.0},{{2},3.0}},0.75,swp)*ValidatedTaylorModel<F>({{{0},3.0},{{1},2.0},{{2},-4.0}},0.5,swp),
                      ValidatedTaylorModel<F>({{{0},3.0},{{1},-4.0},{{2},1.0},{{3},14.0},{{4},-12.0}}, 10.125,swp));

    ARIADNE_TEST_SAME(1-2*x0+3*(x0^2)+r*3/4, ValidatedTaylorModel<F>({{{0,0},1.0},{{1,0},-2.0},{{2,0},3.0}}, 0.75,swp));
    ARIADNE_TEST_CONSTRUCT(ValidatedTaylorModel<F>,t,(1-2*x0+3*(x0^2)+r*3/4));
    ARIADNE_TEST_SAME(t,ValidatedTaylorModel<F>({{{0,0},1.0},{{1,0},-2.0},{{2,0},3.0}},0.75,swp));

    // Reciprocal of a constant
    ARIADNE_TEST_SAME(rec(o*4),o/4);

    ARIADNE_TEST_SAME((1-2*x0+3*x1)*(3+2*x0-4*x1),3-4*x0+5*x1-4*(x0^2)+14*(x0*x1)-12*(x1^2));
    ARIADNE_TEST_SAME(sqr(t),t*t);
    ARIADNE_TEST_SAME(pow(t,0),o);
    ARIADNE_TEST_SAME(pow(t,1),t);
    ARIADNE_TEST_SAME(pow(t,2),t*t);
    ARIADNE_TEST_SAME(pow(t,3),t*t*t);

    F inf_ = F::inf(pr);
    ValidatedTaylorModel<F> tm_inf(Expansion<MultiIndex,FloatType>(2),+inf_,swp);
    ValidatedTaylorModel<F> tm_zero_times_inf=0*tm_inf;
    if(is_nan(tm_zero_times_inf.error().raw())) {
        ARIADNE_TEST_WARN("Multiplying 0+/-inf by 0 yields 0+/-NaN");
    } else if(tm_zero_times_inf.error().raw()==+inf_) {
        ARIADNE_TEST_WARN("Multiplying 0+/-inf by 0 yields 0+/-inf");
    } else if(tm_zero_times_inf.error().raw()==0.0) {
        ARIADNE_TEST_PRINT("Multiplying 0+/-inf by 0 yields 0+/-0");
    }
}

template<class F> Void TestTaylorModel<F>::test_range()
{
    typedef Interval<FloatValue<PR>> ExactRangeType;

    // Test range of cubic, which should be exact
    ValidatedTaylorModel<F> t1 = x0*x0*x0+x0;
    ARIADNE_TEST_SAME(t1.range(),ExactRangeType(-2,+2));

    //x^2+x = (x+1/2)^2-1/4; Range [-0.25,+2.0]
    // Test range of quadratic, which could be exact, but need not be
    ValidatedTaylorModel<F> t2 = x0*x0+x0;
    ARIADNE_TEST_BINARY_PREDICATE(refines,t2.range(),ExactRangeType(-2,+2));
    ARIADNE_TEST_BINARY_PREDICATE(refines,ExactRangeType(-0.25,+2),t2.range());
    if(cast_exact_interval(t2.range())!=ExactRangeType(-0.25,2.0)) {
        ARIADNE_TEST_WARN("ValidatedTaylorModel<F>::range() not exact for quadratic functions."); }
}


template<class F> Void TestTaylorModel<F>::test_functions()
{
    ValidatedTaylorModel<F> x=ValidatedTaylorModel<F>::coordinate(1,0,swp);
    ValidatedTaylorModel<F> hx=x/2;
    ValidatedTaylorModel<F> ophx=1+x/2;
    FloatBounds<PR> e(-1,+1);
    FloatError<PR>::set_output_places(8);

    ARIADNE_TEST_PRINT(exp(x));
    ARIADNE_TEST_PRINT(sin(x));
    ARIADNE_TEST_PRINT(cos(x));
    ARIADNE_TEST_PRINT(x.tolerance());

    FloatValue<PR> LAXITY(128);
    // Threshold based on sweeper characteristics
    static const FloatValue<PR> threshold=FloatValue<PR>(x.tolerance());
    static const FloatBounds<PR> tolerance=threshold*FloatBounds<PR>(-1,+1)*LAXITY;

    // exp, sin and cos have error bound e^N/N!*(1+1/N), where e is bound for |x|, N is the first term omitted
    // ErrorTag bound for rec is e^(N-1); log is e^(N-1)/N; sqrt is ???, where e is bound for |x-1|

    //Functions based on their natural defining points with variable dependence 1.0
    ValidatedTaylorModel<F> expected_exp_x
        = 1+x+(x^2)/2+(x^3)/6+(x^4)/24+(x^5)/120+(x^6)/720+(x^7)/5040+(x^8)/40320+(x^9)/362880+(x^10)/3628800+(x^11)/39916800+(x^12)/479001600+e/5782233600;
    ValidatedTaylorModel<F> expected_sin_x
        = x-(x^3)/6+(x^5)/120-(x^7)/5040+(x^9)/362880-(x^11)/39916800+e/5782233600;
    ValidatedTaylorModel<F> expected_cos_x
        = 1-(x^2)/2+(x^4)/24-(x^6)/720+(x^8)/40320-(x^10)/3628800+(x^12)/479001600+e/81366405120;
    ARIADNE_TEST_REFINES(exp(x),expected_exp_x+tolerance);
    ARIADNE_TEST_REFINES(sin(x),expected_sin_x+tolerance);
    ARIADNE_TEST_REFINES(cos(x),expected_cos_x+tolerance);

    //Functions based on their natural defining points with variable dependence 0.5
    //ValidatedTaylorModel<F> expected_exp_hx = 1+hx+(hx^2)/2+(hx^3)/6+(hx^4)/24+(hx^5)/120+(hx^6)/720+(hx^7)/5040+e/128/4410+tol;
    ValidatedTaylorModel<F> expected_exp_hx
        = 1+x/2+(x^2)/8+(x^3)/48+(x^4)/384+(x^5)/3840+(x^6)/46080+(x^7)/645120+(x^8)/10321920+(x^9)/185794560+(x^10)/3715891200+e/2048/36590400;
    ValidatedTaylorModel<F> expected_sin_hx
        = x/2-(x^3)/48+(x^5)/3840-(x^7)/645120+(x^9)/185794560+e/2048/36590400;
    ValidatedTaylorModel<F> expected_cos_hx
        = 1-(x^2)/8+(x^4)/384-(x^6)/46080+(x^8)/10321920-(x^10)/3715891200+e/4096/5748019200*13;
    // Uncomment the last line so that expected_cos_hx contains error term base on x^8 term rather than x^7 term.
    //ValidatedTaylorModel<F> expected_cos_hx = 1-(hx^2)/2+(hx^4)/24-(hx^6)/720+e/256/35840+tol;

    ARIADNE_TEST_REFINES(exp(hx),expected_exp_hx+tolerance);
    ARIADNE_TEST_REFINES(sin(hx),expected_sin_hx+tolerance);
    ARIADNE_TEST_REFINES(cos(hx),expected_cos_hx+tolerance);

    // Series for sqrt(1+x) is
    // 1, 1/2, -1/8, 1/16, -5/128, 7/256, -21/1024, 33/2048, -273/16384, 35/2048, -357/20480
    // Series for sqrt(1+x/2) is
    // 1, 1/4, -1/32, 1/128, -5/2048, 7/8192, -21/65536, 33/262144, 429/8388608, 715/33554432, 2431/268435456
    ValidatedTaylorModel<F> expected_rec_ophx
        = 1-x/2+(x^2)/4-(x^3)/8+(x^4)/16-(x^5)/32+(x^6)/64-(x^7)/128+(x^8)/256-(x^9)/512+(x^10)/1024+e/1024;
    ValidatedTaylorModel<F> expected_rec_ophx_coarse_error
        = 1-x/2+(x^2)/4-(x^3)/8+(x^4)/16-(x^5)/32+(x^6)/64-(x^7)/128+(x^8)/256-(x^9)/512+(x^10)/1024+e/1024*2;
    ValidatedTaylorModel<F> expected_sqrt_ophx
        = 1+x/4-(x^2)/32+(x^3)/128-(x^4)*5/2048+(x^5)*7/8192-(x^6)*21/65536+(x^7)*33/262144-(x^8)*429/8388608+(x^9)*715/33554432+e*2431/268435456;
    ValidatedTaylorModel<F> expected_sqrt_ophx_coarse_error
        = 1+x/4-(x^2)/32+(x^3)/128-(x^4)*5/2048+(x^5)*7/8192-(x^6)*21/65536+(x^7)*33/262144-(x^8)*429/8388608+(x^9)*715/33554432+e*2431/268435456*2;
    ValidatedTaylorModel<F> expected_log_ophx
        = x/2-(x^2)/4/2+(x^3)/8/3-(x^4)/16/4+(x^5)/32/5-(x^6)/64/6+(x^7)/128/7-(x^8)/256/8+(x^9)/512/9-(x^10)/1024/10+e/1024/10;
    ValidatedTaylorModel<F> expected_log_ophx_coarse_error
        = x/2-(x^2)/4/2+(x^3)/8/3-(x^4)/16/4+(x^5)/32/5-(x^6)/64/6+(x^7)/128/7-(x^8)/256/8+(x^9)/512/9-(x^10)/1024/10+e/1024/10;
    ARIADNE_TEST_REFINES(rec(ophx),expected_rec_ophx_coarse_error+tolerance);
    ARIADNE_TEST_REFINES(sqrt(ophx),expected_sqrt_ophx_coarse_error+tolerance);
    ARIADNE_TEST_REFINES(log(ophx),expected_log_ophx_coarse_error+tolerance);

    ARIADNE_TEST_REFINES(sqrt(ophx)*sqrt(ophx),ophx+tolerance*2);

    // Doubling formulae
    ARIADNE_TEST_REFINES(exp(2*x),(expected_exp_x^2)+tolerance);
    ARIADNE_TEST_REFINES(sin(2*x),2*expected_sin_x*expected_cos_x+tolerance);
    ARIADNE_TEST_REFINES(cos(2*x),2*(expected_cos_x^2)-1+tolerance);

    ARIADNE_TEST_REFINES(3*rec(3*ophx),expected_rec_ophx_coarse_error+tolerance);
    ARIADNE_TEST_REFINES(sqrt(9*ophx)/3,expected_sqrt_ophx_coarse_error+tolerance);
    ARIADNE_TEST_REFINES(log(3*ophx),expected_log_ophx_coarse_error+log(FloatValue<PR>(3))+tolerance);

    Nat rn=3; FloatValue<PR> c(2); FloatValue<PR> frn(rn);
    ARIADNE_TEST_PRINT(c);
    ARIADNE_TEST_PRINT(frn);
    ARIADNE_TEST_REFINES(exp(hx),expected_exp_hx+tolerance);
    ARIADNE_TEST_REFINES(exp(c+frn*hx),exp(c)*pow(expected_exp_hx,rn)+tolerance);
    ARIADNE_TEST_REFINES(sin(c+hx),sin(c)*expected_cos_hx+cos(c)*expected_sin_hx+tolerance);
    ARIADNE_TEST_REFINES(cos(c+hx),cos(c)*expected_cos_hx-sin(c)*expected_sin_hx+tolerance);

    c=FloatValue<PR>(10);
    ARIADNE_TEST_REFINES(sin(c+hx),sin(c)*expected_cos_hx+cos(c)*expected_sin_hx+tolerance);
    ARIADNE_TEST_REFINES(cos(c+hx),cos(c)*expected_cos_hx-sin(c)*expected_sin_hx+tolerance);

    // Test exponential based at log2; exp(log(2)+x/2)=2*exp(x/2)
    FloatValue<PR> log2_apprx(0.693147);
    ARIADNE_TEST_REFINES(exp(log2_apprx+x/2),
                                  2*(1+hx+(hx^2)/2+(hx^3)/6+(hx^4)/24+(hx^5)/120+(hx^6)/720)+e*6/100000);
}


template<class F> Void TestTaylorModel<F>::test_rescale()
{
}

template<class F> Void TestTaylorModel<F>::test_restrict()
{
}

template<class F> Void TestTaylorModel<F>::test_refinement()
{
    ValidatedTaylorModel<F> x=ValidatedTaylorModel<F>::coordinate(2,0,swp);
    ValidatedTaylorModel<F> y=ValidatedTaylorModel<F>::coordinate(2,1,swp);
    ValidatedTaylorModel<F> e=ValidatedTaylorModel<F>::unit_ball(2,swp);

    ARIADNE_TEST_REFINES(1+x+(x^2)/2+(x^3)/4,1+x+(x^2)/2+e/4);
    ARIADNE_TEST_BINARY_PREDICATE(not refines,1+x+(x^2)/2+(x^3)/6+e/pow(2.0_exact,31),1+x+(x^2)/2+e/6);

    // Test refinement with no roundoff errors
    ARIADNE_TEST_SAME(refinement(1-x*3/4+(x^3)*3+(x^4)*13/4+e/2,1+(x^2)/4+(x^3)*2+(x^4)*3+e),1-x*5/8+(x^3)*11/4+(x^4)*13/4+e/2);

    // Test refinement with roundoff errors
    ARIADNE_TEST_SAME(refinement(FloatValue<PR>(2.0/3)*x+e/2,FloatValue<PR>(6.0/5)*x+e/4),
        FloatValue<PR>(1.05833333333333335)*x+FloatValue<PR>(0.108333333333333393)*e);

    // Code below computes expected values for second test
    // FloatDP xv=2./3; FloatDP xe=1./2; FloatDP yv=6./5; FloatDP ye=1./4;
    // FloatDP rl=sub(down,yv,ye); FloatDP ru=add(up,xv,xe); FloatDP rv=add(near,rl,ru)/2; FloatDP re=sub(up,ru,rl)/2;
    // std::cerr << std::setprecision(18) << "xv="<<xv<<" yv="<<yv<<" rl="<<rl<<" ru="<<ru<<" rv="<<rv<<" re="<<re<<"\n";


}

template<class F> Void TestTaylorModel<F>::test_split()
{
    ValidatedTaylorModel<F> x=ValidatedTaylorModel<F>::coordinate(2,0,swp);
    ValidatedTaylorModel<F> y=ValidatedTaylorModel<F>::coordinate(2,1,swp);
    ValidatedTaylorModel<F> z=ValidatedTaylorModel<F>::zero(2,swp);
    ValidatedTaylorModel<F> t=1+3*x+2*y-5*x*x-7*x*y+11*y*y;
    ValidatedTaylorModel<F> es1=-1.75+4*x+5.5*y-1.25*x*x-3.5*x*y+11*y*y;
    ValidatedTaylorModel<F> es2=1+1.5*x+2*y-1.25*x*x-3.5*x*y+11*y*y;
    ValidatedTaylorModel<F> es3=1.25-1*x-1.5*y-1.25*x*x-3.5*x*y+11*y*y;

    ARIADNE_TEST_PRINT(t);
    ARIADNE_TEST_SAME(split(t,0,SplitPart::LOWER),es1);
    ARIADNE_TEST_SAME(split(t,0,SplitPart::MIDDLE),es2);
    ARIADNE_TEST_SAME(split(t,0,SplitPart::UPPER),es3);
}


template<class F> Void TestTaylorModel<F>::test_antiderivative()
{
    FloatBounds<PR> unit_interval(-1,+1);
    ValidatedTaylorModel<F> tm=ValidatedTaylorModel<F>::constant(2,1,swp);
    ValidatedTaylorModel<F> atm=antiderivative(tm,1);

    ARIADNE_TEST_SAME(antiderivative(2*o,0),2*x0);
    ARIADNE_TEST_SAME(antiderivative(2*o,1),2*x1);
    ARIADNE_TEST_SAME(antiderivative(3*x0,0),(x0^2)*3/2);
    ARIADNE_TEST_SAME(antiderivative(2*o+3*x0,0),x0*2+(x0^2)*3/2);
    ARIADNE_TEST_SAME(antiderivative((x0^2)*15/2,0),(x0^3)*5/2);
    ARIADNE_TEST_SAME(antiderivative((x0^2)*(x1^4)*15/2,0),(x0^3)*(x1^4)*5/2);
    ARIADNE_TEST_SAME(antiderivative((x0^2)*(x1^4)*15/2,1),(x0^2)*(x1^5)*3/2);

    if constexpr (IsSame<F,FloatDP>::value) {
        ValidatedTaylorModel<F> x=ValidatedTaylorModel<F>::coordinate(1,0,swp);
        ValidatedTaylorModel<F> e=ValidatedTaylorModel<F>::zero(1,swp)+FloatBounds<PR>(-1,+1);
        ARIADNE_TEST_SAME(antiderivative(2*x*x,0),0.66666666666666663*x*x*x+5.5511151231257827021e-17*e);
        ARIADNE_TEST_SAME(antiderivative(2*x*x+e,0),0.66666666666666663*x*x*x+1.0000000000000002*e);
        ARIADNE_TEST_SAME(antiderivative(2*(x^2),0),FloatValue<PR>(0.66666666666666663)*(x^3)+FloatValue<PR>(5.5511151231257827021e-17)*e);
        ARIADNE_TEST_SAME(antiderivative(2*(x^2)+e,0),FloatValue<PR>(0.66666666666666663)*(x^3)+FloatValue<PR>(1.0000000000000002)*e);

        // Regression test
        ValidatedTaylorModel<F> t1({ {{0,0},1.}, {{1,0},2.}, {{0,1},3.}, {{2,0},4.}, {{1,1},5.}, {{0,2},6.} }, 0., swp);
        ValidatedTaylorModel<F> at1({ {{1,0},1.}, {{2,0},1.}, {{1,1},3.}, {{3,0},1.33333333333333333}, {{2,1},2.5}, {{1,2},6.} }, 1.1102230246251565404e-16, swp);
        ARIADNE_TEST_SAME(antiderivative(t1,0),at1);
    }
}


template<class F> Void TestTaylorModel<F>::test_compose()
{
    ARIADNE_TEST_SAME(compose(2-x0*x0-x1/4,{2-x0*x0-x1/4,x0}),-2-(x0^4)-(x1^2)/16+4*(x0^2)+x1-(x0^2)*x1/2-x0/4);

}

template<class F> Void TestTaylorModel<F>::test_recondition()
{
    ARIADNE_TEST_CONSTRUCT( ValidatedTaylorModel<F>, tm1, ({ {{0,0},2.0}, {{1,0},3.0}, {{0,1},5.0} }, 0.5, swp) );
    ARIADNE_TEST_CONSTRUCT( ValidatedTaylorModel<F>, tm2, ({ {{0,0,1},0.5}, {{0,0,0},2.0}, {{1,0,0},3.0}, {{0,1,0},5.0} }, 0.0, swp) );
    ARIADNE_TEST_CONSTRUCT( ValidatedTaylorModel<F>, tm3, ({ {{0},2.0}, {{1},5.0} }, 3.5, swp) );
    ARIADNE_TEST_SAME(embed_error(tm1),tm2);
    ARIADNE_TEST_SAME(discard_variables(tm2,{2}),tm1);
    ARIADNE_TEST_SAME(discard_variables(tm1,{0}),tm3);
}




Int main() {
    ThresholdSweeper<FloatDP> sweeper_dp(dp,1e-8);
    TestTaylorModel<FloatDP>(sweeper_dp).test();

    MultiplePrecision mp(128);
    ThresholdSweeper<FloatMP> sweeper_mp(mp,std::pow(2.0,-64));
    ARIADNE_TEST_PRINT(sweeper_mp.precision());
    ARIADNE_TEST_PRINT(sweeper_mp);
    TestTaylorModel<FloatMP>(sweeper_mp).test();

    RelativeThresholdSweeper<FloatMP> relative_sweeper_mp(mp,std::pow(2.0,-64));
    TestTaylorModel<FloatMP>(relative_sweeper_mp).test();

    return ARIADNE_TEST_FAILURES;
}



