/***************************************************************************
 *            test_polynomial.cpp
 *
 *  Copyright  2009-20  Pieter Collins
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
#include "config.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/expansion.hpp"
#include "algebra/expansion.inl.hpp"
#include "algebra/algebra.hpp"
#include "function/polynomial.hpp"

#include "../test.hpp"

using namespace std;
using namespace Ariadne;


class TestPolynomial
{
    typedef MultiIndex MI;
    typedef MultivariatePolynomial<FloatDP> P;
  public:
    Void test();
  private:
    Void test_concept();
    Void test_cleanup();
    Void test_constructors();
    Void test_indexing();
    Void test_arithmetic();
    Void test_variables();
    Void test_find();
};


Void TestPolynomial::test()
{
    ARIADNE_TEST_CALL(test_cleanup());
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_indexing());
    ARIADNE_TEST_CALL(test_variables());
    ARIADNE_TEST_CALL(test_find());
}


Void TestPolynomial::test_concept()
{
    FloatDP x(0,dp);
    Vector<FloatDP> v(3,dp);
    MultiIndex a(3);
    MultivariatePolynomial<FloatDP> p(3,dp);
    const MultivariatePolynomial<FloatDP> cp(3,dp);
    Vector< MultivariatePolynomial<FloatDP> > pv(2,cp);

    p=MultivariatePolynomial<FloatDP>(3,dp);
    p=MultivariatePolynomial<FloatDP>(cp);

    p=MultivariatePolynomial<FloatDP>({ {{0,0,0},1.0_x}, {{1,0,0},2.0_x}, {{0,0,0},3.0_x}, {{0,0,1},5.0_x} },dp);

    //p=MultivariatePolynomial<FloatDP>::variable(3u,0u);
    //p=MultivariatePolynomial<FloatDP>::variables(3u)[0u];

    p=x;

    p.reserve(2u);
    p.insert(a,x);

    x=cp[a];
    p[a]=1.0_x;

    cp.argument_size();

    p.erase(p.begin());

    p.check();
    p.cleanup();
    p.clear();

    evaluate(p,v);
    compose(p,pv);

}

Void TestPolynomial::test_cleanup()
{
/*
    {
        MultiIndex a(3);
        MultiIndex b(3); ++b;
        std::vector<ValueType> v;

        for(Nat i=0; i!=20; ++i) {
            if(i%2) { v.push_back(ValueType(a,1/(1.+i)));  ++b; ++b; a=b; ++b; } else { v.push_back(ValueType(b,1/(1.+i)));}
        }
        std::sort(v.begin(),v.end());
    }
    {
        std::vector<Int> v;
        const Int element_size=3;
        const Int index_size=1;
        const Int data_size=2;
        for(Nat i=0; i!=20; ++i) {
            Int word= ( i%2 ? (3*i-2)/2 : (3*i+2)/2 );
            double value=1/(2.+i);
            v.resize(v.size()+3);
            v[v.size()-3]=word;
            reinterpret_cast<double&>(v[v.size()-2])=value;
        }
        for(Nat i=0; i!=v.size()/3; ++i) {
            Int word; double value;
            word=v[3*i];
            value=reinterpret_cast<double&>(v[3*i+1]);
        }

        typedef Expansion<FloatDP>::Iterator Iterator;
        Iterator iter1(3,&*v.begin());
        Iterator iter2(3,&*v.end());
        std::sort(iter1,iter2);

        for(Nat i=0; i!=v.size()/3; ++i) {
            Int word; double value;
            word=v[3*i];
            value=reinterpret_cast<double&>(v[3*i+1]);
        }

    }
*/

    // Test to see if the cleanup/sort operations work.
    // Since these are used in the constructors, we can't use the main constructors to test this
    MultiIndex a(3);
    MultiIndex b(3); ++b;
    MultivariatePolynomial<FloatDP> p(3,dp);
    FloatDP c(1,dp);
    for(Nat i=0; i!=2; ++i) {
        if(i%2) { p.expansion().append(a,c); ++b; ++b; a=b; ++b; } else { p.expansion().append(b,c);}
        c=hlf(c);
    }
    ARIADNE_TEST_PRINT(p);
    ARIADNE_TEST_EXECUTE(p.cleanup());
    ARIADNE_TEST_PRINT(p);

}

Void TestPolynomial::test_constructors()
{
    // Empty polynomial
    ARIADNE_TEST_CONSTRUCT(UnivariatePolynomial<FloatDP>,q1,(SizeOne(),dp));
    // Dense polynomial
    ARIADNE_TEST_CONSTRUCT(UnivariatePolynomial<FloatDP>,q2,({ {0,0.0_x}, {1,0.0_x},{2,5.0_x},{3,2.0_x} },dp));
    ARIADNE_TEST_EQUAL(q2[1],0.0_x);
    ARIADNE_TEST_EQUAL(q2[2],5.0_x);
    ARIADNE_TEST_EQUAL(q2[3],2.0_x);

    // Empty polynomial
    ARIADNE_TEST_CONSTRUCT(MultivariatePolynomial<FloatDP>,p1,(3,dp));
    // Dense polynomial
    ARIADNE_TEST_CONSTRUCT(MultivariatePolynomial<FloatDP>,p2,({ {{0,0,0},0.0_x}, {{1,0,0},0.0_x},{{0,1,0},0.0_x},{{0,0,1},0.0_x}, {{2,0,0},5.0_x},{{1,1,0},2.0_x},{{1,0,1},0.0_x},{{0,2,0},0.0_x},{{0,1,2},3.0_x},{{0,0,2},0.0_x} },dp));
    ARIADNE_TEST_EQUAL(p2[MultiIndex({2,0,0})],5.0_x);
    // Sparse polynomial with unordered indiced
    ARIADNE_TEST_CONSTRUCT(MultivariatePolynomial<FloatDP>,p3,({ {{1,2},5.0_x}, {{0,0},2.0_x}, {{1,0},3.0_x}, {{3,0},7.0_x}, {{0,1},11.0_x} },dp));
    ARIADNE_TEST_EQUAL(p3[MultiIndex({1,2})],5.0_x);
    ARIADNE_TEST_EQUAL(p3[MultiIndex({0,0})],2.0_x);

    // Unordered indices
    ARIADNE_TEST_EQUAL(P({ {{1,2},5.0_x}, {{0,0},2.0_x}, {{1,0},3.0_x}, {{3,0},7.0_x}, {{0,1},11.0_x} },dp),P({ {{0,0},2.0_x}, {{1,0},3.0_x}, {{0,1},11.0_x}, {{3,0},7.0_x}, {{1,2},5.0_x} },dp));
    // Repeated indices
    ARIADNE_TEST_EQUAL(P({ {{1,2},5.0_x}, {{0,0},2.0_x}, {{1,0},3.0_x}, {{1,0},7.0_x}, {{1,2},11.0_x} },dp),P({ {{0,0},2.0_x}, {{1,0},10.0_x}, {{1,2},16.0_x} },dp));

}

Void TestPolynomial::test_indexing()
{
    MultivariatePolynomial<FloatDP> p({ {{0,0,0},2.0_x},  {{1,0,0},3.0_x}, {{1,0,1},5.0_x}, {{2,1,0},7.0_x} },dp);
    const MultivariatePolynomial<FloatDP>& pc=p;
    ARIADNE_TEST_EQUAL(p[MultiIndex({1,0,0})],3.0_x);

    p[MultiIndex({1,0,0})]-=0.5_x;
    ARIADNE_TEST_EQUAL(p[MultiIndex({1,0,0})],2.5_x);

    p[MultiIndex({1,1,0})]=11.0_x;
    ARIADNE_TEST_EQUAL(p[MultiIndex({1,1,0})],11.0_x);

    MultivariatePolynomial<FloatDP> q(3,dp);
    q[MultiIndex({0,0,0})]=2.0_x;
    q[MultiIndex({0,1,0})]=3.0_x;
    ARIADNE_TEST_EQUALS(q.number_of_terms(),2);
    ARIADNE_TEST_EQUALS(q[MultiIndex({0,0,0})],2.0_x);
    ARIADNE_TEST_EQUALS(q[MultiIndex({0,1,0})],3.0_x);
    q[MultiIndex({1,0,0})]=5.0_x;
    ARIADNE_TEST_EQUALS(q[MultiIndex({1,0,0})],5.0_x);
    ARIADNE_TEST_EQUALS(q[MultiIndex({0,1,0})],3.0_x);
    q[MultiIndex({0,0,1})]=7.0_x;
    ARIADNE_TEST_EQUALS(q[MultiIndex({0,0,1})],7.0_x);

    // Test insert at beginning
    p.clear();
    ARIADNE_TEST_PRINT(p.expansion());
    p[MultiIndex({0,1,0})]=2.0_x;
    p[MultiIndex({0,0,1})]=3.0_x;
    p[MultiIndex({2,1,0})]=5.0_x;
    ARIADNE_TEST_PRINT(p.expansion());
    ARIADNE_TEST_EQUALS(pc[MultiIndex({0,0,0})],0.0_x);
    ARIADNE_TEST_EQUALS(p.number_of_terms(),3);
    ARIADNE_TEST_PRINT(p.expansion());
    ARIADNE_TEST_EXECUTE(p[MultiIndex({0,0,0})]=7);
    ARIADNE_TEST_PRINT(p.expansion());
    ARIADNE_TEST_EQUALS(p.number_of_terms(),4);
    p.expansion().graded_sort();
    ARIADNE_TEST_PRINT(p.expansion());

    MultivariatePolynomial<FloatDP>::ConstIterator iter=p.begin();
    ARIADNE_TEST_EQUALS(iter->index(),MultiIndex({0,0,0}));
    ARIADNE_TEST_EQUALS(iter->coefficient(),7.0_x);
    ++iter;
    ARIADNE_TEST_EQUALS(iter->index(),MultiIndex({0,1,0}));
    ARIADNE_TEST_EQUALS(iter->coefficient(),2.0_x);
    ++iter;
    ARIADNE_TEST_EQUALS(iter->index(),MultiIndex({0,0,1}));
    ARIADNE_TEST_EQUALS(iter->coefficient(),3.0_x);
    ++iter;
    ARIADNE_TEST_EQUALS(iter->index(),MultiIndex({2,1,0}));
    ARIADNE_TEST_EQUALS(iter->coefficient(),5.0_x);
}

Void TestPolynomial::test_arithmetic()
{
    ARIADNE_TEST_EQUAL(P(3,dp)+P(3,dp),P(3,dp));
    ARIADNE_TEST_EQUAL(P(3,dp)+P({ {{2,1,0},2.0_x} },dp),P({ {{2,1,0},2.0_x} },dp));
    ARIADNE_TEST_EQUAL(P(3,dp)+P({ {{2,1,0},2.0_x}, {{0,1,0},3.0_x}, {{1,1,0},5.0_x} },dp), P({ {{0,1,0},3.0_x}, {{1,1,0},5.0_x}, {{2,1,0},2.0_x} },dp));

    MultivariatePolynomial<FloatDP> x0(3,dp); x0[MultiIndex({1,0,0})]=1.0_x;
    MultivariatePolynomial<FloatDP> x1(3,dp); x1[MultiIndex({0,1,0})]=1.0_x;
    MultivariatePolynomial<FloatDP> x2=MultivariatePolynomial<FloatDP>::coordinate(3,2,dp);
    UnivariatePolynomial<FloatDP> y=UnivariatePolynomial<FloatDP>::coordinate(SizeOne(),IndexZero(),dp);
    y=UnivariatePolynomial<FloatDP>::coordinate(dp);

    FloatDP w(3,dp);
    Vector<FloatDP> v({3,5,2},dp);

    ARIADNE_TEST_EQUALS(evaluate(2*x0*x0-1,v),2*v[0]*v[0]-1);

    ARIADNE_TEST_EQUALS(evaluate(2*y*y-1,w),2*w*w-1);
    ARIADNE_TEST_EQUALS(evaluate(8*y*y*(y*y-1)+1,w),8*w*w*(w*w-1)+1);

    ARIADNE_TEST_EQUALS(compose(2*y*y-1,x0),2*x0*x0-1);
}

Void TestPolynomial::test_variables()
{
    Vector< MultivariatePolynomial<FloatDP> > x=MultivariatePolynomial<FloatDP>::variables(3,dp);
    Array< Vector<FloatDP> > e=Vector<FloatDP>::basis(2,dp);

    MultivariatePolynomial<FloatDP> p1=x[1]*3.0_x;
    MultivariatePolynomial<FloatDP> p2=p1+x[0]; p2=x[1]*3; p2=x[0]+0;
    MultivariatePolynomial<FloatDP> p3=x[0]*p2; p3=x[0]*(x[1]*3.0_x+x[0]);
    MultivariatePolynomial<FloatDP> p4=x[1]*x[2];
    MultivariatePolynomial<FloatDP> p5=p3+p4;
    MultivariatePolynomial<FloatDP> p=x[0]*(x[1]*3.0_x+x[0])+x[1]*x[2];

    ARIADNE_TEST_EQUAL(x[0], MultivariatePolynomial<FloatDP>({ {{1,0,0},1.0_x} },dp));
    ARIADNE_TEST_EQUAL(x[1], MultivariatePolynomial<FloatDP>({ {{0,1,0},1.0_x} },dp));
    ARIADNE_TEST_EQUAL(x[2], MultivariatePolynomial<FloatDP>({ {{0,0,1},1.0_x} },dp));
    ARIADNE_TEST_EQUAL(x[0]+x[1], MultivariatePolynomial<FloatDP>({ {{1,0,0},1.0_x}, {{0,1,0},1.0_x} },dp));
    ARIADNE_TEST_EQUAL(x[0]*x[1], MultivariatePolynomial<FloatDP>({ {{1,1,0},1.0_x} },dp));
    ARIADNE_TEST_EVALUATE(x[0]*(x[1]*3.0_x+x[0])+x[1]*x[2]);
    ARIADNE_TEST_EQUAL((x[0]*(x[1]*3.0_x+x[0])+x[1]*x[2]), MultivariatePolynomial<FloatDP>({ {{1,1,0},3.0_x}, {{2,0,0},1.0_x}, {{0,1,1},1.0_x} },dp));
    ARIADNE_TEST_EQUAL((e[1]*(x[0]*(x[1]*3.0_x+x[0])+x[1]*x[2]))[1], MultivariatePolynomial<FloatDP>({ {{1,1,0},3.0_x}, {{2,0,0},1.0_x}, {{0,1,1},1.0_x} },dp));
    ARIADNE_TEST_PRINT((e[1]*(x[0]*(x[1]*3.0_x+x[0])+x[1]*x[2]))[0]);
    ARIADNE_TEST_EQUAL((e[1]*(x[0]*(x[1]*3.0_x+x[0])+x[1]*x[2]))[0], MultivariatePolynomial<FloatDP>(3,dp));
    ARIADNE_TEST_EQUAL((e[1]*(x[0]*(x[1]*3.0_x+x[0])+x[1]*x[2]))[0], MultivariatePolynomial<FloatDP>({ {{3,0,0},0.0_x} },dp));

}

Void TestPolynomial::test_find()
{
    MultivariatePolynomial<FloatDP> p({ {{1,2},5.0_x}, {{0,0},2.0_x}, {{1,0},3.0_x}, {{3,0},7.0_x}, {{0,1},11.0_x} },dp);
    MultiIndex a(2);
    a[0]=1; a[1]=2;
    ARIADNE_TEST_PRINT(p);
    ARIADNE_TEST_PRINT(p.find(a)-p.begin());
    ARIADNE_TEST_COMPARE(p.find(a),!=,p.end());
    ARIADNE_TEST_EQUAL(p.find(a)->index(),a);
    ARIADNE_TEST_EQUAL(p.find(a)->coefficient(),5.0_x);
    a[1]=1;
    ARIADNE_TEST_EQUAL(p.find(a),p.end());
}

Int main() {
    TestPolynomial().test();
    return ARIADNE_TEST_FAILURES;
}
