/***************************************************************************
 *            test_polynomial.cc
 *
 *  Copyright 2009--17  Pieter Collins
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
#include "config.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/multi_index.h"
#include "algebra/algebra.h"
#include "function/polynomial.h"

#include "test.h"
using namespace std;
using namespace Ariadne;


class TestPolynomial
{
    typedef MultiIndex MI;
    typedef Polynomial<Float64> P;
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
    Float64 x=0;
    Vector<Float64> v(3);
    MultiIndex a(3);
    Polynomial<Float64> p(3);
    const Polynomial<Float64> cp(3);
    Vector< Polynomial<Float64> > pv(2);

    p=Polynomial<Float64>();
    p=Polynomial<Float64>(3);
    p=Polynomial<Float64>(cp);

    p=Polynomial<Float64>({ {{0,0,0},1}, {{1,0,0},2}, {{0,0,0},3}, {{0,0,1},5.0} });

    //p=Polynomial<Float64>::variable(3u,0u);
    //p=Polynomial<Float64>::variables(3u)[0u];

    p=x;

    p.reserve(2u);
    p.insert(a,x);

    x=cp[a];
    p[a]=1.0;

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

        typedef Expansion<Float64>::Iterator Iterator;
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
    Polynomial<Float64> p(3);
    for(Nat i=0; i!=2; ++i) {
        if(i%2) { p.expansion().append(a,1/(1.+i)); ++b; ++b; a=b; ++b; } else { p.expansion().append(b,1/(1.+i));}
    }
    ARIADNE_TEST_PRINT(p);
    ARIADNE_TEST_EXECUTE(p.cleanup());
    ARIADNE_TEST_PRINT(p);

}

Void TestPolynomial::test_constructors()
{
    // Empty polynomial
    ARIADNE_TEST_CONSTRUCT(Polynomial<Float64>,p1,(3));
    // Dense polynomial
    ARIADNE_TEST_CONSTRUCT(Polynomial<Float64>,p2,({ {{0,0,0},0.}, {{1,0,0},0.},{{0,1,0},0.},{{0,0,1},0.}, {{2,0,0},5.},{{1,1,0},2.},{{1,0,1},0.},{{0,2,0},0.},{{0,1,2},3.},{{0,0,2},0.} }));
    ARIADNE_TEST_EQUAL(p2[MultiIndex({2,0,0})],5.0);
    // Sparse polynomial with unordered indiced
    ARIADNE_TEST_CONSTRUCT(Polynomial<Float64>,p3,({ {{1,2},5.0}, {{0,0},2.0}, {{1,0},3.0}, {{3,0},7.0}, {{0,1},11.0} }));
    ARIADNE_TEST_EQUAL(p3[MultiIndex({1,2})],5.0);
    ARIADNE_TEST_EQUAL(p3[MultiIndex({0,0})],2.0);

    // Unordered indices
    ARIADNE_TEST_EQUAL(P({ {{1,2},5.0}, {{0,0},2.0}, {{1,0},3.0}, {{3,0},7.0}, {{0,1},11.0} }),P({ {{0,0},2.0}, {{1,0},3.0}, {{0,1},11.0}, {{3,0},7.0}, {{1,2},5.0} }));
    // Repeated indices
    ARIADNE_TEST_EQUAL(P({ {{1,2},5.0}, {{0,0},2.0}, {{1,0},3.0}, {{1,0},7.0}, {{1,2},11.0} }),P({ {{0,0},2.0}, {{1,0},10.0}, {{1,2},16.0} }));

}

Void TestPolynomial::test_indexing()
{
    Polynomial<Float64> p({ {{0,0,0},2.0},  {{1,0,0},3.0}, {{1,0,1},5.0}, {{2,1,0},7.0} });
    const Polynomial<Float64>& pc=p;
    ARIADNE_TEST_EQUAL(p[MultiIndex({1,0,0})],3.0);

    p[MultiIndex({1,0,0})]-=0.5;
    ARIADNE_TEST_EQUAL(p[MultiIndex({1,0,0})],2.5);

    p[MultiIndex({1,1,0})]=11.0;
    ARIADNE_TEST_EQUAL(p[MultiIndex({1,1,0})],11.0);

    Polynomial<Float64> q(3);
    q[MultiIndex({0,0,0})]=2.0;
    q[MultiIndex({0,1,0})]=3.0;
    ARIADNE_TEST_EQUALS(q.number_of_terms(),2);
    ARIADNE_TEST_EQUALS(q[MultiIndex({0,0,0})],2.0);
    ARIADNE_TEST_EQUALS(q[MultiIndex({0,1,0})],3.0);
    q[MultiIndex({1,0,0})]=5.0;
    ARIADNE_TEST_EQUALS(q[MultiIndex({1,0,0})],5.0);
    ARIADNE_TEST_EQUALS(q[MultiIndex({0,1,0})],3.0);
    q[MultiIndex({0,0,1})]=7.0;
    ARIADNE_TEST_EQUALS(q[MultiIndex({0,0,1})],7.0);

    // Test insert at beginning
    p.clear();
    ARIADNE_TEST_PRINT(p.expansion());
    p[MultiIndex({0,1,0})]=2.0;
    p[MultiIndex({0,0,1})]=3.0;
    p[MultiIndex({2,1,0})]=5.0;
    ARIADNE_TEST_PRINT(p.expansion());
    ARIADNE_TEST_EQUALS(pc[MultiIndex({0,0,0})],0.0);
    ARIADNE_TEST_EQUALS(p.number_of_terms(),3);
    ARIADNE_TEST_PRINT(p.expansion());
    ARIADNE_TEST_EXECUTE(p[MultiIndex({0,0,0})]=7);
    ARIADNE_TEST_PRINT(p.expansion());
    ARIADNE_TEST_EQUALS(p.number_of_terms(),4);
    p.expansion().graded_sort();
    ARIADNE_TEST_PRINT(p.expansion());

    Polynomial<Float64>::ConstIterator iter=p.begin();
    ARIADNE_TEST_EQUALS(iter->key(),MultiIndex({0,0,0}));
    ARIADNE_TEST_EQUALS(iter->data(),7.0);
    ++iter;
    ARIADNE_TEST_EQUALS(iter->key(),MultiIndex({0,1,0}));
    ARIADNE_TEST_EQUALS(iter->data(),2.0);
    ++iter;
    ARIADNE_TEST_EQUALS(iter->key(),MultiIndex({0,0,1}));
    ARIADNE_TEST_EQUALS(iter->data(),3.0);
    ++iter;
    ARIADNE_TEST_EQUALS(iter->key(),MultiIndex({2,1,0}));
    ARIADNE_TEST_EQUALS(iter->data(),5.0);
}

Void TestPolynomial::test_arithmetic()
{
    typedef Polynomial<Float64> P;
    ARIADNE_TEST_EQUAL(P(3)+P(3),P(3));
    ARIADNE_TEST_EQUAL(P(3)+P({ {{2,1,0},2.0} }),P({ {{2,1,0},2.0} }));
    ARIADNE_TEST_EQUAL(P(3)+P({ {{2,1,0},2.0}, {{0,1,0},3.0}, {{1,1,0},5.0} }), P({ {{0,1,0},3.0}, {{1,1,0},5.0}, {{2,1,0},2.0} }));

    Polynomial<Float64> x0(3); x0[MultiIndex({1,0,0})]=1.0;
    Polynomial<Float64> x1(3); x1[MultiIndex({0,1,0})]=1.0;
}

Void TestPolynomial::test_variables()
{
    Vector< Polynomial<Float64> > x=Polynomial<Float64>::variables(3);
    Array< Vector<Float64> > e=Vector<Float64>::basis(2);

    Polynomial<Float64> p1=x[1]*3.0;
    Polynomial<Float64> p2=p1+x[0]; p2=x[1]*3,0+x[0];
    Polynomial<Float64> p3=x[0]*p2; p3=x[0]*(x[1]*3.0+x[0]);
    Polynomial<Float64> p4=x[1]*x[2];
    Polynomial<Float64> p5=p3+p4;
    Polynomial<Float64> p=x[0]*(x[1]*3.0+x[0])+x[1]*x[2];

    ARIADNE_TEST_EQUAL(x[0], Polynomial<Float64>({ {{1,0,0},1.0} }));
    ARIADNE_TEST_EQUAL(x[1], Polynomial<Float64>({ {{0,1,0},1.0} }));
    ARIADNE_TEST_EQUAL(x[2], Polynomial<Float64>({ {{0,0,1},1.0} }));
    ARIADNE_TEST_EQUAL(x[0]+x[1], Polynomial<Float64>({ {{1,0,0},1.0}, {{0,1,0},1.0} }));
    ARIADNE_TEST_EQUAL(x[0]*x[1], Polynomial<Float64>({ {{1,1,0},1.0} }));
    ARIADNE_TEST_EVALUATE(x[0]*(x[1]*3.0+x[0])+x[1]*x[2]);
    ARIADNE_TEST_EQUAL((x[0]*(x[1]*3.0+x[0])+x[1]*x[2]), Polynomial<Float64>({ {{1,1,0},3.0}, {{2,0,0},1.0}, {{0,1,1},1.0} }));
    ARIADNE_TEST_EQUAL((e[1]*(x[0]*(x[1]*3.0+x[0])+x[1]*x[2]))[1], Polynomial<Float64>({ {{1,1,0},3.0}, {{2,0,0},1.0}, {{0,1,1},1.0} }));
    ARIADNE_TEST_PRINT((e[1]*(x[0]*(x[1]*3.0+x[0])+x[1]*x[2]))[0]);
    ARIADNE_TEST_EQUAL((e[1]*(x[0]*(x[1]*3.0+x[0])+x[1]*x[2]))[0], Polynomial<Float64>(3));
    ARIADNE_TEST_EQUAL((e[1]*(x[0]*(x[1]*3.0+x[0])+x[1]*x[2]))[0], Polynomial<Float64>({ {{3,0,0},0.0} }));

}

Void TestPolynomial::test_find()
{
    Polynomial<Float64> p({ {{1,2},5.0}, {{0,0},2.0}, {{1,0},3.0}, {{3,0},7.0}, {{0,1},11.0} });
    MultiIndex a(2);
    a[0]=1; a[1]=2;
    ARIADNE_TEST_PRINT(p);
    ARIADNE_TEST_PRINT(p.find(a)-p.begin());
    ARIADNE_TEST_COMPARE(p.find(a),!=,p.end());
    ARIADNE_TEST_EQUAL(p.find(a)->key(),a);
    ARIADNE_TEST_EQUAL(p.find(a)->data(),5.0);
    a[1]=1;
    ARIADNE_TEST_EQUAL(p.find(a),p.end());
}

Int main() {
    TestPolynomial().test();
    return ARIADNE_TEST_FAILURES;
}
