/***************************************************************************
 *            test_graded.cpp
 *
 *  Copyright  2018-20  Pieter Collins
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

#include "numeric/numeric.hpp"
#include "algebra/graded.hpp"

#include "../test.hpp"

namespace Ariadne {
template<class T> Bool same(List<T> const& lst1, List<T> const& lst2) {
    if(lst1.size()!=lst2.size()) { return false; }
    for(SizeType i=0; i!=lst1.size(); ++i) {
        if(!same(lst1[i],lst2[i])) { return false; }
    }
    return true;
    auto iter1=lst1.begin();
    auto iter2=lst2.begin();
    for(; iter1!=lst1.end(); ++iter1, ++iter2) {
        if(not same(*iter1,*iter2)) { return false; }
    }
    return true;
}
} // namespace Ariadne

using namespace Ariadne;

template<class X>
class TestGraded
{
    static Graded<X> variable(X const& x0, DegreeType deg) {
        assert(deg>=1); X z=nul(x0); Graded<X> r(x0); r.append(z+1);
        for(DegreeType d=2; d<=deg; ++d) { r.append(z); }
        return r; }

    template<class OP> Graded<X> apply(OP op, Graded<X> const& a) {
        Graded<X> r;
        while(r.size()!=a.size()) { compute(r,op,a); }
        return r;
    }

    template<class OP> Graded<X> apply(OP op, Graded<X> const& a, Int n) {
        Graded<X> r;
        while(r.size()!=a.size()) { compute(r,op,a,n); }
        return r;
    }

  public:
    void test() {
        ARIADNE_TEST_CALL(test_class());
        ARIADNE_TEST_CALL(test_sqr());
        ARIADNE_TEST_CALL(test_rec());
        ARIADNE_TEST_CALL(test_pow());
        ARIADNE_TEST_CALL(test_sqrt());
        ARIADNE_TEST_CALL(test_exp());
        ARIADNE_TEST_CALL(test_log());
        ARIADNE_TEST_CALL(test_sin());
        ARIADNE_TEST_CALL(test_cos());
    }
  private:
    void test_class() {
    }

    void test_sqr() {
        ARIADNE_TEST_EQUALS(apply(Sqr(),variable(2.0,5u)), (Graded<X>{4.0,4.0,1.0,0.0,0.0,0.0}));
        ARIADNE_TEST_EQUALS(apply(Sqr(),Graded<X>{3.0,2.0,1.0,0.0,0.0,0.0}), (Graded<X>{9.0,12.0,10.0,4.0,1.0,0.0}) );
    }

    void test_rec() {
        ARIADNE_TEST_EQUALS(apply(Rec(),variable(2.0,5u)), (Graded<X>{0.5,-0.25,0.125,-0.0625,0.03125,-0.015625}) );
        ARIADNE_TEST_EQUALS(apply(Rec(),variable(1.0,5u)), (Graded<X>{1.0,-1.0,1.0,-1.0,1.0,-1.0}) );
        ARIADNE_TEST_EQUALS(apply(Rec(),variable(0.5,5u)), (Graded<X>{2.0,-4.0,8.0,-16.0,32.0,-64.0}) );
        ARIADNE_TEST_EQUALS(apply(Rec(),variable(-1.0,5u)), (Graded<X>{-1.0,-1.0,-1.0,-1.0,-1.0,-1.0}) );
    }

    void test_pow() {
        ARIADNE_TEST_EQUAL (apply(Pow(),variable(3.0,5u),2), apply(Sqr(),variable(3.0,5u)) );
        ARIADNE_TEST_EQUAL (apply(Pow(),variable(3.0,5u),4), apply(Sqr(),apply(Sqr(),variable(3.0,5u))) );
        ARIADNE_TEST_EQUAL (apply(Pow(),variable(3.0,5u),4), apply(Sqr(),apply(Pow(),variable(3.0,5u),2)) );
        ARIADNE_TEST_EQUAL (apply(Pow(),variable(3.0,5u),4), apply(Pow(),apply(Sqr(),variable(3.0,5u)),2) );

        ARIADNE_TEST_EQUAL (apply(Pow(),variable(2.0,7u),6), (apply(Pow(),apply(Sqr(),variable(2,7u)),3)) )
        ARIADNE_TEST_EQUAL (apply(Pow(),variable(2.0,7u),6), (apply(Sqr(),apply(Pow(),variable(2,7u),3))) )

        ARIADNE_TEST_EQUALS(apply(Pow(),variable(3.0,5u),0), (Graded<X>{1.0,0.0,0.0,0.0,0.0,0.0}) );
        ARIADNE_TEST_EQUALS(apply(Pow(),variable(3.0,5u),1), (Graded<X>{3.0,1.0,0.0,0.0,0.0,0.0}) );
        ARIADNE_TEST_EQUALS(apply(Pow(),variable(3.0,5u),2), (Graded<X>{9.0,6.0,1.0,0.0,0.0,0.0}) );
        ARIADNE_TEST_EQUALS(apply(Pow(),variable(1.0,5u),5), (Graded<X>{1.0,5.0,10.0,10.0,5.0,1.0}) );
        ARIADNE_TEST_EQUALS(apply(Pow(),variable(3.0,5u),5), (Graded<X>{243.0,405.0,270.0,90.0,15.0,1.0}) );
        ARIADNE_TEST_EQUALS(apply(Pow(),variable(-3.0,5u),5), (Graded<X>{-243.0,405.0,-270.0,90.0,-15.0,1.0}) );
        ARIADNE_TEST_EQUAL (apply(Pow(),variable(-3.0,7u),2), (apply(Sqr(),variable(-3,7u))) )
        ARIADNE_TEST_EQUAL (apply(Pow(),variable(2.0,7u),6), (apply(Sqr(),apply(Pow(),variable(2,7u),3))) )
        ARIADNE_TEST_EQUALS(apply(Pow(),Graded<X>{2.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},3), (Graded<X>{8.0,0.0,12.0,0.0,6.0,0.0,1.0,0.0}) );

        ARIADNE_TEST_EQUAL (apply(Pow(),Graded<X>{4.0,4.0,1.0, 0,0,0,0,0},3), (Graded<X>{64.0,192.0,240.0, 160.0,60.0,12.0,1.0,0.0}));

        ARIADNE_TEST_EQUAL (apply(Pow(),variable(2.0,5u),-1), apply(Rec(),variable(2.0,5u)) );
        ARIADNE_TEST_EQUAL (apply(Pow(),variable(2.0,5u),-6), apply(Rec(),apply(Pow(),variable(2.0,5u),6)) );

    }

    void test_sqrt() {
        ARIADNE_TEST_EQUALS( apply(Sqrt(),variable(4.0,5u)),
                             (Graded<X>{2,0.25,-0.015625,0.001953125,-0.00030517578125,0.00005340576171875}) );
        ARIADNE_TEST_EQUALS(apply(Sqrt(),variable(1.0,5u)), (Graded<X>{1,0.5,-0.125,0.0625,-0.0390625,0.02734375}) );
        ARIADNE_TEST_EQUALS(apply(Sqrt(),variable(0.25,5u)), (Graded<X>{0.5,1.0,-1.0,2.0,-5.0,14.0}) );
    }

    void test_exp() {
        static const X exp1=exp(X(1));
        ARIADNE_TEST_EQUALS( apply(Exp(),variable(0.0,5u)),
                             (Graded<X>{1.0,1.0,1.0/2,1.0/6,1.0/24,1.0/120}) );
        ARIADNE_TEST_EQUALS( apply(Exp(),variable(1.0,5u)),
                             (Graded<X>{exp1,exp1/1,exp1/2,exp1/6,exp1/24,exp1/120}) );
    }

    void test_log() {
        static const X log2=log(2);
        ARIADNE_TEST_EQUALS( apply(Log(),variable(1.0,5u)),
                             (Graded<X>{0.0,1.0,-1.0/2,1.0/3,-1.0/4,1.0/5}) );
        ARIADNE_TEST_EQUALS( apply(Log(),variable(2.0,5u)),
                             (Graded<X>{log2,1.0/2,-1.0/8,1.0/24,-1.0/64,1.0/160}) );
    }

    void test_sin() {
        const X one=1;
        const X sin1=sin(one);
        const X cos1=cos(one);

        ARIADNE_TEST_EQUALS( apply(Sin(),variable(0.0,5u)),
                             (Graded<X>{0.0,1.0,0.0,-1.0/6,0.0,1.0/120}) );
        ARIADNE_TEST_EQUALS ( apply(Sin(),variable(1.0,5u)),
                             (Graded<X>{sin1,cos1,-sin1/2,-cos1/6,sin1/24,cos1/120}) );
    }

    void test_cos() {
        const X one=1;
        const X sin1=sin(one);
        const X cos1=cos(one);
        ARIADNE_TEST_EQUALS( apply(Cos(),variable(0.0,5u)),
                             (Graded<X>{1.0,0.0,-1.0/2,0.0,1.0/24,0.0}) );
        // Only go to degree 4 here due to roundoff error
        ARIADNE_TEST_EQUALS ( apply(Cos(),variable(1.0,4u)),
                             (Graded<X>{cos1,-sin1,-cos1/2,sin1/6,cos1/24}) );
    }

};


Int main() {
    TestGraded<FloatDP>().test();

    return ARIADNE_TEST_FAILURES;
}
