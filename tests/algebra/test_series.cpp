/***************************************************************************
 *            test_series.cpp
 *
 *  Copyright  20-15  Pieter Collins
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
#include "algebra/series.hpp"
#include "function/taylor_series.hpp"

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

class TestSeries
{
  public:
    void test() {
        ARIADNE_TEST_CALL(test_class());
        ARIADNE_TEST_CALL(test_rec());
        ARIADNE_TEST_CALL(test_sqrt());
        ARIADNE_TEST_CALL(test_exp());
        ARIADNE_TEST_CALL(test_log());
        ARIADNE_TEST_CALL(test_sin());
        ARIADNE_TEST_CALL(test_cos());

        ARIADNE_TEST_CALL(test_taylor_series());
        ARIADNE_TEST_CALL(test_analytic_function());
    }
  private:
    void test_class() {
        ARIADNE_TEST_CONSTRUCT( Series<FloatDP>, series, (Rec(),1.0) );
        ARIADNE_TEST_EQUALS(series[0],1.0);
        ARIADNE_TEST_EQUALS(series[10],1.0);
        ARIADNE_TEST_EQUALS(series[1],-1.0);
        ARIADNE_TEST_EQUALS(series[127],-1.0);
        ARIADNE_TEST_EQUALS(series[255],-1.0);
//        ARIADNE_TEST_EQUALS(series[32767],1.0);
        std::cerr<<std::setprecision(18);
        FloatDPValue::set_output_places(18);
        FloatDPBounds::set_output_places(18);
    }
    void test_rec() {
        ARIADNE_TEST_EQUALS( Series<FloatDP>(Rec(),2.0).coefficients(5), (List<FloatDP>{0.5,-0.25,0.125,-0.0625,0.03125,-0.015625}) );
        ARIADNE_TEST_EQUALS(Series<FloatDP>(Rec(),1.0).coefficients(5), (List<FloatDP>{1.0,-1.0,1.0,-1.0,1.0,-1.0}) );
        ARIADNE_TEST_EQUALS(Series<FloatDP>(Rec(),0.5).coefficients(5), (List<FloatDP>{2.0,-4.0,8.0,-16.0,32.0,-64.0}) );
        ARIADNE_TEST_EQUALS(Series<FloatDP>(Rec(),-1.0).coefficients(5), (List<FloatDP>{-1.0,-1.0,-1.0,-1.0,-1.0,-1.0}) );
    }

    void test_sqrt() {
        ARIADNE_TEST_EQUALS( Series<FloatDP>(Sqrt(),4.0).coefficients(5),
                             (List<FloatDP>{2,0.25,-0.015625,0.001953125,-0.00030517578125,0.00005340576171875}) );
        ARIADNE_TEST_EQUALS(Series<FloatDP>(Sqrt(),1.0).coefficients(5), (List<FloatDP>{1,0.5,-0.125,0.0625,-0.0390625,0.02734375}) );
        ARIADNE_TEST_EQUALS(Series<FloatDP>(Sqrt(),0.25).coefficients(5), (List<FloatDP>{0.5,1.0,-1.0,2.0,-5.0,14.0}) );
    }

    void test_exp() {
        static const FloatDP exp1=exp(FloatDP(1));
        ARIADNE_TEST_EQUALS( Series<FloatDP>(Exp(),0.0).coefficients(5),
                             (List<FloatDP>{1.0,1.0,1.0/2,1.0/6,1.0/24,1.0/120}) );
        ARIADNE_TEST_EQUALS( Series<FloatDP>(Exp(),1.0).coefficients(5),
                             (List<FloatDP>{exp1,exp1/1,exp1/2,exp1/6,exp1/24,exp1/120}) );
    }

    void test_log() {
        static const FloatDP log2=log(2);
        ARIADNE_TEST_EQUALS( Series<FloatDP>(Log(),1.0).coefficients(5),
                             (List<FloatDP>{0.0,1.0,-1.0/2,1.0/3,-1.0/4,1.0/5}) );
        ARIADNE_TEST_EQUALS( Series<FloatDP>(Log(),2.0).coefficients(5),
                             (List<FloatDP>{log2,1.0/2,-1.0/8,1.0/24,-1.0/64,1.0/160}) );
    }

    void test_sin() {
        const FloatDP one=1;
        const FloatDP sin1=sin(one);
        const FloatDP cos1=cos(one);

        ARIADNE_TEST_EQUALS( Series<FloatDP>(Sin(),0.0).coefficients(5),
                             (List<FloatDP>{0.0,1.0,0.0,-1.0/6,0.0,1.0/120}) );
        ARIADNE_TEST_EQUAL ( Series<FloatDP>(Sin(),1.0).coefficients(5),
                             (List<FloatDP>{sin1,cos1,-sin1/2,-cos1/6,sin1/24,cos1/120}) );
    }

    void test_cos() {
        const FloatDP one=1;
        const FloatDP sin1=sin(one);
        const FloatDP cos1=cos(one);
        ARIADNE_TEST_EQUALS( Series<FloatDP>(Cos(),0.0).coefficients(5),
                             (List<FloatDP>{1.0,0.0,-1.0/2,0.0,1.0/24,0.0}) );
        ARIADNE_TEST_EQUAL ( Series<FloatDP>(Cos(),1.0).coefficients(4),
                             (List<FloatDP>{cos1,-sin1,-cos1/2,sin1/6,cos1/24}) );
    }

    void test_taylor_series() {
        TaylorSeries<FloatDPBounds> exp_series(Exp(),ExactIntervalType(-1,+1),FloatDPValue(0),8u);
        ARIADNE_TEST_PRINT(Series<FloatDPBounds>(Exp(),FloatDPBounds(-1,+1)).coefficients(8u));
        ARIADNE_TEST_PRINT(exp_series);
    }

    void test_analytic_function() {
        Rec rec_operator;
        ARIADNE_TEST_CONSTRUCT( AnalyticFunction, rec_function, (rec_operator) );
        ARIADNE_TEST_SAME( rec_function.series(0.5_approx).coefficients(5), Series<FloatDPApproximation>(Rec(),0.5_approx).coefficients(5) );
    }
};


Int main() {
    TestSeries().test();

    return ARIADNE_TEST_FAILURES;
}
