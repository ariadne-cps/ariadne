/***************************************************************************
 *            test_series.cc
 *
 *  Copyright 20-15  Pieter Collins
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

#include "numeric/numeric.h"
#include "algebra/series.h"
#include "function/taylor_series.h"

#include "test.h"

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
        ARIADNE_TEST_CONSTRUCT( Series<Float>, series, (Rec(),1.0) );
        ARIADNE_TEST_EQUALS(series[0],1.0);
        ARIADNE_TEST_EQUALS(series[10],1.0);
        ARIADNE_TEST_EQUALS(series[1],-1.0);
        ARIADNE_TEST_EQUALS(series[127],-1.0);
        ARIADNE_TEST_EQUALS(series[255],-1.0);
//        ARIADNE_TEST_EQUALS(series[32767],1.0);
        std::cerr<<std::setprecision(18);
        ExactFloat::set_output_precision(18);
        ValidatedFloat::set_output_precision(18);
    }
    void test_rec() {
        ARIADNE_TEST_EQUALS( Series<Float>(Rec(),2.0).coefficients(5), (List<Float>{0.5,-0.25,0.125,-0.0625,0.03125,-0.015625}) );
        ARIADNE_TEST_EQUALS(Series<Float>(Rec(),1.0).coefficients(5), (List<Float>{1.0,-1.0,1.0,-1.0,1.0,-1.0}) );
        ARIADNE_TEST_EQUALS(Series<Float>(Rec(),0.5).coefficients(5), (List<Float>{2.0,-4.0,8.0,-16.0,32.0,-64.0}) );
        ARIADNE_TEST_EQUALS(Series<Float>(Rec(),-1.0).coefficients(5), (List<Float>{-1.0,-1.0,-1.0,-1.0,-1.0,-1.0}) );
    }

    void test_sqrt() {
        ARIADNE_TEST_EQUALS( Series<Float>(Sqrt(),4.0).coefficients(5),
                             (List<Float>{2,0.25,-0.015625,0.001953125,-0.00030517578125,0.00005340576171875}) );
        ARIADNE_TEST_EQUALS(Series<Float>(Sqrt(),1.0).coefficients(5), (List<Float>{1,0.5,-0.125,0.0625,-0.0390625,0.02734375}) );
        ARIADNE_TEST_EQUALS(Series<Float>(Sqrt(),0.25).coefficients(5), (List<Float>{0.5,1.0,-1.0,2.0,-5.0,14.0}) );
    }

    void test_exp() {
        static const Float exp1=exp(Float(1));
        ARIADNE_TEST_EQUALS( Series<Float>(Exp(),0.0).coefficients(5),
                             (List<Float>{1.0,1.0,1.0/2,1.0/6,1.0/24,1.0/120}) );
        ARIADNE_TEST_EQUALS( Series<Float>(Exp(),1.0).coefficients(5),
                             (List<Float>{exp1,exp1/1,exp1/2,exp1/6,exp1/24,exp1/120}) );
    }

    void test_log() {
        static const Float log2=log(2);
        ARIADNE_TEST_EQUALS( Series<Float>(Log(),1.0).coefficients(5),
                             (List<Float>{0.0,1.0,-1.0/2,1.0/3,-1.0/4,1.0/5}) );
        ARIADNE_TEST_EQUALS( Series<Float>(Log(),2.0).coefficients(5),
                             (List<Float>{log2,1.0/2,-1.0/8,1.0/24,-1.0/64,1.0/160}) );
    }

    void test_sin() {
        const Float one=1;
        const Float sin1=sin(one);
        const Float cos1=cos(one);

        ARIADNE_TEST_EQUALS( Series<Float>(Sin(),0.0).coefficients(5),
                             (List<Float>{0.0,1.0,0.0,-1.0/6,0.0,1.0/120}) );
        ARIADNE_TEST_EQUAL ( Series<Float>(Sin(),1.0).coefficients(5),
                             (List<Float>{sin1,cos1,-sin1/2,-cos1/6,sin1/24,cos1/120}) );
    }

    void test_cos() {
        const Float one=1;
        const Float sin1=sin(one);
        const Float cos1=cos(one);
        ARIADNE_TEST_EQUALS( Series<Float>(Cos(),0.0).coefficients(5),
                             (List<Float>{1.0,0.0,-1.0/2,0.0,1.0/24,0.0}) );
        ARIADNE_TEST_EQUAL ( Series<Float>(Cos(),1.0).coefficients(4),
                             (List<Float>{cos1,-sin1,-cos1/2,sin1/6,cos1/24}) );
    }

    void test_taylor_series() {
        TaylorSeries<ValidatedFloat> exp_series(Exp(),ExactInterval(-1,+1),ExactFloat(0),8u);
        ARIADNE_TEST_PRINT(Series<ValidatedNumber>(Exp(),ValidatedNumber(-1,+1)).coefficients(8u));
        ARIADNE_TEST_PRINT(exp_series);
    }

    void test_analytic_function() {
        Rec rec_operator;
        ARIADNE_TEST_CONSTRUCT( AnalyticFunction, rec_function, (rec_operator) );
        ARIADNE_TEST_EQUAL( rec_function.series(0.5).coefficients(5), Series<ApproximateFloat>(Rec(),0.5).coefficients(5) );
    }
};


Int main() {
    TestSeries().test();

    return ARIADNE_TEST_FAILURES;
}
