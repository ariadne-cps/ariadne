/***************************************************************************
 *            test_matrix.cc
 *
 *  Copyright  2006  Pieter Collins, Alberto Casagrande
 *  Email Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
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

#define NO_CBLAS

#include <iostream>
#include <fstream>

#include "test.h"

#include "config.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"


using namespace std;
using namespace Ariadne;


class TestMatrix {
  public:
    void test();
  private:
    void test_concept();
    void test_misc();
};

void
TestMatrix::test()
{
    ARIADNE_TEST_CALL(test_concept());
    ARIADNE_TEST_CALL(test_misc());
}

void
TestMatrix::test_concept()
{
    Float fx(1);
    Interval ix(1);
    Vector<Float> fv;
    Vector<Interval> iv;
    Matrix<Float> fA;
    Matrix<Interval> iA;

    // Redo vector concept tests with new header
    fv=fv+fv; iv=fv+fv; iv=fv+iv; iv=iv+fv; iv=iv+iv;
    fv=fv-fv; iv=fv-fv; iv=fv-iv; iv=iv-fv; iv=iv-iv;
    fv=fx*fv; iv=fx*fv; iv=fx*iv; iv=ix*fv; iv=ix*iv;
    fv=fv*fx; iv=fv*fx; iv=fv*ix; iv=iv*fx; iv=iv*ix;
    fv=fv/fx; iv=fv/fx; iv=fv/ix; iv=iv/fx; iv=iv/ix;

    fA=fA+fA; iA=fA+fA; iA=fA+iA; iA=iA+fA; iA=iA+iA;
    fA=fA-fA; iA=fA-fA; iA=fA-iA; iA=iA-fA; iA=iA-iA;

    fA=fx*fA; iA=fx*fA; iA=fx*iA; iA=ix*fA; iA=ix*iA;
    fA=fA*fx; iA=fA*fx; iA=fA*ix; iA=iA*fx; iA=iA*ix;
    //fv=fA*fv; iv=fA*fv; iv=fA*iv; iv=iA*fv; iv=iA*iv;
    //fA=fA*fA; iA=fA*fA; iA=fA*iA; iA=iA*fA; iA=iA*iA;

    fv=fA*fv; iv=fA*fv; iv=fA*iv; iv=iA*fv; iv=iA*iv;
    fA=fA*fA; iA=fA*fA; iA=fA*iA; iA=iA*fA; iA=iA*iA;
}


void
TestMatrix::test_misc()
{
    Float x=2.25;
    Interval ix(1.5,2.25);
    Float Aptr[9]={-1.0,3.0,1.0, -1.0,1.0,2.0, 2.0,1.0,1.0};
    Interval iAptr[4]={-1.0,3.0, -1.0,1.0};

    Matrix<Float> A0;
    ARIADNE_TEST_PRINT(A0);
    Matrix<Float> A1(3,2);
    ARIADNE_TEST_PRINT(A1);
    Matrix<Float> A2(3,3,Aptr);
    ARIADNE_TEST_PRINT(A2);
    Matrix<Float> A3("[-1.0,3.0,1.0; -1.0,1.0,2.0; 2.0,1.0,1.0]");
    ARIADNE_TEST_PRINT(A3);

    for(size_t i=0; i!=A2.row_size(); ++i) {
        for(size_t j=0; j!=A2.column_size(); ++j) {
            ARIADNE_TEST_EQUAL(A2[i][j],A2.get(i,j));
            ARIADNE_TEST_EQUALS(A2[i][j],Aptr[i*A2.column_size()+j]);
        }
    }

    ARIADNE_TEST_EQUAL(A2,A3);

    A1[0][0]=1.0;
    ARIADNE_TEST_EQUALS(A1[0][0],1.0);
    A1.set(0,1,3.0);
    ARIADNE_TEST_EQUALS(A1[0][1],3.0);

    A1.at(1,0)=2.0;
    ARIADNE_TEST_EQUALS(A1[1][0],2.0);

    A0=Matrix<Float>::zero(2,3);
    ARIADNE_TEST_PRINT(A0);
    A1=Matrix<Float>::identity(4);
    ARIADNE_TEST_PRINT(A1);

    ARIADNE_TEST_EQUALS(+FloatMatrix(2,2,1.,2.,3.,4.),FloatMatrix(2,2,1.,2.,3.,4.));
    ARIADNE_TEST_EQUALS(-FloatMatrix(2,2,1.,2.,3.,4.),FloatMatrix(2,2,-1.,-2.,-3.,-4.));
    ARIADNE_TEST_EQUALS(FloatMatrix(2,2,1.,2.,3.,4.)+FloatMatrix(2,2,5.,7.,8.,6.),FloatMatrix(2,2,6.,9.,11.,10.));
    ARIADNE_TEST_EQUALS(FloatMatrix(2,2,1.,2.,3.,4.)-FloatMatrix(2,2,5.,7.,8.,6.),FloatMatrix(2,2,-4.,-5.,-5.,-2.));
    ARIADNE_TEST_EQUALS(FloatMatrix(2,2,1.,2.,3.,4.)*FloatMatrix(2,2,5.,7.,8.,6.),FloatMatrix(2,2,21.,19.,47.,45.));

}

int main() {
    TestMatrix().test();
    return ARIADNE_TEST_FAILURES;
}

