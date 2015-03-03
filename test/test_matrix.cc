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

#include "config.h"

#include "test.h"

#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"

namespace Ariadne {
typedef Matrix<ApproximateFloat64> ApproximateFloatMatrix;
}

using namespace std;
using namespace Ariadne;


class TestMatrix {
  public:
    Void test();
  private:
    Void test_concept();
    Void test_misc();
};

Void
TestMatrix::test()
{
    ARIADNE_TEST_CALL(test_concept());
    ARIADNE_TEST_CALL(test_misc());
}

Void
TestMatrix::test_concept()
{
    ApproximateFloat64 fx(1);
    BoundedFloat64 ix(1);
    ExactFloat64 ex(1);
    Vector<ApproximateFloat64> fv;
    Vector<BoundedFloat64> iv;
    Vector<ExactFloat64> ev;
    Matrix<ApproximateFloat64> fA;
    Matrix<BoundedFloat64> iA;
    Matrix<ExactFloat64> eA;

    fv=fv+fv; iv=ev+ev; iv=ev+iv; iv=iv+ev; iv=iv+iv;
    fv=fv-fv; iv=ev-ev; iv=ev-iv; iv=iv-ev; iv=iv-iv;
    fv=fx*fv; iv=ex*ev; iv=ex*iv; iv=ix*ev; iv=ix*iv;
    fv=fv*fx; iv=ev*ex; iv=ev*ix; iv=iv*ex; iv=iv*ix;
    fv=fv/fx; iv=ev/ex; iv=ev/ix; iv=iv/ex; iv=iv/ix;

    fA=fA+fA; iA=eA+eA; iA=eA+iA; iA=iA+eA; iA=iA+iA;
    fA=fA-fA; iA=eA-eA; iA=eA-iA; iA=iA-eA; iA=iA-iA;

    fA=fx*fA; iA=ex*eA; iA=ex*iA; iA=ix*eA; iA=ix*iA;
    fA=fA*fx; iA=eA*ex; iA=eA*ix; iA=iA*ex; iA=iA*ix;
    //fv=fA*fv; iv=eA*ev; iv=eA*iv; iv=iA*ev; iv=iA*iv;
    //fA=fA*fA; iA=eA*eA; iA=eA*iA; iA=iA*eA; iA=iA*iA;

    fv=fA*fv; iv=eA*ev; iv=eA*iv; iv=iA*ev; iv=iA*iv;
    fA=fA*fA; iA=eA*eA; iA=eA*iA; iA=iA*eA; iA=iA*iA;
}


Void
TestMatrix::test_misc()
{
    ApproximateFloat64 x=2.25;
    BoundedFloat64 ix(1.5,2.25);
    ApproximateFloat64 Aptr[9]={-1.0,3.0,1.0, -1.0,1.0,2.0, 2.0,1.0,1.0};
    BoundedFloat64 iAptr[4]={-1,3, -1,1};

    Matrix<ApproximateFloat64> A0;
    ARIADNE_TEST_PRINT(A0);
    Matrix<ApproximateFloat64> A1(3,2);
    ARIADNE_TEST_PRINT(A1);
    Matrix<ApproximateFloat64> A2(3,3,Aptr);
    ARIADNE_TEST_PRINT(A2);
    Matrix<ApproximateFloat64> A3({{-1.0,3.0,1.0}, {-1.0,1.0,2.0}, {2.0,1.0,1.0}});
    ARIADNE_TEST_PRINT(A3);

    for(SizeType i=0; i!=A2.row_size(); ++i) {
        for(SizeType j=0; j!=A2.column_size(); ++j) {
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

    A0=Matrix<ApproximateFloat64>::zero(2,3);
    ARIADNE_TEST_PRINT(A0);
    A1=Matrix<ApproximateFloat64>::identity(4);
    ARIADNE_TEST_PRINT(A1);

    typedef ApproximateFloatMatrix MatrixType;
    ARIADNE_TEST_EQUALS(+ApproximateFloatMatrix({{1.,2.},{3.,4.}}),ApproximateFloatMatrix({{1.,2.},{3.,4.}}));
    ARIADNE_TEST_EQUALS(-ApproximateFloatMatrix({{1.,2.},{3.,4.}}),ApproximateFloatMatrix({{-1.,-2.},{-3.,-4.}}));
    ARIADNE_TEST_EQUALS(ApproximateFloatMatrix({{1.,2.},{3.,4.}})+ApproximateFloatMatrix({{5.,7.},{8.,6.}}),ApproximateFloatMatrix({{6.,9.},{11.,10.}}));
    ARIADNE_TEST_EQUALS(ApproximateFloatMatrix({{1.,2.},{3.,4.}})-ApproximateFloatMatrix({{5.,7.},{8.,6.}}),ApproximateFloatMatrix({{-4.,-5.},{-5.,-2.}}));
    ARIADNE_TEST_EQUALS(ApproximateFloatMatrix({{1.,2.},{3.,4.}})*ApproximateFloatMatrix({{5.,7.},{8.,6.}}),ApproximateFloatMatrix({{21.,19.},{47.,45.}}));

    // Transpose operations
    ARIADNE_TEST_EQUALS(ApproximateFloatMatrix(transpose(ApproximateFloatMatrix({{1.,2.,3.},{4.,5.,6.}}))),ApproximateFloatMatrix({{1.,4.},{2.,5.},{3.,6.}}));
    ARIADNE_TEST_EQUALS(transpose(ApproximateFloatMatrix({{1.,2.},{3.,4.}}))*ApproximateFloatMatrix({{5.,7.},{8.,6.}}),ApproximateFloatMatrix({{1.,3.},{2.,4.}})*ApproximateFloatMatrix({{5.,7.},{8.,6.}}));
    ARIADNE_TEST_EQUALS(ApproximateFloatMatrix({{1.,2.},{3.,4.}})*transpose(ApproximateFloatMatrix({{5.,7.},{8.,6.}})),ApproximateFloatMatrix({{1.,2.},{3.,4.}})*ApproximateFloatMatrix({{5.,8.},{7.,6.}}));
    ARIADNE_TEST_EQUALS(transpose(ApproximateFloatMatrix({{1.,2.,3.},{4.,5.,6.}}))*ApproximateFloatVector({5.,7.}),ApproximateFloatMatrix({{1.,4.},{2.,5.},{3.,6.}})*ApproximateFloatVector({5.,7.}));
}

Int main() {
    TestMatrix().test();
    return ARIADNE_TEST_FAILURES;
}

