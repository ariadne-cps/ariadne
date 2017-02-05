/***************************************************************************
 *            test_matrix.cpp
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

#include "test.hpp"

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/covector.hpp"

namespace Ariadne {
typedef Matrix<Float64Approximation> FloatApproximationMatrix;
}

using namespace std;
using namespace Ariadne;


class TestMatrix {
    Precision64 pr;
  public:
    Void test();
  private:
    Void test_concept();
    Void test_project();
    Void test_misc();
};

Void
TestMatrix::test()
{
    ARIADNE_TEST_CALL(test_concept());
    ARIADNE_TEST_CALL(test_project());
    ARIADNE_TEST_CALL(test_misc());
}

Void
TestMatrix::test_concept()
{
    Float64Approximation fx(1);
    Float64Bounds ix(1);
    Float64Value ex(1);
    Vector<Float64Approximation> fv;
    Vector<Float64Bounds> iv;
    Vector<Float64Value> ev;
    Matrix<Float64Approximation> fA;
    Matrix<Float64Bounds> iA;
    Matrix<Float64Value> eA;

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
TestMatrix::test_project()
{
    typedef Float64 X;
    ARIADNE_TEST_CONSTRUCT(Matrix<X>,A,({{11,12,13,14,15},{21,22,23,24,25},{31,32,33,34,35}},pr));
    ARIADNE_TEST_CONSTRUCT(Matrix<X>,B,({{120,130,140},{220,230,240}},pr));

    ARIADNE_TEST_ASSIGN_CONSTRUCT(Vector<X>,v,A[range(0,2)][1]);
    ARIADNE_TEST_EQUAL(v[1],A[1][1]);
    ARIADNE_TEST_ASSIGN_CONSTRUCT(Covector<X>,u,A[2][range(1,4)]);
    ARIADNE_TEST_EQUAL(u[2],A[2][3]);

    ARIADNE_TEST_EQUALS(A[range(0,2)][range(1,4)].row_size(),2);
    ARIADNE_TEST_EQUALS(A[range(0,2)][range(1,4)].column_size(),3);
    ARIADNE_TEST_CONSTRUCT(MatrixRange<Matrix<X>>,AR,(A[range(0,2)][range(1,4)]));
    ARIADNE_TEST_CONSTRUCT(Matrix<X>,ACR,(A[range(0,2)][range(1,4)]));

    ARIADNE_TEST_EQUAL(AR[1][2],A[1][3]);
    ARIADNE_TEST_EQUAL(ACR[1][2],A[1][3]);
    ARIADNE_TEST_EXECUTE(AR=B);
    ARIADNE_TEST_PRINT(A);
    ARIADNE_TEST_EQUALS(A[range(0,2)][range(1,4)],B);
//    ARIADNE_TEST_EQUALS(AR,B);
//    ARIADNE_TEST_EQUALS(ACR,B);

}


Void
TestMatrix::test_misc()
{
    Float64Approximation x={2.25,pr};
    Float64Bounds ix={1.5_exact,2.25_exact};
    Array<Float64Approximation> Aary={{-1.0,3.0,1.0, -1.0,1.0,2.0, 2.0,1.0,1.0},pr};
    Array<Float64Bounds> iAary={{-1.0_x,3.0_x, -1.0_x,1.0_x},pr};
    Float64Approximation* Aptr=Aary.begin();
    Float64Bounds* iAptr=iAary.begin();

    Matrix<Float64Approximation> A0;
    ARIADNE_TEST_PRINT(A0);
    Matrix<Float64Approximation> A1(3,2,pr);
    ARIADNE_TEST_PRINT(A1);
    Matrix<Float64Approximation> A2(3,3,Aary.begin());
    ARIADNE_TEST_PRINT(A2);
    Matrix<Float64Approximation> A3({{-1.0,3.0,1.0}, {-1.0,1.0,2.0}, {2.0,1.0,1.0}},pr);
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

    A0=Matrix<Float64Approximation>::zero(2,3);
    ARIADNE_TEST_PRINT(A0);
    A1=Matrix<Float64Approximation>::identity(4);
    ARIADNE_TEST_PRINT(A1);

    typedef FloatApproximationMatrix MatrixType;
    ARIADNE_TEST_EQUALS(+FloatApproximationMatrix({{1.,2.},{3.,4.}},pr),FloatApproximationMatrix({{1.,2.},{3.,4.}},pr));
    ARIADNE_TEST_EQUALS(-FloatApproximationMatrix({{1.,2.},{3.,4.}},pr),FloatApproximationMatrix({{-1.,-2.},{-3.,-4.}},pr));
    ARIADNE_TEST_EQUALS(FloatApproximationMatrix({{1.,2.},{3.,4.}},pr)+FloatApproximationMatrix({{5.,7.},{8.,6.}},pr),FloatApproximationMatrix({{6.,9.},{11.,10.}},pr));
    ARIADNE_TEST_EQUALS(FloatApproximationMatrix({{1.,2.},{3.,4.}},pr)-FloatApproximationMatrix({{5.,7.},{8.,6.}},pr),FloatApproximationMatrix({{-4.,-5.},{-5.,-2.}},pr));
    ARIADNE_TEST_EQUALS(FloatApproximationMatrix({{1.,2.},{3.,4.}},pr)*FloatApproximationMatrix({{5.,7.},{8.,6.}},pr),FloatApproximationMatrix({{21.,19.},{47.,45.}},pr));

    // Transpose operations
    ARIADNE_TEST_EQUALS(FloatApproximationMatrix(transpose(FloatApproximationMatrix({{1.,2.,3.},{4.,5.,6.}},pr))),FloatApproximationMatrix({{1.,4.},{2.,5.},{3.,6.}},pr));
    ARIADNE_TEST_EQUALS(transpose(FloatApproximationMatrix({{1.,2.},{3.,4.}},pr))*FloatApproximationMatrix({{5.,7.},{8.,6.}},pr),FloatApproximationMatrix({{1.,3.},{2.,4.}},pr)*FloatApproximationMatrix({{5.,7.},{8.,6.}},pr));
    ARIADNE_TEST_EQUALS(FloatApproximationMatrix({{1.,2.},{3.,4.}},pr)*transpose(FloatApproximationMatrix({{5.,7.},{8.,6.}},pr)),FloatApproximationMatrix({{1.,2.},{3.,4.}},pr)*FloatApproximationMatrix({{5.,8.},{7.,6.}},pr));
    ARIADNE_TEST_EQUALS(transpose(FloatApproximationMatrix({{1.,2.,3.},{4.,5.,6.}},pr))*FloatApproximationVector({5.,7.},pr),FloatApproximationMatrix({{1.,4.},{2.,5.},{3.,6.}},pr)*FloatApproximationVector({5.,7.},pr));
}

Int main() {
    TestMatrix().test();
    return ARIADNE_TEST_FAILURES;
}

