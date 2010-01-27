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
}

void 
TestMatrix::test_concept() 
{
    Float fx;
    Interval ix;
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

    fv=prod(fA,fv); iv=prod(fA,fv); iv=prod(fA,iv); iv=prod(iA,fv); iv=prod(iA,iv); 
    fA=prod(fA,fA); iA=prod(fA,fA); iA=prod(fA,iA); iA=prod(iA,fA); iA=prod(iA,iA); 
    
    // Test variadic constructor and comma operator
    fA = Matrix<Float>(2,2);
    fA[0][0] = 1.0; fA[0][1] = 2.0; 
    fA[1][0] = 3.3; fA[1][1] = 4.4;
    Matrix<Float> fAv(2,2, 1.0, 2.0, 3.3, 4.4);
    ARIADNE_TEST_EQUAL(fA,fAv);
    fAv = Matrix<Float>(2,2);
    ARIADNE_TEST_COMPARE(fA,!=,fAv);
    fAv = 1.0, 2.0, 3.3, 4.4;
    ARIADNE_TEST_EQUAL(fA,fAv);
       
    iA = Matrix<Interval>(2,2);
    iA[0][0] = 1.0; iA[0][1] = 2.0; 
    iA[1][0] = 3.3; iA[1][1] = 4.4;
    Matrix<Interval> iAv(2,2, 1.0, 2.0, 3.3, 4.4);
    ARIADNE_TEST_EQUAL(iA,iAv);
    iAv = Matrix<Interval>(2,2);
    ARIADNE_TEST_COMPARE(iA,!=,iAv);
    iAv = Interval(1.0,1.0), Interval(2.0,2.0), Interval(3.3,3.3), Interval(4.4,4.4);
    ARIADNE_TEST_EQUAL(iA,iAv);
       
#ifdef HAVE_RATIONAL
    Matrix<Rational> rA(2,2);
    rA[0][0] = 1.0; rA[0][1] = 2.0; 
    rA[1][0] = 3.3; rA[1][1] = 4.4;
    Matrix<Rational> rAv (2,2, 1.0, 2.0, 3.3, 4.4);
    ARIADNE_TEST_EQUAL(rA,rAv);
    rAv = Matrix<Rational>(2,2);
    ARIADNE_TEST_COMPARE(rA,!=,rAv);
    rAv = 1.0, 2.0, 3.3, 4.4;
    ARIADNE_TEST_EQUAL(rA,rAv);
#endif // HAVE_RATIONAL
         
}


void 
TestMatrix::test_misc()
{
    /*
      Float x=2.25;
      Interval ix(1.5,2.25);
      Float Aptr[9]={-1.0,3.0,1.0, -1.0,1.0,2.0, 2.0,1.0,1.0};
      Interval iAptr[4]={-1.0,3.0, -1.0,1.0};
  
      Matrix<Float> A0;
      cout << "A0=" << A0 << endl;
      Matrix<Float> A1(3,2);
      cout << "A1=" << A1 << endl;
      Matrix<Float> A2(3,3,Aptr,3,1);
      cout << "A2=" << A2 << endl;
      Matrix<Float> A3=make_matrix<Float>("[-1.0,3.0,1.0; -1.0,1.0,2.0; 2.0,1.0,1.0]");
      cout << "A3=" << A3 << endl;

      cout << "A1(0,0)= " << A1(0,0) << endl;
      cout << "A2(0,0)= " << A2(0,0) << endl;
      assert(A2==A3);
  
      A1[0][0]=1.0;
      cout << "A1[0][0]= " << A1[0][0] << endl;
      assert(A1[0][0]==1.0);
      cout << "A2[0][0]= " << A2[0][0] << endl;
      assert(A2==A3);
  
      A1(0,0)=2.0;
      cout << "A1(0,0)= " << A1(0,0) << endl;
      assert(A1(0,0)==2.0);
      assert(A1[0][0]==A1(0,0));
  
      A0=Matrix<Float>::zero(2,3);
      cout << "A0= " << A0 << endl;
      A1=Matrix<Float>::identity(4);
      cout << "A1= " << A1 << endl;
      cout << endl;
  
      A1=make_matrix<Float>("[2,1,1;1,1,-1;1,-1,3]");
      Matrix<Float> fA0;

      fA0=-A1;
      cout << A0 << " = -" << A1 << endl;
      fA0=A1+A2;
      cout << A0 << " = " << A1 << " + " << A2 << endl;
      fA0=A1-A2;
      cout << A0 << " = " << A1 << " - " << A2 << endl;
      fA0=x*A2;
      cout << A0 << " = " << x << " * " << A2 << endl;
      fA0=A1*x;
      cout << A0 << " = " << A1 << " * " << x << endl;
      fA0=Matrix<Float>(A1)/x;
      cout << A0 << " = " << A1 << " / " << x << endl;
      fA0=Matrix<Float>(A1)*A2;
      cout << A0 << " = " << A1 << " * " << A2 << endl;
      cout << endl;

  
      Matrix<Interval> iA0;
      cout << "iA1= " << iA0 << endl;
      Matrix<Interval> iA1(3,2);
      cout << "iA2= " << iA1 << endl;
      Matrix<Interval> iA2(2,2,iAptr,3,1);
      cout << "iA3= " << iA2 << endl;
      Matrix<Interval> iA3(A1);
      cout << "iA4= " << iA3 << endl;
      Matrix<Interval> iA4=make_matrix<Interval>("[[1.875,2.125],[0.75,1.25];[-1.25,-0.875],[0.5,1.75]]");
      cout << "iA5= " << iA4 << endl;
  
      iA0=Matrix<Interval>::zero(2,3);
      cout << "iA0= " << iA0 << endl;
      iA1=Matrix<Interval>::identity(4);
      cout << "iA1= " << iA1 << endl;

      cout << "norm(iA4)=" << norm(iA4) << endl;
      cout << "norm(iA4).upper()=" << norm(iA4).upper() << endl;
    */  
}

int main() {
    TestMatrix().test();
    return ARIADNE_TEST_FAILURES;
}

