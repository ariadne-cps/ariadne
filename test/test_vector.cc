/***************************************************************************
 *            test_vector.cc
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

#include <iostream>
#include <fstream>
#include <cassert>

#include "config.h"
#include "numeric/numeric.h"
#include "numeric/float.h"
#include "geometry/interval.h"
#include "algebra/vector.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

class TestVector {
  public:
    Void test();
  private:
    Void test_concept();
    Void test_constructors();
    Void test_comparisons();
    Void test_arithmetic();
    Void test_misc();
};

Void
TestVector::test()
{
    ARIADNE_TEST_CALL(test_concept());
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_comparisons());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_misc());
}

Void
TestVector::test_concept()
{
    ApproximateFloat ax(1);
    ValidatedFloat ix(1);
    ExactFloat ex(1);
    Vector<ApproximateFloat> av;
    Vector<ValidatedFloat> iv;
    Vector<ExactFloat> ev;

    iv=Vector<ValidatedFloat>(ev);

    av=av+av;
    iv=ev+ev;
    iv=ev+iv;
    iv=iv+ev;
    iv=iv+iv;
    av=av-av; iv=ev-ev; iv=ev-iv; iv=iv-ev; iv=iv-iv;

    av=ax*av; iv=ex*ev; iv=ex*iv; iv=ix*ev; iv=ix*iv;
    av=av*ax; iv=ev*ex; iv=ev*ix; iv=iv*ex; iv=iv*ix;
    av=av/ax; iv=ev/ex; iv=ev/ix; iv=iv/ex; iv=iv/ix;
}


Void
TestVector::test_constructors()
{
    ARIADNE_TEST_CONSTRUCT( RawFloatVector, v0, (3) );
    ARIADNE_TEST_ASSERT( v0.size()==3 );
    ARIADNE_TEST_ASSERT( v0[0]==0.0 );
    ARIADNE_TEST_ASSERT( v0[1]==0.0 );
    ARIADNE_TEST_ASSERT( v0[2]==0.0 );

    ARIADNE_TEST_CONSTRUCT( RawFloatVector, v1, ({3.25,-0.75,0.0,1.375}) );
    ARIADNE_TEST_ASSERT( v1.size()==4 );
    ARIADNE_TEST_ASSERT( v1[0]==3.25 );
    ARIADNE_TEST_ASSERT( v1[1]==-0.75 );
    ARIADNE_TEST_ASSERT( v1[2]==0.0 );
    ARIADNE_TEST_ASSERT( v1[3]==1.375 );
}


Void
TestVector::test_comparisons()
{
    RawFloatVector v(2); v[0]=1.25; v[1]=-1.375;
    ARIADNE_TEST_PRINT( v );
    ARIADNE_TEST_COMPARE( v , ==, RawFloatVector({1.25,-1.375}) );
    ARIADNE_TEST_COMPARE( v , !=, RawFloatVector({1.25,-1.625}) );
}


Void
TestVector::test_arithmetic()
{
    ARIADNE_TEST_EQUAL( + RawFloatVector({2.0,-3.0,5.0}) , RawFloatVector({2.0,-3.0,5.0}) );
    ARIADNE_TEST_EQUAL( - RawFloatVector({2.0,-3.0,5.0}) , RawFloatVector({-2.0,+3.0,-5.0}) );
    ARIADNE_TEST_EQUAL( RawFloatVector({2.0,-3.0,5.0}) + RawFloatVector({1.25,2.75,-3.5}), RawFloatVector({3.25,-0.25,1.5}) );
    ARIADNE_TEST_EQUAL( RawFloatVector({2.0,-3.0,5.0}) - RawFloatVector({1.25,2.75,-3.5}), RawFloatVector({0.75,-5.75,8.5}) );
    ARIADNE_TEST_EQUAL( Float(4.0) * RawFloatVector({2.0,-3.0,5.0}), RawFloatVector({8.0,-12.0,20.0}) );
    ARIADNE_TEST_EQUAL( RawFloatVector({2.0,-3.0,5.0}) * Float(4.0), RawFloatVector({8.0,-12.0,20.0}) );
    ARIADNE_TEST_EQUAL( RawFloatVector({2.0,-3.0,5.0}) / Float(4.0), RawFloatVector({0.5,-0.75,1.25}) );

    ARIADNE_TEST_EQUAL( sup_norm(RawFloatVector({2.0,-3.0,1.0})), Float(3.0) )
}


Void
TestVector::test_misc()
{

    Int n=3;
    ApproximateFloat vptr[3]={-4.0,3.0,1.0};
    ApproximateFloat x=1.5;

    Vector<ApproximateFloat> v0;
    cout << "v0.size()=" << v0.size() << endl;
    cout << "v0=" << flush; cout << v0 << endl;
    ARIADNE_TEST_NOTIFY("Constructor Vector<X>(SizeType, const X*) is unsafe and has been removed.");
    Array<ApproximateFloat> a1(vptr,vptr+3);
    Vector<ApproximateFloat> v1(a1);
    cout << "v1=" << v1 << endl;
    Vector<ApproximateFloat> v2=Vector<ApproximateFloat>({2.375,4.25,-1.25});
    cout << "v2=" << v2 << endl;
    cout << "norm(v1)=" << norm(v1) << "  norm(v2)=" << norm(v2) << endl;
    assert(norm(v1)==4);
    assert(norm(v2)==4.25);

    Vector<ApproximateFloat> v3(1);
    cout << "v3=" << v3 << endl;
    Vector<ApproximateFloat> v4=v2;
    cout << "v4=" << v4 << endl;
    Vector<ApproximateFloat> v5={-4.0,3.0,1.0};
    cout << "v5=" << v5 << endl;
    assert(v1==v5);
    cout << endl;

    Vector<ApproximateFloat> vf0;
    v1=Vector<ApproximateFloat>({0.25,-1.5});
    v2=Vector<ApproximateFloat>({-0.5,2.25});
    vf0=-v1;
    cout << vf0 << " = -" << v1 << endl;
    vf0=Vector<ApproximateFloat>(v1)+v2;
    cout << vf0 << " = " << v1 << " + " << v2 << endl;
    vf0=Vector<ApproximateFloat>(v1)-v2;
    cout << vf0 << " = " << v1 << " - " << v2 << endl;
    vf0=x*Vector<ApproximateFloat>(v2);
    cout << vf0 << " = " << x << " * " << v2 << endl;
    vf0=Vector<ApproximateFloat>(v1)*x;
    cout << vf0 << " = " << v1 << " * " << x << endl;
    vf0=Vector<ApproximateFloat>(v1)/x;
    cout << vf0 << " = " << v1 << " / " << x << endl;
    cout << endl;

    Vector< ValidatedFloat > iv1=Vector<ValidatedFloat>({ValidatedFloat{0.984375,1.015625},{2.25,2.375},{4.0,4.375},{-0.03125,0.015625}});
    cout << "iv1=" << iv1 << endl;
    cout << "norm(iv1)=" << norm(iv1) << endl;
    cout << "norm(iv1).upper()=" << norm(iv1).upper() << endl;

    Vector< ValidatedFloat > iv2=Vector<ValidatedFloat>({{-1,1},{-1,1}});
    cout << "iv2=" << iv2 << endl;
    Vector< ValidatedFloat > iv3(3);
    cout << "iv3=" << iv3 << endl;
    iv3=Vector<ValidatedFloat>({{4.25,4.25},{2.375,2.375}});
    cout << "iv3=" << iv3 << endl;
    ValidatedFloat ix=ValidatedFloat(-2,1);

    Vector< ValidatedFloat > iv0;
    cout << "iv0=" << iv0 << endl;
    iv1=iv0;
    cout << "iv1=" << iv1 << endl;
    iv1=iv2;
    cout << "iv1=" << iv1 << endl;
    cout << endl;

    ValidatedFloat ix2=iv2[0];
    ValidatedFloat ix3=iv3[0];
    ValidatedFloat ix1=ix2+ix3;
    ix1=ix2+ix3;

    cout << "iv2=" << iv2 << ", iv3=" << iv3 << endl;
    iv1=iv2+iv3;
    cout << iv1 << " = " << iv2 << " + " << iv3 << endl;
    iv1=iv2-iv3;
    cout << iv1 << " = " << iv2 << " - " << iv3 << endl;
    iv1=ix*iv3;
    cout << iv1 << " = " << ix << " * " << iv3 << endl;
    iv1=iv2*ix;
    cout << iv1 << " = " << iv2 << " * " << ix << endl;
    ix=ValidatedFloat(1,2);
    iv1=iv2/ix;
    cout << iv1 << " = " << iv2 << " / " << ix << endl;
    cout << endl;

    Vector<ExactFloat> ev1(reinterpret_cast<Vector<ExactFloat>const&>(v1));
    ExactFloat ex(reinterpret_cast<ExactFloat const&>(x));
    iv0=iv1+ev1;
    cout << iv0 << " = " << iv1 << " + " << ev1 << endl;
    iv0=ev1+iv1;
    cout << iv0 << " = " << ev1 << " + " << iv1 << endl;
    iv0=iv1-ev1;
    cout << iv0 << " = " << iv1 << " - " << ev1 << endl;
    iv0=ev1-iv1;
    cout << iv0 << " = " << ev1 << " - " << iv1 << endl;
    iv0=ex*iv1;
    cout << iv0 << " = " << ex << " * " << iv1 << endl;
    iv0=ix*ev1;
    cout << iv0 << " = " << ix << " * " << ev1 << endl;
    iv0=iv1*ex;
    cout << iv0 << " = " << iv1 << " * " << ex << endl;
    iv0=ev1*ix;
    cout << iv0 << " = " << ev1 << " * " << ix << endl;
    iv0=iv1/ex;
    cout << iv0 << " = " << iv1 << " / " << ex << endl;
    iv0=ev1/ix;
    cout << iv0 << " = " << ev1 << " / " << ix << endl;

    iv0=ev1;
    iv0/=ix;
    iv0=Vector<ValidatedFloat>({2,1});
    iv1=Vector<ValidatedFloat>({0,1});
    /*
      ARIADNE_TEST_ASSERT( (iv0+=Vector<ValidatedFloat>("[0,1]")) == Vector<ValidatedFloat>("[2,2]") );
      ARIADNE_TEST_ASSERT( (iv0-=Vector<ValidatedFloat>("[0,1]")) == Vector<ValidatedFloat>("[2,1]") );
      ARIADNE_TEST_ASSERT( (iv0*=2) == Vector<ValidatedFloat>("[4,2]") );
      ARIADNE_TEST_ASSERT( (iv0/=4) == Vector<ValidatedFloat>("[1,0.5]") );
    */

    /*
      cout << "test_vector_slice" << endl;
      v1=Vector<ApproximateFloat>("[-1.25,0.75,-0.5,-4.25,2.375]");
      cout << v1 << endl;
      VectorSlice<ApproximateFloat> vs1(2,v1.begin()+2,2);
      cout << vs1 << endl;
      VectorSlice<ApproximateFloat> vs2(2,v1.begin(),3);
      cout << vs2 << endl;

      iv1=vs1+vs2;
      cout << iv1 << endl;
      iv1=vs1-vs2;
      cout << iv1 << endl;
      iv1=x*vs1;
      cout << iv1 << endl;
      iv1=vs1*x;
      cout << iv1 << endl;
      iv1=vs1/x;
      cout << iv1 << endl;
    */

}


Int main() {
    TestVector().test();

    return ARIADNE_TEST_FAILURES;
}

