/***************************************************************************
 *            test_vector.cpp
 *
 *  Copyright  2006-20  Pieter Collins, Alberto Casagrande
 *  Email Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
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
#include <fstream>
#include <cassert>

#include "config.hpp"
#include "numeric/numeric.hpp"
#include "numeric/float.hpp"
#include "algebra/vector.hpp"

#include "../test.hpp"

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
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_comparisons());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_misc());
}

Void
TestVector::test_concept()
{
    FloatDPApproximation ax(1);
    FloatDPBounds ix(1);
    FloatDPValue ex(1);
    Vector<FloatDPApproximation> av;
    Vector<FloatDPBounds> iv;
    Vector<FloatDPValue> ev;

    iv=Vector<FloatDPBounds>(ev);

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
    ARIADNE_TEST_CONSTRUCT( RawFloatDPVector, v0, (3) );
    ARIADNE_TEST_ASSERT( v0.size()==3 );
    ARIADNE_TEST_ASSERT( v0[0]==0.0 );
    ARIADNE_TEST_ASSERT( v0[1]==0.0 );
    ARIADNE_TEST_ASSERT( v0[2]==0.0 );

    ARIADNE_TEST_CONSTRUCT( RawFloatDPVector, v1, ({3.25,-0.75,0.0,1.375}) );
    ARIADNE_TEST_ASSERT( v1.size()==4 );
    ARIADNE_TEST_ASSERT( v1[0]==3.25 );
    ARIADNE_TEST_ASSERT( v1[1]==-0.75 );
    ARIADNE_TEST_ASSERT( v1[2]==0.0 );
    ARIADNE_TEST_ASSERT( v1[3]==1.375 );
}


Void
TestVector::test_comparisons()
{
    RawFloatDPVector v(2); v[0]=1.25; v[1]=-1.375;
    ARIADNE_TEST_PRINT( v );
    ARIADNE_TEST_COMPARE( v , ==, RawFloatDPVector({1.25,-1.375}) );
    ARIADNE_TEST_COMPARE( v , !=, RawFloatDPVector({1.25,-1.625}) );
}


Void
TestVector::test_arithmetic()
{
    ARIADNE_TEST_EQUAL( + RawFloatDPVector({2.0,-3.0,5.0}) , RawFloatDPVector({2.0,-3.0,5.0}) );
    ARIADNE_TEST_EQUAL( - RawFloatDPVector({2.0,-3.0,5.0}) , RawFloatDPVector({-2.0,+3.0,-5.0}) );
    ARIADNE_TEST_EQUAL( RawFloatDPVector({2.0,-3.0,5.0}) + RawFloatDPVector({1.25,2.75,-3.5}), RawFloatDPVector({3.25,-0.25,1.5}) );
    ARIADNE_TEST_EQUAL( RawFloatDPVector({2.0,-3.0,5.0}) - RawFloatDPVector({1.25,2.75,-3.5}), RawFloatDPVector({0.75,-5.75,8.5}) );
    ARIADNE_TEST_EQUAL( FloatDP(4.0) * RawFloatDPVector({2.0,-3.0,5.0}), RawFloatDPVector({8.0,-12.0,20.0}) );
    ARIADNE_TEST_EQUAL( RawFloatDPVector({2.0,-3.0,5.0}) * FloatDP(4.0), RawFloatDPVector({8.0,-12.0,20.0}) );
    ARIADNE_TEST_EQUAL( RawFloatDPVector({2.0,-3.0,5.0}) / FloatDP(4.0), RawFloatDPVector({0.5,-0.75,1.25}) );

    ARIADNE_TEST_EQUAL( sup_norm(RawFloatDPVector({2.0,-3.0,1.0})), FloatDP(3.0) )
}


Void
TestVector::test_misc()
{
    DoublePrecision pr;
    Array<FloatDPApproximation> vary={{-4.0,3.0,1.0},dp};
    FloatDPApproximation x={1.5,dp};

    Vector<FloatDPApproximation> v0;
    cout << "v0.size()=" << v0.size() << endl;
    cout << "v0=" << flush; cout << v0 << endl;
    ARIADNE_TEST_NOTIFY("Constructor Vector<X>(SizeType, const X*) is unsafe and has been removed.");
    Array<FloatDPApproximation> a1(vary.begin(),vary.end());
    Vector<FloatDPApproximation> v1(a1);
    cout << "v1=" << v1 << endl;
    Vector<FloatDPApproximation> v2=Vector<FloatDPApproximation>({2.375,4.25,-1.25},dp);
    cout << "v2=" << v2 << endl;
    cout << "norm(v1)=" << norm(v1) << "  norm(v2)=" << norm(v2) << endl;
    ARIADNE_TEST_EQUAL(norm(v1).raw(),4);
    ARIADNE_TEST_EQUAL(norm(v2).raw(),4.25);

    Vector<FloatDPApproximation> v3(1);
    cout << "v3=" << v3 << endl;
    Vector<FloatDPApproximation> v4=v2;
    cout << "v4=" << v4 << endl;
    Vector<FloatDPApproximation> v5={{-4.0,3.0,1.0},dp};
    cout << "v5=" << v5 << endl;
    ARIADNE_TEST_EQUAL(v1,v5);
    cout << endl;

    Vector<FloatDPApproximation> vf0;
    v1=Vector<FloatDPApproximation>({0.25,-1.5},dp);
    v2=Vector<FloatDPApproximation>({-0.5,2.25},dp);
    vf0=-v1;
    cout << vf0 << " = -" << v1 << endl;
    vf0=Vector<FloatDPApproximation>(v1)+v2;
    cout << vf0 << " = " << v1 << " + " << v2 << endl;
    vf0=Vector<FloatDPApproximation>(v1)-v2;
    cout << vf0 << " = " << v1 << " - " << v2 << endl;
    vf0=x*Vector<FloatDPApproximation>(v2);
    cout << vf0 << " = " << x << " * " << v2 << endl;
    vf0=Vector<FloatDPApproximation>(v1)*x;
    cout << vf0 << " = " << v1 << " * " << x << endl;
    vf0=Vector<FloatDPApproximation>(v1)/x;
    cout << vf0 << " = " << v1 << " / " << x << endl;
    cout << endl;

    Vector< FloatDPBounds > iv1=Vector<FloatDPBounds>({FloatDPBounds{0.984375_x,1.015625_x,pr},{2.25_x,2.375_x,pr},{4.0_x,4.375_x,pr},{-0.03125_x,0.015625_x,pr}});
    cout << "iv1=" << iv1 << endl;
    cout << "norm(iv1)=" << norm(iv1) << endl;
    cout << "norm(iv1).upper()=" << norm(iv1).upper() << endl;

    Vector< FloatDPBounds > iv2=Vector<FloatDPBounds>({{-1.0_x,1.0_x,pr},{-1.0_x,1.0_x,pr}});
    cout << "iv2=" << iv2 << endl;
    Vector< FloatDPBounds > iv3(3);
    cout << "iv3=" << iv3 << endl;
    iv3=Vector<FloatDPBounds>({{4.25_x,4.25_x,pr},{2.375_x,2.375_x,pr}});
    cout << "iv3=" << iv3 << endl;
    FloatDPBounds ix=FloatDPBounds(-2,1,pr);

    Vector< FloatDPBounds > iv0;
    cout << "iv0=" << iv0 << endl;
    iv1=iv0;
    cout << "iv1=" << iv1 << endl;
    iv1=iv2;
    cout << "iv1=" << iv1 << endl;
    cout << endl;

    FloatDPBounds ix2=iv2[0];
    FloatDPBounds ix3=iv3[0];
    FloatDPBounds ix1=ix2+ix3;
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
    ix=FloatDPBounds(1,2);
    iv1=iv2/ix;
    cout << iv1 << " = " << iv2 << " / " << ix << endl;
    cout << endl;

    Vector<FloatDPValue> ev1(reinterpret_cast<Vector<FloatDPValue>const&>(v1));
    FloatDPValue ex(reinterpret_cast<FloatDPValue const&>(x));
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
    iv0=Vector<FloatDPBounds>({2,1},dp);
    iv1=Vector<FloatDPBounds>({0,1},dp);
    /*
      ARIADNE_TEST_ASSERT( (iv0+=Vector<FloatDPBounds>("[0,1]")) == Vector<FloatDPBounds>("[2,2]") );
      ARIADNE_TEST_ASSERT( (iv0-=Vector<FloatDPBounds>("[0,1]")) == Vector<FloatDPBounds>("[2,1]") );
      ARIADNE_TEST_ASSERT( (iv0*=2) == Vector<FloatDPBounds>("[4,2]") );
      ARIADNE_TEST_ASSERT( (iv0/=4) == Vector<FloatDPBounds>("[1,0.5]") );
    */

    /*
      cout << "test_vector_slice" << endl;
      v1=Vector<FloatDPApproximation>("[-1.25,0.75,-0.5,-4.25,2.375]");
      cout << v1 << endl;
      VectorSlice<FloatDPApproximation> vs1(2,v1.begin()+2,2);
      cout << vs1 << endl;
      VectorSlice<FloatDPApproximation> vs2(2,v1.begin(),3);
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

