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
#include "numeric.h"
#include "vector.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

class TestVector {
  public:
    void test();
  private:
    void test_concept();
    void test_misc();
};

void
TestVector::test()
{
    ARIADNE_TEST_CALL(test_concept());
}

void
TestVector::test_concept()
{
    Float fx;
    Interval ix;
    Vector<Float> fv;
    Vector<Interval> iv;

    fv=fv+fv; iv=fv+fv; iv=fv+iv; iv=iv+fv; iv=iv+iv;
    fv=fv-fv; iv=fv-fv; iv=fv-iv; iv=iv-fv; iv=iv-iv;

    fv=fx*fv; iv=fx*fv; iv=fx*iv; iv=ix*fv; iv=ix*iv;
    fv=fv*fx; iv=fv*fx; iv=fv*ix; iv=iv*fx; iv=iv*ix;
    fv=fv/fx; iv=fv/fx; iv=fv/ix; iv=iv/fx; iv=iv/ix;
}


void
TestVector::test_misc()
{

    int n=3;
    Float vptr[3]={-4.0,3.0,1.0};
    Float x=1.5;

    Vector<Float> v0;
    cout << "v0.size()=" << v0.size() << endl;
    cout << "v0=" << flush; cout << v0 << endl;
    Vector<Float> v1(n,vptr);
    cout << "v1=" << v1 << endl;
    Vector<Float> v2=Vector<Float>("[2.375,4.25,-1.25]");
    cout << "v2=" << v2 << endl;
    cout << "norm(v1)=" << norm(v1) << "  norm(v2)=" << norm(v2) << endl;
    assert(norm(v1)==4);
    assert(norm(v2)==4.25);

    Vector<Float> v3(1);
    cout << "v3=" << v3 << endl;
    Vector<Float> v4=v2;
    cout << "v4=" << v4 << endl;
    cout << endl;

    Vector<Float> vf0;
    v1=Vector<Float>("[0.25,-1.5]");
    v2=Vector<Float>("[-0.5,2.25]");
    vf0=-v1;
    cout << vf0 << " = -" << v1 << endl;
    vf0=Vector<Float>(v1)+v2;
    cout << vf0 << " = " << v1 << " + " << v2 << endl;
    vf0=Vector<Float>(v1)-v2;
    cout << vf0 << " = " << v1 << " - " << v2 << endl;
    vf0=x*Vector<Float>(v2);
    cout << vf0 << " = " << x << " * " << v2 << endl;
    vf0=Vector<Float>(v1)*x;
    cout << vf0 << " = " << v1 << " * " << x << endl;
    vf0=Vector<Float>(v1)/x;
    cout << vf0 << " = " << v1 << " / " << x << endl;
    cout << endl;

    Vector< Interval > iv1=Vector<Interval>("[[0.984375,1.015625],[2.25,2.375],[4.0,4.375],[-0.03125,0.015625]]");
    cout << "iv1=" << iv1 << endl;
    cout << "norm(iv1)=" << norm(iv1) << endl;
    cout << "norm(iv1).upper()=" << norm(iv1).upper() << endl;

    Vector< Interval > iv2=Vector<Interval>("[[-1,1],[-1,1]]");
    cout << "iv2=" << iv2 << endl;
    Vector< Interval > iv3(3);
    cout << "iv3=" << iv3 << endl;
    iv3=Vector<Interval>("[[4.25,4.25],[2.375,2.375]]");
    cout << "iv3=" << iv3 << endl;
    Interval ix=Interval(-2,1);

    Vector< Interval > iv0;
    cout << "iv0=" << iv0 << endl;
    iv1=iv0;
    cout << "iv1=" << iv1 << endl;
    iv1=iv2;
    cout << "iv1=" << iv1 << endl;
    cout << endl;

    Interval ix2=iv2(0);
    Interval ix3=iv3(0);
    Interval ix1=ix2+ix3;
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
    ix=Interval(1,2);
    iv1=iv2/ix;
    cout << iv1 << " = " << iv2 << " / " << ix << endl;
    cout << endl;

    v1==Vector<Float>("[-1.25,0.75]");
    iv0=iv1+v1;
    cout << iv0 << " = " << iv1 << " + " << v1 << endl;
    iv0=v1+iv1;
    cout << iv0 << " = " << v1 << " + " << iv1 << endl;
    iv0=iv1-v1;
    cout << iv0 << " = " << iv1 << " - " << v1 << endl;
    iv0=v1-iv1;
    cout << iv0 << " = " << v1 << " - " << iv1 << endl;
    iv0=x*iv1;
    cout << iv0 << " = " << x << " * " << iv1 << endl;
    iv0=ix*v1;
    cout << iv0 << " = " << ix << " * " << v1 << endl;
    iv0=iv1*x;
    cout << iv0 << " = " << iv1 << " * " << x << endl;
    iv0=v1*ix;
    cout << iv0 << " = " << v1 << " * " << ix << endl;
    iv0=iv1/x;
    cout << iv0 << " = " << iv1 << " / " << x << endl;
    iv0=v1/ix;
    cout << iv0 << " = " << v1 << " / " << ix << endl;

    iv0=v1;
    iv0/=ix;
    iv0=Vector<Interval>("[2,1]");
    iv1=Vector<Interval>("[0,1]");
    /*
      ARIADNE_TEST_ASSERT( (iv0+=Vector<Interval>("[0,1]")) == Vector<Interval>("[2,2]") );
      ARIADNE_TEST_ASSERT( (iv0-=Vector<Interval>("[0,1]")) == Vector<Interval>("[2,1]") );
      ARIADNE_TEST_ASSERT( (iv0*=2) == Vector<Interval>("[4,2]") );
      ARIADNE_TEST_ASSERT( (iv0/=4) == Vector<Interval>("[1,0.5]") );
    */

    /*
      cout << "test_vector_slice" << endl;
      v1=Vector<Float>("[-1.25,0.75,-0.5,-4.25,2.375]");
      cout << v1 << endl;
      VectorSlice<Float> vs1(2,v1.begin()+2,2);
      cout << vs1 << endl;
      VectorSlice<Float> vs2(2,v1.begin(),3);
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


int main() {
    TestVector().test();

    return ARIADNE_TEST_FAILURES;
}

