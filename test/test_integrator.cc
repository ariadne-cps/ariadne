/***************************************************************************
 *            test_integrator.cc
 *
 *  Copyright  2006-8  Pieter Collins
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

#include <fstream>

#include "test_float.h"

#include "ariadne.h"
#include "base/tuple.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "system/affine_vector_field.h"
#include "evaluation/standard_bounder.h"
#include "evaluation/standard_integrator.h"
#include "output/epsstream.h"
#include "output/logging.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Output;
using namespace std;

template<class R> 
class TestIntegrator
{
 public:
  void test() const;
};

int main() {
  TestIntegrator<Flt>().test();
  return ARIADNE_TEST_FAILURES;
}


template<class R> 
void TestIntegrator<R>::test() const
{
  cout << __PRETTY_FUNCTION__ << endl;
  typedef Interval<R> I;

  // R maximum_step_size=0.125;
  R maximum_step_size=0.0625;
  StandardBounder<R> bounder(maximum_step_size);
  StandardIntegrator< Zonotope<R> > standard;

  Box<R> bb=Box<R>("[0.25,1.25]x[0.00,1.00]");
  Box<R> r=Box<R>("[0.96,1.04]x[0.46,0.54]");
  cout << "r=" << r << endl;
  Zonotope<R> z(r);
  cout << "z=" << z << endl;
  Matrix<R> A=Matrix<R>("[-0.25,-1;+1,-0.25]");
  //Matrix<R> A=Matrix<R>("[-0.5,-1;+1,-0.5]");
  cout << "A=" << A << endl;
  Vector<R> b=Vector<R>("[0,0]");
  cout << "b=" << b << endl;
  AffineVectorField<R> avf=AffineVectorField<R>(A,b);
  cout << "avf=" << avf << endl;

  Rational qh(0.125);
  Interval<R> ih(0.125);
  Rational h=qh;
  cout << "h=" << h << endl;
  cout << "z.generators().norm()=" << norm(z.generators()) << endl;


  Rational t0,t1,t2,t3,t4;
  Rational h0,h1,h2,h3,h4;
  Box<R> bb0,bb1,bb2,bb3,bb4;
  Zonotope<R> z0,z1,z2,z3,z4,zr1,zr2,zr3,zr4;
  Zonotope<R> c0z, c1z,afz;
  t0=0;
  z0=z;

  make_lpair(h,bb0)=standard.flow_bounds(avf,z0.bounding_box(),qh);
  t1=t0+h;
  z1=standard.integration_step(avf,z0,h,bb0);
  make_lpair(h,bb1)=standard.flow_bounds(avf,z1.bounding_box(),qh);
  t2=t1+h;
  z2=standard.integration_step(avf,z1,h,bb1);
  make_lpair(h,bb2)=standard.flow_bounds(avf,z2.bounding_box(),qh);
  t3=t2+h;
  z3=standard.integration_step(avf,z2,h,bb2);
  make_lpair(h,bb3)=standard.flow_bounds(avf,z3.bounding_box(),qh);
  t4=t3+h;
  z4=standard.integration_step(avf,z3,h,bb3);
  cout << "t0=" << t0 << " z0=" << z0 << "\n"
       << "t1=" << t1 << " z1=" << z1 << "\n"
       << "t2=" << t2 << " z2=" << z2 << "\n"
       << "t3=" << t3 << " z3=" << z3 << "\n"
       << "t4=" << t4 << " z4=" << z4 << "\n"
       << endl;
  zr1=standard.reachability_step(avf,z0,h,bb0);
  zr2=standard.reachability_step(avf,z1,h,bb1);
  zr3=standard.reachability_step(avf,z2,h,bb2);
  zr4=standard.reachability_step(avf,z3,h,bb3);
  cout << "zr1=" << zr1 << "\nzr2=" << zr2 << "\n"
       << "zr3=" << zr3 << "\nzr4=" << zr4 << "\n" << endl;
  c0z=z4;


  epsfstream eps;
  eps.open("test_integrator-nonlinear.eps",bb);
  eps << fill_colour(green)
      << zr1 << zr2
      << zr3 << zr4;
  eps << fill_colour(blue)
      << z1<< z2
      << z3 << z4;
  eps << fill_colour(yellow)
      << z0;
  eps.close();


  Point<I> pt0 = z0.centre();
  make_lpair(h,bb0)=standard.flow_bounds(avf,Box<R>(pt0),qh);

  Box<R> rbb0=bounder.refine_flow_bounds(avf,Box<R>(pt0),bb0,qh);
  Box<R> rrbb0=bounder.refine_flow_bounds(avf,Box<R>(pt0),rbb0,qh);

  eps.open("test_integrator-affine.eps",bb);
  eps << fill_colour(green)
      << zr1 << zr2
      << zr3 << zr4;
  eps << fill_colour(blue)
      << z1 << z2
      << z3 << z4;
  eps << fill_colour(yellow)
      << z0;
  eps.close();
  
}
