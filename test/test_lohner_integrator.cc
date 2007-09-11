/***************************************************************************
 *            test_lohner_integrator.cc
 *
 *  Copyright  2006  Pieter Collins
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

#include <fstream>

#include "test_float.h"

#include "ariadne.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "system/affine_vector_field.h"
#include "evaluation/bounder_plugin.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"
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

template<class R> int test_lohner_integrator();

int main() {
  test_lohner_integrator<Float>();
  return 0;
}


template<class R> 
int 
test_lohner_integrator()
{
  // set_integrator_verbosity(8);

  cout << __PRETTY_FUNCTION__ << endl;
  typedef Interval<R> I;

  BounderPlugin<R> bounder;
  LohnerIntegrator<R> lohner=LohnerIntegrator<R>();
  C1LohnerIntegrator<R> c1lohner=C1LohnerIntegrator<R>();
  AffineIntegrator<R> affine=AffineIntegrator<R>();

  Rectangle<R> bb=Rectangle<R>("[0.25,1.25]x[0.00,1.00]");
  Rectangle<R> r=Rectangle<R>("[0.96,1.04]x[0.46,0.54]");
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
  Interval<R> h(0.125);
  cout << "h=" << h << endl;
  cout << "z.generators().norm()=" << norm(z.generators()) << endl;

  const VectorFieldInterface<R>& vf=avf;

  Rectangle<R> bb0,bb1,bb2,bb3,bb4;
  Zonotope<I,I> z0,z1,z2,z3,z4,zr1,zr2,zr3,zr4;
  Zonotope<I,I> c0z, c1z,afz;
  z0=z;

  bb0=bounder.estimate_flow_bounds(avf,z0.bounding_box(),qh);
  z1=lohner.integration_step(avf,z0,h,bb0);
  bb1=bounder.estimate_flow_bounds(avf,z1.bounding_box(),qh);
  z2=lohner.integration_step(avf,z1,h,bb1);
  bb2=bounder.estimate_flow_bounds(avf,z2.bounding_box(),qh);
  z3=lohner.integration_step(avf,z2,h,bb2);
  bb3=bounder.estimate_flow_bounds(avf,z3.bounding_box(),qh);
  z4=lohner.integration_step(avf,z3,h,bb3);
  cout << "z0=" << z0 << "\n"
       << "z1=" << z1 << "\nz4=" << z2 << "\n"
       << "z3=" << z3 << "\nz4=" << z4 << endl;
  zr1=lohner.reachability_step(avf,z0,h,bb0);
  zr2=lohner.reachability_step(avf,z1,h,bb1);
  zr3=lohner.reachability_step(avf,z2,h,bb2);
  zr4=lohner.reachability_step(avf,z3,h,bb3);
  cout << "zr1=" << zr1 << "\nzr2=" << zr2 << "\n"
       << "zr3=" << zr3 << "\nzr4=" << zr4 << "\n" << endl;
  c0z=z4;

  epsfstream eps;
  eps.open("test_lohner_integrator-1.eps",bb);
  eps << fill_colour(green)
      << over_approximation(zr1) << over_approximation(zr2)
      << over_approximation(zr3) << over_approximation(zr4);
  eps << fill_colour(blue)
      << over_approximation(z1) << over_approximation(z2)
      << over_approximation(z3) << over_approximation(z4);
  eps << fill_colour(yellow)
      << over_approximation(z0);
  eps.close();


  cout << "\nC1LohnerIntegrator\n";
  bb0=bounder.estimate_flow_bounds(avf,z0.bounding_box(),qh);
  z1=c1lohner.integration_step(avf,z0,h,bb0);
  bb1=bounder.estimate_flow_bounds(avf,z1.bounding_box(),qh);
  z2=c1lohner.integration_step(avf,z1,h,bb1);
  bb2=bounder.estimate_flow_bounds(avf,z2.bounding_box(),qh);
  z3=c1lohner.integration_step(avf,z2,h,bb2);
  bb3=bounder.estimate_flow_bounds(avf,z3.bounding_box(),qh);
  z4=c1lohner.integration_step(avf,z3,h,bb3);
  cout << "z0=" << z0 << "\n"
       << "z1=" << z1 << "\nz4=" << z2 << "\n"
       << "z3=" << z3 << "\nz4=" << z4 << endl;
  zr1=c1lohner.reachability_step(avf,z0,h,bb0);
  zr2=c1lohner.reachability_step(avf,z1,h,bb1);
  zr3=c1lohner.reachability_step(avf,z2,h,bb2);
  zr4=c1lohner.reachability_step(avf,z3,h,bb3);
  cout << "zr1=" << zr1 << "\nzr2=" << zr2 << "\n"
       << "zr3=" << zr3 << "\nzr4=" << zr4 << "\n" << endl;
  c1z=z4;
  
  Point<I> pt0 = z0.centre();
  bb0=bounder.estimate_flow_bounds(avf,Rectangle<R>(pt0),qh);
  Rectangle<R> rbb0=bounder.refine_flow_bounds(avf,Rectangle<R>(pt0),bb0,qh);
  Rectangle<R> rrbb0=bounder.refine_flow_bounds(avf,Rectangle<R>(pt0),rbb0,qh);
  Point<I> pt1 = c1lohner.flow_step(avf,pt0,h,bb0);
  Point<I> rpt1 = c1lohner.flow_step(avf,pt0,h,rbb0);
  Point<I> rrpt1 = c1lohner.flow_step(avf,pt0,h,rrbb0);
  Matrix<I> mx1 = c1lohner.flow_step_jacobian(avf,pt0,h,bb0);
  Matrix<I> rmx1 = c1lohner.flow_step_jacobian(avf,pt0,h,rbb0);
  cout << "pt0=" << pt0 << "\n"
       << "bb0=" << bb0 << ", rbb0=" << rbb0 << ", rrbb0=" << rrbb0 << "\n"
       << "pt1=" << pt1 << ", rpt1=" << rpt1 << ", rrpt1=" << rrpt1 << "\n"
       << "mx1=" << mx1 << ", rmx1=" << rmx1 << "\n" << endl;
  
  eps.open("test_lohner_integrator-2.eps",bb);
  eps << fill_colour(green)
      << over_approximation(zr1) << over_approximation(zr2)
      << over_approximation(zr3) << over_approximation(zr4);
  eps << fill_colour(blue)
      << over_approximation(z1) << over_approximation(z2)
      << over_approximation(z3) << over_approximation(z4);
  eps << fill_colour(yellow)
      << over_approximation(z0);
  eps.close();
  
  cout << "\nAffineIntegrator\n";
  z1=affine.integration_step(avf,z0,h);
  z2=affine.integration_step(avf,z1,h);
  z3=affine.integration_step(avf,z2,h);
  z4=affine.integration_step(avf,z3,h);
  cout << "z0=" << z0 << "\n"
       << "z1=" << z1 << "\nz4=" << z2 << "\n"
       << "z3=" << z3 << "\nz4=" << z4 << endl;
  zr1=affine.reachability_step(avf,z0,h);
  zr2=affine.reachability_step(avf,z1,h);
  zr3=affine.reachability_step(avf,z2,h);
  zr4=affine.reachability_step(avf,z3,h);
  cout << "zr1=" << zr1 << "\nzr2=" << zr2 << "\n"
       << "zr3=" << zr3 << "\nzr4=" << zr4 << "\n" << endl;
  afz=z4;
  
  eps.open("test_lohner_integrator-3.eps",bb);
  eps << fill_colour(green)
      << over_approximation(zr1) << over_approximation(zr2)
      << over_approximation(zr3) << over_approximation(zr4);
  eps << fill_colour(blue)
      << over_approximation(z1) << over_approximation(z2)
      << over_approximation(z3) << over_approximation(z4);
  eps << fill_colour(yellow)
      << over_approximation(z0);
  eps.close();
  
  cout << endl;
  cout << "c0z="<<c0z<<"\n"
       << "c1z="<<c1z<<"\n"
       << "afz="<<afz<<"\n";
  cout << "c0z.bounding_box()="<<c0z.bounding_box()<<"\n"
       << "c1z.bounding_box()="<<c1z.bounding_box()<<"\n"
       << "afz.bounding_box()="<<afz.bounding_box()<<"\n";

  return 0;
}
