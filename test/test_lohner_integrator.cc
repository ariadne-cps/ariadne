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

#include "real_typedef.h"

#include "ariadne.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "system/affine_vector_field.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"
#include "output/epsfstream.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_lohner_integrator();

int main() {
  test_lohner_integrator<Real>();
  return 0;
}


template<class R> 
int 
test_lohner_integrator()
{
  cout << __PRETTY_FUNCTION__ << endl;
  
  LohnerIntegrator<R> lohner=LohnerIntegrator<R>(0.125,0.5,0.0625);
  Rectangle<R> bb=Rectangle<R>("[0.50,1.25]x[0.25,1.00]");
  Rectangle<R> r=Rectangle<R>("[0.96,1.04]x[0.46,0.54]");
  cout << "r=" << r << endl;
  Zonotope<R> z(r);
  cout << "z=" << z << endl;
  Matrix<R> A=Matrix<R>("[-0.25,-1;+1,-0.25]");
  cout << "A=" << A << endl;
  Vector<R> b=Vector<R>("[0,0]");
  cout << "b=" << b << endl;
  AffineVectorField<R> avf=AffineVectorField<R>(A,b);
  cout << "avf=" << avf << endl;

  time_type h(0.125);
  cout << "h=" << h << endl;
  cout << "z.generators().norm()=" << norm(z.generators()) << endl;

  Zonotope< Interval<R> > z0=z;
  Zonotope< Interval<R> > z1=lohner.integration_step(avf,z0,h);
  Zonotope< Interval<R> > z2=lohner.integration_step(avf,z1,h);
  cout << "z0=" << z0 << "\np1=" << z1 << "\nz2=" << z2 << endl;
  Zonotope< Interval<R> > zr1=lohner.reachability_step(avf,z0,h);
  Zonotope< Interval<R> > zr2=lohner.reachability_step(avf,z1,h);
  cout << "zr1=" << zr1 << "\nzr2=" << zr2 << endl;
  
  epsfstream eps;
  eps.open("test_lohner_integrator.eps",bb);
  eps.set_fill_colour("green");
  eps << over_approximation(zr1) << over_approximation(zr2);
  eps.set_fill_colour("blue");
  eps << over_approximation(zr1) << 2;
  eps.set_fill_colour("yellow");
  eps << over_approximation(z0);
  eps.close();
  
  return 0;
}
