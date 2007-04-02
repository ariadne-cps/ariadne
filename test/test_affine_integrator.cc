/***************************************************************************
 *            test_affine_integrator.cc
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
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"
#include "output/epsfstream.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_affine_integrator();

int main() {
  test_affine_integrator<Float>();
  return 0;
}


template<class R> 
int 
test_affine_integrator()
{
  cout << __PRETTY_FUNCTION__ << endl;
  typedef Interval<R> I;

  {
    Matrix<R> T("[2]");
    Matrix<R> I("[1]");
    Vector<R> u("[1]");
    Matrix<R> A("[-2,-1;1,-2]");
    Vector<R> b("[0.125,0.25]");
    time_type h=0.125;
    uint k=1;
    std::cout << gexp(T,u,time_type(0.5),0u) << endl;
    std::cout << gexp(I,u,1,1u) << endl;
    std::cout << gexp(I,u,1,2u) << endl;
    std::cout << gexp(A,b,h,k) << endl;
  }
  
  Matrix<R> A("[-2,-1;1,-2]");
  Vector<R> b("[0.125,0.25]");
  time_type h=0.125;
  time_type hh=h/2;
  time_type th=h*2;
  AffineIntegrator<R> affine(0.125,0.5,0.25);
  AffineVectorField<R> avf(A,b);
  Rectangle<R> bb("[-4,0]x[-2,2]");
  Rectangle<R> r("[-3.125,-2.875]x[-0.125,0.125]");
  Zonotope<I> iz; iz=r;
  
  Zonotope<I> iz1=affine.integration_step(avf,iz,h);
  Zonotope<I> iz2=affine.integration_step(avf,iz1,h);
  Zonotope<I> iz3=affine.integration_step(avf,iz2,h);
  Zonotope<I> iz4=affine.integration_step(avf,iz3,h);
  Zonotope<I> hiz1=affine.integration_step(avf,iz,hh);
  Zonotope<I> hiz2=affine.integration_step(avf,iz1,hh);
  Zonotope<I> hiz3=affine.integration_step(avf,iz2,hh);
  Zonotope<I> hiz4=affine.integration_step(avf,iz3,hh);
  Zonotope<I> riz1=affine.reachability_step(avf,iz,h);
  Zonotope<I> riz2=affine.reachability_step(avf,iz1,h);
  Zonotope<I> riz3=affine.reachability_step(avf,iz2,h);
  Zonotope<I> riz4=affine.reachability_step(avf,iz3,h);
  
  if(h!=0.125) { cout << "h changed from 0.125 to " << h << endl; }
  
  epsfstream eps;
  eps.open("test_affine_integrator.eps",bb);
  eps.set_fill_colour("red");
  eps << over_approximation(riz1) << over_approximation(riz2) << over_approximation(riz3) << over_approximation(riz4);
  eps.set_fill_colour("green");
  eps << approximation(riz1) << approximation(riz2) << approximation(riz3) << approximation(riz4);
  eps.set_fill_colour("magenta");
  eps << approximation(hiz1) << approximation(hiz2) << approximation(hiz3) << approximation(hiz4);
  eps.set_fill_colour("blue");
  eps << approximation(iz1) << approximation(iz2) << approximation(iz3) << approximation(iz4);
  eps.set_fill_colour("yellow");
  eps << over_approximation(iz);
  eps.close();
  
  cout << endl;
  
  return 0;
}
