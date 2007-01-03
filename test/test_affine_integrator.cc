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

template<class R> int test_affine_integrator();

int main() {
  test_affine_integrator<Real>();
  return 0;
}


template<class R> 
int 
test_affine_integrator()
{
  cout << __PRETTY_FUNCTION__ << endl;
  
  {
    Matrix<R> T("[2]");
    Matrix<R> I("[1]");
    Vector<R> u("[1]");
    Matrix<R> A("[-2,-1;1,-2]");
    Vector<R> b("[0.125,0.25]");
    time_type h=0.125;
    uint k=1;
    R err(0.03125);
    std::cout << gexp(T,u,time_type(0.5),0u,err) << endl;
    std::cout << gexp(I,u,1,1u,err) << endl;
    std::cout << gexp(I,u,1,2u,err) << endl;
    std::cout << gexp(A,b,h,k,err) << endl;
  }
  
  Matrix<R> A("[-2,-1;1,-2]");
  Vector<R> b("[0.125,0.25]");
  time_type h=0.125;
  AffineIntegrator<R> affine(0.125,0.5,0.25);
  AffineVectorField<R> avf(A,b);
  Rectangle<R> bb("[-4,4]x[-4,4]");
  Rectangle<R> r("[-3.125,-2.875]x[-0.125,0.125]");
  Zonotope<R> z=r;
  
  Zonotope<R> iz1=affine.integration_step(avf,z,h);
  Zonotope<R> iz2=affine.integration_step(avf,iz1,h);
  Zonotope<R> rz1=affine.reachability_step(avf,z,h);
  Zonotope<R> rz2=affine.reachability_step(avf,iz1,h);
  
  if(h!=0.125) { cout << "h changed from 0.125 to " << h << endl; }
  
  epsfstream eps("test_affine_integrator.eps",bb);
  eps << rz1 << rz2;
  eps.set_fill_colour("blue");
  eps << iz1 << iz2;
  eps.set_fill_colour("yellow");
  eps << z;
  eps.close();
  
  cout << endl;
  
  return 0;
}
