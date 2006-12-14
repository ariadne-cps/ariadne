/***************************************************************************
 *            test_integration_step.cc
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
#include "numeric/numerical_types.h"
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
template<class R> int test_lohner_integrator();

int main() {
  test_lohner_integrator<Real>();
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
  
  epsfstream eps("test_reachability_step.eps",bb);
  eps << rz1 << rz2;
  eps.set_fill_colour("blue");
  eps << iz1 << iz2;
  eps.set_fill_colour("yellow");
  eps << z;
  eps.close();
  
  cout << endl;
  
  return 0;
}


template<class R> 
int 
test_lohner_integrator()
{
  cout << __PRETTY_FUNCTION__ << endl;
  
  LohnerIntegrator<R> lohner=LohnerIntegrator<R>(0.125,0.5,0.0625);
  Rectangle<R> r=Rectangle<R>("[0.96,1.04]x[0.46,0.54]");
  cout << "r=" << r << endl;
  Parallelotope<R> p=Parallelotope<R>(r);
  cout << "p=" << p << endl;
  Matrix<R> A=Matrix<R>("[-0.25,-1;+1,-0.25]");
  cout << "A=" << A << endl;
  Vector<R> b=Vector<R>("[0,0]");
  cout << "b=" << b << endl;
  AffineVectorField<R> avf=AffineVectorField<R>(A,b);
  cout << "avf=" << avf << endl;

  Real x0=0;
  Real x1=0.4;
  Interval<Real> ivl1(0.4);
  Interval<Real> ivl0;
  
  Matrix<Real> fA(2,2);
  fA(0,0)=0.4;
  fA(1,1)=0.4;
  R z=0;
  Interval<Real> ivlm(z,z);
  Interval<Real> ivls(z,z);
  ivls+=abs(fA(0,0));
  cout << fA(0,0) << "  " << ivls << "  " << ivlm << "\n";
  ivls+=abs(fA(0,1));
  ivlm=Numeric::max(ivlm,ivls);
  cout << fA(0,1) << "  " << ivls << "  " << ivlm << "\n";
  ivls=Interval<Real>(0);
  ivls+=abs(fA(1,0));
  cout << fA(1,0) << "  " << ivls << "  " << ivlm << "\n";
  ivls+=abs(fA(1,1));
  ivlm=Numeric::max(ivlm,ivls);
  cout << fA(1,1) << "  " << ivls << "  " << ivlm << "\n";
  
  time_type h(0.125);
  cout << "h=" << h << endl;
  cout << "p.generators().norm()=" << norm(p.generators()) << endl;

  Parallelotope<R> next_paral=lohner.integration_step(avf,p,h);
  cout << next_paral << "\n\n";
  
  return 0;
}
