/***************************************************************************
 *            test_affine_integrator.cc
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
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "system/affine_vector_field.h"
#include "evaluation/affine_integrator.h"
#include "output/epsstream.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

template<class R> 
class TestAffineIntegrator
{
 public:
  void test() const;
};

int main() {
  TestAffineIntegrator<Flt>().test();
  return ARIADNE_TEST_FAILURES;
}


template<class R> 
void TestAffineIntegrator<R>::test() const
{
  cout << __PRETTY_FUNCTION__ << endl;
  typedef Interval<R> I;

  {
    Matrix<I> T("[2]");
    Matrix<I> Id("[1]");
    Vector<I> u("[1]");
    Matrix<I> A("[-2,-1;1,-2]");
    Vector<I> b("[0.125,0.25]");
    time_type h=0.125;
    uint k=1;
    std::cout << gexp(T,u,Interval<R>(0.5),0u) << endl;
    std::cout << gexp(Id,u,Interval<R>(1),1u) << endl;
    std::cout << gexp(Id,u,Interval<R>(1),2u) << endl;
    std::cout << gexp(A,b,Interval<R>(h),k) << endl;
  }
  
  Matrix<R> A("[-2,-1;1,-2]");
  Vector<R> b("[0.125,0.25]");
  time_type h=0.125;
  time_type hh=h/2;
  time_type th=h*2;
  AffineIntegrator< Zonotope<R> > affine;
  AffineVectorField<R> avf(A,b);
  Box<R> bb("[-4,0]x[-2,2]");
  Box<R> r("[-3.125,-2.875]x[-0.125,0.125]");
  cout << "r0=" << r << endl;
  Zonotope<R> iz; iz=r;
  cout << "iz=" << iz << endl;
  Zonotope<R> iz0=iz;
  cout << "iz0=" << iz0 << endl;
  
  Zonotope<R> iz1=affine.integration_step(avf,iz0,h);
  cout << "iz1=" << iz1 << endl;
  Zonotope<R> iz2=affine.integration_step(avf,iz1,h);
  cout << "iz2=" << iz2 << endl;
  Zonotope<R> iz3=affine.integration_step(avf,iz2,h);
  cout << "iz3=" << iz3 << endl;
  Zonotope<R> iz4=affine.integration_step(avf,iz3,h);
  cout << "iz4=" << iz4 << endl;
  Zonotope<R> hiz1=affine.integration_step(avf,iz0,hh);
  Zonotope<R> hiz2=affine.integration_step(avf,iz1,hh);
  Zonotope<R> hiz3=affine.integration_step(avf,iz2,hh);
  Zonotope<R> hiz4=affine.integration_step(avf,iz3,hh);
  cout << "hiz4=" << hiz4 << endl;
  Zonotope<R> riz1=affine.reachability_step(avf,iz0,h);
  cout << "riz1=" << riz1 << endl;
  Zonotope<R> riz2=affine.reachability_step(avf,iz1,h);
  Zonotope<R> riz3=affine.reachability_step(avf,iz2,h);
  Zonotope<R> riz4=affine.reachability_step(avf,iz3,h);
  cout << "riz4=" << riz4 << endl;
  
  if(h!=0.125) { cout << "h changed from 0.125 to " << h << endl; }
  
  epsfstream eps;
  eps.open("test_affine_integrator.eps",bb);
  eps << fill_colour(red) << riz1 << riz2 << riz3 << riz4;
  eps << fill_colour(green) << riz1 << riz2 << riz3 << riz4;
  eps << fill_colour(magenta) << hiz1 << hiz2 << hiz3 << hiz4;
  eps << fill_colour(blue) << iz1 << iz2 << iz3 << iz4;
  eps << fill_colour(yellow) << iz;
  eps.close();
  
  cout << endl;
  
}
