/***************************************************************************
 *            test_reducer.cc
 *
 *  Copyright  2008  Pieter Collins
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

#include <iostream>

#include "test_float.h"

#include "base/pointer.h"
#include "linear_algebra/matrix.h"
#include "geometry/point.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "evaluation/identity_reducer.h"
#include "evaluation/orthogonal_reducer.h"
#include "evaluation/cascade_reducer.h"
#include "output/epsstream.h"
#include "output/logging.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Output;
using namespace std;

template<class R> 
class TestReducer
{
 public:
  void test() const;
};

template<class R> 
void TestReducer<R>::test() const
{
    
  set_evaluation_verbosity(0);
  typedef Interval<R> I;
  typedef Zonotope<R> ZBS;
  typedef Polytope<R> PBS;

  IdentityReducer< Zonotope<R> > identity_reducer;
  OrthogonalReducer< Zonotope<R> > orthogonal_reducer;
  CascadeReducer< Zonotope<R> > cascade_reducer(2);

  Zonotope<R> z1(Point<R>("[0,0]"),Matrix<R>("[2,1,1,1;1,1,0,3]"));
  Zonotope<R> z2(Point<R>("[0,0]"),Matrix<R>("[2,1,1,1,2,0;1,1,0.5,3,0,2]"));
 
  Zonotope<R> oz1=orthogonal_reducer.over_approximate(z1);
  cout << "z1=" << z1 << " oz1=" << oz1 << endl;

  Zonotope<R> cz2=cascade_reducer.over_approximate(z2);
  cout << "z2=" << z2 << " cz2=" << cz2 << endl;
  
  epsfstream eps;
  eps.open("test_reducer-orthogonal.eps",Box<R>("[-4,4]x[-4,4]"));
  eps << oz1 << z1;
  eps.close();
  eps.open("test_reducer-cascade.eps",Box<R>("[-4,4]x[-4,4]"));
  eps << cz2 << z2; 
  eps.close();
  
}

int main() {
  TestReducer<Flt>().test();
  return ARIADNE_TEST_FAILURES;
}
