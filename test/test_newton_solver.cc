/***************************************************************************
 *            test_interval_newton.cc
 *
 *  Copyright 2006  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cassert>


#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "evaluation/exceptions.h"
#include "evaluation/newton.h"

#include "models/henon.h"
#include "test.h"
#include "test_float.h"

using namespace std;
using namespace Ariadne;
using Ariadne::Models::HenonMap;

template<class R> int test_newton();

int main() {
  test_newton<Flt>();
  return 0;
}

template<class R>
int 
test_newton()
{
  R e=1e-10;
  uint n=64;
  IntervalNewtonSolver<R> interval_newton(e,n);

  HenonMap<R> h(Point<R>("(1.5,0.375)"));
  //Box<R> r("[-2.125,-2]x[-2.125,-2]");
  Box<R> r("[-2.25,-2]x[-2.25,-2]");
 
  Point<R> m("(-2,-2)");
  cout << "h.evaluate" << m << "=" << h(m) << endl;
  cout << "h.jacobian" << m << "=" << h.jacobian(m) << endl;
 
  Point< Interval<R> > pt=r;
  

  double afptd[2]={-2.0920128158902654,-2.0920128158902654};
  Point<R> afpt(2,afptd);
  cout << "h" << afpt << "=" << h(afpt) << endl;

  cout << "e=" << e << " n=" << n << endl;
  
  cout << "h=" << h << "\nr=" << r << "\npt=" << pt << endl;
  
  Point< Interval<R> > fpt;
  try {
    fpt=interval_newton.fixed_point(h,pt);
  }
  catch(EvaluationException e) {
    cout << "No solution found" << endl;
    throw e;
  }
  assert(radius(fpt)<e);
  cout << std::setprecision(20);
  cout << fpt << "  " << radius(fpt) << endl;

  assert(encloses(fpt,afpt));
 


  return 0;
}
