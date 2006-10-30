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


#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "evaluation/newton.h"

#include "models/henon.h"
#include "test.h"
#include "real_typedef.h"

using namespace std;
using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;

template<class R> int test_newton();

int main() {
  test_newton<Real>();
  return 0;
}

template<class R>
int 
test_newton()
{
  

  HenonMap<R> h(1.5,0.375);
  Rectangle<R> r("[-2.125,-2]x[-2.125,-2]");
  R e=1e-10;
  
  cout << h << r << e << endl;
  
  Rectangle<R> fr;
  try {
    fr=interval_newton(DifferenceMap<R>(h),r,e);
  }
  catch(EvaluationException e) {
    cout << "No solution found" << endl;
    throw e;
  }
  assert(fr.radius()<e);
  cout << std::setprecision(20);
  cout << fr << "  " << fr.radius() << endl;

  Point<R> fp("(-2.0920128158902654,-2.0920128158902654)");
  assert(fr.contains(fp));

  Point< Interval<R> > ifp=h(fp);
  


  return 0;
}
