/***************************************************************************
 *            test_list_set.cc
 *
 *  24 June 2005
 *  Copyright  2005  Pieter Collins
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

#include <iostream>
#include <fstream>
#include <string>

#include "ariadne.h"
#include "real_typedef.h"
#include "base/exceptions.h"
#include "base/utility.h"
#include "numeric/numerical_types.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/list_set.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

int main() {
  
  Point<Real> s0("(0,0)");
  Point<Real> s1("(1,1)");
  Point<Real> s2("(1.375,1.375)");
  Point<Real> s3("(1.5,1.5)");
  
  cout << s0 << " " << s1 << " " << s2 << " " << s3 << endl;
  
  typedef ListSet<Real, Rectangle> ListSet;
  typedef Rectangle<Real> Rectangle;
  typedef Point<Real> Point;

  Rectangle r0(s0,s1);
  Rectangle r1(s1,s2);
  Rectangle r2(s2,s3);
  
  cout << r0 << " " << r1 << " " << r2 << endl;
  
  ListSet ds1;
  ds1.inplace_union(r0);
  ds1.inplace_union(r1);
  ds1.inplace_union(r2);

  cout << ds1 << endl;
  
  string input("[ [0,1]x[0,1], [1,1.375]x[1,1.375], [1.375,1.5]x[1.375,1.5] ]");
  stringstream is(input);
  
  ListSet ds2;
  is >> ds2;
  cout << ds2 << endl;

  assert(ds1[0]==ds2[0] && ds1[1]==ds2[1] && ds1[2]==ds2[2]);

  /* Test input format */
  try {
    string input("[ ] "
                 "[ [0,2] ] "
                 "[ [0,1], [3/4,4/3], [1,3/2] ] "
                 "{ lower_corner=[0,1], upper_corner=[1,4/3] } " );
    stringstream is(input);
  }
  catch(invalid_input& e) {
    cout << "  invalid_input: " << e.what() << "\n";
    return 1;
  }
  catch(std::invalid_argument& e) {
    cout << "  std::invalid_argument: " << e.what() << "\n";
    return 1;
  }
  catch(...) {
    cout << "  Unknown error\n";
    return 1;
  }

  return 0;
}
