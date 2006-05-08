/***************************************************************************
 *            test_polyhedron.cc
 *
 *  2 May 2005
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
#include <string>

#include "ariadne.h"
#include "base/exception.h"
#include "base/utility.h"
#include "numeric/numerical_types.h"
#include "geometry/point.h"
#include "geometry/polyhedron.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

template class Polyhedron<Rational>;
//template class Polyhedron<double>;

int main() {
  cout << "test_polyhedron: " << flush;

  typedef Rational Real;
  typedef Polyhedron<Real> Polyhedron;
  typedef Point<Real> Point;
  typedef Interval<Real> Interval;
    
  Point s1(4,Rational(1));
  Point s2(4,Rational(4,3));
/*
  Polyhedron p1(s1,s2);
  Polyhedron p2(s2,s1);
  Polyhedron p3;
  Polyhedron p4;
*/
    
  Polyhedron p;
  std::cout << p;
  
  /* Test input format */
  try {
    string input("[ ]  [ [0,2] ]  [ [0,1], [3/4,4/3], [1,3/2] ]  { lower_corner=[0,1], upper_corner=[1,4/3] }");
    stringstream is(input);
  } 
  catch(invalid_input& e) {
    cout << "FAILED\n";
    cout << "  invalid_input: " << e.what() << "\n";
    return 1;
  }
  catch(std::invalid_argument& e) {
    cout << "FAILED\n";
    cout << "  std::invalid_argument: " << e.what() << "\n";
    return 1;
  }
  catch(...) {
    cout << "FAILED\n";
    cout << "  Unknown error\n"; 
    return 1;
  }

  cout << "INCOMPLETE\n";

  return 0;
}
