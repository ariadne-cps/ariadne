/***************************************************************************
 *            test_epsfstream.cc
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
#include "base/exception.h"
#include "base/utility.h"
#include "numeric/numerical_types.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/list_set.h"
#include "utility/epsfstream.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace Ariadne::Postscript;
using namespace std;

int main() {
  cout << "test_epsfstream: " << flush;
  ofstream clog("test_epsfstream.log");
  
  Rectangle<Real> bbox(2), r1, r2;
  Rectangle<Real> r3;
  Rectangle<Real> r4;
  string input("[-0.125,1.125]x[-0.25, 3.25] "
               "[ 0.0,1.0]x[0.0,2.0] "
               "[ 0.5,1.0]x[1.0,3.0] "
               "[ 0,0.3333333]x[2.3333,3] "
               "[ 0.06125,0.125]x[0.5,2.75] "
               );
  stringstream iss(input);

  iss >> bbox >> r1 >> r2 >> r3 >> r4;
  
  clog << bbox << "\n";
  clog << r1 << " " << r2 << " "<< r3 << " " << r4 << std::endl;
  
  epsfstream<Real> eps("../test_epsfstream.eps",bbox);

  eps << r1 << r2 << r3 << r4;


  try {
    string input("[ ] "
                "[ [0,2] ] "
                "[ [0,1], [3/4,4/3], [1,3/2] ] "
                "{ lower_corner=[0,1], upper_corner=[1,4/3] } " );
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

  clog.close();
  cout << "PASS\n";

  return 0;
}
