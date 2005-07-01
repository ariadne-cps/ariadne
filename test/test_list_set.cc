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
#include <string>

#include "ariadne.h"
#include "exception.h"
#include "utility.h"
#include "numerical_type.h"
#include "state.h"
#include "rectangle.h"
#include "list_set.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

template class ListSet<Rational, Rectangle>;
template class ListSet<double, Rectangle>;

int main() {
  cout << "test_list_set: " << flush;

    State<Rational> s0(2,Rational(0));
    State<Rational> s1(2,Rational(1));
    State<Rational> s2(2,Rational(4,3));
    State<Rational> s3(2,Rational(3,2));

    typedef ListSet<Rational, Rectangle> ListSet;
    typedef Rectangle<Rational> Rectangle;
    typedef State<Rational> State;

    Rectangle r0(s0,s1);
    Rectangle r1(s1,s2);
    Rectangle r2(s2,s3);
    ListSet ds1,ds2;

    ds1.inplace_union(r0);
    ds1.inplace_union(r1);
    ds1.inplace_union(r2);

    string input("[ [[0,1],[0,1]], [[1,4/3],[1,4/3]], [[4/3,3/2],[4/3,3/2]] ]");
    stringstream is(input);
    is >> ds2;

    test_assert(ds1[0]==ds2[0] && ds1[1]==ds2[1] && ds1[2]==ds2[2], "stream input");

    /* Test input format */
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

    cout << "PASS\n";

    return 0;
}
