/***************************************************************************
 *            test_rectangle.cc
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
#include "exception.h"
#include "utility.h"
#include "numerical_type.h"
#include "geometry_state.h"
#include "rectangle.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

int main() {
    typedef AriadneRectangle< AriadneState<Rational> > Rectangle;
    typedef Rectangle::State State;
    typedef Rectangle::Real Real;
    typedef interval<Real> Interval;
    
    State s1(4,Rational(1));
    State s2(4,Rational(4,3));
    Rectangle r1(s1,s2);
    Rectangle r2(s2,s1);
    Rectangle r3;
    Rectangle r4;
    
    /* Test input format */
    try {
	string input("[ ]  [ [0,2] ]  [ [0,1], [3/4,4/3], [1,3/2] ]  { lower_corner=[0,1], upper_corner=[1,4/3] }");
	stringstream is(input);
	is >> r1;
	is >> r2;
	is >> r3;
	is >> r4;
   } 
    catch(invalid_input& e) {
	cout << "test_rectangle: FAILED\n";
	cout << "  invalid_input: " << e.what() << "\n";
	return 1;
    }
    catch(std::invalid_argument& e) {
	cout << "test_rectangle: FAILED\n";
	cout << "  std::invalid_argument: " << e.what() << "\n";
	return 1;
    }
    catch(...) {
	cout << "test_rectangle: FAILED\n";
	cout << "  Unknown error\n"; 
	return 1;
    }
	
    cout << "test_rectangle: PASS\n";

    return 0;
}
