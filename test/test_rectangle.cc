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
#include "base/exception.h"
#include "base/utility.h"
#include "geometry/rectangle.h"
#include "geometry/list_set.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

template class Rectangle< double >;
template class Rectangle< Rational >;

int main() {
    cout << "test_rectangle: " << flush;

    typedef Rectangle< Rational > ARectangle;
    typedef ARectangle::Point Point;
    typedef Interval<Rational> Interval;
    
    Point s1(2,Rational(1));
    Point s2(2,Rational(3,2));
    Point s3(2,Rational(4,3));
    Point s4(2,Rational(2));
    ARectangle r0;
    ARectangle r1(s1,s2);
    ARectangle r2(s3,s4);
    ARectangle r3(s3,s2);
    ARectangle r4,r5,r6;
    
    string istr = "[ [0,1],[0,1] ] "
        "[[-1/2,3/2],[-1/3,1/2]] "
        "[[-1/4,2/3],[1/3,3/2]] "
        "[[3/5,6/5],[2/5,7/5]] "
        "[[3/5,6/5],[2/5,1]] "
        "[[0,1],[0,1/2]] ";
    stringstream iss(istr);
    iss >> r1 >> r2 >> r3 >> r4 >> r5 >> r6;
    ARectangle r7=r1;
    
    test_assert(r1==r7,"equality");
    
    ListSet<Rational,Rectangle> cover1,cover2;
    cover1.push_back(r2);
    cover1.push_back(r3);
    cover1.push_back(r4);

    cover2.push_back(r2);
    cover2.push_back(r3);
    cover2.push_back(r5);

    test_assert(!r1.empty(),"empty");
    test_assert(r0.empty(),"empty");

    test_assert(!disjoint(r1,r1),"disjoint");
    test_assert(interiors_intersect (r1,r1),"intersects_interior");
    test_assert(subset(r1,r1),"is_subset_of");
    test_assert(!subset_of_interior(r1,r1),"subset_of_interior");
    test_assert(subset_of_open_cover(r1,cover1),"subset_of_open_cover");
    test_assert(!subset_of_open_cover(r1,cover2),"subset_of_open_cover");
       
    r4=regular_intersection(r1,r2);
    test_assert(r4==r6,"regular_intersection");

    try {
	string input("[ ]  [ [0,2] ]  [ [0,1], [3/4,4/3], [1,3/2] ] "
		     "{ lower_corner=[0,1], upper_corner=[1,4/3] }");
	stringstream is(input);
	is >> r1;
	is >> r2;
	is >> r3;
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
