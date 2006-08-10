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
#include <fstream>
#include <string>

#include "ariadne.h"
#include "real_typedef.h"
#include "base/exception.h"
#include "base/utility.h"
#include "geometry/rectangle.h"
#include "geometry/list_set.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

template class Rectangle< Real >;
template class Point< Real >;

int main() {
    cout << "test_rectangle: " << flush;
    ofstream clog("test_rectangle.log");
  
    Point<Real> s1("(1,1)");
    Point<Real> s2("(1.5,1.5)");
    Point<Real> s3("(1.375,1.375)");
    Point<Real> s4("(2,2)");
    Point<Real> s5("(0.75,0.125)");
    Rectangle<Real> r0(2);
    Rectangle<Real> r1(s1,s2);
    Rectangle<Real> r2(s3,s4);
    Rectangle<Real> r3(s3,s2);
    Rectangle<Real> r4,r5,r6;
    
    string istr = "[0,1]x[0,1] "
        "[-0.5,1.5]x[-0.375,0.5] "
        "[-0.25,0.625]x[0.375,1.5] "
        "[0.5625,1.125]x[0.4375,1.375] "
        "[0.,1.1875]x[0.4375,1] "
        "[0,1]x[0,0.5] ";
    stringstream iss(istr);
    iss >> r1 >> r2 >> r3 >> r4 >> r5 >> r6;
    Rectangle<Real> r7=r1;
    
    clog << "r1=" << r1 << ", r2=" << r2 << ", r3=" << r3 << ", r4=" << r4 << "\n"
         << "r5=" << r5 << ", r6=" << r6 << ", r7=" << r7 << endl;
    
    test_assert(r1==r7,"equality");
    
    ListSet<Real,Rectangle> cover1,cover2;
    cover1.push_back(r2);
    cover1.push_back(r3);
    cover1.push_back(r4);

    cover2.push_back(r2);
    cover2.push_back(r3);
    cover2.push_back(r5);

    clog << "r0=" << r0 << ", r0.dimension()=" << r0.dimension() << ", r0.empty()=" << r0.empty() << endl;
    clog << "r1=" << r1 << ", r1.dimension()=" << r1.dimension() << ", r1.empty()=" << r1.empty() << endl;
    test_assert(!r1.empty(),"empty");
    test_assert(r0.empty(),"empty");

    test_assert(!disjoint(r1,r1),"disjoint");
    test_assert(interiors_intersect (r1,r1),"intersects_interior");
    test_assert(subset(r1,r1),"is_subset_of");
    clog << "r1=" << r1 << ", s1=" << s1 << ", s3=" << s3 << ", s5=" << s5 << endl;
    test_assert(r1.contains(s1),"contains");
    test_assert(!r1.contains(s3),"contains");
    test_assert(r1.contains(s5),"contains");
    test_assert(!r1.interior_contains(s1),"interior_contains");
    test_assert(!r1.interior_contains(s3),"interior_contains");
    test_assert(r1.interior_contains(s5),"interior_contains");
    //test_assert(!subset_of_interior(r1,r1),"subset_of_interior");
    //test_assert(subset_of_open_cover(r1,cover1),"subset_of_open_cover");
    //test_assert(!subset_of_open_cover(r1,cover2),"subset_of_open_cover");
       
    r4=regular_intersection(r1,r2);
    clog << "r1=" << r1 << ",r2=" << r2 << ",r4=" << r4 << ",r6=" << r6 << endl;
    test_assert(r4==r6,"regular_intersection");

    clog << "r1=" << r1;
    r1[1]=Interval<Real>(0.0,0.5);
    clog << " r1=" << r1;
    test_assert(r1==r6,"operator[]");
    r1[1].lower()=-2.25;
    clog << " r1=" << r1 << endl;
    
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
