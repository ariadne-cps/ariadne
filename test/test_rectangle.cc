/***************************************************************************
 *            test_rectangle.cc
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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

#include "real_typedef.h"

#include "ariadne.h"
#include "base/exception.h"
#include "base/utility.h"
#include "geometry/rectangle.h"
#include "geometry/list_set.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

template<class R> int test_rectangle();

int main() {
  test_rectangle<Real>();
}

template<class R>
int
test_rectangle()
{
 
    Point<R> s1("(1,1)");
    Point<R> s2("(1.5,1.5)");
    Point<R> s3("(1.375,1.375)");
    Point<R> s4("(2,2)");
    Point<R> s5("(0.75,0.125)");
    Rectangle<R> r0(2);
    Rectangle<R> r1(s1,s2);
    Rectangle<R> r2(s3,s4);
    Rectangle<R> r3(s3,s2);
    Rectangle<R> r4,r5,r6;
    
    string istr = "[0,1]x[0,1] "
        "[-0.5,1.5]x[-0.375,0.5] "
        "[-0.25,0.625]x[0.375,1.5] "
        "[0.5625,1.125]x[0.4375,1.375] "
        "[0.,1.1875]x[0.4375,1] "
        "[0,1]x[0,0.5] ";
    stringstream iss(istr);
    iss >> r1 >> r2 >> r3 >> r4 >> r5 >> r6;
    Rectangle<R> r7=r1;
    
    cout << "r1=" << r1 << ", r2=" << r2 << ", r3=" << r3 << ", r4=" << r4 << "\n"
         << "r5=" << r5 << ", r6=" << r6 << ", r7=" << r7 << endl;
    
    assert(r1==r7);
    assert(indeterminate(equal(r1,r7)));
    cout << "centre(r2)=" << r2.centre() << endl;
    assert(r2.centre()==Point<R>("(0.5,0.0625)"));
    
    cout << "r4=" << r4 << endl;
    cout << "r4.vertices()=" << flush;
    RectangleVerticesIterator<R> r4e(r4,true);
    for(typename Rectangle<R>::vertices_iterator v=r4.vertices_begin();
        v!=r4.vertices_end(); ++v)
    {
      cout << *v << " " << flush;
    }
    cout << endl;
    
    ListSet<R,Rectangle> cover1,cover2;
    cover1.push_back(r2);
    cover1.push_back(r3);
    cover1.push_back(r4);

    cover2.push_back(r2);
    cover2.push_back(r3);
    cover2.push_back(r5);

    cout << "r0=" << r0 << ", r0.dimension()=" << r0.dimension()
         << ", r0.empty()=" << r0.empty() << endl;
    cout << "r1=" << r1 << ", r1.dimension()=" << r1.dimension()
         << ", r1.empty()=" << r1.empty() << endl;
    assert(!r1.empty());
    assert(r0.empty());

    assert(!disjoint(r1,r1));
    assert(indeterminate(subset(r1,r1)));
    cout << "r1=" << r1 << ", s1=" << s1 << ", s3=" << s3 << ", s5=" << s5 << endl;
    assert(indeterminate(r1.contains(s1)));
    assert(!r1.contains(s3));
    assert(r1.contains(s5));
    //assert(!subset_of_interior(r1,r1));
    //assert(subset_of_open_cover(r1,cover1));
    //assert(!subset_of_open_cover(r1,cover2));
       
    r4=open_intersection(r1,r2);
    cout << "r1=" << r1 << ",r2=" << r2 << ",r4=" << r4 << ",r6=" << r6 << endl;
    assert(r4==r6);

    cout << "r1=" << r1;
    r1[1]=Interval<R>(0.0,0.5);
    cout << " r1=" << r1;
    assert(r1==r6);
    //r1[1].lower()=-2.25;
    cout << " r1=" << r1 << endl;
    
    try {
        string input("[ ]  [ [0,2] ]  [ [0,1], [3/4,4/3], [1,3/2] ] "
                     "{ lower_corner=[0,1], upper_corner=[1,4/3] }");
        stringstream is(input);
        is >> r1;
        is >> r2;
        is >> r3;
   } 
    catch(invalid_input& e) {
        cout << "  invalid_input: " << e.what() << "\n";
        return 1;
    }
    catch(std::invalid_argument& e) {
        cout << "  std::invalid_argument: " << e.what() << "\n";
        return 1;
    }
    catch(std::exception& e) {
        cout << "  Unknown error\n";
        throw e;
    }
        


    return 0;
}
