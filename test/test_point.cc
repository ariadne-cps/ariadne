/***************************************************************************
 *            test_point.cc
 *
 *  2 May 2005
 *  Copyright  2005  Pieter Collins, Alberto Casagrande
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

#include "real_typedef.h"
#include "geometry/point.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;


int main() {
    cout << "test_point: " << flush;
    ofstream clog("test_point.log");
  
    Point<Real> s1(3);
    Point<Real> s2(4);
    Point<Real> s3(2);
    s3[0]=Real(2,3);
    s3[1]=Real(2,3);
    Point<Real> s4(s1);
    Point<Real> s5;

    s1[1] = Real(0.75);
    s5=s1;

    test_assert(s1==s1 && s1!=s2 && s1!=s3 && s1!=s4 && s1==s5,
                "equality tests");

    clog << s1 << " " << s2 << " " << s3 << " " << s4 << " " << s5 << endl;
    /* Test output format */
    string str1("(1, 1.5)");
    string str2;
    istringstream is(str1);
    is >> s1;
    stringstream ss(str2);
    ss << s1;
    getline(ss,str2);
    
    clog << str1 << " " << str2 << endl;
    
    test_assert(str1 == str2,"stream output test");
    
    /* Test input format */
    try {
      string input("(0,0.75,0) (1,1,1,1) (0.666,0.666) \n");
      stringstream is(input);
      
      is >> s1 >> s2 >> s3;
    } 
    catch(std::invalid_argument& e) {
      cerr << "std::invalid_argument " << e.what() << "\n";
      return 1;
    }
    catch(std::runtime_error& e) {
      cerr << "std::runtime_error " << e.what() << "\n";
      return 1;
    }
    catch(...) {
      cerr << "Unknown error\n";
      return 1;
    }
        
    cout << "PASS\n";

    return 0;
}
