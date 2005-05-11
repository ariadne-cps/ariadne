/***************************************************************************
 *            test_state.cc
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

#include "state.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;
using namespace std;

int main() {
    cout << "test_state: " << flush;

    State<Rational> s1(3);
    State<Rational> s2(4);
    State<Rational> s3(2,Rational(2,3));
    State<Rational> s4(s1);
    State<Rational> s5;

    s1[1] = Rational(0.75);
    s5=s1;

    test_assert(s1==s1 && s1!=s2 && s1!=s3 && s1!=s4 && s1==s5,"equality tests");

    /* Test output format */
    string str1("[1, 3/2]");
    string str2;
    istringstream is(str1);
    is >> s1;
    stringstream ss(str2);
    ss << s1;
    getline(ss,str2);
    
    test_assert(str1 == str2,"stream output test");
    
    /* Test input format */
    try {
	string input("[0,3/4,0] [1,1,1,1] [2/3,2/3] \n");
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
