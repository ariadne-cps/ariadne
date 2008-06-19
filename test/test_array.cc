/***************************************************************************
 *            test_array.cc
 *
 *  Copyright  2005-8  Davide Bresolin, Alberto Casagrande, Pieter Collins
 *  bresolin@sci.univr.it, casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
#include <vector>

#include "ariadne.h"

#include "base/stlio.h"
#include "base/array.h"
#include "numeric/rational.h"

#include "test/test.h"

using namespace Ariadne;
using namespace std;

template class array<bool>;
template class array<double>;
template class array<double,3>;
template class array<Rational>;
template class array<Rational,4>;

double input[8]={ 0.1, 1.3, 2.5, 3.7, 2.2, 3.4, 4.6, 5.8 };

class TestArray {
  
	array<double> a0,a1,a2,a3,a4;
  array<double> aa0,aa1,aa2,aa3,aa4;
  array<bool> b1;
  array<float, 3> c1;
 
	public:
  void test_resize() {
    a0.resize(4);
    a1.resize(4);
    ARIADNE_TEST_ASSERT(a0.size()==4)
    ARIADNE_TEST_ASSERT(a1.size()==4)		
	}
	
	void test_fill() {
    a0.resize(4);
    a0.fill(input);
		ARIADNE_TEST_ASSERT(a0[0]==input[0]);
		ARIADNE_TEST_ASSERT(a0[1]==input[1]);
		ARIADNE_TEST_ASSERT(a0[2]==input[2]);
		ARIADNE_TEST_ASSERT(a0[3]==input[3]);	  
	}
	
	void test_assign() {
    a2.assign(input+4,input+8);
		ARIADNE_TEST_ASSERT(a2[0]==input[4]);
		ARIADNE_TEST_ASSERT(a2[1]==input[5]);
		ARIADNE_TEST_ASSERT(a2[2]==input[6]);
		ARIADNE_TEST_ASSERT(a2[3]==input[7]);
  }
	
	void test_value() {
    a0[2]=2.3;
    ARIADNE_TEST_ASSERT(a0[2]==2.3);	
	}

  void test_push_pop() {
	  array_vector<double> fsav(4);

    fsav.push_back(a0);
    fsav.push_back(a1);
    fsav.push_back(a2);
    fsav.pop_back();
    cout << fsav << endl;

    ARIADNE_TEST_ASSERT(fsav.size()==2);
    ARIADNE_TEST_ASSERT(fsav.length()==8);
	}

	void test_reference() {
		array_reference<array_vector<double>::element_iterator> aref=a1;
		ARIADNE_TEST_ASSERT(aref==a1);

		array_reference<array_vector<double>::element_const_iterator> acref=aref;
		ARIADNE_TEST_ASSERT(acref==a1);

		array<double> atmp=a1;
		ARIADNE_TEST_ASSERT(atmp==a1);	
	}
	
	void test_iterator() {
	  array_vector<double> fsav(4);
    fsav.push_back(a0);
    fsav.push_back(a1);
	  
    array_vector<double>::const_iterator fsavi=fsav.begin();
    ARIADNE_TEST_ASSERT(*fsavi==a0);
    ARIADNE_TEST_ASSERT(fsavi!=fsav.end());
    ++fsavi;
    ARIADNE_TEST_ASSERT(*fsavi==a1);
    fsavi+=1;
    ARIADNE_TEST_ASSERT(fsavi==fsav.end());

    a1.resize(3);
    a1.fill(1.1);
    std::vector< array<double> > va;
    va.push_back(a0);
    va.push_back(a1);
    va.push_back(a3);
    va.pop_back();
    va.push_back(a2);

    ARIADNE_TEST_ASSERT(va.size()==3);
    ARIADNE_TEST_ASSERT(va[1]==a1);

    std::vector< array<double> >::const_iterator vai=va.begin();
		ARIADNE_TEST_ASSERT(*vai==a0);
		ARIADNE_TEST_ASSERT(vai!=va.end());
		++vai;
		ARIADNE_TEST_ASSERT(*vai==a1);
		++vai;
		++vai;
		ARIADNE_TEST_ASSERT(vai==va.end());
	}

  void test()
  {
    ARIADNE_TEST_CALL(test_resize());
    ARIADNE_TEST_CALL(test_fill());
    ARIADNE_TEST_CALL(test_assign());
    ARIADNE_TEST_CALL(test_value());
    ARIADNE_TEST_CALL(test_push_pop());
    ARIADNE_TEST_CALL(test_reference());
    ARIADNE_TEST_CALL(test_iterator());
  }
	
};

int main() {
	TestArray().test();
	return ARIADNE_TEST_FAILURES;
}
