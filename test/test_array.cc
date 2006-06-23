/***************************************************************************
 *            test_array.cc
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
#include <vector>

#include "ariadne.h"

#include "base/array.h"
#include "numeric/numerical_types.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

template class array<bool>;
template class array<double>;
template class array<double,3>;
template class array<Rational>;
template class array<Rational,4>;

int main() {
  cout << "test_array: " << flush;
  ofstream clog("test_array.log");

  array<double> a0,a1,a2,a3,a4;
  array<double> aa0,aa1,aa2,aa3,aa4;
  array<bool> b1;
  array<float, 3> c1;
  
  double input[8]={ 0.1, 1.3, 2.5, 3.7, 2.2, 3.4, 4.6, 5.8 };

  a0.resize(4);
  a0.fill(input);
  a1.resize(4);
  a2.assign(input+4,input+8);
  
  clog << "a0=" << a0 << " a1=" << a1 << " a2=" << a2 << endl;
  test_assert(a0[2]==2.5,"array assignment");

  array_vector<double> fsav(4);

  fsav.push_back(a0);
  fsav.push_back(a1);
  fsav.push_back(a2);
  fsav.pop_back();
  clog << fsav << endl;

  test_assert(fsav.size()==2,"array_vector.size");
  test_assert(fsav.length()==8,"array_vector.length");

  array_reference<array_vector<double>::element_iterator> aref=fsav[1];
  test_assert(aref==a1,"array reference equality");
  array_reference<array_vector<double>::element_const_iterator> acref=fsav[1];
  test_assert(acref==a1,"array const reference equality");
  array<double> atmp=fsav[1];
  test_assert(fsav[1]==a1,"array_vector.operator[]");

  array_vector<double>::const_iterator fsavi=fsav.begin();
  test_assert(*fsavi==a0,"array_vector::iterator");
  test_assert(fsavi!=fsav.end(),"array_vector::iterator::operator==");
  ++fsavi;
  test_assert(*fsavi==a1,"array_vector::iterator");
  fsavi+=1;
  test_assert(fsavi==fsav.end(),"array_vector::iterator::end()");

  a1.resize(3);
  a1.fill(1.1);
  std::vector< array<double> > va;
  va.push_back(a0);
  va.push_back(a1);
  va.push_back(a3);
  va.pop_back();
  va.push_back(a2);

  test_assert(va.size()==3,"vector<array>.size");
  test_assert(va[1]==a1,"vector<array>.operator[]");

  std::vector< array<double> >::const_iterator vai=va.begin();
  test_assert(*vai==a0,"vector<array>::iterator::operator*");
  test_assert(vai!=va.end(),"vector<array>::iterator::operator==");
  ++vai;
  test_assert(*vai==a1,"vector<array>::iterator::operator++");
  ++vai;
  ++vai;
  test_assert(vai==va.end(),"vector<array>::iterator::end()");

  clog << fsav;
  
  clog.close();
  cout << "PASS\n";

  return 0;
}
