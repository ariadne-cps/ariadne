/***************************************************************************
 *            test_array.cc
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
#include <vector>

#include "ariadne.h"

#include "base/stlio.h"
#include "base/array.h"
#include "numeric/rational.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Base;
using namespace Ariadne::Numeric;
using namespace std;

template class array<bool>;
template class array<double>;
template class array<double,3>;
template class array<Rational>;
template class array<Rational,4>;

int main() {
  array<double> a0,a1,a2,a3,a4;
  array<double> aa0,aa1,aa2,aa3,aa4;
  array<bool> b1;
  array<float, 3> c1;
  
  double input[8]={ 0.1, 1.3, 2.5, 3.7, 2.2, 3.4, 4.6, 5.8 };

  a0.resize(4);
  a0.fill(input);
  a1.resize(4);
  a2.assign(input+4,input+8);
  
  cout << "a0=" << a0 << " a1=" << a1 << " a2=" << a2 << endl;
  assert(a0[2]==2.5);
  a0[2]=2.3;
  assert(a0[2]==2.3);
  
  array_vector<double> fsav(4);

  fsav.push_back(a0);
  fsav.push_back(a1);
  fsav.push_back(a2);
  fsav.pop_back();
  cout << fsav << endl;

  assert(fsav.size()==2);
  assert(fsav.length()==8);

  array_reference<array_vector<double>::element_iterator> aref=fsav[1];
  assert(aref==a1);
  array_reference<array_vector<double>::element_const_iterator> acref=fsav[1];
  assert(acref==a1);
  array<double> atmp=fsav[1];
  assert(fsav[1]==a1);

  array_vector<double>::const_iterator fsavi=fsav.begin();
  assert(*fsavi==a0);
  assert(fsavi!=fsav.end());
  ++fsavi;
  assert(*fsavi==a1);
  fsavi+=1;
  assert(fsavi==fsav.end());

  a1.resize(3);
  a1.fill(1.1);
  std::vector< array<double> > va;
  va.push_back(a0);
  va.push_back(a1);
  va.push_back(a3);
  va.pop_back();
  va.push_back(a2);

  assert(va.size()==3);
  assert(va[1]==a1);
  
  std::vector< array<double> >::const_iterator vai=va.begin();
  assert(*vai==a0);
  assert(vai!=va.end());
  ++vai;
  assert(*vai==a1);
  ++vai;
  ++vai;
  assert(vai==va.end());

  cout << fsav;
  
  return 0;
}
