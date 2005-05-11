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
#include "array.h"
#include "ariadne.h"
#include "numerical_type.h"

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

  array<double> a0,a1,a2,a3,a4;
  array<double> aa0,aa1,aa2,aa3,aa4;
  array<bool> b1;
  array<float, 3> c1;
  
  double input[8]={ 0.1, 1.3, 2.5, 3.7, 2.2, 3.4, 4.6, 5.8 };

  a0.resize(4);
  a0.fill(input);
  a1.resize(4);
  a2.assign(input+4,input+8);
  
  test_assert(a0[2]==2.5,"array assignment");

  array_block_vector<double> fsav(4);

  fsav.push_back(a0);
  fsav.push_back(a1);
  fsav.push_back(a2);
  fsav.pop_back();

  test_assert(fsav.size()==2,"array_vector_fixed_size.size");
  test_assert(fsav.length()==8,"array_vector_fixed_size.length");
  test_assert(fsav[1]==a1,"array_vector_fixed_size.operator[]");

  array_block_vector<double>::const_iterator fsavi=fsav.begin();
  test_assert(*fsavi==a0,"array_block_vector::iterator");
  test_assert(fsavi!=fsav.end(),"array_block_vector::iterator::operator==");
  ++fsavi;
  test_assert(*fsavi==a1,"array_block_vector::iterator");
  fsavi+=1;
  test_assert(fsavi==fsav.end(),"array_block_vector::iterator::end()");

  a1.resize(3);
  a1.fill(1.1);
  array_vector<double> av;
  av.push_back(a0);
  av.push_back(a1);
  av.push_back(a3);
  av.pop_back();
  av.push_back(a2);

  test_assert(av.size()==3,"array_vector.size");
  test_assert(av.length()==11,"array_vector.length");
  test_assert(av[1]==a1,"array_vector.operator[]");

  array_vector<double>::const_iterator avi=av.begin();
  test_assert(*avi==a0,"array_vector::iterator::operator*");
  test_assert(avi!=av.end(),"array_vector::iterator::operator==");
  ++avi;
  test_assert(*avi==a1,"array_vector::iterator::operator++");
  ++avi;
  ++avi;
  test_assert(avi==av.end(),"array_vector::iterator::end()");

  cout << "PASS\n";

  return 0;
}
