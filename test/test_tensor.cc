/***************************************************************************
 *            test_tensor.cc
 *
 *  Copyright  2006  Pieter Collins, Alberto Casagrande
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
#include <fstream>
#include <cassert>

#include "test/test_float.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/tensor.h"

using namespace std;
using namespace Ariadne;

int test_tensor_index();
template<class R> int test_tensor();
 
int main() {
  test_tensor_index();
  test_tensor<Rational>();
  cerr << "INCOMPLETE ";
  return 0;
}

int
test_tensor_index()
{
  cout << "test_tensor_index()" << endl;
  
  bin(5,2);
  bin(7,3); 
  bin(8,4); 
  bin(1,2);

  MultiIndex mi(4);
  mi.increment_index(0);
  mi.increment_index(0);
  mi.increment_index(2);
  mi.increment_index(2);
  mi.increment_index(3);
  mi.position();

  mi.decrement_index(0);
  mi.increment_index(1);
  mi.position();
  
  mi.decrement_index(3);
  mi.increment_index(2);
  mi.position();
 
  mi.decrement_index(0);
  mi.increment_index(1);
  mi.position();
   
  Index i(4);
  i.push_back(0);
  i.push_back(2);
  i.push_back(1);
  i.push_back(2);
  i.push_back(0);
  
  
  mi=MultiIndex(5);
  mi.set_index(0,4);
  MultiIndex mie=mi; mie.increment_index(0);
  MultiIndexIterator miie(mie);
  uint zi=0;
  cout << endl << endl;
  for(MultiIndexIterator mii(mi); mii!=miie; ++mii) {
    assert(mii->position()==zi++);
    cout << *mii << " " << mii->position() << flush;
    cout << " " << mii->number() << endl;
  }  
  cout << "\n";
  
  return 0;
  
}

template<class R>
int
test_tensor()
{
  //DerivativeTensor<R> t(2,3,2);
  //std::cerr << t.number_of_elements() << endl;
  //std::cerr << t << endl;
  std::cout << "test_tensor<" << name<R>() << ">()" << endl;
  
  R data[6]={21,13,8,5,3,2};
  SymmetricTensor<R> S(3,2,data);
  cout << S << endl;
  
  DerivativeTensor<R> A(2,2,1);
  DerivativeTensor<R> v(2,2,0);
  
  {
    size_type i(0);
    MultiIndex j(2);
    j.set_index(0,1);
    cout << "A(" << i << "," << j << ")=" << flush;
    cout << A(i,j) << endl;
  }
  
  MultiIndex ms(3);
  ms.set_index(0,2);
  S(ms)=13.0;
  ms.decrement_index(0);
  ms.increment_index(1);
  S(ms)=13.0;
  ms.decrement_index(1);
  ms.increment_index(2);
  
  
  
  MultiIndex m(2);
  m.set_index(0,1);
  A(0,m)=8.0;
  A(1,m)=5.0;
  m.set_index(0,0);
  m.set_index(1,1);
  A(0,m)=3.0;
  A(1,m)=2.0;
 
  m.set_index(0,0);
  m.set_index(1,0);
  v(0,m)=11;
  v(1,m)=13;
  
  cout << A << endl;
  cout << v << endl;
  cout << A*v << endl;
  
  A=DerivativeTensor<R>(1,2,2);
  m.set_index(0,2);
  m.set_index(1,0);
  A(0,m)=5.0;
  m.decrement_index(0);
  m.increment_index(1);
  A(0,m)=3.0;
  m.decrement_index(0);
  m.increment_index(1);
  A(0,m)=2.0;
 
  m.set_index(0,0);
  m.set_index(1,0);
  v(0,m)=1;
  v(1,m)=1;
  

  cout << A << endl;
  cout << v << endl;
  cout << A*v << endl;
  cout << (A*v)*v << endl;
  
 return 0;
}
