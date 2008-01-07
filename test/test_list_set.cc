/***************************************************************************
 *            test_list_set.cc
 *
 *  24 June 2005
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
#include <sstream>
#include <string>

#include "ariadne.h"
#include "test_float.h"
#include "base/utility.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace std;
  
template<class R> int test_list_set();
  
int main() {
  test_list_set<Flt>();
  cerr << "INCOMPLETE ";
}

template<class R>
int test_list_set() 
{
  
  Point<R> s0("(0,0)");
  Point<R> s1("(1,1)");
  Point<R> s2("(1.375,1.375)");
  Point<R> s3("(1.5,1.5)");
  cout << s0 << " " << s1 << " " << s2 << " " << s3 << endl;
  
  Box<R> r0(s0,s1);
  Box<R> r1(s1,s2);
  Box<R> r2(s2,s3);
  Box<R> r3(s3,s3);
  cout << r0 << " " << r1 << " " << r2 << endl;
  
  Zonotope<R> z0(r0);
  Zonotope<R> z1(r1);
  
  ListSet< Box<R> > ds1;
  ds1.adjoin(r0);
  ds1.adjoin(r1);
  ds1.adjoin(r2);
  ds1.adjoin(r3);
  cout << ds1 << endl;
  assert(ds1.size()==4);
  assert(ds1[0]==r0);
  assert(ds1[1]==r1);
  assert(ds1[2]==r2);
  assert(ds1[3]==r3);
  
  cout << ds1 << endl;
  
  string input("[ [0,1]x[0,1], [1,1.375]x[1,1.375], [1.375,1.5]x[1.375,1.5] ]");
  stringstream is(input);
  
  ListSet< Box<R> > ds2;
  is >> ds2;
  cout << ds2 << endl;
  assert(ds2.size()==3);
  assert(ds1[0]==ds2[0] && ds1[1]==ds2[1] && ds1[2]==ds2[2]);
  

  ListSet< Zonotope<R> > zds;
  zds.adjoin(z0);
  assert(zds.size()==1);
  zds.adjoin(z1);
  assert(zds.size()==2);
  zds.adjoin(r0);
  assert(zds.size()==3);
  cout << zds << endl;

  return 0;
}
