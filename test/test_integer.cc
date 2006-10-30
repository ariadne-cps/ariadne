/***************************************************************************
 *            test_integer.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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

#include <cassert>
#include <fstream>
#include <iomanip>

#include <gmpxx.h>

#include "numeric/integer.h"

using namespace std;
using namespace Ariadne::Numeric;

template<class R> void test_integer();

int test_integer();

int main() {
  test_integer();
  return 0;
}

int test_integer()
{
  cout << setprecision(20);
  mpf_set_default_prec (8);

  // Default constructor
  Integer i0;
  
  // Construct from an int
  Integer i1(0);
  assert(i0==i1);
  
  // Copy constructor
  Integer i2(i0);
  assert(i0==i2);
  
  // Copy assignment
  i1=i2;
  assert(i1==i2);

  // Arithmetic
  assert(-Integer(-2)==Integer(2));
  assert(Integer(2)+Integer(-5)==Integer(-3));
  assert(Integer(2)-Integer(-5)==Integer(7));
  assert(Integer(2)*Integer(-5)==Integer(-10));
  assert(pow(Integer(5),3u)==Integer(125));

  assert(factorial(Integer(5))==Integer(120));
  assert(gcd(Integer(140),Integer(75))==Integer(5));
  assert(lcm(Integer(140),Integer(75))==Integer(2100));


  return 0;
}
