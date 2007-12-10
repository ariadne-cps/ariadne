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
#include "test/test.h"

using namespace std;
using namespace Ariadne::Numeric;



class TestInteger {
 public:
  TestInteger() {
    cout << setprecision(20);
    mpf_set_default_prec (8);
  }    

  void test_constructors() {
    // Default constructor
    Integer i0;
    ARIADNE_TEST_ASSERT(i0==0);

    // Construct from an int
    Integer i1(0);
    ARIADNE_TEST_ASSERT(i0==i1);
  
    // Copy constructor
    Integer i2(i0);
    ARIADNE_TEST_ASSERT(i0==i2);
    
    // Copy assignment
    i1=i2;
    ARIADNE_TEST_ASSERT(i1==i2);
  }

  void test_cmp() {
  }

  void test_neg() {
    ARIADNE_TEST_ASSERT(-Integer(-2)==Integer(2));
  }

  void test_add() {
    ARIADNE_TEST_ASSERT(Integer(2)+Integer(-5)==Integer(-3));
  }

  void test_sub() {
  ARIADNE_TEST_ASSERT(Integer(2)-Integer(-5)==Integer(7));
 }

  void test_mul() {
  ARIADNE_TEST_ASSERT(Integer(2)*Integer(-5)==Integer(-10));
  }

  void test_pow() {
    ARIADNE_TEST_ASSERT(pow(Integer(5),3u)==Integer(125));
  }

  void test_fac() {
    Integer z=13;
    ARIADNE_TEST_ASSERT(fac(0)==1);
    ARIADNE_TEST_ASSERT(fac<Integer>(0)==1);
    ARIADNE_TEST_ASSERT(Integer(fac(Integer(0)))==1);
    ARIADNE_TEST_ASSERT(fac(12)==479001600);
    ARIADNE_TEST_ASSERT(Integer(fac(Integer(12)))==479001600);
    ARIADNE_TEST_FAIL(fac(13));
    ARIADNE_TEST_ASSERT(fac<Integer>(13)==Integer("6227020800"));
    ARIADNE_TEST_ASSERT(fac(Integer(13))==Integer("6227020800"));
  }

  void test_bin() {
    Integer r; int n=32; int k=12;
    bin_(r,n,k);
    ARIADNE_TEST_EQUAL(bin(7,3),35); 
    ARIADNE_TEST_EQUAL(bin(20u,17u),1140u); 
    ARIADNE_TEST_EQUAL(bin<Integer>(20,17),1140u); 
    ARIADNE_TEST_EQUAL(Integer(bin(20,17)),1140u); 
    ARIADNE_TEST_THROW(bin(42,23),OverflowException);
    ARIADNE_TEST_EQUAL(bin<Integer>(42,23),Integer("446775310800"));
    ARIADNE_TEST_EQUAL(Integer(bin(Integer(42),23)),Integer("446775310800")); 
  }

  void test_gcd() {
    ARIADNE_TEST_ASSERT(gcd(Integer(140),Integer(75))==Integer(5));
  }

  void test_lcm() {
    ARIADNE_TEST_ASSERT(lcm(Integer(140),Integer(75))==Integer(2100));
  }

  void test()
  {
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_cmp());
    ARIADNE_TEST_CALL(test_neg());
    ARIADNE_TEST_CALL(test_add());
    ARIADNE_TEST_CALL(test_sub());
    ARIADNE_TEST_CALL(test_mul());
    ARIADNE_TEST_CALL(test_pow());
    ARIADNE_TEST_CALL(test_fac());
    ARIADNE_TEST_CALL(test_bin());
    ARIADNE_TEST_CALL(test_gcd());
    ARIADNE_TEST_CALL(test_lcm());
  }

};


int main() {
  TestInteger().test();
  //  return ARIADNE_TEST_FAILURES;
}

