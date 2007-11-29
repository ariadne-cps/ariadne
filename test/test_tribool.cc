/***************************************************************************
 *            test_tribool.cc
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

#include <iomanip>


#include "test/test.h"
#include "base/tribool.h"

using namespace std;
using namespace Ariadne::Base;

class TestTribool
{
  tribool t,f,i;
 public:
  TestTribool() : t(true), f(false), i(indeterminate) { std::cout << std::boolalpha; }
  
  void test_not() {
    ARIADNE_TEST_EVALUATE(!t);
    ARIADNE_TEST_EVALUATE(!f);
    ARIADNE_TEST_EVALUATE(!i);
  }

  void test_equal() {
    ARIADNE_TEST_ASSERT(bool(t==true));
    ARIADNE_TEST_ASSERT(!bool(t==indeterminate));
    ARIADNE_TEST_ASSERT(!bool(t==false));
    ARIADNE_TEST_ASSERT(!bool(i==true));
    ARIADNE_TEST_ASSERT(!bool(i==indeterminate));
    ARIADNE_TEST_ASSERT(!bool(i==false));
    ARIADNE_TEST_ASSERT(!bool(f==true));
    ARIADNE_TEST_ASSERT(!bool(f==indeterminate));
    ARIADNE_TEST_ASSERT(bool(f==false));
  }
  
  void test_indeterminate() {
    ARIADNE_TEST_EVALUATE(i==indeterminate);
    ARIADNE_TEST_ASSERT(!bool(i==indeterminate));
    ARIADNE_TEST_ASSERT(!bool(i));
    ARIADNE_TEST_ASSERT(!indeterminate(t));
    ARIADNE_TEST_ASSERT(indeterminate(i));
    ARIADNE_TEST_ASSERT(!indeterminate(f));
  }
  
  void test_possibly() {
    ARIADNE_TEST_ASSERT(possibly(t));
    ARIADNE_TEST_ASSERT(possibly(i));
    ARIADNE_TEST_ASSERT(!possibly(f));
    ARIADNE_TEST_ASSERT(possibly(!t)==false);
    ARIADNE_TEST_ASSERT(!possibly(!t));
    ARIADNE_TEST_ASSERT(possibly(!i));
    ARIADNE_TEST_ASSERT(possibly(!f));
  }
  
  void test_definitely() {
    ARIADNE_TEST_ASSERT(definitely(t));
    ARIADNE_TEST_ASSERT(!definitely(i));
    ARIADNE_TEST_ASSERT(!definitely(f));
    ARIADNE_TEST_ASSERT(definitely(!t)==false);
    ARIADNE_TEST_ASSERT(!definitely(!t));
    ARIADNE_TEST_ASSERT(!definitely(!i));
    ARIADNE_TEST_ASSERT(definitely(!f));
  }
  
  int test() { 
    test_not();
    test_equal();
    test_indeterminate();
    test_possibly();
    test_definitely();
    return 0;
  }
};

int main() {
  return TestTribool().test();
}

