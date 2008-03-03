/***************************************************************************
 *            test_rounding.cc
 *
 *  Copyright  2008  Pieter Collins
 *
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

#include <iostream>
#include <iomanip>
#include "numeric/rounding.h"
#include "numeric/double.h"

#include "test.h"

using namespace std;
using namespace Ariadne;
using namespace Ariadne::Numeric;

class TestRounding 
{
 public:
  void test();
};


int main() {

  cout << setprecision(20);

  TestRounding().test();
  return ARIADNE_TEST_FAILURES;
}


void
TestRounding::test() 
{
  double rn, rl,ru;
  div_(rl,1.0,3.0,round_down);
  div_(ru,1.0,3.0,round_up);
  ARIADNE_TEST_PRINT(rl);
  ARIADNE_TEST_PRINT(ru);
  ARIADNE_TEST_ASSERT(rl<ru);
  div_(rl,2.0,5.0,round_down);
  div_(ru,2.0,5.0,round_up);
  ARIADNE_TEST_PRINT(rl);
  ARIADNE_TEST_PRINT(ru);
  ARIADNE_TEST_ASSERT(rl<ru);

  div_(rl,3.0,7,round_down);
  div_(ru,3.0,7,round_up);
  ARIADNE_TEST_PRINT(rl);
  ARIADNE_TEST_PRINT(ru);
  ARIADNE_TEST_ASSERT(rl<ru);

  div_(rl,3,7.0,round_down);
  div_(ru,3,7.0,round_up);
  ARIADNE_TEST_PRINT(rl);
  ARIADNE_TEST_PRINT(ru);
  ARIADNE_TEST_ASSERT(rl<ru);

  sqrt_(rl,2.0,round_down);
  sqrt_(ru,2.0,round_up);
  ARIADNE_TEST_PRINT(rl);
  ARIADNE_TEST_PRINT(ru);
  ARIADNE_TEST_ASSERT(rl<ru);
  ARIADNE_TEST_ASSERT(1.41<rl);
  ARIADNE_TEST_ASSERT(ru<1.42);

  exp_(rl,1.0,round_down);
  exp_(ru,1.0,round_up);
  ARIADNE_TEST_PRINT(rl);
  ARIADNE_TEST_PRINT(ru);
  ARIADNE_TEST_ASSERT(rl<ru);
  ARIADNE_TEST_ASSERT(2.71<rl);
  //ARIADNE_TEST_ASSERT(ru<2.72);

  /*
  int bad_errors=0;
  int errors=0;
  for(double x=0.0001; x<=2.00; x+=0.0001) {
    cos_(rn,x,hardware_round_near); 
    cos_(rl,x,hardware_round_down); 
    cos_(ru,x,hardware_round_up); 
    cout << x << " " << rn << " " << rl << " " << ru << "\n";
    assert(fabs(ru/rl-1.0)<0.0001);
    if(fabs(ru/rl-1.0)>=0.0001) { ++bad_errors; }
    if(rl>=ru) { ++errors; }
    //assert(rl<=rn);
    //assert(rn<=ru);
  }
  cout << "BAD ERRORS: " << bad_errors << "\n";
  cout << "ERRORS: " << errors << "\n";
  */
}
