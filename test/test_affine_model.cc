/***************************************************************************
 *            test_affine_model.cc
 *
 *  Copyright  2007  Pieter Collins
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

#include <fstream>

#include "test_float.h"

#include "ariadne.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "function/affine_model.h"

#include "test.h"

using namespace Ariadne;
using namespace std;


template<class R>
class TestAffineModel
{
 public:
  typedef Interval<R> I;
  
  AffineModel<R> am1;
  AffineModel<R> am2;
  AffineModel<R> am3;

  TestAffineModel()
    : am1(Vector<I>("[[-1,1],[-1,1],[-1,1]]"),Vector<R>("[0,0,0]"),Vector<I>("[[-0.125,0.125],[0.875,1.125]]"),Matrix<R>("[1,2,3;4,5,6]"))
    , am2(Vector<I>("[[-1,1],[-1,1],[-1,1]]"),Vector<R>("[0,0,0]"),Vector<I>("[[-0.125,0.125],[0.875,1.125]]"),Matrix<I>("[[1,1.1],2,[3,3.1];4,[5,5.5],6]"))
    , am3(Vector<I>("[[-8,8],[-24,24]]"),Vector<R>("[0,0]"),Vector<I>("[[-0.125,0.125],[0.875,1.125]]"),Matrix<I>("[[1,1.1],2;[3,3.1],4]"))
  { }

  void test_copy_constructor() {
    AffineModel<R> am(am1);
    cout << am << endl;
    am=am2;
    cout << am << endl;
  }

  void test_range() {
    ARIADNE_ASSERT(refines(am1.range(),Vector<I>("[[-7,7],[-16,17]]")))
  }

  void test_reduce() {
    cout << reduce(am1,0) << endl;
  }

  void test_add() {
    cout << AffineModel<R>(am1+am2) << endl;
  }

  void test_compose() {
    cout << "am3="<<am3<<"\nam2="<<am2<<endl;
    cout << "am2.range()="<<am2.range()<<endl;
    cout << compose(am3,am2) << endl;
  }

  void test_write() {
    cout << am2 << endl;
  }


  void test() {
    ARIADNE_TEST_CALL(test_copy_constructor());
    ARIADNE_TEST_CALL(test_range());
    //ARIADNE_TEST_CALL(test_reduce());
    ARIADNE_TEST_CALL(test_compose());
    ARIADNE_TEST_CALL(test_write());
    //ARIADNE_TEST_CALL(test_add());
  }
};
  

int main() {
  TestAffineModel<Flt>().test();
  return ARIADNE_TEST_FAILURES;
}
          
