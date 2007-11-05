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
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Function;
using namespace Ariadne::Geometry;
using namespace Ariadne::Output;
using namespace std;


template<class R>
class TestAffineModel
{
 public:
  typedef Interval<R> I;
  
  AffineModel<R> c0model;
  AffineModel<R> c1model;
  AffineModel<R> c1model2;
  // AffineModel<R> exmodel;

  TestAffineModel()
    : c0model(Vector<I>("[[-0.125,0.125],[0.875,1.125]]"),Matrix<R>("[1,2,3;4,5,6]"))
    , c1model(Vector<I>("[[-0.125,0.125],[0.875,1.125]]"),Matrix<I>("[[1,1.1],2,[3,3.1];4,[5,5.5],6]"))
    , c1model2(Vector<I>("[[-0.125,0.125],[0.875,1.125]]"),Matrix<I>("[[1,1.1],2;[4,4.1],[5,5.5]]"))
        { }

  void test_copy_constructor() {
    AffineModel<R> am(c0model);
    cout << am << endl;
    am=c1model;
    cout << am << endl;
  }

  void test_reduce() {
    cout << reduce(c1model,0) << endl;
  }

  void test_add() {
    cout << c0model+c0model << endl;
    cout << c1model+c1model << endl;
  }

  void test_compose() {
    cout << compose(c1model2,c0model) << endl;
  }

  void test_write() {
    cout << c0model << endl;
    cout << c1model << endl;
  }


  void test() {
    test_copy_constructor();
    test_reduce();
    test_compose();
    test_write();
    test_add();
  }
};
  

int main() {
  TestAffineModel<Flt>().test();
  return 0;
}
