/***************************************************************************
 *            test_formula.cpp
 *
 *  Copyright  2010-2018  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>

#include "config.hpp"

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "function/formula.hpp"

#include "../test.hpp"

#include "function/function.hpp"

using namespace std;
using namespace Ariadne;

class TestFormula
{
    EffectiveFormula two, x, y;

  public:
    TestFormula()
      : two(EffectiveFormula::constant(2.0)),
        x(EffectiveFormula::coordinate(0)),
        y(EffectiveFormula::coordinate(1)) { }
  private:
    Void test_construct() {
        EffectiveFormula f1(sqrt(pow(x,2)+pow(y,2))/two);
        ARIADNE_TEST_PRINT(f1);

        Vector<EffectiveFormula> f2({sqrt(pow(x,2)+pow(y,2)), atan(y/x)});
        ARIADNE_TEST_PRINT(f2);
    }

    Void test_identical() {
        ARIADNE_TEST_ASSERT(identical(x,x));
        ARIADNE_TEST_ASSERT(identical(ExactFormula::constant(2.0),ExactFormula::constant(2.0)));
        ARIADNE_TEST_ASSERT(not identical(EffectiveFormula::constant(2.0),EffectiveFormula::constant(2.0)));
        ARIADNE_TEST_ASSERT(identical(x*y,y*x));
        ARIADNE_TEST_ASSERT(identical(x+y,y+x));
        ARIADNE_TEST_ASSERT(identical(cos(x),cos(x)));
        ARIADNE_TEST_ASSERT(identical(pow(x,2),pow(x,2)));
    }

    Void test_is_constant_in() {
        ARIADNE_TEST_ASSERT(is_constant_in(x,{1}));
        ARIADNE_TEST_ASSERT(not is_constant_in(x,{0}));
        ARIADNE_TEST_ASSERT(is_constant_in(x*y,{2,3}));
        ARIADNE_TEST_ASSERT(not is_constant_in(cos(x),{0}));
        ARIADNE_TEST_ASSERT(is_constant_in(pow(x,2),{1}));
        ARIADNE_TEST_ASSERT(is_constant_in(two,{0,1}));
        ARIADNE_TEST_ASSERT(not is_constant_in(3*y,{0,1}));
    }

    Void test_is_affine_in() {
        ARIADNE_TEST_ASSERT(is_affine_in(x+y,{0}));
        ARIADNE_TEST_ASSERT(is_affine_in(x+y,{0,1}));
        ARIADNE_TEST_ASSERT(is_affine_in(x+y,{0,1,2}));
        ARIADNE_TEST_ASSERT(is_affine_in(x*y+y,{0}));
        ARIADNE_TEST_ASSERT(is_affine_in(x*y+y,{1}));
        ARIADNE_TEST_ASSERT(not is_affine_in(x*y+y,{0,1}));
        ARIADNE_TEST_ASSERT(is_affine_in(two,{0,1}));
        ARIADNE_TEST_ASSERT(is_affine_in(pow(x,2),{1}));
        ARIADNE_TEST_ASSERT(not is_affine_in(pow(x,2),{0}));
        ARIADNE_TEST_ASSERT(not is_affine_in(x/y,{1}));
        ARIADNE_TEST_ASSERT(is_affine_in(x/y,{0}));
        ARIADNE_TEST_ASSERT(not is_affine_in(sin(x),{0}));
        ARIADNE_TEST_ASSERT(is_affine_in(y*sin(x),{1}));
        ARIADNE_TEST_ASSERT(is_affine_in(Vector<EffectiveFormula>({y*2+x,x+y}),{0,1}));
    }

  public:
    Void test() {
        ARIADNE_TEST_CALL(test_construct());
        ARIADNE_TEST_CALL(test_identical());
        ARIADNE_TEST_CALL(test_is_constant_in());
        ARIADNE_TEST_CALL(test_is_affine_in());
    }
};


Int main() {
    TestFormula().test();
    return ARIADNE_TEST_FAILURES;
}

