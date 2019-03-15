/***************************************************************************
 *            test_enclosure.cpp
 *
 *  Copyright  2019  Luca Geretti
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
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "config.hpp"
#include "function/taylor_function.hpp"
#include "function/constraint.hpp"
#include "geometry/enclosure.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

class TestEnclosure
{
  public:
    Void test();
  private:
    Void test_recombine();
};


Void
TestEnclosure::test()
{
    ARIADNE_TEST_CALL(test_recombine());
}


Void
TestEnclosure::test_recombine()
{
    TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(DoublePrecision(),1e-8));
    EnclosureConfiguration config(function_factory);

    Enclosure e1(ExactBoxType{{1,3},{0,1}},config);
    Enclosure e2(ExactBoxType{{2,4},{0,2}},config);

    ARIADNE_TEST_PRINT(recombine({e1,e2}));
}


Int main() {

    TestEnclosure().test();

    return ARIADNE_TEST_FAILURES;
}

