/***************************************************************************
 *            test_function_patch.cpp
 *
 *  Copyright  2024  Pieter Collins
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

#include "config.hpp"
#include "algebra/algebra.hpp"
#include "function/taylor_model.hpp"
#include "function/taylor_function.hpp"
#include "function/function_patch.hpp"
#include "function/domain.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/expression_patch.hpp"
#include "symbolic/expression_set.hpp"
#include "symbolic/space.hpp"

#include "../test.hpp"

using namespace std;
using namespace Ariadne;


class TestFunctionPatch
{
  public:
    Void test();
  private:
};

Void TestFunctionPatch::test() {
}

Int main() {
    TestFunctionPatch().test();

    return ARIADNE_TEST_FAILURES;
}
