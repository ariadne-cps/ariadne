/***************************************************************************
 *            test_bounder.cpp
 *
 *  Copyright  2018  Luca Geretti
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "config.hpp"

#include "solvers/bounder.hpp"
#include "function/formula.hpp"
#include "function/taylor_model.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

class TestBounder
{
  private:
    std::unique_ptr<BounderInterface> bounder_ptr;
  public:
    TestBounder(const BounderInterface& b)
        : bounder_ptr(b.clone())
    {

    }

    Int test() {
        ARIADNE_TEST_PRINT(*bounder_ptr);
        ARIADNE_TEST_CALL(test_without_parameters());
        ARIADNE_TEST_CALL(test_with_parameters());
        return 0;
    }

    Void test_without_parameters() {

    }

    Void test_with_parameters() {

    }
};

Int main(Int argc, const char* argv[]) {

    ARIADNE_PRINT_TEST_CASE_TITLE("EulerBounder");
    TestBounder(EulerBounder()).test();
    ARIADNE_PRINT_TEST_CASE_TITLE("HeunBounder");
    TestBounder(HeunBounder()).test();
    ARIADNE_PRINT_TEST_CASE_TITLE("RalstonBounder");
    TestBounder(RalstonBounder()).test();
    ARIADNE_PRINT_TEST_CASE_TITLE("RungeKutta4Bounder");
    TestBounder(RungeKutta4Bounder()).test();

    return ARIADNE_TEST_FAILURES;
}
