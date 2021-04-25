/***************************************************************************
 *            test_array.cpp
 *
 *  Copyright  2009-21  Luca Geretti
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

#include "config.hpp"
#include "utility/array.hpp"
#include "utility/container.hpp"

#include "../test.hpp"

using namespace Ariadne;

struct TestConvertibleTo {
    TestConvertibleTo(int a_) : a(a_) { }
    int a;
};

struct TestClass {
    TestClass(int a_) : a(a_) { }
    explicit TestClass(TestConvertibleTo const& c) : TestClass(c.a) { }
    int a;
};

class TestArray {
  public:

    Void test_convert() {
        Array<TestClass> tca = {TestClass(1), TestClass(2)};
        Array<TestConvertibleTo> tcta = {TestConvertibleTo(1), TestConvertibleTo(2)};
        ARIADNE_TEST_EXECUTE(Array<TestClass> tcac(tcta));
    }

    Void test_print() {
        Array<int> a1;
        ARIADNE_TEST_PRINT(a1);
        Array<int> a2 = {1, 2};
        ARIADNE_TEST_PRINT(a2);
    }

    Void test() {
        ARIADNE_TEST_CALL(test_convert());
        ARIADNE_TEST_CALL(test_print());
    }

};

Int main() {
    TestArray().test();
    return ARIADNE_TEST_FAILURES;
}
