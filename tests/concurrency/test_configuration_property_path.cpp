/***************************************************************************
 *            test_configuration_property_path.cpp
 *
 *  Copyright  2008-20  Luca Geretti
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

#include "symbolic/identifier.hpp"
#include "concurrency/configuration_property_path.hpp"
#include "../test.hpp"

using namespace Ariadne;

class TestConfigurationPropertyPath {
  public:

    void test_construction() {
        ConfigurationPropertyPath p;
        ARIADNE_TEST_EQUALS(p.repr(),"./");
        ConfigurationPropertyPath p2(p);
        ARIADNE_TEST_EQUALS(p2.repr(),"./");
    }

    void test_append() {
        ConfigurationPropertyPath p;
        p.append("child1");
        ARIADNE_TEST_EQUALS(p.repr(),"./child1/");
        p.append("child2");
        ARIADNE_TEST_EQUALS(p.repr(),"./child1/child2/");
    }

    void test_prepend() {
        ConfigurationPropertyPath p;
        p.prepend("child2");
        ARIADNE_TEST_EQUALS(p.repr(),"./child2/");
        p.prepend("child1");
        ARIADNE_TEST_EQUALS(p.repr(),"./child1/child2/");
    }

    void test_subpath() {
        ConfigurationPropertyPath p;
        p.append("child");
        ARIADNE_TEST_EQUALS(p.repr(),"./child/");
        auto ps = p.subpath();
        ARIADNE_TEST_EQUALS(ps.repr(),"./");
    }

    void test_copy() {
        ConfigurationPropertyPath p1;
        p1.append("child1");
        auto p2 = p1;
        p2.append("child2");
        ARIADNE_TEST_EQUALS(p1.repr(),"./child1/");
        ARIADNE_TEST_EQUALS(p2.repr(),"./child1/child2/");
    }

    void test() {
        ARIADNE_TEST_CALL(test_construction());
        ARIADNE_TEST_CALL(test_append());
        ARIADNE_TEST_CALL(test_subpath());
        ARIADNE_TEST_CALL(test_copy());
    }
};

int main() {
    TestConfigurationPropertyPath().test();
    return ARIADNE_TEST_FAILURES;
}
