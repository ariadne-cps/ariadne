/***************************************************************************
 *            test_configuration.cpp
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

#include "solvers/configuration_interface.hpp"
#include "../test.hpp"

using namespace Ariadne;

class A;

template<> class Configuration<A> : public ConfigurationInterface {
public:
    Configuration() : _int_parameter(0), _bool_parameter(false) { }

    int int_parameter() const { return _int_parameter; }
    bool bool_parameter() const { return _bool_parameter; }
    void set_int_parameter(int v) { _int_parameter = v; }
    void set_bool_parameter(bool v) { _bool_parameter = v; }

    OutputStream& _write(OutputStream& os) const override { os << "int_parameter = " << _int_parameter << ", bool_parameter = " << _bool_parameter; return os; }

private:
    int _int_parameter;
    bool _bool_parameter;
};

class A : public Configurable<A>, public WritableInterface {
  public:
    A() : Configurable<A>(Configuration<A>()) { }
    OutputStream& _write(OutputStream& os) const override { os << "(A's configuration:" << configuration() << ")"; return os; }
};

class TestConfiguration {
  public:

    void test_simple_configuration() {
        A a;
        ARIADNE_TEST_PRINT(a);
    }

    void test_nested_configuration() {

    }

    void test() {
        ARIADNE_TEST_CALL(test_simple_configuration());
        ARIADNE_TEST_CALL(test_nested_configuration());
    }
};

int main() {

    TestConfiguration().test();
    return ARIADNE_TEST_FAILURES;
}
