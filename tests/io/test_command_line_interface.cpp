/***************************************************************************
 *            test_command_line_interface.cpp
 *
 *  Copyright  2008-21  Luca Geretti
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
#include "io/command_line_interface.hpp"
#include "../test.hpp"

using namespace Ariadne;

class TestCommandLineInterface {
  public:

    Void test() {
        ARIADNE_TEST_CALL(test_empty_argument_stream())
        ARIADNE_TEST_CALL(test_nonempty_argument_stream())
        ARIADNE_TEST_CALL(test_cli_instantiation())
    }

    void test_empty_argument_stream() {
        const char* argv[] = {nullptr};
        auto stream = ArgumentStream(1,argv);
        ARIADNE_TEST_ASSERT(stream.empty())
        ARIADNE_TEST_EQUALS(stream.size(),0)
        ARIADNE_TEST_FAIL(stream.peek())
        ARIADNE_TEST_FAIL(stream.pop())
    }

    void test_nonempty_argument_stream() {
        const char* argv[] = {nullptr, "a", "b"};
        auto stream = ArgumentStream(3,argv);
        ARIADNE_TEST_ASSERT(not stream.empty())
        ARIADNE_TEST_EQUALS(stream.size(),2)
        ARIADNE_TEST_EQUALS(stream.peek(),"a")
        auto val1 = stream.pop();
        ARIADNE_TEST_EQUALS(val1,"a")
        auto val2 = stream.pop();
        ARIADNE_TEST_EQUALS(val2,"b")
        ARIADNE_TEST_ASSERT(stream.empty())
        ARIADNE_TEST_FAIL(stream.peek())
        ARIADNE_TEST_FAIL(stream.pop())
    }

    void test_cli_instantiation() {
        ARIADNE_TEST_EXECUTE(CommandLineInterface::instance())
    }

};

Int main(Int argc, const char* argv[]) {
    TestCommandLineInterface().test();
    return 0;
}
