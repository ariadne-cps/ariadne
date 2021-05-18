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

#include "io/command_line_interface.hpp"
#include "io/logging.hpp"
#include "concurrency/task_manager.hpp"
#include "../test.hpp"

using namespace Ariadne;

class TestCommandLineInterface {
  public:

    Void test() {
        ARIADNE_TEST_CALL(test_empty_argument_stream())
        ARIADNE_TEST_CALL(test_nonempty_argument_stream())
        ARIADNE_TEST_CALL(test_cli_instantiation())
        ARIADNE_TEST_CALL(test_verbosity_parsing())
        ARIADNE_TEST_CALL(test_concurrency_parsing())
        ARIADNE_TEST_CALL(test_multiple_argument_parsing())
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

    void test_verbosity_parsing() {
        const char* argv[] = {nullptr, "-v", "5"};
        Bool success = CommandLineInterface::instance().acquire(3,argv);
        ARIADNE_TEST_ASSERT(success)
        ARIADNE_TEST_EQUALS(Logger::instance().configuration().verbosity(),5)
        const char* argv2[] = {nullptr, "--verbosity", "0"};
        Bool success2 = CommandLineInterface::instance().acquire(3,argv2);
        ARIADNE_TEST_ASSERT(success2)
        ARIADNE_TEST_EQUALS(Logger::instance().configuration().verbosity(),0)
        const char* argv3[] = {nullptr, "-v", "-2"};
        Bool success3 = CommandLineInterface::instance().acquire(3,argv3);
        ARIADNE_TEST_ASSERT(not success3)
        const char* argv4[] = {nullptr, "-v", "q"};
        Bool success4 = CommandLineInterface::instance().acquire(3,argv4);
        ARIADNE_TEST_ASSERT(not success4)
        const char* argv5[] = {nullptr, "-v"};
        Bool success5 = CommandLineInterface::instance().acquire(2,argv5);
        ARIADNE_TEST_ASSERT(not success5)
    }

    void test_concurrency_parsing() {
        const char* argv[] = {nullptr, "-c", "5"};
        Bool success = CommandLineInterface::instance().acquire(3,argv);
        ARIADNE_TEST_ASSERT(success)
        ARIADNE_TEST_EQUALS(TaskManager::instance().concurrency(),5)
        const char* argv2[] = {nullptr, "--concurrency", "0"};
        Bool success2 = CommandLineInterface::instance().acquire(3,argv2);
        ARIADNE_TEST_ASSERT(success2)
        ARIADNE_TEST_EQUALS(TaskManager::instance().concurrency(),0)
        const char* argv3[] = {nullptr, "-c", "-2"};
        Bool success3 = CommandLineInterface::instance().acquire(3,argv3);
        ARIADNE_TEST_ASSERT(not success3)
        const char* argv4[] = {nullptr, "-c", "q"};
        Bool success4 = CommandLineInterface::instance().acquire(3,argv4);
        ARIADNE_TEST_ASSERT(not success4)
        const char* argv5[] = {nullptr, "-c"};
        Bool success5 = CommandLineInterface::instance().acquire(2,argv5);
        ARIADNE_TEST_ASSERT(not success5)
        const char* argv6[] = {nullptr, "-c", "max"};
        Bool success6 = CommandLineInterface::instance().acquire(3,argv6);
        ARIADNE_TEST_ASSERT(success6)
        ARIADNE_TEST_EQUALS(TaskManager::instance().concurrency(),TaskManager::instance().maximum_concurrency())
    }

    void test_multiple_argument_parsing() {
        const char* argv[] = {nullptr, "-c", "2", "--verbosity", "4"};
        Bool success = CommandLineInterface::instance().acquire(5,argv);
        ARIADNE_TEST_ASSERT(success)
        ARIADNE_TEST_EQUALS(TaskManager::instance().concurrency(),2)
        ARIADNE_TEST_EQUALS(Logger::instance().configuration().verbosity(),4)
    }

};

Int main(Int argc, const char* argv[]) {
    TestCommandLineInterface().test();
    return 0;
}

