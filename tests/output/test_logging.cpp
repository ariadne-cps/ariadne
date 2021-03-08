/***************************************************************************
 *            test_logging.cpp
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

#include "config.hpp"
#include "../test.hpp"

#include "output/logging.hpp"

using namespace Ariadne;

void sample_function() {
    ARIADNE_LOG_SCOPE_CREATE;

    ARIADNE_LOG_PRINTLN("val=inf, x0=2.0^3*1.32424242432423[2,3], y>[0.1:0.2] (z={0:1})");
}

class TestLogging {
  public:

    TestLogging() {
        Logger::use_immediate_scheduler();
        Logger::configuration().set_prints_level_on_change_only(false);
    }

    Int test() {
        ARIADNE_TEST_CALL(test_shown_single_print());
        ARIADNE_TEST_CALL(test_hidden_single_print());
        ARIADNE_TEST_CALL(test_shown_call_function_with_entrance_and_exit());
        ARIADNE_TEST_CALL(test_hide_call_function_with_entrance_and_exit());
        ARIADNE_TEST_CALL(test_indents_based_on_level());
        ARIADNE_TEST_CALL(test_dark_theme());
        ARIADNE_TEST_CALL(test_light_theme());
        ARIADNE_TEST_CALL(test_handles_multiline_output());
        ARIADNE_TEST_CALL(test_discards_newlines_and_indentation());
        return 0;
    }

    Void test_shown_single_print() {
        Logger::configuration().set_verbosity(1);
        ARIADNE_LOG_PRINTLN("This is a call on level 1");
    }

    Void test_hidden_single_print() {
        Logger::configuration().set_verbosity(0);
        ARIADNE_LOG_PRINTLN("This is a hidden call on level 1");
    }

    Void test_shown_call_function_with_entrance_and_exit() {
        Logger::configuration().set_verbosity(2);
        Logger::configuration().set_prints_scope_entrance(true);
        Logger::configuration().set_prints_scope_exit(true);
        ARIADNE_LOG_PRINTLN("This is a call on level 1");
        ARIADNE_LOG_RUN_AT(0,sample_function());
        ARIADNE_LOG_PRINTLN("This is again a call on level 1");
    }

    Void test_hide_call_function_with_entrance_and_exit() {
        Logger::configuration().set_verbosity(2);
        Logger::configuration().set_prints_scope_entrance(true);
        Logger::configuration().set_prints_scope_exit(true);
        ARIADNE_LOG_PRINTLN("This is a call on level 1");
        ARIADNE_LOG_RUN_AT(1,sample_function());
        ARIADNE_LOG_PRINTLN("This is again a call on level 1");
    }

    Void test_indents_based_on_level() {
        Logger::configuration().set_verbosity(2);
        Logger::configuration().set_indents_based_on_level(true);
        ARIADNE_LOG_PRINTLN("Call at level 1");
        ARIADNE_LOG_PRINTLN_AT(1,"Call at level 2");
        Logger::configuration().set_indents_based_on_level(false);
        ARIADNE_LOG_PRINTLN("Call at level 1");
        ARIADNE_LOG_PRINTLN_AT(1,"Call at level 2");
    }

    Void test_handles_multiline_output() {
        Logger::configuration().set_verbosity(2);
        Logger::configuration().set_handles_multiline_output(true);
        ARIADNE_LOG_PRINTLN("This is a very long string for this test that will most definitely be longer than just one line, at least if the number of columns in the terminal is not excessively large, but this should suffice I believe. Just to be sure, let's add some more characters to the line and the result should get into a second line.");
        Logger::configuration().set_handles_multiline_output(false);
        ARIADNE_LOG_PRINTLN("This is a very long string for this test that will most definitely be longer than just one line, at least if the number of columns in the terminal is not excessively large, but this should suffice I believe. Just to be sure, let's add some more characters to the line and the result should get into a second line.");
    }

    Void test_discards_newlines_and_indentation() {
        Logger::configuration().set_verbosity(2);
        Logger::configuration().set_discards_newlines_and_indentation(true);
        ARIADNE_LOG_PRINTLN("This text should just be in a single line \n       with no extra whitespaces.");
        Logger::configuration().set_discards_newlines_and_indentation(false);
        ARIADNE_LOG_PRINTLN("This text should be in two lines\n          where the second one starts several whitespaces in.");
    }

    Void test_dark_theme() {
        Logger::configuration().set_verbosity(2);
        Logger::configuration().set_theme(TT_THEME_DARK);
        ARIADNE_LOG_PRINTLN("This is a call on level 1");
        ARIADNE_LOG_RUN_AT(0,sample_function());
        ARIADNE_LOG_PRINTLN("This is again a call on level 1");
    }

    Void test_light_theme() {
        Logger::configuration().set_verbosity(2);
        Logger::configuration().set_theme(TT_THEME_LIGHT);
        ARIADNE_LOG_PRINTLN("This is a call on level 1");
        ARIADNE_LOG_RUN_AT(0,sample_function());
        ARIADNE_LOG_PRINTLN("This is again a call on level 1");
    }

};

Int main(Int argc, const char* argv[]) {

    TestLogging().test();

    return 0;
}

