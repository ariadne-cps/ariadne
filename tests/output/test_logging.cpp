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
        ARIADNE_TEST_CALL(test_dark_theme());
        ARIADNE_TEST_CALL(test_light_theme());
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

