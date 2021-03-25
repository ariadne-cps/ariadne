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
#include "output/progress_indicator.hpp"
#include "output/logging.hpp"

using namespace Ariadne;

void sample_function() {
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("val=inf, x0=2.0^3*1.32424242432423[2,3], y>[0.1:0.2] (z={0:1})");
}

Void print_something1() {
    ARIADNE_LOG_PRINTLN("This is a call from thread id " << std::this_thread::get_id());
}

Void print_something2() {
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("This is a call from thread id " << std::this_thread::get_id());
}

class TestLogging {
  public:

    TestLogging() {
        Logger::instance().configuration().set_prints_level_on_change_only(false);
    }

    Int test() {
        ARIADNE_TEST_CALL(test_print_configuration());
        ARIADNE_TEST_CALL(test_shown_single_print());
        ARIADNE_TEST_CALL(test_hidden_single_print());
        ARIADNE_TEST_CALL(test_use_blocking_scheduler());
        ARIADNE_TEST_CALL(test_use_nonblocking_scheduler());
        ARIADNE_TEST_CALL(test_shown_call_function_with_entrance_and_exit());
        ARIADNE_TEST_CALL(test_hide_call_function_with_entrance_and_exit());
        ARIADNE_TEST_CALL(test_indents_based_on_level());
        ARIADNE_TEST_CALL(test_hold_line());
        ARIADNE_TEST_CALL(test_dark_theme());
        ARIADNE_TEST_CALL(test_light_theme());
        ARIADNE_TEST_CALL(test_handles_multiline_output());
        ARIADNE_TEST_CALL(test_discards_newlines_and_indentation());
        ARIADNE_TEST_CALL(test_multiple_threads());
        ARIADNE_TEST_CALL(test_redirect());
        return 0;
    }

    Void test_print_configuration() {
        Logger::instance().use_immediate_scheduler();
        ARIADNE_LOG_SET_VERBOSITY(1);
        ARIADNE_LOG_PRINTLN(Logger::instance().configuration());
    }

    Void test_shown_single_print() {
        Logger::instance().use_immediate_scheduler();
        ARIADNE_LOG_SET_VERBOSITY(1);
        ARIADNE_LOG_PRINTLN("This is a call on level 1");
    }

    Void test_hidden_single_print() {
        Logger::instance().use_immediate_scheduler();
        ARIADNE_LOG_SET_VERBOSITY(0);
        ARIADNE_LOG_PRINTLN("This is a hidden call on level 1");
    }

    Void test_use_blocking_scheduler() {
        Logger::instance().use_blocking_scheduler();
        ARIADNE_LOG_SET_VERBOSITY(1);
        ARIADNE_LOG_PRINTLN("This is a call");
        ARIADNE_LOG_PRINTLN("This is another call");
    }

    Void test_use_nonblocking_scheduler() {
        Logger::instance().use_nonblocking_scheduler();
        ARIADNE_LOG_SET_VERBOSITY(1);
        ARIADNE_LOG_PRINTLN("This is a call");
        ARIADNE_LOG_PRINTLN("This is another call");
    }

    Void test_shown_call_function_with_entrance_and_exit() {
        Logger::instance().use_immediate_scheduler();
        ARIADNE_LOG_SET_VERBOSITY(2);
        Logger::instance().configuration().set_prints_scope_entrance(true);
        Logger::instance().configuration().set_prints_scope_exit(true);
        ARIADNE_LOG_PRINTLN("This is a call on level 1");
        ARIADNE_LOG_RUN_AT(0,sample_function());
        ARIADNE_LOG_PRINTLN("This is again a call on level 1");
    }

    Void test_hide_call_function_with_entrance_and_exit() {
        Logger::instance().use_immediate_scheduler();
        ARIADNE_LOG_SET_VERBOSITY(2);
        Logger::instance().configuration().set_prints_scope_entrance(true);
        Logger::instance().configuration().set_prints_scope_exit(true);
        ARIADNE_LOG_PRINTLN("This is a call on level 1");
        ARIADNE_LOG_RUN_AT(1,sample_function());
        ARIADNE_LOG_PRINTLN("This is again a call on level 1");
    }

    Void test_indents_based_on_level() {
        Logger::instance().use_immediate_scheduler();
        ARIADNE_LOG_SET_VERBOSITY(2);
        Logger::instance().configuration().set_indents_based_on_level(true);
        ARIADNE_LOG_PRINTLN("Call at level 1");
        ARIADNE_LOG_PRINTLN_AT(1,"Call at level 2");
        Logger::instance().configuration().set_indents_based_on_level(false);
        ARIADNE_LOG_PRINTLN("Call at level 1");
        ARIADNE_LOG_PRINTLN_AT(1,"Call at level 2");
    }

    Void test_handles_multiline_output() {
        Logger::instance().use_immediate_scheduler();
        ARIADNE_LOG_SET_VERBOSITY(2);
        Logger::instance().configuration().set_handles_multiline_output(true);
        ARIADNE_LOG_PRINTLN("This is a very long string for this test that will most definitely be longer than just one line, at least if the number of columns in the terminal is not excessively large, but this should suffice I believe. Just to be sure, let's add some more characters to the line and the result should get into a second line.");
        Logger::instance().configuration().set_handles_multiline_output(false);
        ARIADNE_LOG_PRINTLN("This is a very long string for this test that will most definitely be longer than just one line, at least if the number of columns in the terminal is not excessively large, but this should suffice I believe. Just to be sure, let's add some more characters to the line and the result should get into a second line.");
    }

    Void test_discards_newlines_and_indentation() {
        Logger::instance().use_immediate_scheduler();
        ARIADNE_LOG_SET_VERBOSITY(2);
        Logger::instance().configuration().set_discards_newlines_and_indentation(true);
        ARIADNE_LOG_PRINTLN("This text should just be in a single line \n       with no extra whitespaces.");
        Logger::instance().configuration().set_discards_newlines_and_indentation(false);
        ARIADNE_LOG_PRINTLN("This text should be in two lines\n          where the second one starts several whitespaces in.");
    }

    Void test_hold_line() {
        Logger::instance().use_immediate_scheduler();
        ARIADNE_LOG_SCOPE_CREATE;
        ProgressIndicator indicator(10.0);
        for (unsigned int i=0; i<10; ++i) {
            indicator.update_current(i);
            std::this_thread::sleep_for(std::chrono::milliseconds(20));
            ARIADNE_LOG_SCOPE_PRINTHOLD("[" << indicator.symbol() << "] " << indicator.percentage() << "%");
        }
    }

    Void test_dark_theme() {
        Logger::instance().use_immediate_scheduler();
        ARIADNE_LOG_SET_VERBOSITY(2);
        Logger::instance().configuration().set_theme(TT_THEME_DARK);
        std::clog << TT_THEME_DARK << std::endl;
        ARIADNE_LOG_PRINTLN("This is a call on level 1");
        ARIADNE_LOG_RUN_AT(0,sample_function());
        ARIADNE_LOG_PRINTLN("This is again a call on level 1");
        Logger::instance().configuration().set_theme(TT_THEME_NONE);
    }

    Void test_light_theme() {
        Logger::instance().use_immediate_scheduler();
        ARIADNE_LOG_SET_VERBOSITY(2);
        Logger::instance().configuration().set_theme(TT_THEME_LIGHT);
        std::clog << TT_THEME_LIGHT << std::endl;
        ARIADNE_LOG_PRINTLN("This is a call on level 1");
        ARIADNE_LOG_RUN_AT(0,sample_function());
        ARIADNE_LOG_PRINTLN("This is again a call on level 1");
        Logger::instance().configuration().set_theme(TT_THEME_NONE);
    }

    Void test_redirect() {
        Logger::instance().use_immediate_scheduler();
        ARIADNE_LOG_SET_VERBOSITY(1);
        ARIADNE_LOG_PRINTLN("This is call 1");
        Logger::instance().redirect_to_file("log.txt");
        ARIADNE_LOG_PRINTLN("This is call 2");
        ARIADNE_LOG_PRINTLN("This is call 3");
        Logger::instance().redirect_to_console();
        ARIADNE_LOG_PRINTLN("This is call 4");
        ARIADNE_LOG_PRINTLN("This is call 5");

        std::string line;
        std::ifstream file("log.txt");
        unsigned int count = 0;
        if(file.is_open()) {
            while(!file.eof()) {
                getline(file,line);
                count++;
            }
            file.close();
        }
        ARIADNE_TEST_EQUALS(count,3);
    }

    Void test_multiple_threads() {
        Logger::instance().use_nonblocking_scheduler();
        ARIADNE_LOG_SET_VERBOSITY(3);
        Logger::instance().configuration().set_thread_name_printing_policy(ThreadNamePrintingPolicy::BEFORE);
        ARIADNE_LOG_PRINTLN("Printing on the main thread without other threads");
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        std::thread thread1([]() { print_something1(); }), thread2([]() { print_something2(); });
        Logger::instance().register_thread(thread1.get_id(),"thr1");
        Logger::instance().register_thread(thread2.get_id(),"thr2");
        ARIADNE_LOG_PRINTLN("Printing again on the main thread, but with other threads");
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        thread1.join();
        thread2.join();
        Logger::instance().unregister_thread(thread1.get_id());
        Logger::instance().unregister_thread(thread2.get_id());
    }
};

Int main(Int argc, const char* argv[]) {

    TestLogging().test();

    return 0;
}

