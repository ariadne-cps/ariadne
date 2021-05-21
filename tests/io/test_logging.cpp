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

#include <sys/ioctl.h>
#include <unistd.h>
#include "config.hpp"
#include "concurrency/thread.hpp"
#include "concurrency/buffered_thread.hpp"
#include "io/progress_indicator.hpp"
#include "io/logging.hpp"
#include "../test.hpp"

using namespace Ariadne;

void sample_function() {
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("val=inf, x0=2.0^3*1.32424242432423[2,3], y>[0.1:0.2] (z={0:1}), 1, x0, x11, true@1.");
}

Void print_something1() {
    ARIADNE_LOG_PRINTLN("This is a call from thread id " << std::this_thread::get_id() << " named '" << Logger::instance().current_thread_name() << "'");
}

Void print_something2() {
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("This is a call from thread id " << std::this_thread::get_id() << " named '" << Logger::instance().current_thread_name() << "'");
}

class TestLogging {
  public:

    TestLogging() {
        Logger::instance().configuration().set_prints_level_on_change_only(false);
    }

    Void test() {
        ARIADNE_TEST_CALL(test_print_configuration())
        ARIADNE_TEST_CALL(test_shown_single_print())
        ARIADNE_TEST_CALL(test_hidden_single_print())
        ARIADNE_TEST_CALL(test_use_blocking_scheduler())
        ARIADNE_TEST_CALL(test_use_nonblocking_scheduler())
        ARIADNE_TEST_CALL(test_shown_call_function_with_entrance_and_exit())
        ARIADNE_TEST_CALL(test_hide_call_function_with_entrance_and_exit())
        ARIADNE_TEST_CALL(test_indents_based_on_level())
        ARIADNE_TEST_CALL(test_hold_line())
        ARIADNE_TEST_CALL(test_hold_long_line())
        ARIADNE_TEST_CALL(test_light_theme())
        ARIADNE_TEST_CALL(test_dark_theme())
        ARIADNE_TEST_CALL(test_theme_custom_keyword())
        ARIADNE_TEST_CALL(test_theme_bgcolor_bold_underline())
        ARIADNE_TEST_CALL(test_handles_multiline_output())
        ARIADNE_TEST_CALL(test_discards_newlines_and_indentation())
        ARIADNE_TEST_CALL(test_redirect())
        ARIADNE_TEST_CALL(test_multiple_threads_with_blocking_scheduler())
        ARIADNE_TEST_CALL(test_multiple_threads_with_nonblocking_scheduler())
    }

    Void test_print_configuration() {
        Logger::instance().use_immediate_scheduler();
        Logger::instance().configuration().set_verbosity(1);
        ARIADNE_LOG_PRINTLN(Logger::instance().configuration());
    }

    Void test_shown_single_print() {
        Logger::instance().use_immediate_scheduler();
        Logger::instance().configuration().set_verbosity(1);
        ARIADNE_LOG_PRINTLN("This is a call on level 1");
    }

    Void test_hidden_single_print() {
        Logger::instance().use_immediate_scheduler();
        Logger::instance().configuration().set_verbosity(0);
        ARIADNE_LOG_PRINTLN("This is a hidden call on level 1");
    }

    Void test_use_blocking_scheduler() {
        Logger::instance().use_blocking_scheduler();
        Logger::instance().configuration().set_verbosity(1);
        ARIADNE_LOG_PRINTLN("This is a call");
        ARIADNE_LOG_PRINTLN("This is another call");
    }

    Void test_use_nonblocking_scheduler() {
        Logger::instance().use_nonblocking_scheduler();
        Logger::instance().configuration().set_verbosity(1);
        ARIADNE_LOG_PRINTLN("This is a call");
        ARIADNE_LOG_PRINTLN("This is another call");
    }

    Void test_shown_call_function_with_entrance_and_exit() {
        Logger::instance().use_immediate_scheduler();
        Logger::instance().configuration().set_verbosity(2);
        Logger::instance().configuration().set_prints_scope_entrance(true);
        Logger::instance().configuration().set_prints_scope_exit(true);
        ARIADNE_LOG_PRINTLN("This is a call on level 1");
        ARIADNE_LOG_RUN_AT(0,sample_function());
        ARIADNE_LOG_PRINTLN("This is again a call on level 1");
    }

    Void test_hide_call_function_with_entrance_and_exit() {
        Logger::instance().use_immediate_scheduler();
        Logger::instance().configuration().set_verbosity(2);
        Logger::instance().configuration().set_prints_scope_entrance(false);
        Logger::instance().configuration().set_prints_scope_exit(false);
        ARIADNE_LOG_PRINTLN("This is a call on level 1");
        ARIADNE_LOG_RUN_AT(1,sample_function());
        ARIADNE_LOG_PRINTLN("This is again a call on level 1");
    }

    Void test_indents_based_on_level() {
        Logger::instance().use_immediate_scheduler();
        Logger::instance().configuration().set_verbosity(2);
        Logger::instance().configuration().set_indents_based_on_level(true);
        ARIADNE_LOG_PRINTLN("Call at level 1");
        ARIADNE_LOG_PRINTLN_AT(1,"Call at level 2");
        Logger::instance().configuration().set_indents_based_on_level(false);
        ARIADNE_LOG_PRINTLN("Call at level 1");
        ARIADNE_LOG_PRINTLN_AT(1,"Call at level 2");
    }

    Void test_handles_multiline_output() {
        Logger::instance().use_immediate_scheduler();
        Logger::instance().configuration().set_verbosity(2);
        Logger::instance().configuration().set_handles_multiline_output(true);
        ARIADNE_LOG_PRINTLN("This is a very long string for this test that will most definitely be longer than just one line, at least if the number of columns in the terminal is not excessively large, but this should suffice I believe. Just to be sure, let's add some more characters to the line and the result should get into a second line.");
        Logger::instance().configuration().set_handles_multiline_output(false);
        ARIADNE_LOG_PRINTLN("This is a very long string for this test that will most definitely be longer than just one line, at least if the number of columns in the terminal is not excessively large, but this should suffice I believe. Just to be sure, let's add some more characters to the line and the result should get into a second line.");
    }

    Void test_discards_newlines_and_indentation() {
        Logger::instance().use_immediate_scheduler();
        Logger::instance().configuration().set_verbosity(2);
        Logger::instance().configuration().set_discards_newlines_and_indentation(true);
        ARIADNE_LOG_PRINTLN("This text should just be in a single line \n       with no extra whitespaces.");
        Logger::instance().configuration().set_discards_newlines_and_indentation(false);
        ARIADNE_LOG_PRINTLN("This text should be in two lines\n          where the second one starts several whitespaces in.");
    }

    Void _hold_short_line() {
        ARIADNE_LOG_SCOPE_CREATE;
        ProgressIndicator indicator(10.0);
        for (unsigned int i=0; i<10; ++i) {
            indicator.update_current(i);
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            ARIADNE_LOG_SCOPE_PRINTHOLD("[" << indicator.symbol() << "] " << indicator.percentage() << "%");
        }
    }

    Void test_hold_line() {
        Logger::instance().use_immediate_scheduler();
        _hold_short_line();
        Logger::instance().use_blocking_scheduler();
        _hold_short_line();
        Logger::instance().use_nonblocking_scheduler();
        _hold_short_line();
    }

    Void test_hold_long_line() {
        Logger::instance().use_immediate_scheduler();
        ARIADNE_LOG_SCOPE_CREATE;

        const unsigned int DEFAULT_COLUMNS = 80;
        struct winsize ws;
        ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws);
        SizeType num_cols = (ws.ws_col > 0 ? ws.ws_col : DEFAULT_COLUMNS);

        std::string exactly_str(num_cols-4,'x'); // exactly the length required to fill the columns (given a prefix of 4 chars)
        std::string larger_str(num_cols,'x'); // larger enough

        for (unsigned int i=0; i<10; ++i) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            ARIADNE_LOG_SCOPE_PRINTHOLD(exactly_str);
        }

        for (unsigned int i=0; i<10; ++i) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            ARIADNE_LOG_SCOPE_PRINTHOLD(larger_str);
        }
    }

    Void test_light_theme() {
        Logger::instance().use_immediate_scheduler();
        Logger::instance().configuration().set_verbosity(2);
        Logger::instance().configuration().set_theme(TT_THEME_LIGHT);
        std::clog << TT_THEME_LIGHT << std::endl;
        ARIADNE_LOG_PRINTLN("This is a call on level 1");
        ARIADNE_LOG_RUN_AT(0,sample_function());
        ARIADNE_LOG_PRINTLN("This is again a call on level 1");
        Logger::instance().configuration().set_theme(TT_THEME_NONE);
    }

    Void test_dark_theme() {
        Logger::instance().use_immediate_scheduler();
        Logger::instance().configuration().set_verbosity(2);
        Logger::instance().configuration().set_theme(TT_THEME_DARK);
        std::clog << TT_THEME_DARK << std::endl;
        ARIADNE_LOG_PRINTLN("This is a call on level 1");
        ARIADNE_LOG_RUN_AT(0,sample_function());
        ARIADNE_LOG_PRINTLN("This is again a call on level 1");
    }

    Void test_theme_custom_keyword() {
        Logger::instance().use_immediate_scheduler();
        Logger::instance().configuration().set_verbosity(1);
        Logger::instance().configuration().add_custom_keyword("first");
        Logger::instance().configuration().add_custom_keyword("second",TT_STYLE_CREAM);
        ARIADNE_LOG_PRINTLN("This is a first custom keyword and a custom styled second one.");
    }

    Void test_theme_bgcolor_bold_underline() {
        Logger::instance().use_immediate_scheduler();
        std::list<TerminalTextStyle> styles;
        styles.push_back(TerminalTextStyle(0,0,true,false));
        styles.push_back(TerminalTextStyle(0,0,false,true));
        styles.push_back(TerminalTextStyle(0,0,true,true));
        styles.push_back(TerminalTextStyle(0,88,false,false));
        styles.push_back(TerminalTextStyle(0,88,true,false));
        styles.push_back(TerminalTextStyle(0,88,false,true));
        styles.push_back(TerminalTextStyle(0,88,true,true));
        std::ostringstream ss;
        for (auto style: styles) {
            ss << style() << "x" << TerminalTextStyle::RESET << " ";
        }
        std::clog << ss.str() << std::endl;
    }

    Void test_redirect() {
        Logger::instance().use_immediate_scheduler();
        Logger::instance().configuration().set_verbosity(1);
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

    Void test_multiple_threads_with_blocking_scheduler() {
        Logger::instance().use_blocking_scheduler();
        Logger::instance().configuration().set_verbosity(3);
        Logger::instance().configuration().set_thread_name_printing_policy(ThreadNamePrintingPolicy::BEFORE);
        ARIADNE_LOG_PRINTLN("Printing on the " << Logger::instance().current_thread_name() << " thread without other threads");
        ARIADNE_TEST_EQUALS(Logger::instance().cached_last_printed_thread_name().compare("main"),0);
        BufferedThread thread1("thr1"), thread2("thr2");
        thread1.enqueue([] { print_something1(); });
        thread2.enqueue([] { print_something2(); });
        ARIADNE_LOG_PRINTLN("Printing again on the main thread, but with other threads");
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        ARIADNE_TEST_PRINT(Logger::instance().cached_last_printed_thread_name());
        ARIADNE_TEST_ASSERT(Logger::instance().cached_last_printed_thread_name().compare("thr1") == 0 or
                                    Logger::instance().cached_last_printed_thread_name().compare("thr2") == 0);
    }

    Void test_multiple_threads_with_nonblocking_scheduler() {
        Logger::instance().use_nonblocking_scheduler();
        Logger::instance().configuration().set_theme(TT_THEME_DARK);
        Logger::instance().configuration().set_verbosity(3);
        Logger::instance().configuration().set_thread_name_printing_policy(ThreadNamePrintingPolicy::BEFORE);

        ARIADNE_LOG_PRINTLN("Printing on the " << Logger::instance().current_thread_name() << " thread without other threads");
        BufferedThread thread1, thread2, thread3, thread4, thread5, thread6;
        thread1.enqueue([] { print_something1(); });
        thread2.enqueue([] { print_something1(); });
        thread3.enqueue([] { print_something1(); });
        thread4.enqueue([] { print_something1(); });
        thread5.enqueue([] { print_something1(); });
        thread6.enqueue([] { print_something1(); });
        ARIADNE_LOG_PRINTLN("Printing again on the main thread, but with other threads");
    }
};

Int main(Int argc, const char* argv[]) {

    TestLogging().test();

    return ARIADNE_TEST_FAILURES;
}

