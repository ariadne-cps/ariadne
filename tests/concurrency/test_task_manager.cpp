/***************************************************************************
 *            test_task_manager.cpp
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

#include <thread>
#include "concurrency/task_manager.hpp"
#include "../test.hpp"

using namespace Ariadne;
using namespace std::chrono_literals;

class TestWorkloadAdvancement {
  public:

    void test_set_concurrency() {
        auto max_concurrency = TaskManager::instance().maximum_concurrency();
        TaskManager::instance().set_concurrency(max_concurrency);
        ARIADNE_TEST_EQUALS(TaskManager::instance().concurrency(), max_concurrency)
        TaskManager::instance().set_maximum_concurrency();
        ARIADNE_TEST_EQUALS(TaskManager::instance().concurrency(), max_concurrency)
        ARIADNE_TEST_FAIL(TaskManager::instance().set_concurrency(1 + max_concurrency))
    }

    void test_run_task_with_one_thread() {
        TaskManager::instance().set_concurrency(1);
        int a = 10;
        auto result = TaskManager::instance().enqueue([&a]{ return a * a; }).get();
        ARIADNE_TEST_EQUALS(result,100)
    }

    void test_run_task_with_multiple_threads() {
        TaskManager::instance().set_concurrency(TaskManager::instance().maximum_concurrency());
        int a = 10;
        auto result = TaskManager::instance().enqueue([&a]{ return a * a; }).get();
        ARIADNE_TEST_EQUALS(result,100)
    }

    void test_run_task_with_no_threads() {
        TaskManager::instance().set_concurrency(0);
        int a = 10;
        auto result = TaskManager::instance().enqueue([&a]{ return a * a; }).get();
        ARIADNE_TEST_EQUALS(result,100)
    }

    void test_change_concurrency_and_log_scheduler() {
        ARIADNE_TEST_EXECUTE(TaskManager::instance().set_concurrency(1))
        ARIADNE_TEST_FAIL(TaskManager::instance().set_logging_immediate_scheduler())
        ARIADNE_TEST_FAIL(TaskManager::instance().set_logging_blocking_scheduler())
        ARIADNE_TEST_FAIL(TaskManager::instance().set_logging_nonblocking_scheduler())
        ARIADNE_TEST_EXECUTE(TaskManager::instance().set_concurrency(0))
        ARIADNE_TEST_EXECUTE(TaskManager::instance().set_logging_immediate_scheduler())
        ARIADNE_TEST_EXECUTE(TaskManager::instance().set_logging_blocking_scheduler())
        ARIADNE_TEST_EXECUTE(TaskManager::instance().set_logging_nonblocking_scheduler())
        ARIADNE_TEST_EXECUTE(TaskManager::instance().set_concurrency(1))
        ARIADNE_TEST_EXECUTE(TaskManager::instance().set_concurrency(0))
    }

    void test() {
        ARIADNE_TEST_CALL(test_set_concurrency())
        ARIADNE_TEST_CALL(test_run_task_with_one_thread())
        ARIADNE_TEST_CALL(test_run_task_with_multiple_threads())
        ARIADNE_TEST_CALL(test_run_task_with_no_threads())
        ARIADNE_TEST_CALL(test_change_concurrency_and_log_scheduler())
    }
};

int main() {
    TestWorkloadAdvancement().test();
    return ARIADNE_TEST_FAILURES;
}
