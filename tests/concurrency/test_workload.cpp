/***************************************************************************
 *            test_workload.cpp
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

#include "concurrency/workload.hpp"
#include "../test.hpp"

using namespace Ariadne;

template<class T> class SynchronisedList : public List<T> {
  public:
    Void append(T const& v) { LockGuard<Mutex> guard(_mux); return List<T>::append(v); }
  private:
    Mutex _mux;
};

using StaticWorkloadType = StaticWorkload<int,std::shared_ptr<std::atomic<int>>>;
using DynamicWorkloadType = DynamicWorkload<int,std::shared_ptr<SynchronisedList<int>>>;

void sum_all(int const& val, std::shared_ptr<std::atomic<int>> result) {
    result->operator+=(val);
}

void print(int const& val) {
    CONCLOG_PRINTLN_VAR(val)
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
}

Void square_and_store(DynamicWorkloadType::Access& wla, int const& val, std::shared_ptr<SynchronisedList<int>> results) {
    int next_val = val*val;
    if (next_val < 46340) {
        wla.append(next_val);
    }
    results->append(next_val);
}

Void progress_acknowledge(int const& val, std::shared_ptr<ProgressIndicator> indicator) {
    indicator->update_current(val);
    indicator->update_final(std::numeric_limits<int>::max());
}

Void throw_exception_immediately(DynamicWorkloadType::Access& wla, int const& val, std::shared_ptr<SynchronisedList<int>> results) {
    throw new std::exception();
}

Void throw_exception_later(DynamicWorkloadType::Access& wla, int const& val, std::shared_ptr<SynchronisedList<int>> results) {
    int next_val = val+1;
    if (next_val > 4) throw new std::exception();
    else wla.append(next_val);
}

class TestWorkload {
  public:

    void test_construct_static() {
        TaskManager::instance().set_concurrency(0);
        auto result = std::make_shared<std::atomic<int>>();
        StaticWorkloadType wl(&sum_all, result);
    }

    void test_construct_dynamic() {
        TaskManager::instance().set_concurrency(0);
        std::shared_ptr<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        DynamicWorkloadType wl(&progress_acknowledge, &square_and_store, result);
    }

    void test_append() {
        TaskManager::instance().set_concurrency(0);
        std::shared_ptr<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        DynamicWorkloadType wl(&progress_acknowledge, &square_and_store, result);
        wl.append(2);
        ARIADNE_TEST_EQUALS(wl.size(),1)
        wl.append({10,20});
        ARIADNE_TEST_EQUALS(wl.size(),3)
    }

    void test_process_nothing() {
        TaskManager::instance().set_maximum_concurrency();
        auto result = std::make_shared<std::atomic<int>>();
        StaticWorkloadType wl(&sum_all, result);
        ARIADNE_TEST_EXECUTE(wl.process())
    }

    void test_serial_processing_static() {
        TaskManager::instance().set_concurrency(0);
        std::shared_ptr<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        result->append(2);
        DynamicWorkloadType wl(&progress_acknowledge, &square_and_store, result);
        wl.append(2);
        wl.process();
        ARIADNE_TEST_PRINT(*result)
        ARIADNE_TEST_EQUALS(result->size(),5)
    }

    void test_serial_processing_dynamic() {
        TaskManager::instance().set_concurrency(0);
        std::shared_ptr<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        result->append(2);
        DynamicWorkloadType wl(&progress_acknowledge, &square_and_store, result);
        wl.append(2);
        wl.process();
        ARIADNE_TEST_PRINT(*result)
        ARIADNE_TEST_EQUALS(result->size(),5)
    }

    void test_concurrent_processing_static() {
        TaskManager::instance().set_maximum_concurrency();
        auto result = std::make_shared<std::atomic<int>>();
        *result = 0;
        StaticWorkloadType wl(&sum_all, result);
        wl.append({2,7,-3,5,8,10,5,8});
        wl.process();
        ARIADNE_TEST_EQUALS(*result,42)
    }

    void test_concurrent_processing_dynamic() {
        TaskManager::instance().set_maximum_concurrency();
        std::shared_ptr<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        result->append(2);
        DynamicWorkloadType wl(&progress_acknowledge, &square_and_store, result);
        wl.append(2);
        wl.process();
        ARIADNE_TEST_PRINT(*result)
        ARIADNE_TEST_EQUALS(result->size(),5)
    }

    void test_print_hold() {
        TaskManager::instance().set_concurrency(0);
        Logger::instance().configuration().set_verbosity(2);
        StaticWorkload<int> wl(&print);
        wl.append({1,2,3,4,5});
        wl.process();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        Logger::instance().configuration().set_verbosity(0);
    }

    void test_throw_serial_exception_immediately() {
        TaskManager::instance().set_concurrency(0);
        std::shared_ptr<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        DynamicWorkloadType wl(&progress_acknowledge, &throw_exception_immediately, result);
        wl.append(2);
        ARIADNE_TEST_FAIL(wl.process())
    }

    void test_throw_serial_exception_later() {
        TaskManager::instance().set_concurrency(0);
        std::shared_ptr<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        DynamicWorkloadType wl(&progress_acknowledge, &throw_exception_later, result);
        wl.append(2);
        ARIADNE_TEST_FAIL(wl.process())
    }

    void test_throw_concurrent_exception_immediately() {
        TaskManager::instance().set_maximum_concurrency();
        std::shared_ptr<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        DynamicWorkloadType wl(&progress_acknowledge, &throw_exception_immediately, result);
        wl.append(2);
        ARIADNE_TEST_FAIL(wl.process())
    }

    void test_throw_concurrent_exception_later() {
        TaskManager::instance().set_maximum_concurrency();
        std::shared_ptr<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        DynamicWorkloadType wl(&progress_acknowledge, &throw_exception_later, result);
        wl.append(2);
        ARIADNE_TEST_FAIL(wl.process())
    }

    void test_multiple_append() {
        TaskManager::instance().set_maximum_concurrency();
        std::shared_ptr<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        DynamicWorkloadType wl(&progress_acknowledge, &square_and_store, result);
        result->append(2);
        result->append(3);
        wl.append({2,3});
        wl.process();
        ARIADNE_TEST_PRINT(*result)
        ARIADNE_TEST_EQUALS(result->size(),10)
    }

    void test_multiple_process() {
        TaskManager::instance().set_maximum_concurrency();
        std::shared_ptr<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        result->append(2);
        DynamicWorkloadType wl(&progress_acknowledge, &square_and_store, result);
        wl.append(2);
        wl.process();
        result->clear();
        result->append(3);
        wl.append(3);
        wl.process();
        ARIADNE_TEST_PRINT(*result)
        ARIADNE_TEST_EQUALS(result->size(),5)
    }

    void test() {
        ARIADNE_TEST_CALL(test_construct_static())
        ARIADNE_TEST_CALL(test_construct_dynamic())
        ARIADNE_TEST_CALL(test_append())
        ARIADNE_TEST_CALL(test_process_nothing())
        ARIADNE_TEST_CALL(test_serial_processing_static())
        ARIADNE_TEST_CALL(test_serial_processing_dynamic())
        ARIADNE_TEST_CALL(test_concurrent_processing_static())
        ARIADNE_TEST_CALL(test_concurrent_processing_dynamic())
        ARIADNE_TEST_CALL(test_print_hold())
        ARIADNE_TEST_CALL(test_throw_serial_exception_immediately())
        ARIADNE_TEST_CALL(test_throw_serial_exception_later())
        ARIADNE_TEST_CALL(test_throw_concurrent_exception_immediately())
        ARIADNE_TEST_CALL(test_throw_concurrent_exception_later())
        ARIADNE_TEST_CALL(test_multiple_append())
        ARIADNE_TEST_CALL(test_multiple_process())
    }

};

int main() {
    TestWorkload().test();
    return ARIADNE_TEST_FAILURES;
}
