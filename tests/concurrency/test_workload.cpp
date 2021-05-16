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

using StaticWorkloadType = StaticWorkload<int,SharedPointer<std::atomic<int>>>;
using DynamicWorkloadType = DynamicWorkload<int,SharedPointer<SynchronisedList<int>>>;

void sum_all(int const& val, SharedPointer<std::atomic<int>> result) {
    result->operator+=(val);
}

Void square_and_store(DynamicWorkloadType::Access& wla, int const& val, SharedPointer<SynchronisedList<int>> results) {
    int next_val = val*val;
    if (next_val < 46340) {
        wla.append(next_val);
    }
    results->append(next_val);
}

Void check_processing(DynamicWorkloadType::Access& wla, int const& val, SharedPointer<SynchronisedList<int>> results) {
    ARIADNE_TEST_EQUALS(wla.advancement().processing(),1)
}

Void throw_exception_immediately(DynamicWorkloadType::Access& wla, int const& val, SharedPointer<SynchronisedList<int>> results) {
    throw new std::exception();
}

Void throw_exception_later(DynamicWorkloadType::Access& wla, int const& val, SharedPointer<SynchronisedList<int>> results) {
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
        SharedPointer<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        DynamicWorkloadType wl(&square_and_store, result);
    }

    void test_invalid_process() {
        TaskManager::instance().set_concurrency(0);
        SharedPointer<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        DynamicWorkloadType wl(&square_and_store, result);
        ARIADNE_TEST_FAIL(wl.process())
    }

    void test_append() {
        TaskManager::instance().set_concurrency(0);
        SharedPointer<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        DynamicWorkloadType wl(&square_and_store, result);
        wl.append(2);
        ARIADNE_TEST_EQUALS(wl.size(),1)
        wl.append({10,20});
        ARIADNE_TEST_EQUALS(wl.size(),3)
    }

    void test_serial_processing() {
        TaskManager::instance().set_concurrency(0);
        SharedPointer<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        result->append(2);
        DynamicWorkloadType wl(&square_and_store, result);
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
        TaskManager::instance().set_concurrency(1);
        SharedPointer<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        result->append(2);
        DynamicWorkloadType wl(&square_and_store, result);
        wl.append(2);
        wl.process();
        ARIADNE_TEST_PRINT(*result)
        ARIADNE_TEST_EQUALS(result->size(),5)
    }

    void test_throw_serial_exception_immediately() {
        TaskManager::instance().set_concurrency(0);
        SharedPointer<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        DynamicWorkloadType wl(&throw_exception_immediately, result);
        wl.append(2);
        ARIADNE_TEST_FAIL(wl.process())
    }

    void test_throw_serial_exception_later() {
        TaskManager::instance().set_concurrency(0);
        SharedPointer<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        DynamicWorkloadType wl(&throw_exception_later, result);
        wl.append(2);
        ARIADNE_TEST_FAIL(wl.process())
    }

    void test_throw_concurrent_exception_immediately() {
        TaskManager::instance().set_concurrency(1);
        SharedPointer<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        DynamicWorkloadType wl(&throw_exception_immediately, result);
        wl.append(2);
        ARIADNE_TEST_FAIL(wl.process())
    }

    void test_throw_concurrent_exception_later() {
        TaskManager::instance().set_concurrency(1);
        SharedPointer<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        DynamicWorkloadType wl(&throw_exception_later, result);
        wl.append(2);
        ARIADNE_TEST_FAIL(wl.process())
    }

    void test_multiple_append() {
        TaskManager::instance().set_concurrency(2);
        SharedPointer<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        DynamicWorkloadType wl(&square_and_store, result);
        result->append(2);
        result->append(3);
        wl.append({2,3});
        wl.process();
        ARIADNE_TEST_PRINT(*result)
        ARIADNE_TEST_EQUALS(result->size(),10)
    }

    void test_multiple_process() {
        TaskManager::instance().set_concurrency(2);
        SharedPointer<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        result->append(2);
        DynamicWorkloadType wl(&square_and_store, result);
        wl.append(2);
        wl.process();
        result->clear();
        result->append(3);
        wl.append(3);
        wl.process();
        ARIADNE_TEST_PRINT(*result)
        ARIADNE_TEST_EQUALS(result->size(),5)
    }

    void test_advancement() {
        TaskManager::instance().set_concurrency(1);
        SharedPointer<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        DynamicWorkloadType wl(&check_processing, result);
        wl.append(0);
        wl.process();
    }

    void test() {
        ARIADNE_TEST_CALL(test_construct_static())
        ARIADNE_TEST_CALL(test_construct_dynamic())
        ARIADNE_TEST_CALL(test_append())
        ARIADNE_TEST_CALL(test_invalid_process())
        ARIADNE_TEST_CALL(test_serial_processing())
        ARIADNE_TEST_CALL(test_concurrent_processing_static())
        ARIADNE_TEST_CALL(test_concurrent_processing_dynamic())
        ARIADNE_TEST_CALL(test_throw_serial_exception_immediately())
        ARIADNE_TEST_CALL(test_throw_serial_exception_later())
        ARIADNE_TEST_CALL(test_throw_concurrent_exception_immediately())
        ARIADNE_TEST_CALL(test_throw_concurrent_exception_later())
        ARIADNE_TEST_CALL(test_multiple_append())
        ARIADNE_TEST_CALL(test_multiple_process())
        ARIADNE_TEST_CALL(test_advancement())
    }

};

int main() {
    TestWorkload().test();
    return ARIADNE_TEST_FAILURES;
}
