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
    Void append(T const& v) {
        LockGuard<Mutex> guard(_mux);
        return List<T>::append(v);
    }
    SizeType size() const {
        return List<T>::size();
    }
  private:
    Mutex mutable _mux;
};

using WorkloadType = Workload<int,SharedPointer<SynchronisedList<int>>>;
using WorkloadAccessType = WorkloadType::StackAccess;

Void square_and_store(WorkloadAccessType& wl, int val, SharedPointer<SynchronisedList<int>> results) {
    val *= val;
    if (val < 46340) {
        wl.push(val);
    }
    results->append(val);
}

class TestWorkload {
  public:

    void test_construction() {
        SharedPointer<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        WorkloadType wl(&square_and_store,result);
    }

    void test_serial_processing() {
        TaskManager::instance().set_concurrency(0);
        SharedPointer<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        result->append(2);
        WorkloadType wl(&square_and_store,result);
        wl.add(2);
        wl.process();
        ARIADNE_TEST_PRINT(*result);
        ARIADNE_TEST_EQUALS(result->size(),5);
    }

    void test_concurrent_processing() {
        TaskManager::instance().set_concurrency(1);
        SharedPointer<SynchronisedList<int>> result = std::make_shared<SynchronisedList<int>>();
        result->append(2);
        WorkloadType wl(&square_and_store,result);
        wl.add(2);
        wl.process();
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        ARIADNE_TEST_PRINT(*result);
        ARIADNE_TEST_EQUALS(result->size(),5);
    }

    void test() {
        ARIADNE_TEST_CALL(test_construction());
        ARIADNE_TEST_CALL(test_serial_processing());
        ARIADNE_TEST_CALL(test_concurrent_processing());
    }

};

int main() {
    TestWorkload().test();
    return ARIADNE_TEST_FAILURES;
}
