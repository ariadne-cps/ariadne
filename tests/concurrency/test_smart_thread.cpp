/***************************************************************************
 *            test_smart_thread.cpp
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

#include "concurrency/smart_thread.hpp"
#include "concurrency/concurrency_manager.hpp"
#include "../test.hpp"

using namespace Ariadne;

class TestSmartThread {
  public:

    void test_create_simplest() const {
        SmartThread thread("thr",[](){});
        ARIADNE_TEST_EXECUTE(thread.id());
        ARIADNE_TEST_ASSERT(thread.has_started());
        ARIADNE_TEST_EQUALS(thread.name(),"thr");
    }

    void test_start_deferred() const {
        SmartThread thread("thr",[](){},false);
        thread.start();
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        ARIADNE_TEST_ASSERT(thread.has_started());
    }

    void test_task() const {
        int a = 0;
        SmartThread thread("thr",[&a](){ a++; });
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        ARIADNE_TEST_EQUALS(a,1);
    }

    void test_entry_exit() const {
        int a = 0;
        {
            SmartThread thread("thr",[&a](){ a*=2; },[&a](){ a++; },[&a](){ a++; });
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            ARIADNE_TEST_EQUALS(a,2);
        }
        ARIADNE_TEST_EQUALS(a,3);
    }

    void test_atomic_multiple_threads() const {
        SizeType n_threads = 10*ConcurrencyManager::instance().maximum_concurrency();
        ARIADNE_TEST_PRINT(n_threads);
        List<SharedPointer<SmartThread>> threads;

        std::atomic<SizeType> a = 0;
        for (SizeType i=0; i<n_threads; ++i)
            threads.append(SharedPointer<SmartThread>(new SmartThread("++"+to_string(i),[&a](){ a++; })));

        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        ARIADNE_TEST_EQUALS(a,n_threads);
    }

    void test() {
        ARIADNE_TEST_CALL(test_create_simplest());
        ARIADNE_TEST_CALL(test_start_deferred());
        ARIADNE_TEST_CALL(test_task());
        ARIADNE_TEST_CALL(test_entry_exit());
        ARIADNE_TEST_CALL(test_atomic_multiple_threads());
    }
};

int main() {
    TestSmartThread().test();
    return ARIADNE_TEST_FAILURES;
}
