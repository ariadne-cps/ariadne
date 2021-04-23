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

    void test_create() const {
        SmartThread thread1("thr");
        ARIADNE_TEST_EXECUTE(thread1.id());
        ARIADNE_TEST_EQUALS(thread1.name(),"thr");
        SmartThread thread2;
        ARIADNE_TEST_EQUALS(to_string(thread2.id()),thread2.name());
    }

    void test_task_return() const {
        SmartThread thread;
        auto result = thread.execute([] { return 42; });
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        ARIADNE_TEST_EQUALS(result.get(),42);
    }

    void test_task_capture() const {
        int a = 0;
        SmartThread thread;
        thread.execute([&a] { a++; });
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        ARIADNE_TEST_EQUALS(a,1);
    }

    void test_task_arguments() const {
        int x = 3;
        int y = 5;
        SmartThread thread;
        auto future = thread.execute([](int x, int y) { return x * y; }, x, y);
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        auto r = future.get();
        ARIADNE_TEST_EQUALS(r,15);
    }

    void test_multiple_tasks() const {
        SmartThread thread;
        int a = 4;
        thread.execute([&a]{ a+=2; return a; });
        auto future = thread.execute([&a]{ a*=7; return a; });
        int r = future.get();
        ARIADNE_TEST_EQUALS(r,42);
    }

    void test_atomic_multiple_threads() const {
        SizeType n_threads = 10*ConcurrencyManager::instance().maximum_concurrency();
        ARIADNE_TEST_PRINT(n_threads);
        List<SharedPointer<SmartThread>> threads;

        std::atomic<SizeType> a = 0;
        for (SizeType i=0; i<n_threads; ++i) {
            threads.append(SharedPointer<SmartThread>(new SmartThread("add" + to_string(i))));
            threads.at(i)->execute([&a] { a++; });
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        ARIADNE_TEST_EQUALS(a,n_threads);
    }

    void test() {
        ARIADNE_TEST_CALL(test_create());
        ARIADNE_TEST_CALL(test_task_return());
        ARIADNE_TEST_CALL(test_task_capture());
        ARIADNE_TEST_CALL(test_task_arguments());
        ARIADNE_TEST_CALL(test_multiple_tasks());
        ARIADNE_TEST_CALL(test_atomic_multiple_threads());
    }

};

int main() {
    TestSmartThread().test();
    return ARIADNE_TEST_FAILURES;
}
