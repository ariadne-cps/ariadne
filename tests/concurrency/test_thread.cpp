/***************************************************************************
 *            test_thread.cpp
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

#include "utility/container.hpp"
#include "concurrency/thread.hpp"
#include "conclog/logging.hpp"
#include "conclog/thread_registry_interface.hpp"
#include "../test.hpp"

using namespace Ariadne;

using namespace std::chrono_literals;

class ThreadRegistry : public ConcLog::ThreadRegistryInterface {
public:
    ThreadRegistry() : _threads_registered(0) { }
    bool has_threads_registered() const override { return _threads_registered > 0; }
    void set_threads_registered(unsigned int threads_registered) { _threads_registered = threads_registered; }
private:
    unsigned int _threads_registered;
};

class TestThread {
  public:

    void test_create() const {
        Thread thread1([]{}, "thr");
        ARIADNE_TEST_EXECUTE(thread1.id())
        ARIADNE_TEST_EQUALS(thread1.name(),"thr")
        Thread thread2([]{});
        ARIADNE_TEST_EQUALS(to_string(thread2.id()),thread2.name())
    }

    void test_destroy_before_completion() const {
        Thread thread([] { std::this_thread::sleep_for(100ms); });
    }

    void test_task() const {
        int a = 0;
        Thread thread([&a] { a++; });
        std::this_thread::sleep_for(10ms);
        ARIADNE_TEST_EQUALS(a,1)
        ARIADNE_TEST_ASSERT(thread.exception() == nullptr)
    }

    void test_exception() const {
        Thread thread([] { throw new std::exception(); });
        std::this_thread::sleep_for(10ms);
        ARIADNE_TEST_ASSERT(thread.exception() != nullptr)
    }

    void test_atomic_multiple_threads() const {
        SizeType n_threads = 10*std::thread::hardware_concurrency();
        ARIADNE_TEST_PRINT(n_threads)
        List<SharedPointer<Thread>> threads;

        std::atomic<SizeType> a = 0;
        for (SizeType i=0; i<n_threads; ++i) {
            threads.append(std::make_shared<Thread>([&a] { a++; }));
        }

        std::this_thread::sleep_for(100ms);
        ARIADNE_TEST_EQUALS(a,n_threads)
        threads.clear();
    }

    void test() {
        ARIADNE_TEST_CALL(test_create())
        ARIADNE_TEST_CALL(test_destroy_before_completion())
        ARIADNE_TEST_CALL(test_task())
        ARIADNE_TEST_CALL(test_exception())
        ARIADNE_TEST_CALL(test_atomic_multiple_threads())
    }

};

int main() {
    ThreadRegistry registry;
    ConcLog::Logger::instance().attach_thread_registry(&registry);
    TestThread().test();
    return ARIADNE_TEST_FAILURES;
}
