/***************************************************************************
 *            test_buffered_thread.cpp
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
#include "concurrency/buffered_thread.hpp"
#include "conclog/logging.hpp"
#include "conclog/thread_registry_interface.hpp"
#include "../test.hpp"

using namespace Ariadne;

class ThreadRegistry : public ConcLog::ThreadRegistryInterface {
public:
    ThreadRegistry() : _threads_registered(0) { }
    bool has_threads_registered() const override { return _threads_registered > 0; }
    void set_threads_registered(unsigned int threads_registered) { _threads_registered = threads_registered; }
private:
    unsigned int _threads_registered;
};

class TestBufferedThread {
  public:

    void test_create() const {
        BufferedThread thread1("thr");
        ARIADNE_TEST_EXECUTE(thread1.id());
        ARIADNE_TEST_EQUALS(thread1.name(),"thr");
        ARIADNE_TEST_EQUALS(thread1.queue_size(),0);
        ARIADNE_TEST_EQUALS(thread1.queue_capacity(),1);
        BufferedThread thread2;
        ARIADNE_TEST_EQUALS(to_string(thread2.id()),thread2.name());
    }

    void test_set_queue_capacity() const {
        BufferedThread thread;
        ARIADNE_TEST_FAIL(thread.set_queue_capacity(0));
        ARIADNE_TEST_EXECUTE(thread.set_queue_capacity(2));
        ARIADNE_TEST_EXECUTE(thread.set_queue_capacity(1));
    }

    void test_destroy_before_completion() const {
        BufferedThread thread;
        thread.enqueue([] { std::this_thread::sleep_for(std::chrono::milliseconds(100)); });
    }

    void test_exception() const {
        BufferedThread thread;
        auto future = thread.enqueue([] { throw new std::exception(); });
        ARIADNE_TEST_FAIL(future.get());
    }

    void test_has_queued_tasks() const {
        BufferedThread thread;
        thread.set_queue_capacity(2);
        thread.enqueue([] { std::this_thread::sleep_for(std::chrono::milliseconds(100)); });
        thread.enqueue([] { std::this_thread::sleep_for(std::chrono::milliseconds(100)); });
        ARIADNE_TEST_ASSERT(thread.queue_size()>0);
        std::this_thread::sleep_for(std::chrono::milliseconds(300));
        ARIADNE_TEST_EQUALS(thread.queue_size(),0);
    }

    void test_set_queue_capacity_down_failure() const {
        BufferedThread thread;
        thread.set_queue_capacity(3);
        VoidFunction fn([]{ std::this_thread::sleep_for(std::chrono::milliseconds(100)); });
        thread.enqueue(fn);
        thread.enqueue(fn);
        thread.enqueue(fn);
        ARIADNE_TEST_FAIL(thread.set_queue_capacity(1));
        std::this_thread::sleep_for(std::chrono::milliseconds(400));
        ARIADNE_TEST_EXECUTE(thread.set_queue_capacity(1));
    }

    void test_task_return() const {
        BufferedThread thread;
        auto result = thread.enqueue([] { return 42; });
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        ARIADNE_TEST_EQUALS(result.get(),42);
    }

    void test_task_capture() const {
        int a = 0;
        BufferedThread thread;
        thread.enqueue([&a] { a++; });
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        ARIADNE_TEST_EQUALS(a,1);
    }

    void test_task_arguments() const {
        int x = 3;
        int y = 5;
        BufferedThread thread;
        auto future = thread.enqueue([](int x, int y) { return x * y; }, x, y);
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        auto r = future.get();
        ARIADNE_TEST_EQUALS(r,15);
    }

    void test_multiple_tasks() const {
        BufferedThread thread;
        int a = 4;
        thread.enqueue([&a] {
            a += 2;
            return a;
        });
        auto future = thread.enqueue([&a] {
            a *= 7;
            return a;
        });
        int r = future.get();
        ARIADNE_TEST_EQUALS(r,42);
    }

    void test_atomic_multiple_threads() const {
        SizeType n_threads = 10*std::thread::hardware_concurrency();
        ARIADNE_TEST_PRINT(n_threads);
        List<SharedPointer<BufferedThread>> threads;

        std::atomic<SizeType> a = 0;
        for (SizeType i=0; i<n_threads; ++i) {
            threads.append(make_shared<BufferedThread>(("add" + to_string(i))));
            threads.at(i)->enqueue([&a] { a++; });
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        ARIADNE_TEST_EQUALS(a,n_threads);
        threads.clear();
    }

    void test() {
        ARIADNE_TEST_CALL(test_create());
        ARIADNE_TEST_CALL(test_set_queue_capacity());
        ARIADNE_TEST_CALL(test_destroy_before_completion());
        ARIADNE_TEST_CALL(test_exception());
        ARIADNE_TEST_CALL(test_has_queued_tasks());
        ARIADNE_TEST_CALL(test_set_queue_capacity_down_failure());
        ARIADNE_TEST_CALL(test_task_return());
        ARIADNE_TEST_CALL(test_task_capture());
        ARIADNE_TEST_CALL(test_task_arguments());
        ARIADNE_TEST_CALL(test_multiple_tasks());
        ARIADNE_TEST_CALL(test_atomic_multiple_threads());
    }

};

int main() {
    ThreadRegistry registry;
    ConcLog::Logger::instance().attach_thread_registry(&registry);
    TestBufferedThread().test();
    return ARIADNE_TEST_FAILURES;
}
