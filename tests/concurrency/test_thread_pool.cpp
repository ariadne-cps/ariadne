/***************************************************************************
 *            test_thread_pool.cpp
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

#include "concurrency/thread_pool.hpp"
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

class TestSmartThreadPool {
  public:

    void test_construct_thread_name() const {
        ARIADNE_TEST_EQUALS(construct_thread_name("name",9,9),"name9");
        ARIADNE_TEST_EQUALS(construct_thread_name("name",9,10),"name09");
        ARIADNE_TEST_EQUALS(construct_thread_name("name",10,11),"name10");
    }

    void test_construct() {
        auto max_concurrency = std::thread::hardware_concurrency();
        ThreadPool pool(max_concurrency);
        ARIADNE_TEST_EQUALS(pool.num_threads(),max_concurrency);
        ARIADNE_TEST_EQUALS(pool.queue_size(),0);
    }

    void test_construct_empty() {
        ThreadPool pool(0);
        ARIADNE_TEST_EQUALS(pool.num_threads(),0);
        VoidFunction fn([]{ std::this_thread::sleep_for(100ms); });
        pool.enqueue(fn);
        ARIADNE_TEST_EQUALS(pool.queue_size(),1);
    }

    void test_construct_with_name() {
        ThreadPool pool(1);
        ARIADNE_TEST_EQUALS(pool.name(),THREAD_POOL_DEFAULT_NAME);
        ThreadPool pool2(1,"name");
        ARIADNE_TEST_EQUALS(pool2.name(),"name");
    }

    void test_execute_single() {
        ThreadPool pool(1);
        ARIADNE_TEST_EQUALS(pool.num_threads(),1);
        VoidFunction fn([]{ std::this_thread::sleep_for(100ms); });
        pool.enqueue(fn);
        std::this_thread::sleep_for(200ms);
        ARIADNE_TEST_EQUALS(pool.queue_size(),0);
    }

    void test_exception() {
        ThreadPool pool(1);
        auto future = pool.enqueue([]{ throw new std::exception(); });
        ARIADNE_TEST_FAIL(future.get());
    }

    void test_destroy_before_completion() {
        ThreadPool pool(1);
        pool.enqueue([]{ std::this_thread::sleep_for(100ms); });
    }

    void test_execute_multiple_sequentially() {
        ThreadPool pool(1);
        ARIADNE_TEST_EQUALS(pool.num_threads(),1);
        ARIADNE_TEST_EQUALS(pool.queue_size(),0);
        VoidFunction fn([]{ std::this_thread::sleep_for(100ms); });
        for (SizeType i=0; i<2; ++i) pool.enqueue(fn);
        ARIADNE_TEST_ASSERT(pool.queue_size() > 0);
        std::this_thread::sleep_for(400ms);
        ARIADNE_TEST_EQUALS(pool.queue_size(),0);
    }

    void test_execute_multiple_concurrently() {
        SizeType num_threads = 2;
        ThreadPool pool(num_threads);
        ARIADNE_TEST_EQUALS(pool.num_threads(),2);
        VoidFunction fn([]{ std::this_thread::sleep_for(100ms); });
        for (SizeType i=0; i<2; ++i) pool.enqueue(fn);
        std::this_thread::sleep_for(std::chrono::milliseconds(400*num_threads));
    }

    void test_execute_multiple_concurrently_sequentially() {
        SizeType num_threads = 2;
        ThreadPool pool(num_threads);
        VoidFunction fn([]{ std::this_thread::sleep_for(100ms); });
        for (SizeType i=0; i<2*num_threads; ++i) pool.enqueue(fn);
        ARIADNE_TEST_ASSERT(pool.queue_size() > 0);
        std::this_thread::sleep_for(std::chrono::milliseconds(400*num_threads));
        ARIADNE_TEST_EQUALS(pool.queue_size(),0);
    }

    void test_process_on_atomic_type() {
        auto max_concurrency = std::thread::hardware_concurrency();
        ThreadPool pool(max_concurrency);
        std::vector<Future<SizeType>> results;
        std::atomic<SizeType> x;

        for (SizeType i = 0; i < 2 * max_concurrency; ++i) {
            results.emplace_back(pool.enqueue([&x] {
                                     SizeType r = ++x;
                                     return r * r;
                                 })
            );
        }
        std::this_thread::sleep_for(100ms);
        ARIADNE_TEST_EQUALS(x,2*max_concurrency);

        SizeType actual_sum = 0, expected_sum = 0;
        for (SizeType i = 0; i < 2 * max_concurrency; ++i) {
            actual_sum += results[i].get();
            expected_sum += (i+1)*(i+1);
        }
        ARIADNE_TEST_EQUAL(actual_sum,expected_sum);
    }

    void test_set_num_threads_up_statically() const {
        ThreadPool pool(0);
        ARIADNE_TEST_EXECUTE(pool.set_num_threads(1));
        ARIADNE_TEST_EQUALS(pool.num_threads(),1);
        ARIADNE_TEST_EXECUTE(pool.set_num_threads(3));
        ARIADNE_TEST_EQUALS(pool.num_threads(),3);
    }

    void test_set_num_threads_same_statically() const {
        ThreadPool pool(3);
        ARIADNE_TEST_EXECUTE(pool.set_num_threads(3));
        ARIADNE_TEST_EQUALS(pool.num_threads(),3);
    }

    void test_set_num_threads_down_statically() const {
        ThreadPool pool(3);
        ARIADNE_TEST_EXECUTE(pool.set_num_threads(1));
        ARIADNE_TEST_EQUAL(pool.num_threads(),1);
    }

    void test_set_num_threads_up_dynamically() const {
        ThreadPool pool(0);
        VoidFunction fn([] { std::this_thread::sleep_for(100ms); });
        pool.enqueue(fn);
        std::this_thread::sleep_for(100ms);
        ARIADNE_TEST_EQUALS(pool.queue_size(),1);
        ARIADNE_TEST_EXECUTE(pool.set_num_threads(1));
        ARIADNE_TEST_EQUALS(pool.num_threads(),1);
        std::this_thread::sleep_for(100ms);
        ARIADNE_TEST_EQUALS(pool.queue_size(),0);
        pool.enqueue(fn);
        pool.enqueue(fn);
        ARIADNE_TEST_EXECUTE(pool.set_num_threads(3));
        ARIADNE_TEST_EQUALS(pool.num_threads(),3);
    }

    void test_set_num_threads_down_dynamically() const {
        ThreadPool pool(3);
        VoidFunction fn([] { std::this_thread::sleep_for(100ms); });
        for (SizeType i=0; i<5; ++i)
            pool.enqueue(fn);
        ARIADNE_TEST_EXECUTE(pool.set_num_threads(2));
        ARIADNE_TEST_EQUAL(pool.num_threads(),2);
        std::this_thread::sleep_for(200ms);
        ARIADNE_TEST_EQUALS(pool.queue_size(),0);
    }

    void test_set_num_threads_to_zero_dynamically() const {
        ThreadPool pool(3);
        VoidFunction fn([] { std::this_thread::sleep_for(100ms); });
        for (SizeType i=0; i<5; ++i)
            pool.enqueue(fn);
        ARIADNE_TEST_EXECUTE(pool.set_num_threads(0));
        ARIADNE_TEST_EQUAL(pool.num_threads(),0);
        std::this_thread::sleep_for(100ms);
        ARIADNE_TEST_ASSERT(pool.queue_size() > 0);
    }

    void test() {
        ARIADNE_TEST_CALL(test_construct_thread_name());
        ARIADNE_TEST_CALL(test_construct());
        ARIADNE_TEST_CALL(test_construct_empty());
        ARIADNE_TEST_CALL(test_construct_with_name());
        ARIADNE_TEST_CALL(test_execute_single());
        ARIADNE_TEST_CALL(test_exception());
        ARIADNE_TEST_CALL(test_destroy_before_completion());
        ARIADNE_TEST_CALL(test_execute_multiple_sequentially());
        ARIADNE_TEST_CALL(test_execute_multiple_concurrently());
        ARIADNE_TEST_CALL(test_execute_multiple_concurrently_sequentially());
        ARIADNE_TEST_CALL(test_process_on_atomic_type());
        ARIADNE_TEST_CALL(test_set_num_threads_up_statically());
        ARIADNE_TEST_CALL(test_set_num_threads_same_statically());
        ARIADNE_TEST_CALL(test_set_num_threads_down_statically());
        ARIADNE_TEST_CALL(test_set_num_threads_up_dynamically());
        ARIADNE_TEST_CALL(test_set_num_threads_down_dynamically());
        ARIADNE_TEST_CALL(test_set_num_threads_to_zero_dynamically());
    }
};

Int main() {
    ThreadRegistry registry;
    ConcLog::Logger::instance().attach_thread_registry(&registry);
    TestSmartThreadPool().test();
    return ARIADNE_TEST_FAILURES;
}
