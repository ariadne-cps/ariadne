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
#include "../test.hpp"

using namespace Ariadne;

class TestSmartThreadPool {
  public:

    void test_construct() {
        auto max_concurrency = std::thread::hardware_concurrency();
        ThreadPool pool(max_concurrency);
        ARIADNE_TEST_EQUALS(pool.num_threads(),max_concurrency);
        ARIADNE_TEST_EQUALS(pool.queue_size(),0);
    }

    void test_construct_empty() {
        ThreadPool pool(0);
        ARIADNE_TEST_EQUALS(pool.num_threads(),0);
        VoidFunction fn([]{ std::this_thread::sleep_for(std::chrono::milliseconds(100)); });
        pool.enqueue(fn);
        ARIADNE_TEST_EQUALS(pool.queue_size(),1);
    }

    void test_execute_single() {
        ThreadPool pool(1);
        ARIADNE_TEST_EQUALS(pool.num_threads(),1);
        VoidFunction fn([]{ std::this_thread::sleep_for(std::chrono::milliseconds(100)); });
        pool.enqueue(fn);
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        ARIADNE_TEST_EQUALS(pool.queue_size(),0);
    }

    void test_destroy_before_completion() {
        ThreadPool pool(1);
        VoidFunction fn([]{ std::this_thread::sleep_for(std::chrono::milliseconds(100)); });
        pool.enqueue(fn);
    }

    void test_execute_multiple_sequentially() {
        ThreadPool pool(1);
        ARIADNE_TEST_EQUALS(pool.num_threads(),1);
        ARIADNE_TEST_EQUALS(pool.queue_size(),0);
        VoidFunction fn([]{ std::this_thread::sleep_for(std::chrono::milliseconds(100)); });
        for (SizeType i=0; i<2; ++i) pool.enqueue(fn);
        ARIADNE_TEST_ASSERT(pool.queue_size() > 0);
        std::this_thread::sleep_for(std::chrono::milliseconds(400));
        ARIADNE_TEST_EQUALS(pool.queue_size(),0);
    }

    void test_execute_multiple_concurrently() {
        SizeType num_threads = 2;
        ThreadPool pool(num_threads);
        ARIADNE_TEST_EQUALS(pool.num_threads(),2);
        VoidFunction fn([]{ std::this_thread::sleep_for(std::chrono::milliseconds(100)); });
        for (SizeType i=0; i<2; ++i) pool.enqueue(fn);
        std::this_thread::sleep_for(std::chrono::milliseconds(400*num_threads));
    }

    void test_execute_multiple_concurrently_sequentially() {
        SizeType num_threads = 2;
        ThreadPool pool(num_threads);
        VoidFunction fn([]{ std::this_thread::sleep_for(std::chrono::milliseconds(100)); });
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
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        ARIADNE_TEST_EQUALS(x,2*max_concurrency);

        SizeType actual_sum = 0, expected_sum = 0;
        for (SizeType i = 0; i < 2 * max_concurrency; ++i) {
            actual_sum += results[i].get();
            expected_sum += (i+1)*(i+1);
        }
        ARIADNE_TEST_EQUAL(actual_sum,expected_sum);
    }

    void test_add_threads() const {
        ThreadPool pool(0);
        VoidFunction fn([] { std::this_thread::sleep_for(std::chrono::milliseconds(100)); });
        pool.enqueue(fn);
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        ARIADNE_TEST_EQUALS(pool.queue_size(),1);
        ARIADNE_TEST_EXECUTE(pool.add_threads(0));
        ARIADNE_TEST_EXECUTE(pool.add_threads(1));
        ARIADNE_TEST_EQUALS(pool.num_threads(),1);
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        ARIADNE_TEST_EQUALS(pool.queue_size(),0);
        pool.enqueue(fn);
        pool.enqueue(fn);
        ARIADNE_TEST_EXECUTE(pool.add_threads(2));
        ARIADNE_TEST_EQUALS(pool.num_threads(),3);
    }

    void test_schedule_stop_threads() const {
        ThreadPool pool(2);
        VoidFunction fn([] { std::this_thread::sleep_for(std::chrono::milliseconds(100)); });
        for (SizeType i=0; i<16; ++i)
            pool.enqueue(fn);
        pool.schedule_stop_threads();
        std::this_thread::sleep_for(std::chrono::milliseconds(300));
        ARIADNE_TEST_ASSERT(pool.num_threads() > 0);
        ARIADNE_TEST_EQUALS(pool.num_stopped_threads(),pool.num_threads());
        ARIADNE_TEST_ASSERT(pool.queue_size() > 0);
        pool.remove_threads();
        ARIADNE_TEST_EQUALS(pool.num_threads(),0);
    }

    void test_improperly_remove_threads() const {
        ThreadPool pool(2);
        ARIADNE_TEST_FAIL(pool.remove_threads());
    }

    void test_resize_pool_down() const {
        ThreadPool pool(3);
        VoidFunction fn([] { std::this_thread::sleep_for(std::chrono::milliseconds(100)); });
        for (SizeType i=0; i<6; ++i)
            pool.enqueue(fn);
        pool.schedule_stop_threads();
        while(true) {
            try {
                pool.remove_threads();
                break;
            } catch(...) {
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
            }
        }
        ARIADNE_TEST_ASSERT(pool.queue_size() > 0);
        pool.add_threads(2);
        std::this_thread::sleep_for(std::chrono::milliseconds(300));
        ARIADNE_TEST_EQUALS(pool.queue_size(),0);
    }

    void test() {
        ARIADNE_TEST_CALL(test_construct());
        ARIADNE_TEST_CALL(test_construct_empty());
        ARIADNE_TEST_CALL(test_execute_single());
        ARIADNE_TEST_CALL(test_destroy_before_completion());
        ARIADNE_TEST_CALL(test_execute_multiple_sequentially());
        ARIADNE_TEST_CALL(test_execute_multiple_concurrently());
        ARIADNE_TEST_CALL(test_execute_multiple_concurrently_sequentially());
        ARIADNE_TEST_CALL(test_process_on_atomic_type());
        ARIADNE_TEST_CALL(test_add_threads());
        ARIADNE_TEST_CALL(test_schedule_stop_threads());
        ARIADNE_TEST_CALL(test_improperly_remove_threads());
        ARIADNE_TEST_CALL(test_resize_pool_down());
    }
};

Int main() {
    TestSmartThreadPool().test();
    return ARIADNE_TEST_FAILURES;
}
