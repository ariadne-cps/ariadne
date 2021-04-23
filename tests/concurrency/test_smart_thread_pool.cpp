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

#include "concurrency/smart_thread_pool.hpp"
#include "../test.hpp"

using namespace Ariadne;

class TestSmartThreadPool {
  public:

    void test_push() {

        auto max_concurrency = std::thread::hardware_concurrency();

        SmartThreadPool pool(max_concurrency);
        std::vector<std::future<SizeType> > results;
        std::atomic<SizeType> x;

        for (SizeType i = 0; i < 2 * max_concurrency; ++i) {
            results.emplace_back(
                    pool.execute([&x] {
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

    void test() {
        ARIADNE_TEST_CALL(test_push());
    }
};

int main() {
    TestSmartThreadPool().test();
    return ARIADNE_TEST_FAILURES;
}
