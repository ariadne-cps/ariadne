/***************************************************************************
 *            test_workload_advancement.cpp
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

#include "concurrency/workload_advancement.hpp"
#include "../test.hpp"

using namespace Ariadne;

class TestWorkloadAdvancement {
  public:

    void test_creation() {
        WorkloadAdvancement wp(5);
        ARIADNE_TEST_EQUALS(wp.completion_rate(),0.0);
        ARIADNE_TEST_EQUALS(wp.waiting(),5);
        ARIADNE_TEST_EQUALS(wp.processing(),0);
        ARIADNE_TEST_EQUALS(wp.completed(),0);
        ARIADNE_TEST_ASSERT(not wp.has_finished());
    }

    void test_advance() {
        WorkloadAdvancement wp(3);
        wp.add_to_waiting();
        ARIADNE_TEST_EQUALS(wp.waiting(),4);
        wp.add_to_processing();
        ARIADNE_TEST_EQUALS(wp.waiting(),3);
        ARIADNE_TEST_EQUALS(wp.processing(),1);
        wp.add_to_completed();
        ARIADNE_TEST_EQUALS(wp.processing(),0);
        ARIADNE_TEST_EQUALS(wp.completed(),1);
        ARIADNE_TEST_EQUALS(wp.completion_rate(),0.25);
    }

    void test_finished() {
        WorkloadAdvancement wp;
        wp.add_to_waiting(2);
        wp.add_to_processing(2);
        wp.add_to_completed(2);
        ARIADNE_TEST_EQUALS(wp.completion_rate(),1.0);
        ARIADNE_TEST_ASSERT(wp.has_finished());
    }

    void test_invalid_transitions() {
        WorkloadAdvancement wp(4);
        ARIADNE_TEST_FAIL(wp.add_to_processing(5));
        ARIADNE_TEST_FAIL(wp.add_to_completed());
        wp.add_to_processing(2);
        ARIADNE_TEST_FAIL(wp.add_to_completed(3));
        wp.add_to_completed(1);
        ARIADNE_TEST_EQUALS(wp.completion_rate(),0.25);
    }

    void test() {
        ARIADNE_TEST_CALL(test_creation());
        ARIADNE_TEST_CALL(test_advance());
        ARIADNE_TEST_CALL(test_finished());
        ARIADNE_TEST_CALL(test_invalid_transitions());
    }

};

int main() {
    TestWorkloadAdvancement().test();
    return ARIADNE_TEST_FAILURES;
}
