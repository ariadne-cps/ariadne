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

using WorkloadType = Workload<int,SharedPointer<List<int>>>;
using WorkloadAccessType = WorkloadType::StackAccess;

Void square_and_store(WorkloadAccessType& wl, int val, SharedPointer<List<int>> results) {
    val *= val;
    results->push_back(val);
    if (val < 46340) {
        wl.push(val);
    }
}

struct WorkloadUser {

    List<int> process(int initial) {
        SharedPointer<List<int>> result = std::make_shared<List<int>>();
        WorkloadType wl([=,this](WorkloadAccessType& wl, int val, SharedPointer<List<int>> results){ _square_and_store(wl, val, results); },result);
        result->push_back(initial);
        wl.add(initial);
        wl.process();
        return *result;
    }

  private:

    Void _square_and_store(WorkloadAccessType& wl, int val, SharedPointer<List<int>> results) {
        square_and_store(wl,val,results);
    }
};

class TestWorkload {
  public:

    void test_construction_from_function() {
        SharedPointer<List<int>> result = std::make_shared<List<int>>();
        WorkloadType wl(&square_and_store,result);
    }

    void test_construction_from_method() {
        ARIADNE_TEST_EXECUTE(WorkloadUser());
    }

    void test_method_processing() {
        WorkloadUser wu;
        auto result = wu.process(2);
        ARIADNE_TEST_EQUALS(result.size(),5);
    }

    void test() {
        ARIADNE_TEST_CALL(test_construction_from_function());
        ARIADNE_TEST_CALL(test_construction_from_method());
        ARIADNE_TEST_CALL(test_method_processing());
    }

};

int main() {
    TestWorkload().test();
    return ARIADNE_TEST_FAILURES;
}
