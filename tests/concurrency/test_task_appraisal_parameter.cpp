/***************************************************************************
 *            test_task_appraisal_parameter.cpp
 *
 *  Copyright  2008-20  Luca Geretti
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

#include "concurrency/task_appraisal_parameter.hpp"
#include "concurrency/task_runner_interface.hpp"
#include "utility/array.hpp"
#include "../test.hpp"

using namespace Ariadne;

class TestRunnable : public TaskRunnable<TestRunnable> { };
typedef TestRunnable R;

namespace Ariadne {
template<> struct TaskInput<R> {
    TaskInput(int i1_, Array<int> i2_) : i1(i1_), i2(i2_) { }
    int i1;
    Array<int> i2;
};
template<> struct TaskOutput<R> {
    TaskOutput(int o_) : o(o_) { }
    int o;
};
}

typedef TaskInput<R> I;
typedef TaskOutput<R> O;

class TestTaskAppraisalParameter {
  private:
    DurationType _duration;
  public:

    TestTaskAppraisalParameter() {
        auto now = std::chrono::high_resolution_clock::now();
        _duration = std::chrono::duration_cast<DurationType>(now-now);
    }

    Void test_scalar_appraisal_creation() {
        ScalarAppraisalParameter<R> p("chosen_step_size",TaskAppraisalParameterOptimisation::MAXIMISE,
        [](I const& input, O const& output, DurationType const& duration) { return output.o + duration.count() + input.i1; });
        auto input = I(2,{1,2});
        auto output = O(7);
        auto cost = p.appraise(input,output,_duration);
        ARIADNE_TEST_PRINT(p);
        ARIADNE_TEST_ASSERT(p.is_scalar());
        ARIADNE_TEST_EQUALS(cost,9);
        ARIADNE_TEST_EQUALS(p.dimension(input),1);
        ARIADNE_TEST_EQUALS(p.optimisation(),TaskAppraisalParameterOptimisation::MAXIMISE);
    }

    Void test_vector_appraisal_creation() {
        VectorAppraisalParameter<R> p("enclosure_widths",TaskAppraisalParameterOptimisation::MINIMISE,
                                            [](I const& input, O const& output, DurationType const& duration, SizeType const& idx) {
                                                return output.o + duration.count() + input.i2[idx]; },
                                            [](I const& input) { return input.i2.size(); });
        auto input = I(2,{1,2});
        auto output = O(7);

        ARIADNE_TEST_PRINT(p);
        ARIADNE_TEST_ASSERT(not p.is_scalar());
        ARIADNE_TEST_EQUALS(p.appraise(input,output,_duration,0),8);
        ARIADNE_TEST_EQUALS(p.appraise(input,output,_duration,1),9);
        ARIADNE_TEST_EQUALS(p.dimension(input),2);
        ARIADNE_TEST_EQUALS(p.optimisation(),TaskAppraisalParameterOptimisation::MINIMISE);
    }

    Void test_task_appraisal_set() {
        ScalarAppraisalParameter<R> p1("chosen_step_size",TaskAppraisalParameterOptimisation::MAXIMISE,
                                      [](I const& input, O const& output, DurationType const& duration) { return output.o + duration.count() + input.i1; });
        VectorAppraisalParameter<R> p2("enclosure_widths",TaskAppraisalParameterOptimisation::MINIMISE,
                                      [](I const& input, O const& output, DurationType const& duration, SizeType const& idx) {
                                          return output.o + duration.count() + input.i2[idx]; },
                                      [](I const& input) { return input.i2.size(); });

        Set<TaskAppraisalParameter<R>> ps = {p1,p2};

        ARIADNE_TEST_PRINT(ps);
    }

    Void test() {
        ARIADNE_TEST_CALL(test_scalar_appraisal_creation());
        ARIADNE_TEST_CALL(test_vector_appraisal_creation());
        ARIADNE_TEST_CALL(test_task_appraisal_set());
    }
};

int main() {
    TestTaskAppraisalParameter().test();
    return ARIADNE_TEST_FAILURES;
}