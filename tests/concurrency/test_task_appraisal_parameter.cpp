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
#include "utility/array.hpp"
#include "../test.hpp"

using namespace Ariadne;

class TestTaskAppraisalParameter {
  private:
    DurationType _duration;
  public:

    TestTaskAppraisalParameter() {
        auto now = std::chrono::high_resolution_clock::now();
        _duration = std::chrono::duration_cast<DurationType>(now-now);
    }

    Void test_scalar_appraisal_creation() {
        ScalarAppraisalParameter<int,int> p("chosen_step_size",TaskAppraisalParameterOptimisation::MAXIMISE,
        [](int const& input, int const& output, DurationType const& duration) { return output + duration.count() + input; });
        auto cost = p.appraise(2,5,_duration);
        ARIADNE_TEST_PRINT(p);
        ARIADNE_TEST_ASSERT(p.is_scalar());
        ARIADNE_TEST_EQUALS(cost,7);
        ARIADNE_TEST_EQUALS(p.dimension(2),1);
        ARIADNE_TEST_EQUALS(p.optimisation(),TaskAppraisalParameterOptimisation::MAXIMISE);
    }

    Void test_vector_appraisal_creation() {
        VectorAppraisalParameter<Array<int>,int> p("enclosure_widths",TaskAppraisalParameterOptimisation::MINIMISE,
                                            [](Array<int> const& input, int const& output, DurationType const& duration,SizeType const& idx) {
                                                return output + duration.count() + input[idx]; },
                                            [](Array<int> const& input) { return input.size(); });
        Array<int> input = {1,-3};
        int output = 5;

        ARIADNE_TEST_PRINT(p);
        ARIADNE_TEST_ASSERT(not p.is_scalar());
        ARIADNE_TEST_EQUALS(p.appraise(input,output,_duration,0),6);
        ARIADNE_TEST_EQUALS(p.appraise(input,output,_duration,1),2);
        ARIADNE_TEST_EQUALS(p.dimension(input),2);
        ARIADNE_TEST_EQUALS(p.optimisation(),TaskAppraisalParameterOptimisation::MINIMISE);
    }

    Void test_task_appraisal_weight() {
        typedef ScalarAppraisalParameter<int,int> AppraisalParameter;
        AppraisalParameter p1("chosen_step_size",TaskAppraisalParameterOptimisation::MAXIMISE,
                                            [](int const& input, int const& output, DurationType const& duration) { return output + duration.count() + input; },1.5);
        ARIADNE_TEST_PRINT(p1);
        ARIADNE_TEST_EQUALS(p1.weight(),1.5);


        ARIADNE_TEST_FAIL(AppraisalParameter("chosen_step_size",TaskAppraisalParameterOptimisation::MAXIMISE,
                                             [](int const& input, int const& output, DurationType const& duration) { return output + duration.count() + input; },-1.0));
    }

    Void test_task_appraisal_set() {
        TaskAppraisalParameter<Array<int>,int> p1 = ScalarAppraisalParameter<Array<int>,int>("chosen_step_size",TaskAppraisalParameterOptimisation::MAXIMISE,
            [](Array<int> const& input, int const& output, DurationType const& duration) { return output + duration.count() + input[0]; });

        TaskAppraisalParameter<Array<int>,int> p2 = VectorAppraisalParameter<Array<int>,int>("enclosure_widths",TaskAppraisalParameterOptimisation::MINIMISE,
            [](Array<int> const& input, int const& output, DurationType const& duration,SizeType const& idx) { return output + duration.count() + input[idx]; },
            [](Array<int> const& input) { return input.size(); });

        Set<TaskAppraisalParameter<Array<int>,int>> ps = {p1,p2};

        ARIADNE_TEST_PRINT(ps);
    }

    Void test() {
        ARIADNE_TEST_CALL(test_scalar_appraisal_creation());
        ARIADNE_TEST_CALL(test_vector_appraisal_creation());
        ARIADNE_TEST_CALL(test_task_appraisal_weight());
        ARIADNE_TEST_CALL(test_task_appraisal_set());
    }
};

int main() {
    TestTaskAppraisalParameter().test();
    return ARIADNE_TEST_FAILURES;
}
