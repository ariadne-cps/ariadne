/***************************************************************************
 *            test_task_parameter.cpp
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

#include "symbolic/expression_set.hpp"
#include "dynamics/task_parameter.hpp"

#include "../test.hpp"

using namespace Ariadne;

class TestInterface { };
class A : public TestInterface { };
class B : public TestInterface { };
class C : public TestInterface { };

class TestTaskParameter {
  public:

    Void test_boolean_task_parameter() {
        TaskParameter p = BooleanTaskParameter("use_subdivisions");
        ARIADNE_TEST_PRINT(p);
    }

    Void test_metric_task_parameter() {
        TaskParameter unbounded = MetricTaskParameter("step_size");
        ARIADNE_TEST_PRINT(unbounded);
        TaskParameter bounded = MetricTaskParameter("sweep_threshold",10);
        ARIADNE_TEST_PRINT(bounded);
    }

    Void test_enumeration_task_parameter() {
        TaskParameter p = EnumerationTaskParameter<TestInterface>("test_class",{A(),B(),C()});
        ARIADNE_TEST_PRINT(p);
        ARIADNE_TEST_ASSERT(p.is_upper_bounded());
        ARIADNE_TEST_ASSERT(not p.is_metric());
        ARIADNE_TEST_EQUAL(p.upper_bound(),2);
    }

    Void test_parameter_space() {
        BooleanTaskParameter boolean_task("boolean_task");
        MetricTaskParameter metric_task("metric_task");
        EnumerationTaskParameter<TestInterface> enum_task("enum_task",{A(),B(),C()});

        TaskParameterSpace space({boolean_task,metric_task,enum_task});
        ARIADNE_TEST_PRINT(space);
        ARIADNE_TEST_EQUAL(space.dimension(),3);
        ARIADNE_TEST_EQUAL(space.index(boolean_task),0);
        ARIADNE_TEST_EQUAL(space.index(enum_task),2);
        ARIADNE_TEST_EQUAL(space.index(metric_task),1);
    }

    Void test_parameter_point_creation() {
        BooleanTaskParameter boolean_task("boolean_task");
        MetricTaskParameter metric_task("metric_task");
        EnumerationTaskParameter<TestInterface> enum_task("enum_task",{A(),B(),C()});

        ParameterBindingsMap bindings = {{boolean_task,1},{metric_task,5},{enum_task,2}};
        ARIADNE_TEST_PRINT(bindings);
        TaskParameterPoint point(bindings);
        ARIADNE_TEST_PRINT(point);
        TaskParameterSpace space = point.space();
        ARIADNE_TEST_EQUAL(space.dimension(),3);
        ARIADNE_TEST_EQUAL(space.index(boolean_task),0);
        ARIADNE_TEST_EQUAL(space.index(enum_task),1);
        ARIADNE_TEST_EQUAL(space.index(metric_task),2);
    }

    Void test() {
        ARIADNE_TEST_CALL(test_boolean_task_parameter());
        ARIADNE_TEST_CALL(test_metric_task_parameter());
        ARIADNE_TEST_CALL(test_enumeration_task_parameter());
        ARIADNE_TEST_CALL(test_parameter_space());
        ARIADNE_TEST_CALL(test_parameter_point_creation());

    }
};

int main() {
    TestTaskParameter().test();
    return ARIADNE_TEST_FAILURES;
}
