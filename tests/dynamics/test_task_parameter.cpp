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
class D : public TestInterface { };

class TestTaskParameter {
  public:

    TestTaskParameter() { srand (time(NULL)); }

    Void test_boolean_task_parameter() {
        TaskParameter p = BooleanTaskParameter("use_subdivisions");
        ARIADNE_TEST_PRINT(p);
    }

    Void test_boolean_task_parameter_shift() {
        BooleanTaskParameter b("boolean_task");
        ARIADNE_TEST_EQUALS(b.shifted_value_from(0),1);
        ARIADNE_TEST_EQUALS(b.shifted_value_from(1),0);
    }

    Void test_metric_task_parameter() {
        TaskParameter metric = MetricTaskParameter("sweep_threshold",10);
        ARIADNE_TEST_PRINT(metric);
    }

    Void test_metric_task_parameter_shift() {
        MetricTaskParameter metric("sweep_threshold",10);
        ARIADNE_TEST_EQUALS(metric.shifted_value_from(0),1);
        auto from_1 = metric.shifted_value_from(1);
        ARIADNE_TEST_ASSERT(from_1 == 0 or from_1 == 2);

        ARIADNE_TEST_EQUALS(metric.shifted_value_from(0),1);
        auto from_ub = metric.shifted_value_from(metric.upper_bound());
        ARIADNE_TEST_EQUALS(from_ub,metric.upper_bound()-1);
    }

    Void test_enumeration_task_parameter() {
        TaskParameter p = EnumerationTaskParameter<TestInterface>("test_class",{A(),B(),C(),D()});
        ARIADNE_TEST_PRINT(p);
        ARIADNE_TEST_ASSERT(not p.is_metric());
        ARIADNE_TEST_EQUALS(p.upper_bound(),3);
    }

    Void test_enumeration_task_parameter_shift() {
        EnumerationTaskParameter<TestInterface> e("enum_task",{A(),B(),C(),D()});
        auto shifted_0 = e.shifted_value_from(0);
        ARIADNE_TEST_ASSERT(shifted_0 == 1 or shifted_0 == 2 or shifted_0 == 3);
        auto shifted_1 = e.shifted_value_from(1);
        ARIADNE_TEST_ASSERT(shifted_1 == 0 or shifted_1 == 2 or shifted_1 == 3);
        auto shifted_2 = e.shifted_value_from(2);
        ARIADNE_TEST_ASSERT(shifted_2 == 0 or shifted_2 == 1 or shifted_2 == 3);
        auto shifted_3 = e.shifted_value_from(3);
        ARIADNE_TEST_ASSERT(shifted_3 == 0 or shifted_3 == 1 or shifted_3 == 2);
    }

    Void test_parameter_space() {
        BooleanTaskParameter boolean_task("boolean_task");
        MetricTaskParameter metric_task("metric_task",6);
        EnumerationTaskParameter<TestInterface> enum_task("enum_task",{A(),B(),C(),D()});

        TaskParameterSpace space({boolean_task,metric_task,enum_task});
        ARIADNE_TEST_PRINT(space);
        ARIADNE_TEST_EQUALS(space.dimension(),3);
        ARIADNE_TEST_EQUALS(space.index(boolean_task),0);
        ARIADNE_TEST_EQUALS(space.index(enum_task),1);
        ARIADNE_TEST_EQUALS(space.index(metric_task),2);
    }

    Void test_parameter_point_creation() {
        BooleanTaskParameter boolean_task("boolean_task");
        MetricTaskParameter metric_task("metric_task",6);
        EnumerationTaskParameter<TestInterface> enum_task("enum_task",{A(),B(),C(),D()});

        ParameterBindingsMap bindings = {{boolean_task,1},{metric_task,5},{enum_task,2}};
        ARIADNE_TEST_PRINT(bindings);
        TaskParameterPoint point(bindings);
        ARIADNE_TEST_PRINT(point);
        TaskParameterSpace space = point.space();
        ARIADNE_TEST_EQUALS(space.dimension(),3);
        ARIADNE_TEST_EQUALS(space.index(boolean_task),0);
        ARIADNE_TEST_EQUALS(space.index(enum_task),1);
        ARIADNE_TEST_EQUALS(space.index(metric_task),2);
    }

    Void test_parameter_point_equality() {
        BooleanTaskParameter b("boolean_task");
        MetricTaskParameter m("metric_task",6);
        EnumerationTaskParameter<TestInterface> e("enum_task",{A(),B(),C(),D()});

        TaskParameterPoint point1({{b,1},{m,5},{e,2}});
        TaskParameterPoint point2({{b,1},{e,2},{m,5}});
        TaskParameterPoint point3({{b,1},{e,3},{m,5}});
        ARIADNE_TEST_EQUAL(point1,point2);
        ARIADNE_TEST_NOT_EQUAL(point1,point3);
    }

    Void test_parameter_point_distance() {
        BooleanTaskParameter b("boolean_task");
        MetricTaskParameter m("metric_task",6);
        EnumerationTaskParameter<TestInterface> e("enum_task",{A(),B(),C(),D()});

        TaskParameterPoint point({{b,1},{m,5},{e,2}});
        ARIADNE_TEST_PRINT(point);
        TaskParameterPoint point2({{b,1},{m,5},{e,2}});
        ARIADNE_TEST_PRINT(point2);
        ARIADNE_TEST_EQUALS(point.distance(point2),0);
        TaskParameterPoint point3({{b,0},{m,5},{e,2}});
        ARIADNE_TEST_PRINT(point3);
        ARIADNE_TEST_EQUALS(point.distance(point3),1);
        TaskParameterPoint point4({{b,1},{m,5},{e,0}});
        ARIADNE_TEST_PRINT(point4);
        ARIADNE_TEST_EQUALS(point.distance(point4),1);
        TaskParameterPoint point5({{b,1},{m,8},{e,2}});
        ARIADNE_TEST_PRINT(point5);
        ARIADNE_TEST_EQUALS(point.distance(point5),3);
        TaskParameterPoint point6({{b,0},{m,4},{e,0}});
        ARIADNE_TEST_PRINT(point6);
        ARIADNE_TEST_EQUALS(point.distance(point6),3);
    }

    Void test_parameter_point_adjacent_shift() {
        BooleanTaskParameter b("boolean_task");
        MetricTaskParameter m("metric_task",6);
        EnumerationTaskParameter<TestInterface> e("enum_task",{A(),B(),C(),D()});

        TaskParameterPoint point1({{b,1},{m,5},{e,2}});
        auto points = point1.make_adjacent_shifted(1);
        auto point2 = *points.begin();
        ARIADNE_TEST_EQUALS(points.size(),1);
        ARIADNE_TEST_PRINT(point1);
        ARIADNE_TEST_PRINT(point2);
        ARIADNE_TEST_NOT_EQUAL(point1,point2);

        ARIADNE_PRINT_TEST_COMMENT("Checking multiple shifted points");
        points = point2.make_adjacent_shifted(3);
        ARIADNE_TEST_EQUALS(points.size(),3);
        for (auto point : points) {
            ARIADNE_TEST_PRINT(point);
            ARIADNE_TEST_NOT_EQUAL(point2,point);
        }
    }

    Void test_parameter_point_random_shift() {
        BooleanTaskParameter b("boolean_task");
        MetricTaskParameter m("metric_task",6);
        EnumerationTaskParameter<TestInterface> e("enum_task",{A(),B(),C(),D()});

        TaskParameterPoint point1({{b,1},{m,5},{e,2}});
        auto points = point1.make_random_shifted(5);
        ARIADNE_TEST_EQUALS(points.size(),5);
        for (auto point : points) {
            ARIADNE_TEST_PRINT(point);
        }
    }

    Void test() {
        ARIADNE_TEST_CALL(test_boolean_task_parameter());
        ARIADNE_TEST_CALL(test_boolean_task_parameter_shift());
        ARIADNE_TEST_CALL(test_metric_task_parameter());
        ARIADNE_TEST_CALL(test_metric_task_parameter_shift());
        ARIADNE_TEST_CALL(test_enumeration_task_parameter());
        ARIADNE_TEST_CALL(test_enumeration_task_parameter_shift());
        ARIADNE_TEST_CALL(test_parameter_space());
        ARIADNE_TEST_CALL(test_parameter_point_creation());
        ARIADNE_TEST_CALL(test_parameter_point_equality());
        ARIADNE_TEST_CALL(test_parameter_point_distance());
        ARIADNE_TEST_CALL(test_parameter_point_adjacent_shift());
        ARIADNE_TEST_CALL(test_parameter_point_random_shift());
    }
};

int main() {
    TestTaskParameter().test();
    return ARIADNE_TEST_FAILURES;
}
