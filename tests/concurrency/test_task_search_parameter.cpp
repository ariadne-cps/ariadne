/***************************************************************************
 *            test_task_search_parameter.cpp
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
#include "concurrency/task_search_point.hpp"
#include "concurrency/task_search_space.hpp"
#include "concurrency/task_appraisal.hpp"

#include "../test.hpp"

using namespace Ariadne;

class TestTaskSearchParameter {
  public:

    Void test_task_parameter_creation() {
        TaskSearchParameter p("use_subdivisions", false, List<int>({0,2}));
        ARIADNE_TEST_PRINT(p);
    }

    Void test_metric_task_parameter_shift() {
        TaskSearchParameter metric("sweep_threshold", true, List<int>({8,9,10,11}));
        ARIADNE_TEST_EQUALS(metric.shifted_value_from(8),9);
        ARIADNE_TEST_EQUALS(metric.shifted_value_from(11),10);
        auto from_1 = metric.shifted_value_from(10);
        ARIADNE_TEST_ASSERT(from_1 == 9 or from_1 == 11);
    }

    Void test_parameter_space() {
        TaskSearchParameter bp("use_subdivisions", false,List<int>({0,1}));
        TaskSearchParameter mp("sweep_threshold", true, List<int>({3,4,5}));
        TaskSearchSpace space({bp, mp});
        ARIADNE_TEST_PRINT(space);
        ARIADNE_TEST_PRINT(space.parameters());
        ARIADNE_TEST_EQUALS(space.dimension(),2);
        ARIADNE_TEST_EQUALS(space.index(bp),1);
        ARIADNE_TEST_EQUALS(space.index(mp),0);
        ARIADNE_TEST_EQUALS(space.total_points(),6);
    }

    Void test_parameter_point_creation() {
        TaskSearchParameter bp("use_subdivisions", false,List<int>({0,1}));
        TaskSearchParameter mp("sweep_threshold", true, List<int>({3,4,5}));
        TaskSearchSpace space({bp, mp});
        ARIADNE_TEST_PRINT(space.initial_point());
        Map<Identifier,int> bindings = {{"use_subdivisions",1},
                                          {"sweep_threshold",5}};
        ARIADNE_TEST_PRINT(bindings);
        TaskSearchPoint point = space.make_point(bindings);
        ARIADNE_TEST_PRINT(point);
        ARIADNE_TEST_PRINT(point.space());
    }

    Void test_parameter_point_equality() {
        TaskSearchParameter bp("use_subdivisions", false,List<int>({0,1}));
        TaskSearchParameter mp("sweep_threshold", true, List<int>({3,4,5}));
        TaskSearchSpace space({bp, mp});

        TaskSearchPoint point1 = space.make_point({{"use_subdivisions", 1}, {"sweep_threshold", 2}});
        TaskSearchPoint point2 = space.make_point({{"use_subdivisions", 1}, {"sweep_threshold", 2}});
        TaskSearchPoint point3 = space.make_point({{"use_subdivisions", 1}, {"sweep_threshold", 3}});
        ARIADNE_TEST_EQUAL(point1,point2);
        ARIADNE_TEST_NOT_EQUAL(point1,point3);
    }

    Void test_parameter_point_distance() {
        TaskSearchParameter bp("use_subdivisions", false,List<int>({0,1}));
        TaskSearchParameter mp("sweep_threshold", true, List<int>({3,4,5}));
        TaskSearchSpace space({bp, mp});

        TaskSearchPoint point1 = space.make_point({{"use_subdivisions", 1}, {"sweep_threshold", 2}});
        TaskSearchPoint point2 = space.make_point({{"use_subdivisions", 1}, {"sweep_threshold", 2}});
        TaskSearchPoint point3 = space.make_point({{"use_subdivisions", 1}, {"sweep_threshold", 3}});
        TaskSearchPoint point4 = space.make_point({{"use_subdivisions", 0}, {"sweep_threshold", 5}});
        ARIADNE_TEST_EQUALS(point1.distance(point2),0);
        ARIADNE_TEST_EQUALS(point1.distance(point3),1);
        ARIADNE_TEST_EQUALS(point3.distance(point4),3);
    }

    Void test_parameter_point_adjacent_shift() {
        TaskSearchParameter bp("use_subdivisions", false,List<int>({0,1}));
        TaskSearchParameter mp("sweep_threshold", true, List<int>({3,4,5}));
        TaskSearchSpace space({bp, mp});

        TaskSearchPoint point1 = space.make_point({{"use_subdivisions", 1}, {"sweep_threshold", 5}});
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
        TaskSearchParameter bp("use_subdivisions", false,List<int>({0,1}));
        TaskSearchParameter mp("sweep_threshold", true, List<int>({3,4,5,6,7}));
        TaskSearchSpace space({bp, mp});

        TaskSearchPoint point1 = space.make_point({{"use_subdivisions", 1}, {"sweep_threshold", 5}});
        auto points = point1.make_random_shifted(3);
        ARIADNE_TEST_EQUALS(points.size(),3);
        for (auto point : points) {
            ARIADNE_TEST_PRINT(point);
        }
    }

    Void test_parameter_point_adjacent_set_shift() {
        TaskSearchParameter bp("use_subdivisions", false,List<int>({0,1}));
        TaskSearchParameter mp("sweep_threshold", true, List<int>({3,4,5,6,7}));
        TaskSearchSpace space({bp, mp});

        TaskSearchPoint point1 = space.make_point({{"use_subdivisions", 1}, {"sweep_threshold", 5}});
        Set<TaskSearchPoint> points = point1.make_random_shifted(3);
        ARIADNE_TEST_PRINT(points);
        auto all_points = make_extended_set_by_shifting(points, 5);
        ARIADNE_TEST_EQUALS(all_points.size(),5);
        ARIADNE_TEST_PRINT(all_points);
    }

    Void test_parameter_point_appraisal() {
        TaskSearchParameter bp("use_subdivisions", false,List<int>({0,1}));
        TaskSearchParameter mp("sweep_threshold", true, List<int>({3,4,5,6,7}));
        TaskSearchSpace space({bp, mp});

        TaskSearchPoint point1 = space.make_point({{"use_subdivisions", 1}, {"sweep_threshold", 2}});
        TaskSearchPoint point2 = space.make_point({{"use_subdivisions", 1}, {"sweep_threshold", 2}});
        TaskSearchPoint point3 = space.make_point({{"use_subdivisions", 1}, {"sweep_threshold", 3}});
        TaskSearchPoint point4 = space.make_point({{"use_subdivisions", 0}, {"sweep_threshold", 4}});
        TaskSearchPointAppraisal a1(point1,3,0,0);
        TaskSearchPointAppraisal a2(point2,2,1,0);
        TaskSearchPointAppraisal a3(point3,4,0,0);
        TaskSearchPointAppraisal a4(point4,4,0,1);
        Set<TaskSearchPointAppraisal> as = {a1, a2, a3, a4};

        ARIADNE_TEST_PRINT(as);
        ARIADNE_TEST_ASSERT(a1 < a2);
        ARIADNE_TEST_ASSERT(a1 < a3);
        ARIADNE_TEST_ASSERT(a3 < a2);
        ARIADNE_TEST_ASSERT(a2 < a4);
        ARIADNE_TEST_ASSERT(a3 < a4);
    }

    Void test() {
        ARIADNE_TEST_CALL(test_metric_task_parameter_shift());
        ARIADNE_TEST_CALL(test_parameter_space());
        ARIADNE_TEST_CALL(test_parameter_point_creation());
        ARIADNE_TEST_CALL(test_parameter_point_equality());
        ARIADNE_TEST_CALL(test_parameter_point_distance());
        ARIADNE_TEST_CALL(test_parameter_point_adjacent_shift());
        ARIADNE_TEST_CALL(test_parameter_point_random_shift());
        ARIADNE_TEST_CALL(test_parameter_point_adjacent_set_shift());
        ARIADNE_TEST_CALL(test_parameter_point_appraisal());
    }
};

int main() {
    TestTaskSearchParameter().test();
    return ARIADNE_TEST_FAILURES;
}
