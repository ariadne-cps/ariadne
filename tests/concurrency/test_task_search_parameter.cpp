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

#include "../test.hpp"

using namespace Ariadne;

struct TestInterfaceBase {
    TestInterfaceBase(Nat val) : _val(val) { }
    bool operator==(TestInterfaceBase const& o) const { return _val==o._val; }
  private:
    Nat _val;
};
struct A : public TestInterfaceBase { A() : TestInterfaceBase(0) { }};
struct B : public TestInterfaceBase { B() : TestInterfaceBase(1) { }};
struct C : public TestInterfaceBase { C() : TestInterfaceBase(2) { }};
struct D : public TestInterfaceBase { D() : TestInterfaceBase(3) { }};

class TestTaskSearchParameter {
  public:

    TestTaskSearchParameter() { srand (time(NULL)); }

    Void test_boolean_task_parameter() {
        TaskSearchParameter p = BooleanSearchParameter(RealVariable("use_subdivisions"), false);
        ARIADNE_TEST_PRINT(p);
    }

    Void test_boolean_task_parameter_shift() {
        BooleanSearchParameter b(RealVariable("use_subdivisions"), false);
        ARIADNE_TEST_EQUALS(b.shifted_value_from(0),1);
        ARIADNE_TEST_EQUALS(b.shifted_value_from(1),0);
    }

    Void test_metric_task_parameter() {
        RealVariable st("sweep_threshold");
        RealExpression expr = st;
        TaskSearchParameter metric = MetricSearchParameter(st, expr, 10, 8);
        ARIADNE_TEST_PRINT(metric);
    }

    Void test_metric_task_parameter_shift() {
        RealVariable st("sweep_threshold");
        MetricSearchParameter metric(st, 10, 8);
        ARIADNE_TEST_EQUALS(metric.shifted_value_from(0),1);
        auto from_1 = metric.shifted_value_from(1);
        ARIADNE_TEST_ASSERT(from_1 == 0 or from_1 == 2);

        ARIADNE_TEST_EQUALS(metric.shifted_value_from(0),1);
        auto from_ub = metric.shifted_value_from(metric.upper_bound());
        ARIADNE_TEST_EQUALS(from_ub,metric.upper_bound()-1);
    }

    Void test_enumeration_task_parameter() {
        RealVariable tc("integrator");
        TaskSearchParameter p = EnumerationSearchParameter<TestInterfaceBase>(tc, {A(), B(), C(), D()}, B());
        auto etp_ptr = dynamic_cast<EnumerationSearchParameter<TestInterfaceBase>*>(p.ptr());
        ARIADNE_TEST_EQUALS(etp_ptr->elements().size(),4);
        ARIADNE_TEST_PRINT(p);
        ARIADNE_TEST_EQUALS(p.upper_bound(),3);
    }

    Void test_enumeration_task_parameter_shift() {
        RealVariable et("integrator");
        EnumerationSearchParameter<TestInterfaceBase> e(et, {A(), B(), C(), D()}, B());
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
        RealVariable b("use_subdivisions"), m("sweep_threshold"), e("integrator");
        BooleanSearchParameter bp(b, false);
        MetricSearchParameter mp(m, 10, 8);
        EnumerationSearchParameter<TestInterfaceBase> ep(e, {A(), B(), C(), D()}, B());
        TaskSearchSpace space({bp, mp, ep}, b * m + e);
        ARIADNE_TEST_PRINT(space);
        ARIADNE_TEST_PRINT(space.parameters());
        ARIADNE_TEST_PRINT(space.time_cost_estimator());
        ARIADNE_TEST_EQUALS(space.dimension(),3);
        ARIADNE_TEST_EQUALS(space.index(bp),2);
        ARIADNE_TEST_EQUALS(space.index(ep),0);
        ARIADNE_TEST_EQUALS(space.index(mp),1);
    }

    Void test_parameter_point_creation() {
        RealVariable b("use_subdivisions"), m("sweep_threshold"), e("integrator");
        TaskSearchSpace space({BooleanSearchParameter(b, false),
                               MetricSearchParameter(m, 10, 8),
                               EnumerationSearchParameter<TestInterfaceBase>(e, {A(), B(), C(), D()}, B())},
                                 b*m+e);
        ARIADNE_TEST_PRINT(space.initial_point());
        Map<RealVariable,Nat> bindings = {{b,1},{m,5},{e,2}};
        ARIADNE_TEST_PRINT(bindings);
        TaskSearchPoint point = space.make_point(bindings);
        ARIADNE_TEST_PRINT(point);
        ARIADNE_TEST_PRINT(point.space());
    }

    Void test_parameter_point_time_cost_estimation() {
        RealVariable mss("maximum_step_size"), st("sweep_threshold"), mto("maximum_temporal_order"), rsp("relative_set_parameters");
        RealVariable lip("lipschitz_step");
        TaskSearchSpace space({MetricSearchParameter(mss, lip * exp(-mss), 10, 0),
                               MetricSearchParameter(st, exp(-st), 10, 8),
                               MetricSearchParameter(mto, 10, 5),
                               MetricSearchParameter(rsp, 6, 3)
                                  },(st*mto+rsp)/mss);

        TaskSearchPoint point = space.make_point({{mss, 2}, {mto, 6}, {rsp, 3}, {st, 8}});
        Map<RealVariable,Real> constant_values;
        constant_values[lip] =0.2_decimal;
        auto estimate = point.time_cost_estimate(constant_values);
        ARIADNE_TEST_PRINT(point.values(constant_values));
        ARIADNE_TEST_PRINT(estimate);
        ARIADNE_TEST_PRINT(estimate.get(DoublePrecision()));
    }

    Void test_parameter_point_equality() {
        RealVariable b("use_subdivisions"), m("sweep_threshold"), e("integrator");
        TaskSearchSpace space({BooleanSearchParameter(b, false),
                               MetricSearchParameter(m, 10, 8),
                               EnumerationSearchParameter<TestInterfaceBase>(e, {A(), B(), C(), D()}, B())},
                                 b*m+e);

        TaskSearchPoint point1 = space.make_point({{b, 1}, {m, 5}, {e, 2}});
        TaskSearchPoint point2 = space.make_point({{b, 1}, {e, 2}, {m, 5}});
        TaskSearchPoint point3 = space.make_point({{b, 1}, {e, 3}, {m, 5}});
        ARIADNE_TEST_EQUAL(point1,point2);
        ARIADNE_TEST_NOT_EQUAL(point1,point3);
    }

    Void test_parameter_point_distance() {
        RealVariable b("use_subdivisions"), m("sweep_threshold"), e("integrator");
        TaskSearchSpace space({BooleanSearchParameter(b, false),
                               MetricSearchParameter(m, 10, 8),
                               EnumerationSearchParameter<TestInterfaceBase>(e, {A(), B(), C(), D()}, B())},
                                 b*m+e);

        TaskSearchPoint point = space.make_point({{b, 1}, {m, 5}, {e, 2}});
        ARIADNE_TEST_PRINT(point);
        TaskSearchPoint point2 = space.make_point({{b, 1}, {m, 5}, {e, 2}});
        ARIADNE_TEST_PRINT(point2);
        ARIADNE_TEST_EQUALS(point.distance(point2),0);
        TaskSearchPoint point3 = space.make_point({{b, 0}, {m, 5}, {e, 2}});
        ARIADNE_TEST_PRINT(point3);
        ARIADNE_TEST_EQUALS(point.distance(point3),1);
        TaskSearchPoint point4 = space.make_point({{b, 1}, {m, 5}, {e, 0}});
        ARIADNE_TEST_PRINT(point4);
        ARIADNE_TEST_EQUALS(point.distance(point4),1);
        TaskSearchPoint point5 = space.make_point({{b, 1}, {m, 8}, {e, 2}});
        ARIADNE_TEST_PRINT(point5);
        ARIADNE_TEST_EQUALS(point.distance(point5),3);
        TaskSearchPoint point6 = space.make_point({{b, 0}, {m, 4}, {e, 0}});
        ARIADNE_TEST_PRINT(point6);
        ARIADNE_TEST_EQUALS(point.distance(point6),3);
    }

    Void test_parameter_point_adjacent_shift() {
        RealVariable b("use_subdivisions"), m("sweep_threshold"), e("integrator");
        TaskSearchSpace space({BooleanSearchParameter(b, false),
                               MetricSearchParameter(m, 10, 8),
                               EnumerationSearchParameter<TestInterfaceBase>(e, {A(), B(), C(), D()}, B())},
                                 b*m+e);

        TaskSearchPoint point1 = space.make_point({{b, 1}, {m, 5}, {e, 2}});
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
        RealVariable b("use_subdivisions"), m("sweep_threshold"), e("integrator");
        TaskSearchSpace space({BooleanSearchParameter(b, false),
                               MetricSearchParameter(m, 10, 8),
                               EnumerationSearchParameter<TestInterfaceBase>(e, {A(), B(), C(), D()}, B())},
                                 b*m+e);

        TaskSearchPoint point1 = space.make_point({{b, 1}, {m, 5}, {e, 2}});
        auto points = point1.make_random_shifted(5);
        ARIADNE_TEST_EQUALS(points.size(),5);
        for (auto point : points) {
            ARIADNE_TEST_PRINT(point);
        }
    }

    Void test_parameter_point_adjacent_set_shift() {
        RealVariable b("use_subdivisions"), m1("sweep_threshold"), m2("maximum_step_size"), e("integrator");
        TaskSearchSpace space({BooleanSearchParameter(b, false),
                               MetricSearchParameter(m1, 10, 8),
                               MetricSearchParameter(m2, 6, 0),
                               EnumerationSearchParameter<TestInterfaceBase>(e, {A(), B(), C(), D()}, B())},
                                 b*m1+m2+e);

        TaskSearchPoint starting_point = space.make_point({{b, 1}, {m1, 5}, {m2, 2}, {e, 2}});
        Set<TaskSearchPoint> points = starting_point.make_random_shifted(8);
        ARIADNE_TEST_PRINT(points);
        auto all_points = make_extended_set_by_shifting(points, 13);
        ARIADNE_TEST_EQUALS(all_points.size(),13);
        ARIADNE_TEST_PRINT(all_points);
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
        ARIADNE_TEST_CALL(test_parameter_point_time_cost_estimation());
        ARIADNE_TEST_CALL(test_parameter_point_equality());
        ARIADNE_TEST_CALL(test_parameter_point_distance());
        ARIADNE_TEST_CALL(test_parameter_point_adjacent_shift());
        ARIADNE_TEST_CALL(test_parameter_point_random_shift());
        ARIADNE_TEST_CALL(test_parameter_point_adjacent_set_shift());
    }
};

int main() {
    TestTaskSearchParameter().test();
    return ARIADNE_TEST_FAILURES;
}
