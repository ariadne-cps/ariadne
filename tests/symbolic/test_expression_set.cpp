/***************************************************************************
 *            test_expression_set.cpp
 *
 *  Copyright  2009-21  Luca Geretti
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

#include <iostream>

#include "config.hpp"
#include "utility/container.hpp"
#include "utility/stlio.hpp"
#include "numeric/numeric.hpp"
#include "symbolic/expression_set.hpp"
#include "symbolic/assignment.hpp"
#include "symbolic/valuation.hpp"
#include "symbolic/space.hpp"
#include "function/formula.hpp"
#include "geometry/function_set.hpp"
#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "output/logging.hpp"

#include "../test.hpp"

using namespace Ariadne;

using RealVariablesPoint = VariablesPoint<Real>;

class TestExpressionSet {
  private:
    RealVariable x, y, z;

  public:
    TestExpressionSet() : x("x"), y("y"), z("z") { }

    Void test() {
        ARIADNE_TEST_CALL(test_labelled_set());
        ARIADNE_TEST_CALL(test_variables_point());
        ARIADNE_TEST_CALL(test_variables_box());
        ARIADNE_TEST_CALL(test_real_expression_constraint_set());
        ARIADNE_TEST_CALL(test_real_expression_bounded_constraint_set());
    }

    Void test_labelled_set() const {
        LabelledSet<RealBox> set(RealSpace({x,y,z}),RealBox({{0.1_dec,0.5_dec},{-2.3_dec,1},{0,3}}));
        ARIADNE_TEST_EQUALS(set.dimension(),3);
        RealVariable u("u");
        ARIADNE_TEST_EXECUTE(set.euclidean_set(RealSpace({y,z})));
        ARIADNE_TEST_FAIL(set.euclidean_set(RealSpace({y,u})));
    }

    Void test_variables_point() const {
        Map<RealVariable,Real> bindings;
        bindings.insert(x,0.1_dec);
        bindings.insert(y,-1);
        RealVariablesPoint p(bindings);
        ARIADNE_TEST_EQUALS(p.variables().size(),2);
        ARIADNE_TEST_EQUALS(p.euclidean_set(RealSpace({x,y})).dimension(),2);
        ARIADNE_TEST_FAIL(p.euclidean_set(RealSpace({x,z})));

        auto labelled_set = LabelledSet<RealPoint>(p);
        ARIADNE_TEST_EQUALS(labelled_set.dimension(),2);
    }

    Void test_variables_box() const {
        RealVariablesBox bx1 = {0<=x<=1.2_dec};
        ARIADNE_TEST_ASSERT(not bx1.is_empty());

        RealVariablesBox bx2(RealSpace({x,y}),RealBox({{0.1_dec,0.5_dec},{-2.3_dec,1}}));
        ARIADNE_TEST_EQUALS(bx2.variables().size(),2);

        auto labelled_set = LabelledSet<RealBox>(bx2);
        ARIADNE_TEST_EQUALS(labelled_set.dimension(),2);
    }

    Void test_real_expression_constraint_set() const {
        RealExpressionConstraintSet set({y >= 0, z + y <= 2, sqr(z) + sqr(y) >= 0});
        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_EQUALS(set.variables().size(), 2);
        ARIADNE_TEST_EQUALS(set.constraints().size(), 3);
        ARIADNE_TEST_FAIL(set.euclidean_set(RealSpace({y})));
        ARIADNE_TEST_FAIL(set.euclidean_set(RealSpace({x})));
        ARIADNE_TEST_EXECUTE(set.euclidean_set(RealSpace({y, z})));
        ARIADNE_TEST_EXECUTE(set.euclidean_set(RealSpace({z, y})));
    }

    Void test_real_expression_bounded_constraint_set() const {
        RealExpressionConstraintSet unbounded_set({y>=0,z+y<=2,sqr(z)+sqr(y)>=0});
        RealVariable u("u");
        ARIADNE_TEST_EXECUTE(RealExpressionBoundedConstraintSet(RealVariablesBox({0.1_dec<=x<=0.5_dec,-2.3_dec<=u<=1}),unbounded_set));
        ARIADNE_TEST_EXECUTE(RealExpressionBoundedConstraintSet(RealVariablesBox({0.1_dec<=x<=0.5_dec,-2.3_dec<=z<=1,-1<=y<=2}),unbounded_set));
        ARIADNE_TEST_EXECUTE(intersection(RealVariablesBox({0.1_dec<=z<=0.5_dec,-2.3_dec<=y<=1}),unbounded_set));
    }
};

Int main(Int argc, const char* argv[])
{
    ARIADNE_LOG_SET_VERBOSITY(get_verbosity(argc,argv));

    TestExpressionSet().test();
    return ARIADNE_TEST_FAILURES;
}
