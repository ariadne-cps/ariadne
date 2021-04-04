/***************************************************************************
 *            test_enclosure.cpp
 *
 *  Copyright  2006-21  Luca Geretti
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

#include <fstream>
#include <iostream>

#include "config.hpp"
#include "utility/tuple.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/constraint.hpp"
#include "function/taylor_function.hpp"
#include "dynamics/enclosure.hpp"
#include "geometry/box.hpp"
#include "geometry/list_set.hpp"
#include "geometry/function_set.hpp"
#include "solvers/integrator.hpp"
#include "symbolic/expression_set.hpp"
#include "output/graphics.hpp"
#include "output/drawer.hpp"
#include "output/logging.hpp"

#include "function/user_function.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

class TestEnclosure
{
  public:
    Void test() const {
        ARIADNE_TEST_CALL(test_construct_without_domain());
        ARIADNE_TEST_CALL(test_construct_with_domain());
        ARIADNE_TEST_CALL(test_is_bounded());
        ARIADNE_TEST_CALL(test_subset_separated());
        ARIADNE_TEST_CALL(test_restriction());
        ARIADNE_TEST_CALL(test_auxiliary_map());
        ARIADNE_TEST_CALL(test_constraints());
        ARIADNE_TEST_CALL(test_draw());
    }

    Void test_construct_without_domain() const {
        Enclosure encl;
        TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp,1e-8));
        EnclosureConfiguration configuration(function_factory);
        Enclosure encl_config(configuration);
    }

    Void test_construct_with_domain() const {
        ExactBoxType dom({{0.0_x,2.0_x},{1.0_x,3.0_x}});
        TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp,1e-8));
        EnclosureConfiguration configuration(function_factory);
        Enclosure encl(dom,configuration);
        ARIADNE_TEST_PRINT(configuration);
        ARIADNE_TEST_PRINT(encl);
        ARIADNE_TEST_EQUALS(encl.centre(),FloatDPValuePoint({FloatDPValue(1.0_x,DoublePrecision()),FloatDPValue(2.0_x,DoublePrecision())}));
        ARIADNE_TEST_EQUALS(encl.radius().raw(),FloatDP(1.0_x,DoublePrecision()));
        ARIADNE_TEST_ASSERT(dom.inside(encl.codomain()));
    }

    Void test_is_bounded() const {
        ExactBoxType bounded_dom({{0.0_x,2.0_x},{1.0_x,3.0_x}});
        EnclosureConfiguration configuration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,1e-8)));
        Enclosure encl(bounded_dom,configuration);
        ARIADNE_TEST_ASSERT(encl.is_bounded());
        ExactBoxType unbounded_dom({{0.0_x,infty},{1.0_x,3.0_x}});
        Enclosure encl2(unbounded_dom,configuration);
        ARIADNE_TEST_ASSERT(is_indeterminate(encl2.is_bounded()));
    }

    Void test_subset_separated() const {
        ExactBoxType box1({{0.0_x,2.0_x},{1.0_x,3.0_x}});
        Enclosure encl(box1,EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,1e-8))));
        ARIADNE_TEST_ASSERT(definitely(encl.subset(box1)));
        ARIADNE_TEST_ASSERT(possibly(not encl.separated(box1)));
        ExactBoxType box2({{2.0_x,3.0_x},{-1.0_x,0.0_x}});
        ARIADNE_TEST_ASSERT(possibly(not encl.subset(box2)));
        ARIADNE_TEST_ASSERT(definitely(encl.separated(box2)));
    }

    Void test_restriction() const {
        ExactBoxType box1({{0.0_x,2.0_x},{1.0_x,3.0_x}});
        Enclosure encl(box1,EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,1e-8))));
        auto codomain1 = encl.codomain();
        ExactBoxType box2({{0.5_x,1.5_x},{1.5_x,2.5_x}});
        auto encl2 = restriction(encl,box2);
        ARIADNE_TEST_ASSERT(definitely(encl2.subset(codomain1)));
    }

    Void test_auxiliary_map() const {
        ExactBoxType dom({{0.0_x,2.0_x},{1.0_x,3.0_x}});
        Enclosure encl(dom,EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,1e-8))));
        auto x0 = EffectiveScalarMultivariateFunction::coordinate(2,0);
        auto x1 = EffectiveScalarMultivariateFunction::coordinate(2,1);
        EffectiveVectorMultivariateFunction auxiliary(1,x0+sqr(x1));
        encl.set_auxiliary_mapping(auxiliary);
        ARIADNE_TEST_PRINT(encl.state_auxiliary_function());
        ARIADNE_TEST_PRINT(encl.state_auxiliary_set());
        for (SizeType i=0; i<4; ++i)
            ARIADNE_TEST_PRINT(encl.get_function(i));
    }

    Void test_constraints() const {
        ExactBoxType dom({{0.0_x,2.0_x},{1.0_x,3.0_x}});
        Enclosure encl(dom,EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,1e-8))));
        auto x0 = EffectiveScalarMultivariateFunction::coordinate(2,0);
        auto x1 = EffectiveScalarMultivariateFunction::coordinate(2,1);
        encl.new_zero_state_constraint(x0+sqr(x1));
        encl.new_positive_state_constraint(x0);
        ARIADNE_TEST_EQUALS(encl.number_of_constraints(),2);
        ARIADNE_TEST_PRINT(encl.constraint_function());
        ARIADNE_TEST_PRINT(encl.constraint_models());
        ARIADNE_TEST_PRINT(encl.constraint(0));
        ARIADNE_TEST_PRINT(encl.constraint(1));
        ARIADNE_TEST_PRINT(encl.constraint_bounds());
        ARIADNE_TEST_ASSERT(definitely(encl.satisfies(x1+x0)));
        ARIADNE_TEST_ASSERT(definitely(not encl.satisfies(-sqr(x1))));
        ARIADNE_TEST_ASSERT(is_indeterminate(encl.satisfies(-sqr(x0))));
    }

    Void _draw(Drawer const& drawer, String suffix) const {
        ARIADNE_PRINT_TEST_COMMENT("Drawing with " + suffix + " method");
        ARIADNE_TEST_PRINT(drawer);
        ExactBoxType dom({{0.0_x,1.0_x},{1.0_x,2.0_x}});
        auto swp = ThresholdSweeper<FloatDP>(dp,1e-8);
        auto fnc = antiderivative(ValidatedVectorMultivariateTaylorFunctionModelDP::identity(dom,swp),0);
        TaylorFunctionFactory function_factory(swp);
        EnclosureConfiguration configuration(function_factory);
        configuration.set_drawer(drawer);
        Enclosure encl(dom,fnc,configuration);
        Figure fig(dom,Projection2d(2,0,1));;
        fig << line_style(true) << encl;
        fig.write(("test_enclosure_"+suffix).c_str());
    }

    Void test_draw() const {
        _draw(BoxDrawer(),"box");
        _draw(GridDrawer(4),"grid");
        _draw(AffineDrawer(1),"affine");
    }
};

Int main()
{
    ARIADNE_TEST_CALL(TestEnclosure().test());
    return ARIADNE_TEST_FAILURES;
}
