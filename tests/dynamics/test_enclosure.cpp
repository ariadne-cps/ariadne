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
        Figure fig;
        fig.set_bounding_box(encl.bounding_box());
        fig << line_style(true) << encl;
        fig.write(("test_enclosure_"+suffix).c_str());
    }

    Void test_draw() const {
        _draw(BoxDrawer(),"Box");
        _draw(GridDrawer(4),"Grid");
        _draw(AffineDrawer(1),"Affine");
    }
};

Int main()
{
    ARIADNE_TEST_CALL(TestEnclosure().test());
    return ARIADNE_TEST_FAILURES;
}
