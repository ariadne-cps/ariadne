/***************************************************************************
 *            test_hybrid_enclosure.cpp
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

#include <iostream>

#include "config.hpp"
#include "../test.hpp"

#include "algebra/algebra.hpp"
#include "function/taylor_function.hpp"
#include "hybrid/hybrid_set.hpp"
#include "hybrid/hybrid_enclosure.hpp"

using namespace std;
using namespace Ariadne;

class TestHybridEnclosure {
  public:
    Void test() {
        ARIADNE_TEST_CALL(test_construct())
        ARIADNE_TEST_CALL(test_space_having_time())
    }

    void test_construct() {
        HybridEnclosure encl;

        TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp,1e-8));
        EnclosureConfiguration config(function_factory);
        HybridEnclosure encl2(config);
        ARIADNE_TEST_PRINT(encl2.configuration())

        DiscreteLocation location;
        RealVariable x("x"), y("y");
        HybridRealBox box(location,{1<=x<=2,0<=y<=1});
        HybridEnclosure encl3(box,config);

        Enclosure cencl(config);
        HybridEnclosure(location,RealSpace({x,y}),cencl);
    }

    void test_space_having_time() {
        TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp,1e-8));
        EnclosureConfiguration config(function_factory);
        RealVariable x("x"), t("t");
        HybridRealBox box(DiscreteLocation(),{1<=t<=2,0<=x<=1});
        HybridEnclosure encl(box,config);
        ARIADNE_TEST_PRINT(encl.state_time_auxiliary_space())
    }
};

Int main() {
    TestHybridEnclosure().test();
    return ARIADNE_TEST_FAILURES;
}

