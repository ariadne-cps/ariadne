/***************************************************************************
 *            test_multifunction.cpp
 *
 *  Copyright  2020  Pieter Collins
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

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <fenv.h>

#include "config.hpp"
#include "algebra/sweeper.hpp"
#include "function/multifunction.hpp"
#include "numeric/numeric.hpp"

#include "symbolic/variables.hpp"
#include "symbolic/space.hpp"
#include "symbolic/expression.hpp"

#include "../test.hpp"

using namespace std;
using namespace Ariadne;

class TestMultifunction
{
  public:
    Void test();
  private:
    Void test_concept();
    Void test_evaluate();
};

Void TestMultifunction::test()
{
    Dyadic::set_default_writer(DecimalWriter());
    ARIADNE_TEST_CALL(test_evaluate());
}

Void TestMultifunction::test_concept()
{
}

Void TestMultifunction::test_evaluate()
{
    RealVariable x0("x0"), x1("x1"), p("p");

    BoxDomainType dom({{-1,3},{-1,3}});
    BoxDomainType prms({{-1,1}});;

    EffectiveVectorMultivariateFunction f=make_function({x0,x1,p},{1.5_dy-x0*x0-0.375_dy*x1+p,x0});

    Sweeper<FloatDP> swp;

    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateMultifunctionPatch,mvfp,(f,prms));
    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateMultifunction,mvf,(mvfp));

    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateMultifunctionModel<DP>,mvfmdp,(dom,f,prms,swp));;
    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateMultifunctionModel<DP>,mvfdp,(mvfmdp));;

    Vector<ValidatedNumber> v(2,1);
    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedLocatedSet<RealVector>,res,mvfp(v));
    ARIADNE_TEST_ASSIGN(res,mvf(v));
    ARIADNE_TEST_ASSIGN(res,mvfmdp(v));
    ARIADNE_TEST_ASSIGN(res,mvfdp(v));
}


Int main() {
    TestMultifunction().test();

    return ARIADNE_TEST_FAILURES;
}

