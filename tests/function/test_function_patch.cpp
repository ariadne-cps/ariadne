/***************************************************************************
 *            test_function_patch.cpp
 *
 *  Copyright  2022  Pieter Collins
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
#include "function/function_patch.hpp"
#include "function/function.hpp"
#include "function/formula.hpp"
#include "function/symbolic_function.hpp"
#include "function/taylor_model.hpp"
#include "function/taylor_function.hpp"
#include "algebra/algebra.hpp"

#include "numeric/numeric.hpp"

#include "../test.hpp"

using namespace std;
using namespace Ariadne;

class TestFunctionPatch
{
  public:
    Void test();
  private:
    Void test_concept();
};

Void TestFunctionPatch::test()
{
    ARIADNE_TEST_CALL(test_concept());
}

Void TestFunctionPatch::test_concept()
{
    auto swp=ThresholdSweeper<FloatDP>(dp,1e-8);
    auto dom=BoxDomainType({{1,5},{2,3}});
    auto fm=ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate(dom,0,swp);
    auto fp=ValidatedScalarMultivariateFunctionPatch(fm);
    auto f=ValidatedScalarMultivariateFunction(fp);

#warning
//    static_assert(AFunction<ValidatedScalarMultivariateFunctionPatch,ValidatedTag,RealScalar(RealVector)>);

    auto x=FloatDPBoundsVector({2,2.75_x},dp);
    f(x);
}

#warning
namespace Ariadne {

template<class P, class SIG> class FunctionArchetype;
template<class P, class SIG> class FunctionPatchArchetype;
template<class P, class SIG> class EntireFunctionArchetype;

template<> class FunctionArchetype<ApproximateTag,Real(Real)> {
  public:
    SizeOne argument_size() const;
    SizeOne result_size() const;

    template<class T> T operator() (T const& a) const;
};

template<> class FunctionPatchArchetype<ApproximateTag,Real(Real)> : public FunctionArchetype<ApproximateTag,Real(Real)> {
  public:
    IntervalDomainType domain() const;
};

template<> class EntireFunctionArchetype<ApproximateTag,Real(Real)> : public FunctionArchetype<ApproximateTag,Real(Real)> {
  public:
    RealDomain domain() const;
};

#warning Find some better way of communicating what we expect from function domains
/*
static_assert(not Convertible<RealDomain,IntervalDomainType>);

static_assert(AFunction<FunctionArchetype<ApproximateTag,Real(Real)>,ApproximateTag,Real(Real)>);
static_assert(AFunction<FunctionPatchArchetype<ApproximateTag,Real(Real)>,ApproximateTag,Real(Real)>);
static_assert(AFunction<EntireFunctionArchetype<ApproximateTag,Real(Real)>,ApproximateTag,Real(Real)>);

static_assert(not IsFunctionPatch<FunctionArchetype<ApproximateTag,Real(Real)>,ApproximateTag,Real(Real)>);
static_assert(IsFunctionPatch<FunctionPatchArchetype<ApproximateTag,Real(Real)>,ApproximateTag,Real(Real)>);
static_assert(not IsFunctionPatch<EntireFunctionArchetype<ApproximateTag,Real(Real)>,ApproximateTag,Real(Real)>);

static_assert(not IsEntireFunction<FunctionArchetype<ApproximateTag,Real(Real)>,ApproximateTag,Real(Real)>);
static_assert(not IsEntireFunction<FunctionPatchArchetype<ApproximateTag,Real(Real)>,ApproximateTag,Real(Real)>);
static_assert(IsEntireFunction<EntireFunctionArchetype<ApproximateTag,Real(Real)>,ApproximateTag,Real(Real)>);
*/

} // namespace Ariadne


int main() {
    TestFunctionPatch().test();

    return ARIADNE_TEST_FAILURES;
}

