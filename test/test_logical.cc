/***************************************************************************
 *            test_logical.cc
 *
 *  Copyright  2015  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "utility/metaprogramming.h"
#include "numeric/paradigm.h"
#include "numeric/logical.h"

#include "test.h"

using namespace std;
using namespace Ariadne;


class TestParadigm
{
  public:
    Void test();
  private:
    Void test_concept();
};

class TestLogical
{
  public:
    Void test();
  private:
    Void test_conversion_to_bool();
};


Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    ARIADNE_TEST_CLASS(TestParadigm,TestParadigm());
    ARIADNE_TEST_CLASS(TestLogical,TestLogical());

    return ARIADNE_TEST_FAILURES;
}


Void
TestParadigm::test()
{
    ARIADNE_TEST_CALL(test_concept());
}

// Test that the type implements all operations of
// the Float64 concept without testing correctness
Void
TestParadigm::test_concept()
{
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<Approximate,Lower>);
    ARIADNE_TEST_STATIC_ASSERT(IsStronger<Lower,Approximate>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsWeaker<EffectiveUpper,Bounded>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsStronger<EffectiveUpper,Bounded>>);

    ARIADNE_TEST_STATIC_ASSERT(IsSame<Upper,Weaker<EffectiveUpper,Bounded>>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<Upper,Weaker<PositiveEffectiveUpper,Bounded>>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<Approximate,ParadigmClass<ParadigmCode::APPROXIMATE>>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<Lower,ParadigmClass<ParadigmCode::LOWER>>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<Upper,ParadigmClass<ParadigmCode::UPPER>>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<Bounded,ParadigmClass<ParadigmCode::BOUNDED>>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<Metric,ParadigmClass<ParadigmCode::METRIC>>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<Validated,ParadigmClass<ParadigmCode::VALIDATED>>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<PositiveUpper,ParadigmClass<ParadigmCode::POSITIVE_UPPER>>);

    ARIADNE_TEST_STATIC_ASSERT(IsSame<ParadigmClass<strengthen(ParadigmCode::BOUNDED)>,Validated>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<ParadigmClass<strengthen(ParadigmCode::METRIC)>,Validated>);

    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<Validated,Validated>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<Validated,Bounded>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<Validated,Metric>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<Bounded,Validated>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<Bounded,Bounded>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<Bounded,Metric>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<Metric,Validated>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<Metric,Bounded>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<Metric,Metric>);

    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<Approximate,Validated>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<Validated,Effective>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<Effective,Exact>);

    ARIADNE_TEST_STATIC_ASSERT(IsSame<Weaken<Effective>,Validated>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<Weaken<Validated>,Approximate>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<Weaken<Approximate>,Void>);
}

Void
TestLogical::test()
{
    ARIADNE_TEST_CALL(test_conversion_to_bool());
}

Void
TestLogical::test_conversion_to_bool()
{
    ARIADNE_TEST_STATIC_ASSERT(IsConvertible<Logical<Exact>,Bool>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<Logical<Effective>,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<Logical<EffectiveUpper>,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<Logical<EffectiveLower>,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<Logical<Validated>,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<Logical<Upper>,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<Logical<Lower>,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<Logical<Approximate>,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(IsConvertible<Boolean,Bool>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<Tribool,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<Sierpinski,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<NegSierpinski,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<Fuzzy,Bool>>);
}
