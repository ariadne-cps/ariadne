/***************************************************************************
 *            test_logical.cpp
 *
 *  Copyright  2015-20  Pieter Collins
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

#include "utility/metaprogramming.hpp"
#include "numeric/paradigm.hpp"
#include "numeric/logical.hpp"
#include "numeric/integer.hpp"
#include "numeric/sequence.hpp"

#include "../test.hpp"

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
    Void test_concept();
    Void test_conversion_to_bool();
    Void test_conversion();
    Void test_disjunction();
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
}

// Test that the type implements all operations of
// the FloatDP concept without testing correctness
Void
TestParadigm::test_concept()
{
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<ApproximateTag,LowerTag>);
    ARIADNE_TEST_STATIC_ASSERT(IsStronger<LowerTag,ApproximateTag>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsWeaker<EffectiveUpperTag,BoundedTag>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsStronger<EffectiveUpperTag,BoundedTag>>);

    ARIADNE_TEST_STATIC_ASSERT(IsSame<UpperTag,Weaker<EffectiveUpperTag,BoundedTag>>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<UpperTag,Weaker<PositiveEffectiveUpperTag,BoundedTag>>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<ApproximateTag,ParadigmClass<ParadigmCode::APPROXIMATE>>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<LowerTag,ParadigmClass<ParadigmCode::LOWER>>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<UpperTag,ParadigmClass<ParadigmCode::UPPER>>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<BoundedTag,ParadigmClass<ParadigmCode::BOUNDED>>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<MetricTag,ParadigmClass<ParadigmCode::METRIC>>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<ValidatedTag,ParadigmClass<ParadigmCode::VALIDATED>>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<PositiveUpperTag,ParadigmClass<ParadigmCode::POSITIVE_UPPER>>);

    ARIADNE_TEST_STATIC_ASSERT(IsSame<ParadigmClass<strengthen(ParadigmCode::BOUNDED)>,ValidatedTag>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<ParadigmClass<strengthen(ParadigmCode::METRIC)>,ValidatedTag>);

    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<ValidatedTag,ValidatedTag>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<ValidatedTag,BoundedTag>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<ValidatedTag,MetricTag>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<BoundedTag,ValidatedTag>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<BoundedTag,BoundedTag>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<BoundedTag,MetricTag>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<MetricTag,ValidatedTag>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<MetricTag,BoundedTag>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<MetricTag,MetricTag>);

    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<ApproximateTag,ValidatedTag>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<ValidatedTag,EffectiveTag>);
    ARIADNE_TEST_STATIC_ASSERT(IsWeaker<EffectiveTag,ExactTag>);

    ARIADNE_TEST_STATIC_ASSERT(IsSame<Weaken<EffectiveTag>,ValidatedTag>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<Weaken<ValidatedTag>,ApproximateTag>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<Weaken<ApproximateTag>,Void>);
}

Void
TestLogical::test()
{
    ARIADNE_TEST_CALL(test_conversion_to_bool());
    ARIADNE_TEST_CALL(test_conversion());
    ARIADNE_TEST_CALL(test_disjunction());
}

Void
TestLogical::test_concept()
{
    // Check to see if we can perform operations on computational and specification logical types
    LogicalType<ExactTag> xl(true);
    LogicalType<EffectiveTag> el(true);
    LogicalType<ValidatedTag> vl(true);
    Effort eff(0);

    vl=el.check(eff);
    vl=check(el,eff);

    vl=indeterminate;
    vl=LogicalType<ValidatedTag>(LogicalValue::LIKELY);

    xl = xl && xl;
    el = xl && el;
    vl = xl && vl;
    el = el && el;
    vl = vl && vl;
}

Void
TestLogical::test_conversion_to_bool()
{
    ARIADNE_TEST_STATIC_ASSERT(IsConvertible<LogicalType<ExactTag>,Bool>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<LogicalType<EffectiveTag>,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<LogicalType<EffectiveUpperTag>,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<LogicalType<EffectiveLowerTag>,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<LogicalType<ValidatedTag>,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<LogicalType<UpperTag>,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<LogicalType<LowerTag>,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<LogicalType<ApproximateTag>,Bool>>);

    ARIADNE_TEST_STATIC_ASSERT(IsConvertible<Boolean,Bool>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<ValidatedKleenean,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<ValidatedSierpinskian,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<ValidatedNegatedSierpinskian,Bool>>);
    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<Fuzzy,Bool>>);
}

Void
TestLogical::test_conversion()
{
    if(IsConvertible<LogicalType<EffectiveTag>,LogicalType<ValidatedTag>>::value) {
        ARIADNE_TEST_NOTIFY("EffectiveTag logical types may be converted to values using default Effort.");
    } else if(IsConvertible<LogicalType<EffectiveTag>,LogicalType<ValidatedTag>>::value) {
        ARIADNE_TEST_NOTIFY("EffectiveTag logical types may be explicitly converted to values using default Effort.");
    } else {
        ARIADNE_TEST_NOTIFY("EffectiveTag logical types cannot be converted to values; the Effort used must be specified.");
    }

    try {
        if(decide(indeterminate)) {
            ARIADNE_TEST_NOTIFY("decide(...) is true on INDETERMINATE value.");
        } else {
            ARIADNE_TEST_NOTIFY("decide(...) is false on INDETERMINATE value.");
        }
    } catch(...) {
        ARIADNE_TEST_NOTIFY("decide(...) is throws error on INDETERMINATE value.");
    }

    ARIADNE_TEST_STATIC_ASSERT(Not<IsConvertible<Indeterminate,Boolean>>);
    ARIADNE_TEST_STATIC_ASSERT(IsConvertible<Indeterminate,Sierpinskian>);
    ARIADNE_TEST_STATIC_ASSERT(IsConvertible<Indeterminate,NegatedSierpinskian>);
    ARIADNE_TEST_STATIC_ASSERT(IsConvertible<Indeterminate,Kleenean>);
    ARIADNE_TEST_STATIC_ASSERT(IsConvertible<Indeterminate,LowerKleenean>);
    ARIADNE_TEST_STATIC_ASSERT(IsConvertible<Indeterminate,UpperKleenean>);
//    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(indeterminate and true),Kleenean>);
//    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(indeterminate and Boolean(true)),Kleenean>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(indeterminate and Sierpinskian(true)),Sierpinskian>);
    ARIADNE_TEST_STATIC_ASSERT(IsSame<decltype(indeterminate and Kleenean(true)),Kleenean>);

    ARIADNE_TEST_CONSTRUCT(LogicalType<ValidatedTag>,vl,(LogicalValue::LIKELY))
    ARIADNE_TEST_EQUAL(definitely(vl),false);
    ARIADNE_TEST_EQUAL(possibly(vl),true);
    ARIADNE_TEST_EQUAL(decide(vl),true);

    ARIADNE_TEST_CONSTRUCT(LogicalType<ValidatedTag>,vi,(LogicalValue::INDETERMINATE))
    ARIADNE_TEST_EQUAL(definitely(vl),false);
    ARIADNE_TEST_EQUAL(possibly(vl),true);
}

Void
TestLogical::test_disjunction()
{
    Sequence<LowerKleenean> seq([](Natural n){return n==2 ? LowerKleenean(true) : LowerKleenean(indeterminate);});
    ARIADNE_TEST_ASSIGN_CONSTRUCT(LowerKleenean, some, disjunction(seq));
    ARIADNE_TEST_ASSERT(possibly(not some.check(2_eff)));
    ARIADNE_TEST_ASSERT(definitely(some.check(3_eff)));
    ARIADNE_TEST_ASSERT(definitely(some.check(4_eff)));

    ARIADNE_TEST_ASSIGN_CONSTRUCT(
        UpperKleenean, all, conjunction(Sequence<UpperKleenean>([](Natural n){return n==2 ? UpperKleenean(false) : UpperKleenean(indeterminate);})));
    ARIADNE_TEST_ASSERT(possibly(all.check(2_eff)));
    ARIADNE_TEST_ASSERT(not possibly(all.check(3_eff)));
    ARIADNE_TEST_ASSERT(definitely(not all.check(4_eff)));

}

