/***************************************************************************
 *            test_logical_sequence.cpp
 *
 *  Copyright  2026  Ariadne contributors
 *
 ****************************************************************************/
#include "foundation/logical.hpp"
#include "numeric/integer.hpp"
#include "numeric/sequence.hpp"
#include "numeric/logical_sequence.hpp"
#include "utility/test.hpp"
using namespace Ariadne;
class TestLogicalSequence { public: Void test(); };
Int main() {
    ARIADNE_TEST_CLASS(TestLogicalSequence,TestLogicalSequence());
    return ARIADNE_TEST_FAILURES;
}
Void TestLogicalSequence::test() {
    Sequence<LowerKleenean> seq([](Natural n){return n==2 ? LowerKleenean(true) : LowerKleenean(indeterminate);});
    ARIADNE_TEST_ASSIGN_CONSTRUCT(LowerKleenean, some, disjunction(seq));
    ARIADNE_TEST_ASSERT(possibly(not some.check(2_eff)));
    ARIADNE_TEST_ASSERT(definitely(some.check(3_eff)));
    ARIADNE_TEST_ASSERT(definitely(some.check(4_eff)));
    Sequence<UpperKleenean> upper_seq([](Natural n){return n==2 ? UpperKleenean(false) : UpperKleenean(indeterminate);});
    ARIADNE_TEST_ASSIGN_CONSTRUCT(UpperKleenean, all, conjunction(upper_seq));
    ARIADNE_TEST_ASSERT(possibly(all.check(2_eff)));
    ARIADNE_TEST_ASSERT(not possibly(all.check(3_eff)));
    ARIADNE_TEST_ASSERT(definitely(not all.check(4_eff)));
}
