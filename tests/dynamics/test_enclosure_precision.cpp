/***************************************************************************
 *            test_enclosure_precision.cpp
 *
 *  Regression test for precision loss in RealBox enclosure construction.
 ****************************************************************************/

#include "config.hpp"
#include "algebra/vector.hpp"
#include "function/taylor_function.hpp"
#include "function/constraint.hpp"
#include "geometry/box.hpp"
#include "dynamics/enclosure.hpp"

#include "../test.hpp"

using namespace Ariadne;

Int main()
{
    RealBox box({{Decimal("1.25"),Decimal("1.55")},
                 {Decimal("2.35"),Decimal("2.45")}});

    TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp,1e-14));
    Enclosure enclosure(box,EnclosureConfiguration(function_factory));

    // The RealBox-to-Enclosure conversion must not quantise the endpoints to
    // Float32.  The old conversion left a remainder of about 9.5e-8 here.
    auto const& state_model =
        dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(
            enclosure.state_function().reference());
    for(SizeType i=0; i!=state_model.size(); ++i) {
        ARIADNE_TEST_COMPARE(state_model[i].error().raw(),<,1e-14_pr);
    }

    return ARIADNE_TEST_FAILURES;
}
