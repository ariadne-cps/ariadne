/***************************************************************************
 *            test_smt_neural_benchmarks.cpp
 *
 *  Fast regression smoke test for the published Barr3 neural barrier.
 *  Performance and scaling measurements live in benchmark_smt_barr3_verification.
 ***************************************************************************/

#include <iostream>
#include <limits>

#include "solvers/smt_solver.hpp"

#include "smt_barr3_full64.hpp"

#include "utility/test.hpp"

using namespace Ariadne;

class TestSmtNeuralBenchmarks {
  public:
    Void test() {
        ARIADNE_TEST_CALL(test_barr3_full64_smoke());
    }

  private:
    Void test_barr3_full64_smoke() {
        std::cout << std::endl
                  << "=== SMT neural regression: published Barr3 2-64-64-1 ==="
                  << std::endl;

        RealVariable bx("barr3_64_x"), by("barr3_64_y");
        RealExpression ebx=bx, eby=by;
        RealSpace space({bx,by});
        RealExpression barrier=TestBarr3Full64::network(ebx,eby);

        auto barrier_function=make_function(space,barrier);
        UpperBoxType origin=ExactBoxType({
            ExactIntervalType(0.0_x,0.0_x),
            ExactIntervalType(0.0_x,0.0_x)
        });
        UpperIntervalType origin_image=apply(barrier_function,origin);
        ARIADNE_TEST_ASSERT(definitely(subset(
            origin_image,
            ExactIntervalType(3.809_x,3.810_x))));

        auto alternatives=normalize_smt_theory_literal(
            make_smt_theory_literal(barrier==0));
        ARIADNE_TEST_EQUAL(alternatives.size(),1u);
        ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);

        List<SmtTheoryPrimitiveLiteral> literals({alternatives[0][0]});
        ExactBoxType barr3_domain({
            ExactIntervalType(-3.0_x,2.5_x),
            ExactIntervalType(-2.0_x,1.0_x)
        });
        SmtSolver solver(SmtSolverConfiguration(
            0.01_x,
            std::numeric_limits<SizeType>::max(),
            std::numeric_limits<SizeType>::max(),
            1u,
            false,
            false,
            false,
            false,
            false,
            false));

        SmtResult solve_result=solver.solve(space,barr3_domain,literals);

        ARIADNE_TEST_ASSERT(solve_result.is_unknown());
        ARIADNE_TEST_EQUAL(
            solve_result.unknown_reason(),SmtUnknownReason::RESOURCE_EXHAUSTED);
        ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
        ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,1u);
        ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,1u);
        ARIADNE_TEST_EQUAL(
            solve_result.statistics().fused_direct_classification_boxes,1u);
        ARIADNE_TEST_EQUAL(
            solve_result.statistics().sensitivity_guided_splits,0u);
    }

};

Int main() {
    TestSmtNeuralBenchmarks().test();
    return ARIADNE_TEST_FAILURES;
}
