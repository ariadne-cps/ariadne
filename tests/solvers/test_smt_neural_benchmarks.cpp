/***************************************************************************
 *            test_smt_neural_benchmarks.cpp
 *
 *  Deterministic neural-network regression and scaling benchmarks for the
 *  epsilon-SMT solver. Keep timing output isolated from the general SMT tests.
 ***************************************************************************/

#include <iostream>
#include <limits>

#include "utility/stopwatch.hpp"
#include "solvers/smt_solver.hpp"

#include "smt_barr3_prefix8.hpp"
#include "smt_barr3_prefix16.hpp"
#include "smt_barr3_prefix32.hpp"

#include "../test.hpp"

using namespace Ariadne;

class TestSmtNeuralBenchmarks {
  public:
    Void test() {
        ARIADNE_TEST_CALL(test_barr3_prefix8());
        ARIADNE_TEST_CALL(test_barr3_prefix16());
        ARIADNE_TEST_CALL(test_barr3_prefix32());
    }

  private:
    Void test_barr3_prefix32() {
        std::cout << std::endl
                  << "=== SMT neural benchmark: Barr3-derived 2-32-32-1 ==="
                  << std::endl;

        RealVariable bx("barr3_32_x"), by("barr3_32_y");
        RealExpression ebx=bx, eby=by;
        RealSpace space({bx,by});
        RealExpression barrier=TestBarr3Prefix32::network(ebx,eby);

        auto barrier_function=make_function(space,barrier);
        UpperBoxType origin=ExactBoxType({
            ExactIntervalType(0.0_x,0.0_x),
            ExactIntervalType(0.0_x,0.0_x)
        });
        UpperIntervalType origin_image=apply(barrier_function,origin);
        ARIADNE_TEST_ASSERT(definitely(subset(
            origin_image,
            ExactIntervalType(1.190_x,1.191_x))));

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
            false));

        Stopwatch<Milliseconds> stopwatch;
        SmtResult solve_result=solver.solve(space,barr3_domain,literals);
        stopwatch.click();

        std::cout << "[timing] solve=" << stopwatch.elapsed_seconds() << " s"
                  << " boxes=" << solve_result.statistics().boxes_processed
                  << " split=" << solve_result.statistics().boxes_split
                  << " sensitivity="
                  << solve_result.statistics().sensitivity_guided_splits
                  << std::endl;

        ARIADNE_TEST_ASSERT(solve_result.is_unknown());
        ARIADNE_TEST_EQUAL(
            solve_result.unknown_reason(),SmtUnknownReason::RESOURCE_EXHAUSTED);
        ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
        ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,1u);
        ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,1u);
        ARIADNE_TEST_ASSERT(
            solve_result.statistics().sensitivity_guided_splits>=1u);
    }

    Void test_barr3_prefix16() {
        std::cout << std::endl
                  << "=== SMT neural benchmark: Barr3-derived 2-16-16-1 ==="
                  << std::endl;

        RealVariable bx("barr3_16_x"), by("barr3_16_y");
        RealExpression ebx=bx, eby=by;
        RealSpace space({bx,by});
        RealExpression barrier=TestBarr3Prefix16::network(ebx,eby);

        auto barrier_function=make_function(space,barrier);
        UpperBoxType origin=ExactBoxType({
            ExactIntervalType(0.0_x,0.0_x),
            ExactIntervalType(0.0_x,0.0_x)
        });
        UpperIntervalType origin_image=apply(barrier_function,origin);
        ARIADNE_TEST_ASSERT(definitely(subset(
            origin_image,
            ExactIntervalType(-0.466_x,-0.465_x))));

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
            false));

        Stopwatch<Milliseconds> stopwatch;
        SmtResult solve_result=solver.solve(space,barr3_domain,literals);
        stopwatch.click();

        std::cout << "[timing] solve=" << stopwatch.elapsed_seconds() << " s"
                  << " boxes=" << solve_result.statistics().boxes_processed
                  << " split=" << solve_result.statistics().boxes_split
                  << " sensitivity="
                  << solve_result.statistics().sensitivity_guided_splits
                  << std::endl;

        ARIADNE_TEST_ASSERT(solve_result.is_unknown());
        ARIADNE_TEST_EQUAL(
            solve_result.unknown_reason(),SmtUnknownReason::RESOURCE_EXHAUSTED);
        ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
        ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,1u);
        ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,1u);
        ARIADNE_TEST_ASSERT(
            solve_result.statistics().sensitivity_guided_splits>=1u);
    }

    Void test_barr3_prefix8() {
        std::cout << std::endl
                  << "=== SMT neural benchmark: Barr3-derived 2-8-8-1 ==="
                  << std::endl;

        RealVariable bx("barr3_x"), by("barr3_y");
        RealExpression ebx=bx, eby=by;
        RealSpace space({bx,by});
        RealExpression barrier=TestBarr3Prefix8::network(ebx,eby);

        auto barrier_function=make_function(space,barrier);
        UpperBoxType origin=ExactBoxType({
            ExactIntervalType(0.0_x,0.0_x),
            ExactIntervalType(0.0_x,0.0_x)
        });
        UpperIntervalType origin_image=apply(barrier_function,origin);
        ARIADNE_TEST_ASSERT(definitely(subset(
            origin_image,
            ExactIntervalType(-0.275_x,-0.274_x))));

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
            false));

        Stopwatch<Milliseconds> stopwatch;
        SmtResult solve_result=solver.solve(space,barr3_domain,literals);
        stopwatch.click();

        std::cout << "[timing] solve=" << stopwatch.elapsed_seconds() << " s"
                  << " boxes=" << solve_result.statistics().boxes_processed
                  << " split=" << solve_result.statistics().boxes_split
                  << " sensitivity="
                  << solve_result.statistics().sensitivity_guided_splits
                  << std::endl;

        ARIADNE_TEST_ASSERT(solve_result.is_unknown());
        ARIADNE_TEST_EQUAL(
            solve_result.unknown_reason(),SmtUnknownReason::RESOURCE_EXHAUSTED);
        ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
        ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,1u);
        ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,1u);
        ARIADNE_TEST_ASSERT(
            solve_result.statistics().sensitivity_guided_splits>=1u);
    }
};

Int main() {
    TestSmtNeuralBenchmarks().test();
    return ARIADNE_TEST_FAILURES;
}
