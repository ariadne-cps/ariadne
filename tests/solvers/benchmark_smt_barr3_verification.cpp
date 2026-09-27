/***************************************************************************
 * benchmark_smt_barr3_verification.cpp
 *
 * Standalone reproduction of the published FOSSIL/dReal Barr3 BarrierAlt
 * counterexample queries using Ariadne's epsilon-SMT solver.
 ***************************************************************************/

#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "utility/stopwatch.hpp"
#include "solvers/smt_solver.hpp"

#include "smt_barr3_full64.hpp"

using namespace Ariadne;

namespace {

SizeType box_limit_from_argument(Int argc,const char* argv[]) {
    if(argc<=1) { return 1u; }
    String argument(argv[1]);
    if(argument=="full") {
        return std::numeric_limits<SizeType>::max();
    }
    char* end=nullptr;
    unsigned long long parsed=std::strtoull(argv[1],&end,10);
    if(end==argv[1] || *end!='\0' || parsed==0u) {
        throw std::runtime_error(
            "Usage: benchmark_smt_barr3_verification [positive-box-limit|full]");
    }
    return static_cast<SizeType>(parsed);
}

Bool sensitivity_from_argument(Int argc,const char* argv[]) {
    if(argc<=2) { return true; }
    String argument(argv[2]);
    if(argument=="sensitivity") { return true; }
    if(argument=="geometric") { return false; }
    throw std::runtime_error(
        "Usage: benchmark_smt_barr3_verification "
        "[positive-box-limit|full] [sensitivity|geometric]");
}

SmtTheoryPrimitiveLiteral primitive(ContinuousPredicate const& predicate) {
    auto alternatives=normalize_smt_theory_literal(
        make_smt_theory_literal(predicate));
    if(alternatives.size()!=1u || alternatives[0].size()!=1u) {
        throw std::runtime_error("Expected one primitive SMT literal");
    }
    return alternatives[0][0];
}

Void print_result(String const& name,SmtResult const& result,double seconds) {
    std::cout << "[" << name << "] status=" << result.status();
    if(result.is_unknown()) {
        std::cout << " reason=" << result.unknown_reason();
    }
    std::cout << " time=" << seconds << " s"
              << " boxes=" << result.statistics().boxes_processed
              << " pruned=" << result.statistics().boxes_pruned
              << " split=" << result.statistics().boxes_split
              << " epsilon-certified=" << result.statistics().epsilon_box_certifications
              << " candidate-searches=" << result.statistics().candidate_witness_searches
              << " candidate-successes=" << result.statistics().candidate_witness_successes
              << " hull-rounds=" << result.statistics().hull_reduction_rounds
              << " shaving-rounds=" << result.statistics().shaving_reduction_rounds
              << " shaving-evals=" << result.statistics().shaving_function_evaluations
              << " monotone-rounds=" << result.statistics().monotone_reduction_rounds
              << " sensitivity-splits=" << result.statistics().sensitivity_guided_splits
              << " sensitivity-overrides="
              << result.statistics().sensitivity_overrides_geometric_splits
              << " compile=" << result.statistics().theory_compile_seconds
              << " reduce=" << result.statistics().reduction_seconds
              << " epsilon=" << result.statistics().epsilon_check_seconds
              << " witness=" << result.statistics().witness_probe_seconds
              << " split-phase=" << result.statistics().split_seconds
              << " sensitivity-builds=" << result.statistics().sensitivity_derivatives_built
              << " sensitivity-evals=" << result.statistics().sensitivity_derivative_evaluations
              << " sensitivity-build-time="
              << result.statistics().sensitivity_derivative_build_seconds
              << " sensitivity-eval-time="
              << result.statistics().sensitivity_derivative_evaluation_seconds
              << " candidate=" << result.statistics().candidate_search_seconds
              << std::endl;
}

SmtResult timed_solve(
    String const& name,
    SmtSolver const& solver,
    RealSpace const& space,
    ExactBoxType const& domain,
    List<SmtTheoryPrimitiveLiteral> const& literals)
{
    std::cout << "[start] " << name
              << " literals=" << literals.size() << std::endl << std::flush;
    Stopwatch<Milliseconds> stopwatch;
    SmtResult result=solver.solve(space,domain,literals);
    stopwatch.click();
    print_result(name,result,stopwatch.elapsed_seconds());
    return result;
}

} // namespace

Int main(Int argc,const char* argv[]) {
    SizeType const box_limit=box_limit_from_argument(argc,argv);
    Bool const sensitivity_enabled=sensitivity_from_argument(argc,argv);

    std::cout << "=== Published Barr3 2-64-64-1 verification ===" << std::endl;
    std::cout << "epsilon=1e-5 box-limit=";
    if(box_limit==std::numeric_limits<SizeType>::max()) {
        std::cout << "unlimited";
    } else {
        std::cout << box_limit;
    }
    std::cout << " split-policy="
              << (sensitivity_enabled ? "sensitivity" : "geometric")
              << std::endl;

    RealVariable x("barr3_x"), y("barr3_y");
    RealExpression ex=x, ey=y;
    RealSpace space({x,y});
    ExactBoxType domain({
        ExactIntervalType(-3.0_x,2.5_x),
        ExactIntervalType(-2.0_x,1.0_x)
    });

    Stopwatch<Milliseconds> construction_stopwatch;
    auto network=TestBarr3Full64::network_and_lie(ex,ey);
    construction_stopwatch.click();
    std::cout << "[construction] network+Lie="
              << construction_stopwatch.elapsed_seconds() << " s" << std::endl;

    // For diagnosis and performance, use direct primitive conjunctions rather
    // than routing the unsafe union through DPLL(T). The decomposition is
    // logically exact: XU is the union of these three components, so the
    // unsafe obligation is UNSAT iff all three component queries are UNSAT.
    //
    // Rectangle membership is encoded directly in each query domain. The
    // spherical component uses its tight bounding box plus the sphere literal.
    ExactBoxType unsafe_sphere_domain({
        ExactIntervalType(-1.4_x,-0.6_x),
        ExactIntervalType(-1.4_x,-0.6_x)
    });
    ExactBoxType unsafe_rectangle_1_domain({
        ExactIntervalType(0.4_x,0.6_x),
        ExactIntervalType(0.1_x,0.5_x)
    });
    ExactBoxType unsafe_rectangle_2_domain({
        ExactIntervalType(0.4_x,0.8_x),
        ExactIntervalType(0.1_x,0.3_x)
    });

    RealExpression sphere_value=(ex+1)*(ex+1)+(ey+1)*(ey+1);
    SmtTheoryPrimitiveLiteral barrier_nonnegative=
        primitive(network.barrier>=0);
    SmtTheoryPrimitiveLiteral sphere_inside=
        primitive(sphere_value<=0.16_x);
    SmtTheoryPrimitiveLiteral lie_violation=
        primitive((network.lie+network.barrier)<0);

    // First isolate the cost of validated box processing. Candidate search and
    // derivative-assisted monotone contraction are deliberately disabled here:
    // the previous version enabled both and a nominal one-box run could spend
    // an unbounded amount of wall time inside the per-box nonlinear candidate
    // optimiser / derivative contractors before the box budget was observed.
    SmtSolver solver(SmtSolverConfiguration(
        0.00001_x,
        std::numeric_limits<SizeType>::max(),
        std::numeric_limits<SizeType>::max(),
        box_limit,
        false,
        false,
        sensitivity_enabled));

    List<SmtTheoryPrimitiveLiteral> sphere_literals({
        sphere_inside,barrier_nonnegative});
    List<SmtTheoryPrimitiveLiteral> barrier_literals({
        barrier_nonnegative});
    List<SmtTheoryPrimitiveLiteral> lie_literals({
        barrier_nonnegative,lie_violation});

    timed_solve(
        "unsafe-sphere",solver,space,unsafe_sphere_domain,sphere_literals);
    timed_solve(
        "unsafe-rectangle-1",solver,space,
        unsafe_rectangle_1_domain,barrier_literals);
    timed_solve(
        "unsafe-rectangle-2",solver,space,
        unsafe_rectangle_2_domain,barrier_literals);
    timed_solve(
        "lie",solver,space,domain,lie_literals);

    return 0;
}
