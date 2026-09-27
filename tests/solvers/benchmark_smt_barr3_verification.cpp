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
              << " monotone-rounds=" << result.statistics().monotone_reduction_rounds
              << " sensitivity-splits=" << result.statistics().sensitivity_guided_splits
              << std::endl;
}

} // namespace

Int main(Int argc,const char* argv[]) {
    SizeType const box_limit=box_limit_from_argument(argc,argv);

    std::cout << "=== Published Barr3 2-64-64-1 verification ===" << std::endl;
    std::cout << "epsilon=1e-5 box-limit=";
    if(box_limit==std::numeric_limits<SizeType>::max()) {
        std::cout << "unlimited";
    } else {
        std::cout << box_limit;
    }
    std::cout << std::endl;

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

    RealExpression sphere_value=(ex+1)*(ex+1)+(ey+1)*(ey+1);
    ContinuousPredicate unsafe_sphere=(sphere_value<=0.16_x);
    ContinuousPredicate unsafe_rectangle_1=
        (ex>=0.4_x)&&(ex<=0.6_x)&&(ey>=0.1_x)&&(ey<=0.5_x);
    ContinuousPredicate unsafe_rectangle_2=
        (ex>=0.4_x)&&(ex<=0.8_x)&&(ey>=0.1_x)&&(ey<=0.3_x);
    ContinuousPredicate unsafe=
        unsafe_sphere||unsafe_rectangle_1||unsafe_rectangle_2;

    // These are exactly the two negated BarrierAlt conditions in the public
    // Barr3 verifier. An UNSAT answer means the corresponding certificate
    // obligation is established; EPSILON_SAT is a candidate counterexample.
    ContinuousPredicate unsafe_counterexample=
        unsafe&&(network.barrier>=0);
    ContinuousPredicate lie_counterexample=
        (network.barrier>=0)&&((network.lie+network.barrier)<0);

    SmtSolver solver(SmtSolverConfiguration(
        0.00001_x,
        std::numeric_limits<SizeType>::max(),
        std::numeric_limits<SizeType>::max(),
        box_limit,
        true,
        true));

    Stopwatch<Milliseconds> unsafe_stopwatch;
    SmtResult unsafe_result=solver.solve(space,domain,unsafe_counterexample);
    unsafe_stopwatch.click();
    print_result("unsafe",unsafe_result,unsafe_stopwatch.elapsed_seconds());

    Stopwatch<Milliseconds> lie_stopwatch;
    SmtResult lie_result=solver.solve(space,domain,lie_counterexample);
    lie_stopwatch.click();
    print_result("lie",lie_result,lie_stopwatch.elapsed_seconds());

    return 0;
}
