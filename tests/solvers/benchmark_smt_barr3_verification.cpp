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
#include "function/procedure.hpp"
#include "function/procedure.tpl.hpp"
#include "function/taylor_function.hpp"
#include "function/affine_model.hpp"
#include "symbolic/expression.hpp"
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

String query_from_argument(Int argc,const char* argv[]) {
    if(argc<=6) { return "all"; }
    String argument(argv[6]);
    if(argument=="all" || argument=="lie" || argument=="lie-only" || argument=="eval" || argument=="taylor" || argument=="affine") { return argument; }
    throw std::runtime_error(
        "Usage: benchmark_smt_barr3_verification "
        "[positive-box-limit|full] [sensitivity|geometric|lookahead] "
        "[witness|no-witness] [shaving|no-shaving] [hull|no-hull] "
        "[all|lie|lie-only|eval|taylor|affine]");
}

Bool monotone_from_argument(Int argc,const char* argv[]) {
    if(argc<=7) { return false; }
    String argument(argv[7]);
    if(argument=="monotone") { return true; }
    if(argument=="no-monotone") { return false; }
    throw std::runtime_error(
        "Usage: benchmark_smt_barr3_verification "
        "[positive-box-limit|full] [sensitivity|geometric|lookahead] "
        "[witness|no-witness] [shaving|no-shaving] [hull|no-hull] "
        "[all|lie|lie-only|eval|taylor|affine] [monotone|no-monotone]");
}

String lie_literal_order_from_argument(Int argc,const char* argv[]) {
    if(argc<=8) { return "barrier-first"; }
    String argument(argv[8]);
    if(argument=="barrier-first" || argument=="lie-first") { return argument; }
    throw std::runtime_error(
        "Usage: benchmark_smt_barr3_verification "
        "[positive-box-limit|full] [sensitivity|geometric|lookahead] "
        "[witness|no-witness] [shaving|no-shaving] [hull|no-hull] "
        "[all|lie|lie-only|eval|taylor|affine] [monotone|no-monotone] "
        "[barrier-first|lie-first]");
}

String child_order_from_argument(Int argc,const char* argv[]) {
    if(argc<=9) { return "lower-first"; }
    String argument(argv[9]);
    if(argument=="lower-first" || argument=="upper-first") { return argument; }
    throw std::runtime_error(
        "Usage: benchmark_smt_barr3_verification "
        "[positive-box-limit|full] [sensitivity|geometric|lookahead] "
        "[witness|no-witness] [shaving|no-shaving] [hull|no-hull] "
        "[all|lie|lie-only|eval|taylor|affine] [monotone|no-monotone] "
        "[barrier-first|lie-first] [lower-first|upper-first]");
}

Bool hull_from_argument(Int argc,const char* argv[]) {
    if(argc<=5) { return true; }
    String argument(argv[5]);
    if(argument=="hull") { return true; }
    if(argument=="no-hull") { return false; }
    throw std::runtime_error(
        "Usage: benchmark_smt_barr3_verification "
        "[positive-box-limit|full] [sensitivity|geometric|lookahead] "
        "[witness|no-witness] [shaving|no-shaving] [hull|no-hull]");
}

Bool shaving_from_argument(Int argc,const char* argv[]) {
    if(argc<=4) { return true; }
    String argument(argv[4]);
    if(argument=="shaving") { return true; }
    if(argument=="no-shaving") { return false; }
    throw std::runtime_error(
        "Usage: benchmark_smt_barr3_verification "
        "[positive-box-limit|full] [sensitivity|geometric|lookahead] "
        "[witness|no-witness] [shaving|no-shaving]");
}

Bool witness_probing_from_argument(Int argc,const char* argv[]) {
    if(argc<=3) { return true; }
    String argument(argv[3]);
    if(argument=="witness") { return true; }
    if(argument=="no-witness") { return false; }
    throw std::runtime_error(
        "Usage: benchmark_smt_barr3_verification "
        "[positive-box-limit|full] [sensitivity|geometric|lookahead] "
        "[witness|no-witness]");
}

String split_policy_from_argument(Int argc,const char* argv[]) {
    if(argc<=2) { return "sensitivity"; }
    String argument(argv[2]);
    if(argument=="sensitivity" || argument=="geometric" || argument=="lookahead") {
        return argument;
    }
    throw std::runtime_error(
        "Usage: benchmark_smt_barr3_verification "
        "[positive-box-limit|full] [sensitivity|geometric|lookahead]");
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
              << " fused-direct=" << result.statistics().fused_direct_classification_boxes
              << " fused-literal-evals="
              << result.statistics().fused_direct_literal_evaluations
              << " candidate-searches=" << result.statistics().candidate_witness_searches
              << " candidate-successes=" << result.statistics().candidate_witness_successes
              << " hull-rounds=" << result.statistics().hull_reduction_rounds
              << " hull-effective=" << result.statistics().hull_effective_reductions
              << " hull-procedure-builds=" << result.statistics().hull_procedure_builds
              << " hull-procedure-build-time="
              << result.statistics().hull_procedure_build_seconds
              << " hull-contract-time=" << result.statistics().hull_contraction_seconds
              << " hull-reject-time="
              << result.statistics().hull_direct_rejection_seconds
              << " hull-temp-time="
              << result.statistics().hull_temporary_allocation_seconds
              << " hull-forward-time="
              << result.statistics().hull_forward_execution_seconds
              << " hull-backward-time="
              << result.statistics().hull_backward_propagation_seconds
              << " shaving-rounds=" << result.statistics().shaving_reduction_rounds
              << " shaving-evals=" << result.statistics().shaving_function_evaluations
              << " monotone-rounds=" << result.statistics().monotone_reduction_rounds
              << " monotone-effective=" << result.statistics().monotone_effective_reductions
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
              << " lookahead-splits="
              << result.statistics().interval_lookahead_guided_splits
              << " lookahead-overrides="
              << result.statistics().interval_lookahead_overrides_geometric_splits
              << " lookahead-evals="
              << result.statistics().interval_lookahead_function_evaluations
              << " lookahead-eval-time="
              << result.statistics().interval_lookahead_evaluation_seconds
              << " candidate=" << result.statistics().candidate_search_seconds
              << std::endl;
}

FloatDPBounds stable_tanh(FloatDPBounds const& x)
{
    return tanh(x);
}

struct DirectBarr3Evaluation {
    FloatDPBounds barrier;
    FloatDPBounds db_dx;
    FloatDPBounds db_dy;
    FloatDPBounds lie;
    FloatDPBounds lie_plus_barrier;
};

DirectBarr3Evaluation direct_barr3_evaluate(Vector<FloatDPBounds> const& x)
{
    auto const& p=TestBarr3Full64::parameters();
    constexpr SizeType width=TestBarr3Full64::width;

    auto c=[&](SizeType i) {
        return FloatDPBounds(ExactDouble(p[i]),dp);
    };

    Vector<FloatDPBounds> h1(width,FloatDPBounds(0,dp));
    Vector<FloatDPBounds> dh1_dx(width,FloatDPBounds(0,dp));
    Vector<FloatDPBounds> dh1_dy(width,FloatDPBounds(0,dp));
    for(SizeType i=0u;i!=width;++i) {
        FloatDPBounds z=c(TestBarr3Full64::b1_offset+i)
            + c(TestBarr3Full64::w1_offset+2u*i)*x[0]
            + c(TestBarr3Full64::w1_offset+2u*i+1u)*x[1];
        h1[i]=stable_tanh(z);
        FloatDPBounds factor=FloatDPBounds(1,dp)-sqr(h1[i]);
        dh1_dx[i]=factor*c(TestBarr3Full64::w1_offset+2u*i);
        dh1_dy[i]=factor*c(TestBarr3Full64::w1_offset+2u*i+1u);
    }

    Vector<FloatDPBounds> h2(width,FloatDPBounds(0,dp));
    Vector<FloatDPBounds> dh2_dx(width,FloatDPBounds(0,dp));
    Vector<FloatDPBounds> dh2_dy(width,FloatDPBounds(0,dp));
    for(SizeType i=0u;i!=width;++i) {
        FloatDPBounds z=c(TestBarr3Full64::b2_offset+i);
        FloatDPBounds dz_dx(0,dp);
        FloatDPBounds dz_dy(0,dp);
        for(SizeType j=0u;j!=width;++j) {
            FloatDPBounds weight=c(TestBarr3Full64::w2_offset+width*i+j);
            z+=weight*h1[j];
            dz_dx+=weight*dh1_dx[j];
            dz_dy+=weight*dh1_dy[j];
        }
        h2[i]=stable_tanh(z);
        FloatDPBounds factor=FloatDPBounds(1,dp)-sqr(h2[i]);
        dh2_dx[i]=factor*dz_dx;
        dh2_dy[i]=factor*dz_dy;
    }

    FloatDPBounds barrier=c(TestBarr3Full64::b3_offset);
    FloatDPBounds db_dx(0,dp);
    FloatDPBounds db_dy(0,dp);
    for(SizeType i=0u;i!=width;++i) {
        FloatDPBounds weight=c(TestBarr3Full64::w3_offset+i);
        barrier+=weight*h2[i];
        db_dx+=weight*dh2_dx[i];
        db_dy+=weight*dh2_dy[i];
    }

    FloatDPBounds dx=x[1];
    FloatDPBounds dy=-x[0]-x[1]+(x[0]*x[0]*x[0])/FloatDPBounds(3,dp);
    FloatDPBounds lie=db_dx*dx+db_dy*dy;
    return {barrier,db_dx,db_dy,lie,lie+barrier};
}

Void profile_direct_barr3(UpperBoxType const& domain, SizeType repetitions)
{
    Vector<FloatDPBounds> x(
        domain.size(),FloatDPBounds(DoublePrecision()));
    for(SizeType i=0u;i!=domain.size();++i) {
        x[i]=FloatDPBounds(
            domain[i].lower_bound().raw(),
            domain[i].upper_bound().raw());
    }

    DirectBarr3Evaluation result=direct_barr3_evaluate(x);
    Stopwatch<Milliseconds> stopwatch;
    for(SizeType i=0u;i!=repetitions;++i) {
        result=direct_barr3_evaluate(x);
    }
    stopwatch.click();

    std::cout << "[direct-profile] repetitions=" << repetitions
              << " evaluate=" << stopwatch.elapsed_seconds()
              << " barrier=" << result.barrier
              << " db/dx=" << result.db_dx
              << " db/dy=" << result.db_dy
              << " lie=" << result.lie
              << " lie+barrier=" << result.lie_plus_barrier
              << std::endl;
}

Void profile_formula_evaluator(
    String const& name,
    RealExpression const& expression,
    RealSpace const& space,
    UpperBoxType const& domain,
    SizeType repetitions)
{
    std::cout << "[formula-profile] " << name
              << " expression-nodes=" << count_nodes(expression)
              << " distinct-node-pointers="
              << count_distinct_node_pointers(expression)
              << std::endl;

    Stopwatch<Milliseconds> build_stopwatch;
    Formula<EffectiveNumber> formula=make_formula(expression,space);
    build_stopwatch.click();

    Stopwatch<Milliseconds> procedure_build_stopwatch;
    EffectiveProcedure expression_procedure(space.dimension(),formula);
    procedure_build_stopwatch.click();

    Vector<UpperIntervalType> upper_arguments=cast_vector(domain);
    UpperIntervalType upper_image=evaluate(formula,upper_arguments);
    Stopwatch<Milliseconds> upper_stopwatch;
    for(SizeType i=0u; i!=repetitions; ++i) {
        upper_image=evaluate(formula,upper_arguments);
    }
    upper_stopwatch.click();

    Vector<FloatDPBounds> bounds_arguments(
        domain.size(),FloatDPBounds(DoublePrecision()));
    for(SizeType i=0u; i!=domain.size(); ++i) {
        bounds_arguments[i]=FloatDPBounds(
            domain[i].lower_bound().raw(),
            domain[i].upper_bound().raw());
    }
    FloatDPBounds bounds_image=evaluate(formula,bounds_arguments);
    Stopwatch<Milliseconds> bounds_stopwatch;
    for(SizeType i=0u; i!=repetitions; ++i) {
        bounds_image=evaluate(formula,bounds_arguments);
    }
    bounds_stopwatch.click();

    FloatDPBounds expression_procedure_image=
        evaluate(expression_procedure,bounds_arguments);
    Stopwatch<Milliseconds> expression_procedure_stopwatch;
    for(SizeType i=0u; i!=repetitions; ++i) {
        expression_procedure_image=evaluate(
            expression_procedure,bounds_arguments);
    }
    expression_procedure_stopwatch.click();

    std::cout << "[formula-profile] " << name
              << " repetitions=" << repetitions
              << " formula-build=" << build_stopwatch.elapsed_seconds()
              << " expression-procedure-build="
              << procedure_build_stopwatch.elapsed_seconds()
              << " expression-procedure-instructions="
              << expression_procedure._instructions.size()
              << " formula-upper-evaluate=" << upper_stopwatch.elapsed_seconds()
              << " formula-bounds-evaluate=" << bounds_stopwatch.elapsed_seconds()
              << " expression-procedure-bounds-evaluate="
              << expression_procedure_stopwatch.elapsed_seconds()
              << " upper-image=" << upper_image
              << " bounds-image=" << bounds_image
              << " expression-procedure-image=" << expression_procedure_image
              << std::endl;
}

Void profile_expression_bounds(
    String const& name,
    RealExpression const& expression,
    RealSpace const& space,
    UpperBoxType const& domain)
{
    Formula<EffectiveNumber> formula=make_formula(expression,space);
    EffectiveProcedure procedure(space.dimension(),formula);
    Vector<FloatDPBounds> arguments(
        domain.size(),FloatDPBounds(DoublePrecision()));
    for(SizeType i=0u;i!=domain.size();++i) {
        arguments[i]=FloatDPBounds(
            domain[i].lower_bound().raw(),
            domain[i].upper_bound().raw());
    }
    FloatDPBounds image=evaluate(procedure,arguments);
    std::cout << "[bounds-diagnostic] " << name
              << " instructions=" << procedure._instructions.size()
              << " image=" << image
              << std::endl;
}

Void profile_evaluator(
    String const& name,
    ValidatedScalarMultivariateFunction const& function,
    UpperBoxType const& domain,
    SizeType repetitions)
{
    UpperIntervalType function_image;
    Stopwatch<Milliseconds> function_stopwatch;
    for(SizeType i=0u; i!=repetitions; ++i) {
        function_image=apply(function,domain);
    }
    function_stopwatch.click();

    Stopwatch<Milliseconds> build_stopwatch;
    ValidatedProcedure procedure(function);
    build_stopwatch.click();

    UpperIntervalType procedure_image;
    Vector<UpperIntervalType> arguments=cast_vector(domain);
    Stopwatch<Milliseconds> procedure_stopwatch;
    for(SizeType i=0u; i!=repetitions; ++i) {
        procedure_image=evaluate(procedure,arguments);
    }
    procedure_stopwatch.click();

    Vector<FloatDPBounds> bounds_arguments(
        domain.size(),FloatDPBounds(DoublePrecision()));
    for(SizeType i=0u; i!=domain.size(); ++i) {
        bounds_arguments[i]=FloatDPBounds(
            domain[i].lower_bound().raw(),
            domain[i].upper_bound().raw());
    }
    FloatDPBounds procedure_bounds_image=evaluate(procedure,bounds_arguments);
    Stopwatch<Milliseconds> procedure_bounds_stopwatch;
    for(SizeType i=0u; i!=repetitions; ++i) {
        procedure_bounds_image=evaluate(procedure,bounds_arguments);
    }
    procedure_bounds_stopwatch.click();

    std::cout << "[eval-profile] " << name
              << " repetitions=" << repetitions
              << " function-apply=" << function_stopwatch.elapsed_seconds()
              << " procedure-build=" << build_stopwatch.elapsed_seconds()
              << " procedure-instructions=" << procedure._instructions.size()
              << " procedure-upper-evaluate=" << procedure_stopwatch.elapsed_seconds()
              << " procedure-bounds-evaluate="
              << procedure_bounds_stopwatch.elapsed_seconds()
              << " function-image=" << function_image
              << " procedure-image=" << procedure_image
              << " procedure-bounds-image=" << procedure_bounds_image
              << std::endl;
}

Void profile_taylor_range(
    String const& name,
    ValidatedScalarMultivariateFunction const& function,
    ExactBoxType const& domain)
{
    UpperBoxType upper_domain(domain);
    Stopwatch<Milliseconds> interval_stopwatch;
    UpperIntervalType interval_image=apply(function,upper_domain);
    interval_stopwatch.click();

    ThresholdSweeper<FloatDP> sweeper(dp,1e-8);
    Stopwatch<Milliseconds> build_stopwatch;
    ValidatedScalarMultivariateTaylorFunctionModelDP model(
        domain,function,sweeper);
    build_stopwatch.click();

    Stopwatch<Milliseconds> range_stopwatch;
    auto model_range=model.range();
    range_stopwatch.click();

    std::cout << "[taylor-profile] " << name
              << " interval-time=" << interval_stopwatch.elapsed_seconds()
              << " interval-image=" << interval_image
              << " model-build=" << build_stopwatch.elapsed_seconds()
              << " model-range-time=" << range_stopwatch.elapsed_seconds()
              << " model-range=" << model_range
              << std::endl;
}

Void profile_affine_range(
    String const& name,
    ValidatedScalarMultivariateFunction const& function,
    ExactBoxType const& domain)
{
    Stopwatch<Milliseconds> build_stopwatch;
    auto model=affine_model(domain,function,dp);
    build_stopwatch.click();

    Stopwatch<Milliseconds> range_stopwatch;
    auto model_range=model.range();
    range_stopwatch.click();

    std::cout << "[affine-profile] " << name
              << " model-build=" << build_stopwatch.elapsed_seconds()
              << " model-range-time=" << range_stopwatch.elapsed_seconds()
              << " model-range=" << model_range
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
    String const split_policy=split_policy_from_argument(argc,argv);
    Bool const sensitivity_enabled=split_policy=="sensitivity";
    Bool const lookahead_enabled=split_policy=="lookahead";
    Bool const witness_probing_enabled=witness_probing_from_argument(argc,argv);
    Bool const shaving_enabled=shaving_from_argument(argc,argv);
    Bool const hull_enabled=hull_from_argument(argc,argv);
    String const query=query_from_argument(argc,argv);
    Bool const monotone_enabled=monotone_from_argument(argc,argv);
    String const lie_literal_order=lie_literal_order_from_argument(argc,argv);
    String const child_order=child_order_from_argument(argc,argv);
    Bool const upper_child_first=child_order=="upper-first";

    std::cout << "=== Published Barr3 2-64-64-1 verification ===" << std::endl;
    std::cout << "epsilon=1e-5 box-limit=";
    if(box_limit==std::numeric_limits<SizeType>::max()) {
        std::cout << "unlimited";
    } else {
        std::cout << box_limit;
    }
    std::cout << " split-policy=" << split_policy
              << " witness-probing="
              << (witness_probing_enabled ? "enabled" : "disabled")
              << " shaving="
              << (shaving_enabled ? "enabled" : "disabled")
              << " hull="
              << (hull_enabled ? "enabled" : "disabled")
              << " query=" << query
              << " monotone="
              << (monotone_enabled ? "enabled" : "disabled")
              << " lie-literal-order=" << lie_literal_order
              << " child-order=" << child_order
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

    if(query=="eval") {
        ValidatedScalarMultivariateFunction barrier_function=
            make_function(space,network.barrier);
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,network.lie+network.barrier);
        profile_evaluator("barrier",barrier_function,UpperBoxType(domain),8u);
        profile_formula_evaluator(
            "barrier",network.barrier,space,UpperBoxType(domain),8u);
        profile_evaluator("lie+barrier",lie_function,UpperBoxType(domain),8u);
        profile_formula_evaluator(
            "lie+barrier",network.lie+network.barrier,
            space,UpperBoxType(domain),8u);
        profile_expression_bounds("db/dx",network.db_dx,space,UpperBoxType(domain));
        profile_expression_bounds("db/dy",network.db_dy,space,UpperBoxType(domain));
        profile_expression_bounds("lie",network.lie,space,UpperBoxType(domain));
        profile_direct_barr3(UpperBoxType(domain),8u);
        return 0;
    }

    if(query=="taylor") {
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,network.lie+network.barrier);
        profile_taylor_range(
            "lie+barrier",lie_function,domain);
        return 0;
    }

    if(query=="affine") {
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,network.lie+network.barrier);
        profile_affine_range(
            "lie+barrier",lie_function,domain);
        return 0;
    }

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
        monotone_enabled,
        sensitivity_enabled,
        witness_probing_enabled,
        shaving_enabled,
        hull_enabled,
        lookahead_enabled,
        upper_child_first));

    List<SmtTheoryPrimitiveLiteral> sphere_literals({
        sphere_inside,barrier_nonnegative});
    List<SmtTheoryPrimitiveLiteral> barrier_literals({
        barrier_nonnegative});
    List<SmtTheoryPrimitiveLiteral> lie_literals;
    if(lie_literal_order=="lie-first") {
        lie_literals.append(lie_violation);
        lie_literals.append(barrier_nonnegative);
    } else {
        lie_literals.append(barrier_nonnegative);
        lie_literals.append(lie_violation);
    }

    if(query=="all") {
        timed_solve(
            "unsafe-sphere",solver,space,unsafe_sphere_domain,sphere_literals);
        timed_solve(
            "unsafe-rectangle-1",solver,space,
            unsafe_rectangle_1_domain,barrier_literals);
        timed_solve(
            "unsafe-rectangle-2",solver,space,
            unsafe_rectangle_2_domain,barrier_literals);
    }
    if(query=="lie-only") {
        List<SmtTheoryPrimitiveLiteral> lie_only_literals({lie_violation});
        timed_solve(
            "lie-only",solver,space,domain,lie_only_literals);
    } else {
        timed_solve(
            "lie",solver,space,domain,lie_literals);
    }

    return 0;
}
