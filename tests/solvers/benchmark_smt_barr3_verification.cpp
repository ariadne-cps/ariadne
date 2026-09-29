/***************************************************************************
 * benchmark_smt_barr3_verification.cpp
 *
 * Standalone reproduction of the published FOSSIL/dReal Barr3 BarrierAlt
 * counterexample queries using Ariadne's epsilon-SMT solver.
 ***************************************************************************/

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

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
    if(argument=="all" || argument=="lie" || argument=="lie-only" || argument=="eval" || argument=="taylor" || argument=="affine" || argument=="mean-value" || argument=="lie-components" || argument=="lie-split-profile" || argument=="lie-dynamics-rewrite" || argument=="lie-correlation-profile" || argument=="lie-gradient-reassociation" || argument=="lie-gradient-sign-profile" || argument=="lie-gradient-split-profile" || argument=="lie-gradient-mean-value-profile" || argument=="lie-gradient-quadrant-profile" || argument=="lie-gradient-composite-profile" || argument=="lie-gradient-symbolic-procedure-profile" || argument=="lie-gradient-shared-procedure-profile" || argument=="lie-gradient-shared-procedure-check" || argument=="lie-gradient-expression-cse-profile" || argument=="lie-gradient-expression-cse-frontier" || argument=="lie-gradient-cse-prune-profile" || argument=="lie-root-cse-range-profile" || argument=="lie-root-cse-search" || argument=="lie-algebraic-form-profile" || argument=="lie-directional-propagation-profile" || argument=="lie-layer2-directional-profile" || argument=="lie-gradient-block-reassociation-profile" || argument=="lie-gradient-width-attribution-profile" || argument=="lie-gradient-correlation-attribution-profile") { return argument; }
    throw std::runtime_error(
        "Usage: benchmark_smt_barr3_verification "
        "[positive-box-limit|full] [sensitivity|geometric|lookahead] "
        "[witness|no-witness] [shaving|no-shaving] [hull|no-hull] "
        "[all|lie|lie-only|eval|taylor|affine|mean-value|lie-components|lie-split-profile|lie-dynamics-rewrite|lie-correlation-profile|lie-gradient-reassociation|lie-gradient-sign-profile|lie-gradient-split-profile|lie-gradient-mean-value-profile|lie-gradient-quadrant-profile|lie-gradient-composite-profile|lie-gradient-symbolic-procedure-profile|lie-gradient-shared-procedure-profile|lie-gradient-shared-procedure-check|lie-gradient-expression-cse-profile|lie-gradient-expression-cse-frontier|lie-gradient-cse-prune-profile|lie-root-cse-range-profile|lie-root-cse-search|lie-algebraic-form-profile|lie-directional-propagation-profile|lie-layer2-directional-profile|lie-gradient-block-reassociation-profile|lie-gradient-width-attribution-profile|lie-gradient-correlation-attribution-profile]");
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
        "[all|lie|lie-only|eval|taylor|affine|mean-value|lie-components|lie-split-profile|lie-dynamics-rewrite|lie-correlation-profile|lie-gradient-reassociation|lie-gradient-sign-profile|lie-gradient-split-profile|lie-gradient-mean-value-profile|lie-gradient-quadrant-profile|lie-gradient-composite-profile|lie-gradient-symbolic-procedure-profile|lie-gradient-shared-procedure-profile|lie-gradient-shared-procedure-check|lie-gradient-expression-cse-profile|lie-gradient-expression-cse-frontier|lie-gradient-cse-prune-profile|lie-root-cse-range-profile|lie-root-cse-search|lie-algebraic-form-profile|lie-directional-propagation-profile|lie-layer2-directional-profile|lie-gradient-block-reassociation-profile|lie-gradient-width-attribution-profile|lie-gradient-correlation-attribution-profile] [monotone|no-monotone]");
}

String lie_literal_order_from_argument(Int argc,const char* argv[]) {
    if(argc<=8) { return "barrier-first"; }
    String argument(argv[8]);
    if(argument=="barrier-first" || argument=="lie-first") { return argument; }
    throw std::runtime_error(
        "Usage: benchmark_smt_barr3_verification "
        "[positive-box-limit|full] [sensitivity|geometric|lookahead] "
        "[witness|no-witness] [shaving|no-shaving] [hull|no-hull] "
        "[all|lie|lie-only|eval|taylor|affine|mean-value|lie-components|lie-split-profile|lie-dynamics-rewrite|lie-correlation-profile|lie-gradient-reassociation|lie-gradient-sign-profile|lie-gradient-split-profile|lie-gradient-mean-value-profile|lie-gradient-quadrant-profile|lie-gradient-composite-profile|lie-gradient-symbolic-procedure-profile|lie-gradient-shared-procedure-profile|lie-gradient-shared-procedure-check|lie-gradient-expression-cse-profile|lie-gradient-expression-cse-frontier|lie-gradient-cse-prune-profile|lie-root-cse-range-profile|lie-root-cse-search|lie-algebraic-form-profile|lie-directional-propagation-profile|lie-layer2-directional-profile|lie-gradient-block-reassociation-profile|lie-gradient-width-attribution-profile|lie-gradient-correlation-attribution-profile] [monotone|no-monotone] "
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
        "[all|lie|lie-only|eval|taylor|affine|mean-value|lie-components|lie-split-profile|lie-dynamics-rewrite|lie-correlation-profile|lie-gradient-reassociation|lie-gradient-sign-profile|lie-gradient-split-profile|lie-gradient-mean-value-profile|lie-gradient-quadrant-profile|lie-gradient-composite-profile|lie-gradient-symbolic-procedure-profile|lie-gradient-shared-procedure-profile|lie-gradient-shared-procedure-check|lie-gradient-expression-cse-profile|lie-gradient-expression-cse-frontier|lie-gradient-cse-prune-profile|lie-root-cse-range-profile|lie-root-cse-search|lie-algebraic-form-profile|lie-directional-propagation-profile|lie-layer2-directional-profile|lie-gradient-block-reassociation-profile|lie-gradient-width-attribution-profile|lie-gradient-correlation-attribution-profile] [monotone|no-monotone] "
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

Void profile_interval_expression(
    String const& name,
    RealExpression const& expression,
    RealSpace const& space,
    UpperBoxType const& domain)
{
    ValidatedScalarMultivariateFunction function=make_function(space,expression);
    Stopwatch<Milliseconds> stopwatch;
    UpperIntervalType image=apply(function,domain);
    stopwatch.click();
    std::cout << "[lie-component] " << name
              << " time=" << stopwatch.elapsed_seconds()
              << " image=" << image
              << std::endl;
}

Void profile_split_candidate(
    String const& name,
    ValidatedScalarMultivariateFunction const& function,
    UpperBoxType const& domain,
    SizeType coordinate)
{
    auto children=domain.split(coordinate);

    Stopwatch<Milliseconds> stopwatch;
    UpperIntervalType first_image=apply(function,children.first);
    UpperIntervalType second_image=apply(function,children.second);
    stopwatch.click();

    std::cout << "[lie-split-profile] " << name
              << " coordinate=" << coordinate
              << " time=" << stopwatch.elapsed_seconds()
              << " first=" << first_image
              << " first-width=" << first_image.width()
              << " second=" << second_image
              << " second-width=" << second_image.width()
              << " width-sum=" << (first_image.width()+second_image.width())
              << std::endl;
}

struct SignProfileCounts {
    SizeType processed = 0u;
    SizeType pruned = 0u;
    SizeType split = 0u;
    SizeType split_db_dy_crosses_zero = 0u;
    SizeType split_dy_crosses_zero = 0u;
    SizeType split_both_cross_zero = 0u;
    SizeType split_both_sign_definite = 0u;
    double split_db_dy_width_sum = 0.0;
    double split_dy_width_sum = 0.0;
    double split_lie_width_sum = 0.0;
};

Bool crosses_zero(UpperIntervalType const& interval) {
    return possibly(interval.lower_bound()<=0)
        && possibly(interval.upper_bound()>=0);
}

Void profile_lie_correlation_frontier(
    SizeType box_limit,
    UpperBoxType const& domain,
    ValidatedScalarMultivariateFunction const& db_dy_function,
    ValidatedScalarMultivariateFunction const& dy_function,
    ValidatedScalarMultivariateFunction const& lie_function)
{
    std::vector<UpperBoxType> pending;
    pending.push_back(domain);
    SignProfileCounts counts;

    Stopwatch<Milliseconds> stopwatch;
    while(not pending.empty() && counts.processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++counts.processed;

        UpperIntervalType lie_image=apply(lie_function,box);
        if(definitely(lie_image.lower_bound()>=0)) {
            ++counts.pruned;
            continue;
        }

        auto children=box.split();
        if(definitely(children.first==children.second)) {
            continue;
        }

        ++counts.split;
        UpperIntervalType db_dy_image=apply(db_dy_function,box);
        UpperIntervalType dy_image=apply(dy_function,box);
        Bool const db_dy_crosses=crosses_zero(db_dy_image);
        Bool const dy_crosses=crosses_zero(dy_image);

        if(db_dy_crosses) { ++counts.split_db_dy_crosses_zero; }
        if(dy_crosses) { ++counts.split_dy_crosses_zero; }
        if(db_dy_crosses && dy_crosses) {
            ++counts.split_both_cross_zero;
        }
        if(not db_dy_crosses && not dy_crosses) {
            ++counts.split_both_sign_definite;
        }

        counts.split_db_dy_width_sum+=db_dy_image.width().raw().get_d();
        counts.split_dy_width_sum+=dy_image.width().raw().get_d();
        counts.split_lie_width_sum+=lie_image.width().raw().get_d();

        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }
    stopwatch.click();

    auto average=[](double sum,SizeType count) {
        return count==0u ? 0.0 : sum/static_cast<double>(count);
    };

    std::cout << "[lie-correlation-profile]"
              << " time=" << stopwatch.elapsed_seconds()
              << " processed=" << counts.processed
              << " pruned=" << counts.pruned
              << " split=" << counts.split
              << " split-db/dy-cross-zero="
              << counts.split_db_dy_crosses_zero
              << " split-dy-cross-zero="
              << counts.split_dy_crosses_zero
              << " split-both-cross-zero="
              << counts.split_both_cross_zero
              << " split-both-sign-definite="
              << counts.split_both_sign_definite
              << " avg-split-db/dy-width="
              << average(counts.split_db_dy_width_sum,counts.split)
              << " avg-split-dy-width="
              << average(counts.split_dy_width_sum,counts.split)
              << " avg-split-lie-width="
              << average(counts.split_lie_width_sum,counts.split)
              << std::endl;
}

struct GradientSignSampleCounts {
    SizeType processed = 0u;
    SizeType pruned = 0u;
    SizeType split = 0u;
    SizeType db_dy_cross_zero = 0u;
    SizeType sampled_mixed_sign = 0u;
    SizeType sampled_positive_only = 0u;
    SizeType sampled_negative_only = 0u;
    SizeType sampled_ambiguous = 0u;
    SizeType point_evaluations = 0u;
};

Int point_sign(
    ValidatedScalarMultivariateFunction const& function,
    UpperBoxType const& point)
{
    UpperIntervalType image=apply(function,point);
    if(definitely(image.lower_bound()>0)) { return 1; }
    if(definitely(image.upper_bound()<0)) { return -1; }
    return 0;
}

Void profile_lie_gradient_sign_frontier(
    SizeType box_limit,
    UpperBoxType const& domain,
    ValidatedScalarMultivariateFunction const& db_dy_function,
    ValidatedScalarMultivariateFunction const& lie_function)
{
    std::vector<UpperBoxType> pending;
    pending.push_back(domain);
    GradientSignSampleCounts counts;

    Stopwatch<Milliseconds> stopwatch;
    while(not pending.empty() && counts.processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++counts.processed;

        UpperIntervalType lie_image=apply(lie_function,box);
        if(definitely(lie_image.lower_bound()>=0)) {
            ++counts.pruned;
            continue;
        }

        auto children=box.split();
        if(definitely(children.first==children.second)) {
            continue;
        }
        ++counts.split;

        UpperIntervalType db_dy_image=apply(db_dy_function,box);
        if(crosses_zero(db_dy_image)) {
            ++counts.db_dy_cross_zero;
            Bool positive=false;
            Bool negative=false;
            Bool ambiguous=false;

            UpperBoxType midpoint_box(box);
            for(SizeType i=0u;i!=box.dimension();++i) {
                auto midpoint=box[i].midpoint();
                midpoint_box[i]=UpperIntervalType(midpoint.raw(),midpoint.raw());
            }
            Int sign=point_sign(db_dy_function,midpoint_box);
            ++counts.point_evaluations;
            positive=positive || sign>0;
            negative=negative || sign<0;
            ambiguous=ambiguous || sign==0;

            SizeType const corner_count=SizeType(1u)<<box.dimension();
            for(SizeType corner=0u;corner!=corner_count;++corner) {
                UpperBoxType point(box);
                for(SizeType i=0u;i!=box.dimension();++i) {
                    auto const endpoint=((corner>>i)&SizeType(1u))==0u
                        ? box[i].lower_bound().raw()
                        : box[i].upper_bound().raw();
                    point[i]=UpperIntervalType(endpoint,endpoint);
                }
                sign=point_sign(db_dy_function,point);
                ++counts.point_evaluations;
                positive=positive || sign>0;
                negative=negative || sign<0;
                ambiguous=ambiguous || sign==0;
            }

            if(positive && negative) {
                ++counts.sampled_mixed_sign;
            } else if(ambiguous) {
                ++counts.sampled_ambiguous;
            } else if(positive) {
                ++counts.sampled_positive_only;
            } else if(negative) {
                ++counts.sampled_negative_only;
            } else {
                ++counts.sampled_ambiguous;
            }
        }

        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }
    stopwatch.click();

    std::cout << "[lie-gradient-sign-profile]"
              << " time=" << stopwatch.elapsed_seconds()
              << " processed=" << counts.processed
              << " pruned=" << counts.pruned
              << " split=" << counts.split
              << " db/dy-cross-zero=" << counts.db_dy_cross_zero
              << " sampled-mixed-sign=" << counts.sampled_mixed_sign
              << " sampled-positive-only=" << counts.sampled_positive_only
              << " sampled-negative-only=" << counts.sampled_negative_only
              << " sampled-ambiguous=" << counts.sampled_ambiguous
              << " point-evals=" << counts.point_evaluations
              << std::endl;
}

struct GradientSplitProfileCounts {
    SizeType processed = 0u;
    SizeType pruned = 0u;
    SizeType split = 0u;
    SizeType db_dy_cross_zero = 0u;
    SizeType x_sign_definite_children = 0u;
    SizeType y_sign_definite_children = 0u;
    SizeType x_both_sign_definite = 0u;
    SizeType y_both_sign_definite = 0u;
    SizeType x_width_score_better = 0u;
    SizeType y_width_score_better = 0u;
    SizeType width_score_ties = 0u;
    double x_width_score_sum = 0.0;
    double y_width_score_sum = 0.0;
    SizeType candidate_evaluations = 0u;
};

Void profile_lie_gradient_split_frontier(
    SizeType box_limit,
    UpperBoxType const& domain,
    ValidatedScalarMultivariateFunction const& db_dy_function,
    ValidatedScalarMultivariateFunction const& lie_function)
{
    std::vector<UpperBoxType> pending;
    pending.push_back(domain);
    GradientSplitProfileCounts counts;

    Stopwatch<Milliseconds> stopwatch;
    while(not pending.empty() && counts.processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++counts.processed;

        UpperIntervalType lie_image=apply(lie_function,box);
        if(definitely(lie_image.lower_bound()>=0)) {
            ++counts.pruned;
            continue;
        }

        auto children=box.split();
        if(definitely(children.first==children.second)) {
            continue;
        }
        ++counts.split;

        UpperIntervalType db_dy_image=apply(db_dy_function,box);
        if(crosses_zero(db_dy_image)) {
            ++counts.db_dy_cross_zero;

            double scores[2]={0.0,0.0};
            SizeType sign_definite_children[2]={0u,0u};
            Bool both_sign_definite[2]={false,false};

            for(SizeType coordinate=0u;coordinate!=2u;++coordinate) {
                auto candidate_children=box.split(coordinate);
                UpperIntervalType first_image=
                    apply(db_dy_function,candidate_children.first);
                UpperIntervalType second_image=
                    apply(db_dy_function,candidate_children.second);
                counts.candidate_evaluations+=2u;

                Bool const first_definite=not crosses_zero(first_image);
                Bool const second_definite=not crosses_zero(second_image);
                sign_definite_children[coordinate]+=
                    static_cast<SizeType>(first_definite)
                    + static_cast<SizeType>(second_definite);
                both_sign_definite[coordinate]=first_definite && second_definite;

                scores[coordinate]=
                    first_image.width().raw().get_d()
                    + second_image.width().raw().get_d();
            }

            counts.x_sign_definite_children+=sign_definite_children[0];
            counts.y_sign_definite_children+=sign_definite_children[1];
            if(both_sign_definite[0]) { ++counts.x_both_sign_definite; }
            if(both_sign_definite[1]) { ++counts.y_both_sign_definite; }
            counts.x_width_score_sum+=scores[0];
            counts.y_width_score_sum+=scores[1];

            if(scores[0]<scores[1]) {
                ++counts.x_width_score_better;
            } else if(scores[1]<scores[0]) {
                ++counts.y_width_score_better;
            } else {
                ++counts.width_score_ties;
            }
        }

        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }
    stopwatch.click();

    auto average=[](double sum,SizeType count) {
        return count==0u ? 0.0 : sum/static_cast<double>(count);
    };

    std::cout << "[lie-gradient-split-profile]"
              << " time=" << stopwatch.elapsed_seconds()
              << " processed=" << counts.processed
              << " pruned=" << counts.pruned
              << " split=" << counts.split
              << " db/dy-cross-zero=" << counts.db_dy_cross_zero
              << " x-sign-definite-children="
              << counts.x_sign_definite_children
              << " y-sign-definite-children="
              << counts.y_sign_definite_children
              << " x-both-sign-definite=" << counts.x_both_sign_definite
              << " y-both-sign-definite=" << counts.y_both_sign_definite
              << " x-width-score-better=" << counts.x_width_score_better
              << " y-width-score-better=" << counts.y_width_score_better
              << " width-score-ties=" << counts.width_score_ties
              << " avg-x-width-score="
              << average(counts.x_width_score_sum,counts.db_dy_cross_zero)
              << " avg-y-width-score="
              << average(counts.y_width_score_sum,counts.db_dy_cross_zero)
              << " candidate-evals=" << counts.candidate_evaluations
              << std::endl;
}

struct GradientMeanValueProfileCounts {
    SizeType processed = 0u;
    SizeType pruned = 0u;
    SizeType split = 0u;
    SizeType db_dy_cross_zero = 0u;
    SizeType mean_value_sign_definite = 0u;
    SizeType intersected_sign_definite = 0u;
    double direct_width_sum = 0.0;
    double mean_value_width_sum = 0.0;
    double intersected_width_sum = 0.0;
    SizeType derivative_evaluations = 0u;
    SizeType midpoint_evaluations = 0u;
};

Void profile_lie_gradient_mean_value_frontier(
    SizeType box_limit,
    UpperBoxType const& domain,
    ValidatedScalarMultivariateFunction const& db_dy_function,
    ValidatedScalarMultivariateFunction const& lie_function)
{
    std::vector<ValidatedScalarMultivariateFunction> derivatives;
    derivatives.reserve(domain.dimension());
    Stopwatch<Milliseconds> build_stopwatch;
    for(SizeType i=0u;i!=domain.dimension();++i) {
        derivatives.push_back(db_dy_function.derivative(i));
    }
    build_stopwatch.click();

    std::vector<UpperBoxType> pending;
    pending.push_back(domain);
    GradientMeanValueProfileCounts counts;

    Stopwatch<Milliseconds> stopwatch;
    while(not pending.empty() && counts.processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++counts.processed;

        UpperIntervalType lie_image=apply(lie_function,box);
        if(definitely(lie_image.lower_bound()>=0)) {
            ++counts.pruned;
            continue;
        }

        auto children=box.split();
        if(definitely(children.first==children.second)) {
            continue;
        }
        ++counts.split;

        UpperIntervalType direct_image=apply(db_dy_function,box);
        if(crosses_zero(direct_image)) {
            ++counts.db_dy_cross_zero;
            counts.direct_width_sum+=direct_image.width().raw().get_d();

            UpperBoxType midpoint_box(box);
            for(SizeType i=0u;i!=box.dimension();++i) {
                auto midpoint=box[i].midpoint();
                midpoint_box[i]=UpperIntervalType(midpoint.raw(),midpoint.raw());
            }

            UpperIntervalType mean_value_image=apply(db_dy_function,midpoint_box);
            ++counts.midpoint_evaluations;
            for(SizeType i=0u;i!=box.dimension();++i) {
                UpperIntervalType derivative_image=apply(derivatives[i],box);
                ++counts.derivative_evaluations;
                auto midpoint=box[i].midpoint();
                UpperIntervalType midpoint_interval(midpoint.raw(),midpoint.raw());
                mean_value_image+=
                    derivative_image*(box[i]-midpoint_interval);
            }

            UpperIntervalType intersected_image=
                intersection(direct_image,mean_value_image);
            counts.mean_value_width_sum+=
                mean_value_image.width().raw().get_d();
            counts.intersected_width_sum+=
                intersected_image.width().raw().get_d();

            if(not crosses_zero(mean_value_image)) {
                ++counts.mean_value_sign_definite;
            }
            if(not crosses_zero(intersected_image)) {
                ++counts.intersected_sign_definite;
            }
        }

        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }
    stopwatch.click();

    auto average=[](double sum,SizeType count) {
        return count==0u ? 0.0 : sum/static_cast<double>(count);
    };

    std::cout << "[lie-gradient-mean-value-profile]"
              << " time=" << stopwatch.elapsed_seconds()
              << " derivative-build-time=" << build_stopwatch.elapsed_seconds()
              << " processed=" << counts.processed
              << " pruned=" << counts.pruned
              << " split=" << counts.split
              << " db/dy-cross-zero=" << counts.db_dy_cross_zero
              << " mean-value-sign-definite="
              << counts.mean_value_sign_definite
              << " intersected-sign-definite="
              << counts.intersected_sign_definite
              << " avg-direct-width="
              << average(counts.direct_width_sum,counts.db_dy_cross_zero)
              << " avg-mean-value-width="
              << average(counts.mean_value_width_sum,counts.db_dy_cross_zero)
              << " avg-intersected-width="
              << average(counts.intersected_width_sum,counts.db_dy_cross_zero)
              << " midpoint-evals=" << counts.midpoint_evaluations
              << " derivative-evals=" << counts.derivative_evaluations
              << std::endl;
}

struct GradientQuadrantProfileCounts {
    SizeType processed = 0u;
    SizeType pruned = 0u;
    SizeType split = 0u;
    SizeType db_dy_cross_zero = 0u;
    SizeType quadrant_hull_sign_definite = 0u;
    SizeType all_quadrants_sign_definite = 0u;
    SizeType mixed_sign_quadrants = 0u;
    SizeType residual_crossing_quadrants = 0u;
    SizeType sign_definite_quadrants = 0u;
    double direct_width_sum = 0.0;
    double quadrant_hull_width_sum = 0.0;
    SizeType candidate_evaluations = 0u;
};

Void profile_lie_gradient_quadrant_frontier(
    SizeType box_limit,
    UpperBoxType const& domain,
    ValidatedScalarMultivariateFunction const& db_dy_function,
    ValidatedScalarMultivariateFunction const& lie_function)
{
    std::vector<UpperBoxType> pending;
    pending.push_back(domain);
    GradientQuadrantProfileCounts counts;

    Stopwatch<Milliseconds> stopwatch;
    while(not pending.empty() && counts.processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++counts.processed;

        UpperIntervalType lie_image=apply(lie_function,box);
        if(definitely(lie_image.lower_bound()>=0)) {
            ++counts.pruned;
            continue;
        }

        auto children=box.split();
        if(definitely(children.first==children.second)) {
            continue;
        }
        ++counts.split;

        UpperIntervalType direct_image=apply(db_dy_function,box);
        if(crosses_zero(direct_image)) {
            ++counts.db_dy_cross_zero;
            counts.direct_width_sum+=direct_image.width().raw().get_d();

            auto x_children=box.split(0u);
            auto lower_quadrants=x_children.first.split(1u);
            auto upper_quadrants=x_children.second.split(1u);
            UpperBoxType quadrants[4]={
                lower_quadrants.first,
                lower_quadrants.second,
                upper_quadrants.first,
                upper_quadrants.second
            };

            UpperIntervalType quadrant_image=apply(db_dy_function,quadrants[0]);
            ++counts.candidate_evaluations;
            UpperIntervalType quadrant_hull=quadrant_image;
            Bool all_definite=not crosses_zero(quadrant_image);
            Bool any_positive=definitely(quadrant_image.lower_bound()>0);
            Bool any_negative=definitely(quadrant_image.upper_bound()<0);
            if(all_definite) { ++counts.sign_definite_quadrants; }

            for(SizeType i=1u;i!=4u;++i) {
                quadrant_image=apply(db_dy_function,quadrants[i]);
                ++counts.candidate_evaluations;
                quadrant_hull=hull(quadrant_hull,quadrant_image);

                Bool const definite=not crosses_zero(quadrant_image);
                all_definite=all_definite && definite;
                any_positive=any_positive
                    || definitely(quadrant_image.lower_bound()>0);
                any_negative=any_negative
                    || definitely(quadrant_image.upper_bound()<0);
                if(definite) { ++counts.sign_definite_quadrants; }
            }

            counts.quadrant_hull_width_sum+=
                quadrant_hull.width().raw().get_d();

            if(not crosses_zero(quadrant_hull)) {
                ++counts.quadrant_hull_sign_definite;
            }
            if(all_definite) {
                ++counts.all_quadrants_sign_definite;
            } else {
                ++counts.residual_crossing_quadrants;
            }
            if(any_positive && any_negative) {
                ++counts.mixed_sign_quadrants;
            }
        }

        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }
    stopwatch.click();

    auto average=[](double sum,SizeType count) {
        return count==0u ? 0.0 : sum/static_cast<double>(count);
    };

    std::cout << "[lie-gradient-quadrant-profile]"
              << " time=" << stopwatch.elapsed_seconds()
              << " processed=" << counts.processed
              << " pruned=" << counts.pruned
              << " split=" << counts.split
              << " db/dy-cross-zero=" << counts.db_dy_cross_zero
              << " quadrant-hull-sign-definite="
              << counts.quadrant_hull_sign_definite
              << " all-quadrants-sign-definite="
              << counts.all_quadrants_sign_definite
              << " mixed-sign-quadrants="
              << counts.mixed_sign_quadrants
              << " residual-crossing-quadrants="
              << counts.residual_crossing_quadrants
              << " sign-definite-quadrants="
              << counts.sign_definite_quadrants
              << " avg-direct-width="
              << average(counts.direct_width_sum,counts.db_dy_cross_zero)
              << " avg-quadrant-hull-width="
              << average(counts.quadrant_hull_width_sum,counts.db_dy_cross_zero)
              << " candidate-evals=" << counts.candidate_evaluations
              << std::endl;
}

struct GradientCompositeProfileCounts {
    SizeType processed = 0u;
    SizeType pruned = 0u;
    SizeType split = 0u;
    SizeType natural_cross_zero = 0u;
    SizeType centered_sign_definite = 0u;
    SizeType monotone_sign_definite = 0u;
    SizeType combined_sign_definite = 0u;
    SizeType monotone_boxes = 0u;
    SizeType monotone_coordinates = 0u;
    SizeType midpoint_evaluations = 0u;
    SizeType monotone_evaluations = 0u;
    SizeType gradient_evaluations = 0u;
    double natural_width_sum = 0.0;
    double centered_width_sum = 0.0;
    double monotone_width_sum = 0.0;
    double combined_width_sum = 0.0;
    double gradient_seconds = 0.0;
    double midpoint_seconds = 0.0;
    double monotone_seconds = 0.0;
};

Vector<FloatDPBounds> bounds_arguments(UpperBoxType const& box)
{
    Vector<FloatDPBounds> arguments(
        box.dimension(),FloatDPBounds(DoublePrecision()));
    for(SizeType i=0u;i!=box.dimension();++i) {
        arguments[i]=FloatDPBounds(
            box[i].lower_bound().raw(),
            box[i].upper_bound().raw());
    }
    return arguments;
}

Void profile_lie_gradient_composite_frontier(
    SizeType box_limit,
    UpperBoxType const& domain,
    ValidatedScalarMultivariateFunction const& db_dy_function,
    ValidatedScalarMultivariateFunction const& lie_function)
{
    Stopwatch<Milliseconds> build_stopwatch;
    ValidatedProcedure procedure(db_dy_function);
    build_stopwatch.click();

    std::vector<UpperBoxType> pending;
    pending.push_back(domain);
    GradientCompositeProfileCounts counts;

    Stopwatch<Milliseconds> stopwatch;
    while(not pending.empty() && counts.processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++counts.processed;

        UpperIntervalType lie_image=apply(lie_function,box);
        if(definitely(lie_image.lower_bound()>=0)) {
            ++counts.pruned;
            continue;
        }

        auto children=box.split();
        if(definitely(children.first==children.second)) {
            continue;
        }
        ++counts.split;

        Vector<FloatDPBounds> arguments=bounds_arguments(box);
        UpperIntervalType natural_image=
            make_interval(evaluate(procedure,arguments));
        if(crosses_zero(natural_image)) {
            ++counts.natural_cross_zero;
            counts.natural_width_sum+=natural_image.width().raw().get_d();

            Stopwatch<Milliseconds> gradient_stopwatch;
            Covector<FloatDPBounds> gradient_image=
                gradient(procedure,arguments);
            gradient_stopwatch.click();
            counts.gradient_seconds+=gradient_stopwatch.elapsed_seconds();
            ++counts.gradient_evaluations;

            Vector<FloatDPBounds> midpoint_arguments(arguments);
            for(SizeType i=0u;i!=box.dimension();++i) {
                auto midpoint=box[i].midpoint();
                midpoint_arguments[i]=FloatDPBounds(
                    midpoint.raw(),midpoint.raw());
            }

            Stopwatch<Milliseconds> midpoint_stopwatch;
            FloatDPBounds centered_bounds=
                evaluate(procedure,midpoint_arguments);
            midpoint_stopwatch.click();
            counts.midpoint_seconds+=midpoint_stopwatch.elapsed_seconds();
            ++counts.midpoint_evaluations;

            for(SizeType i=0u;i!=box.dimension();++i) {
                centered_bounds+=gradient_image[i]
                    *(arguments[i]-midpoint_arguments[i]);
            }
            UpperIntervalType centered_image=make_interval(centered_bounds);
            counts.centered_width_sum+=
                centered_image.width().raw().get_d();
            if(not crosses_zero(centered_image)) {
                ++counts.centered_sign_definite;
            }

            Vector<FloatDPBounds> lower_arguments(arguments);
            Vector<FloatDPBounds> upper_arguments(arguments);
            SizeType monotone_coordinates=0u;
            for(SizeType i=0u;i!=box.dimension();++i) {
                UpperIntervalType derivative_image=
                    make_interval(gradient_image[i]);
                if(definitely(derivative_image.lower_bound()>0)) {
                    auto lower=box[i].lower_bound().raw();
                    auto upper=box[i].upper_bound().raw();
                    lower_arguments[i]=FloatDPBounds(lower,lower);
                    upper_arguments[i]=FloatDPBounds(upper,upper);
                    ++monotone_coordinates;
                } else if(definitely(derivative_image.upper_bound()<0)) {
                    auto lower=box[i].lower_bound().raw();
                    auto upper=box[i].upper_bound().raw();
                    lower_arguments[i]=FloatDPBounds(upper,upper);
                    upper_arguments[i]=FloatDPBounds(lower,lower);
                    ++monotone_coordinates;
                }
            }

            UpperIntervalType monotone_image=natural_image;
            if(monotone_coordinates!=0u) {
                ++counts.monotone_boxes;
                counts.monotone_coordinates+=monotone_coordinates;

                Stopwatch<Milliseconds> monotone_stopwatch;
                FloatDPBounds lower_image=
                    evaluate(procedure,lower_arguments);
                FloatDPBounds upper_image=
                    evaluate(procedure,upper_arguments);
                monotone_stopwatch.click();
                counts.monotone_seconds+=monotone_stopwatch.elapsed_seconds();
                counts.monotone_evaluations+=2u;

                monotone_image=make_interval(
                    FloatDPBounds(lower_image.lower(),upper_image.upper()));
            }
            counts.monotone_width_sum+=
                monotone_image.width().raw().get_d();
            if(not crosses_zero(monotone_image)) {
                ++counts.monotone_sign_definite;
            }

            UpperIntervalType combined_image=
                intersection(natural_image,centered_image);
            combined_image=intersection(combined_image,monotone_image);
            counts.combined_width_sum+=
                combined_image.width().raw().get_d();
            if(not crosses_zero(combined_image)) {
                ++counts.combined_sign_definite;
            }
        }

        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }
    stopwatch.click();

    auto average=[](double sum,SizeType count) {
        return count==0u ? 0.0 : sum/static_cast<double>(count);
    };

    std::cout << "[lie-gradient-composite-profile]"
              << " time=" << stopwatch.elapsed_seconds()
              << " procedure-build-time=" << build_stopwatch.elapsed_seconds()
              << " processed=" << counts.processed
              << " pruned=" << counts.pruned
              << " split=" << counts.split
              << " natural-cross-zero=" << counts.natural_cross_zero
              << " centered-sign-definite="
              << counts.centered_sign_definite
              << " monotone-sign-definite="
              << counts.monotone_sign_definite
              << " combined-sign-definite="
              << counts.combined_sign_definite
              << " monotone-boxes=" << counts.monotone_boxes
              << " monotone-coordinates=" << counts.monotone_coordinates
              << " avg-natural-width="
              << average(counts.natural_width_sum,counts.natural_cross_zero)
              << " avg-centered-width="
              << average(counts.centered_width_sum,counts.natural_cross_zero)
              << " avg-monotone-width="
              << average(counts.monotone_width_sum,counts.natural_cross_zero)
              << " avg-combined-width="
              << average(counts.combined_width_sum,counts.natural_cross_zero)
              << " gradient-evals=" << counts.gradient_evaluations
              << " gradient-time=" << counts.gradient_seconds
              << " midpoint-evals=" << counts.midpoint_evaluations
              << " midpoint-time=" << counts.midpoint_seconds
              << " monotone-evals=" << counts.monotone_evaluations
              << " monotone-time=" << counts.monotone_seconds
              << std::endl;
}

struct GradientSymbolicProcedureProfileCounts {
    SizeType processed = 0u;
    SizeType pruned = 0u;
    SizeType split = 0u;
    SizeType natural_cross_zero = 0u;
    SizeType centered_sign_definite = 0u;
    SizeType monotone_sign_definite = 0u;
    SizeType combined_sign_definite = 0u;
    SizeType monotone_boxes = 0u;
    SizeType monotone_coordinates = 0u;
    SizeType natural_evaluations = 0u;
    SizeType midpoint_evaluations = 0u;
    SizeType derivative_evaluations = 0u;
    SizeType monotone_evaluations = 0u;
    double natural_width_sum = 0.0;
    double centered_width_sum = 0.0;
    double monotone_width_sum = 0.0;
    double combined_width_sum = 0.0;
    double natural_seconds = 0.0;
    double midpoint_seconds = 0.0;
    double derivative_seconds = 0.0;
    double monotone_seconds = 0.0;
};

Void profile_lie_gradient_symbolic_procedure_frontier(
    SizeType box_limit,
    UpperBoxType const& domain,
    ValidatedScalarMultivariateFunction const& db_dy_function,
    ValidatedScalarMultivariateFunction const& lie_function)
{
    Stopwatch<Milliseconds> derivative_build_stopwatch;
    std::vector<ValidatedScalarMultivariateFunction> derivatives;
    derivatives.reserve(domain.dimension());
    for(SizeType i=0u;i!=domain.dimension();++i) {
        derivatives.push_back(db_dy_function.derivative(i));
    }
    derivative_build_stopwatch.click();

    Stopwatch<Milliseconds> procedure_build_stopwatch;
    ValidatedProcedure procedure(db_dy_function);
    std::vector<ValidatedProcedure> derivative_procedures;
    derivative_procedures.reserve(domain.dimension());
    for(SizeType i=0u;i!=domain.dimension();++i) {
        derivative_procedures.emplace_back(derivatives[i]);
    }
    procedure_build_stopwatch.click();

    std::vector<UpperBoxType> pending;
    pending.push_back(domain);
    GradientSymbolicProcedureProfileCounts counts;

    Stopwatch<Milliseconds> stopwatch;
    while(not pending.empty() && counts.processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++counts.processed;

        UpperIntervalType lie_image=apply(lie_function,box);
        if(definitely(lie_image.lower_bound()>=0)) {
            ++counts.pruned;
            continue;
        }

        auto children=box.split();
        if(definitely(children.first==children.second)) {
            continue;
        }
        ++counts.split;

        Vector<FloatDPBounds> arguments=bounds_arguments(box);

        Stopwatch<Milliseconds> natural_stopwatch;
        UpperIntervalType natural_image=
            make_interval(evaluate(procedure,arguments));
        natural_stopwatch.click();
        counts.natural_seconds+=natural_stopwatch.elapsed_seconds();
        ++counts.natural_evaluations;

        if(crosses_zero(natural_image)) {
            ++counts.natural_cross_zero;
            counts.natural_width_sum+=natural_image.width().raw().get_d();

            Vector<FloatDPBounds> midpoint_arguments(arguments);
            for(SizeType i=0u;i!=box.dimension();++i) {
                auto midpoint=box[i].midpoint();
                midpoint_arguments[i]=FloatDPBounds(
                    midpoint.raw(),midpoint.raw());
            }

            Stopwatch<Milliseconds> midpoint_stopwatch;
            FloatDPBounds centered_bounds=
                evaluate(procedure,midpoint_arguments);
            midpoint_stopwatch.click();
            counts.midpoint_seconds+=midpoint_stopwatch.elapsed_seconds();
            ++counts.midpoint_evaluations;

            std::vector<FloatDPBounds> derivative_images;
            derivative_images.reserve(domain.dimension());

            Stopwatch<Milliseconds> derivative_stopwatch;
            for(SizeType i=0u;i!=domain.dimension();++i) {
                derivative_images.push_back(
                    evaluate(derivative_procedures[i],arguments));
                ++counts.derivative_evaluations;
            }
            derivative_stopwatch.click();
            counts.derivative_seconds+=derivative_stopwatch.elapsed_seconds();

            for(SizeType i=0u;i!=box.dimension();++i) {
                centered_bounds+=derivative_images[i]
                    *(arguments[i]-midpoint_arguments[i]);
            }
            UpperIntervalType centered_image=make_interval(centered_bounds);
            counts.centered_width_sum+=
                centered_image.width().raw().get_d();
            if(not crosses_zero(centered_image)) {
                ++counts.centered_sign_definite;
            }

            Vector<FloatDPBounds> lower_arguments(arguments);
            Vector<FloatDPBounds> upper_arguments(arguments);
            SizeType monotone_coordinates=0u;
            for(SizeType i=0u;i!=box.dimension();++i) {
                UpperIntervalType derivative_image=
                    make_interval(derivative_images[i]);
                if(definitely(derivative_image.lower_bound()>0)) {
                    auto lower=box[i].lower_bound().raw();
                    auto upper=box[i].upper_bound().raw();
                    lower_arguments[i]=FloatDPBounds(lower,lower);
                    upper_arguments[i]=FloatDPBounds(upper,upper);
                    ++monotone_coordinates;
                } else if(definitely(derivative_image.upper_bound()<0)) {
                    auto lower=box[i].lower_bound().raw();
                    auto upper=box[i].upper_bound().raw();
                    lower_arguments[i]=FloatDPBounds(upper,upper);
                    upper_arguments[i]=FloatDPBounds(lower,lower);
                    ++monotone_coordinates;
                }
            }

            UpperIntervalType monotone_image=natural_image;
            if(monotone_coordinates!=0u) {
                ++counts.monotone_boxes;
                counts.monotone_coordinates+=monotone_coordinates;

                Stopwatch<Milliseconds> monotone_stopwatch;
                FloatDPBounds lower_image=
                    evaluate(procedure,lower_arguments);
                FloatDPBounds upper_image=
                    evaluate(procedure,upper_arguments);
                monotone_stopwatch.click();
                counts.monotone_seconds+=monotone_stopwatch.elapsed_seconds();
                counts.monotone_evaluations+=2u;

                monotone_image=make_interval(
                    FloatDPBounds(lower_image.lower(),upper_image.upper()));
            }
            counts.monotone_width_sum+=
                monotone_image.width().raw().get_d();
            if(not crosses_zero(monotone_image)) {
                ++counts.monotone_sign_definite;
            }

            UpperIntervalType combined_image=
                intersection(natural_image,centered_image);
            combined_image=intersection(combined_image,monotone_image);
            counts.combined_width_sum+=
                combined_image.width().raw().get_d();
            if(not crosses_zero(combined_image)) {
                ++counts.combined_sign_definite;
            }
        }

        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }
    stopwatch.click();

    auto average=[](double sum,SizeType count) {
        return count==0u ? 0.0 : sum/static_cast<double>(count);
    };

    std::cout << "[lie-gradient-symbolic-procedure-profile]"
              << " time=" << stopwatch.elapsed_seconds()
              << " derivative-build-time="
              << derivative_build_stopwatch.elapsed_seconds()
              << " procedure-build-time="
              << procedure_build_stopwatch.elapsed_seconds()
              << " processed=" << counts.processed
              << " pruned=" << counts.pruned
              << " split=" << counts.split
              << " natural-cross-zero=" << counts.natural_cross_zero
              << " centered-sign-definite="
              << counts.centered_sign_definite
              << " monotone-sign-definite="
              << counts.monotone_sign_definite
              << " combined-sign-definite="
              << counts.combined_sign_definite
              << " monotone-boxes=" << counts.monotone_boxes
              << " monotone-coordinates=" << counts.monotone_coordinates
              << " avg-natural-width="
              << average(counts.natural_width_sum,counts.natural_cross_zero)
              << " avg-centered-width="
              << average(counts.centered_width_sum,counts.natural_cross_zero)
              << " avg-monotone-width="
              << average(counts.monotone_width_sum,counts.natural_cross_zero)
              << " avg-combined-width="
              << average(counts.combined_width_sum,counts.natural_cross_zero)
              << " natural-evals=" << counts.natural_evaluations
              << " natural-time=" << counts.natural_seconds
              << " midpoint-evals=" << counts.midpoint_evaluations
              << " midpoint-time=" << counts.midpoint_seconds
              << " derivative-evals=" << counts.derivative_evaluations
              << " derivative-time=" << counts.derivative_seconds
              << " monotone-evals=" << counts.monotone_evaluations
              << " monotone-time=" << counts.monotone_seconds
              << std::endl;
}

struct GradientSharedProcedureProfileCounts {
    SizeType processed = 0u;
    SizeType pruned = 0u;
    SizeType split = 0u;
    SizeType natural_cross_zero = 0u;
    SizeType centered_sign_definite = 0u;
    SizeType monotone_sign_definite = 0u;
    SizeType combined_sign_definite = 0u;
    SizeType monotone_boxes = 0u;
    SizeType monotone_coordinates = 0u;
    SizeType shared_evaluations = 0u;
    SizeType midpoint_evaluations = 0u;
    SizeType monotone_evaluations = 0u;
    double natural_width_sum = 0.0;
    double centered_width_sum = 0.0;
    double monotone_width_sum = 0.0;
    double combined_width_sum = 0.0;
    double shared_seconds = 0.0;
    double midpoint_seconds = 0.0;
    double monotone_seconds = 0.0;
};

Void profile_lie_gradient_shared_procedure_frontier(
    SizeType box_limit,
    UpperBoxType const& domain,
    ValidatedScalarMultivariateFunction const& db_dy_function,
    ValidatedScalarMultivariateFunction const& lie_function)
{
    Stopwatch<Milliseconds> derivative_build_stopwatch;
    ValidatedScalarMultivariateFunction derivative_x=
        db_dy_function.derivative(0u);
    ValidatedScalarMultivariateFunction derivative_y=
        db_dy_function.derivative(1u);
    derivative_build_stopwatch.click();

    Stopwatch<Milliseconds> procedure_build_stopwatch;
    ValidatedVectorMultivariateFunction shared_function=
        join(db_dy_function,join(derivative_x,derivative_y));
    Vector<ValidatedProcedure> shared_procedure(shared_function);
    ValidatedProcedure scalar_procedure(db_dy_function);
    ValidatedProcedure derivative_x_procedure(derivative_x);
    ValidatedProcedure derivative_y_procedure(derivative_y);
    procedure_build_stopwatch.click();

    SizeType const separate_instruction_count=
        scalar_procedure._instructions.size()
        + derivative_x_procedure._instructions.size()
        + derivative_y_procedure._instructions.size();
    SizeType const shared_instruction_count=
        shared_procedure.temporaries_size();

    std::vector<UpperBoxType> pending;
    pending.push_back(domain);
    GradientSharedProcedureProfileCounts counts;

    Stopwatch<Milliseconds> stopwatch;
    while(not pending.empty() && counts.processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++counts.processed;

        UpperIntervalType lie_image=apply(lie_function,box);
        if(definitely(lie_image.lower_bound()>=0)) {
            ++counts.pruned;
            continue;
        }

        auto children=box.split();
        if(definitely(children.first==children.second)) {
            continue;
        }
        ++counts.split;

        Vector<FloatDPBounds> arguments=bounds_arguments(box);

        Stopwatch<Milliseconds> shared_stopwatch;
        Vector<FloatDPBounds> shared_images=
            evaluate(shared_procedure,arguments);
        shared_stopwatch.click();
        counts.shared_seconds+=shared_stopwatch.elapsed_seconds();
        ++counts.shared_evaluations;

        UpperIntervalType natural_image=make_interval(shared_images[0u]);
        if(crosses_zero(natural_image)) {
            ++counts.natural_cross_zero;
            counts.natural_width_sum+=natural_image.width().raw().get_d();

            Vector<FloatDPBounds> midpoint_arguments(arguments);
            for(SizeType i=0u;i!=box.dimension();++i) {
                auto midpoint=box[i].midpoint();
                midpoint_arguments[i]=FloatDPBounds(
                    midpoint.raw(),midpoint.raw());
            }

            Stopwatch<Milliseconds> midpoint_stopwatch;
            FloatDPBounds centered_bounds=
                evaluate(scalar_procedure,midpoint_arguments);
            midpoint_stopwatch.click();
            counts.midpoint_seconds+=midpoint_stopwatch.elapsed_seconds();
            ++counts.midpoint_evaluations;

            for(SizeType i=0u;i!=box.dimension();++i) {
                centered_bounds+=shared_images[i+1u]
                    *(arguments[i]-midpoint_arguments[i]);
            }
            UpperIntervalType centered_image=make_interval(centered_bounds);
            counts.centered_width_sum+=
                centered_image.width().raw().get_d();
            if(not crosses_zero(centered_image)) {
                ++counts.centered_sign_definite;
            }

            Vector<FloatDPBounds> lower_arguments(arguments);
            Vector<FloatDPBounds> upper_arguments(arguments);
            SizeType monotone_coordinates=0u;
            for(SizeType i=0u;i!=box.dimension();++i) {
                UpperIntervalType derivative_image=
                    make_interval(shared_images[i+1u]);
                if(definitely(derivative_image.lower_bound()>0)) {
                    auto lower=box[i].lower_bound().raw();
                    auto upper=box[i].upper_bound().raw();
                    lower_arguments[i]=FloatDPBounds(lower,lower);
                    upper_arguments[i]=FloatDPBounds(upper,upper);
                    ++monotone_coordinates;
                } else if(definitely(derivative_image.upper_bound()<0)) {
                    auto lower=box[i].lower_bound().raw();
                    auto upper=box[i].upper_bound().raw();
                    lower_arguments[i]=FloatDPBounds(upper,upper);
                    upper_arguments[i]=FloatDPBounds(lower,lower);
                    ++monotone_coordinates;
                }
            }

            UpperIntervalType monotone_image=natural_image;
            if(monotone_coordinates!=0u) {
                ++counts.monotone_boxes;
                counts.monotone_coordinates+=monotone_coordinates;

                Stopwatch<Milliseconds> monotone_stopwatch;
                FloatDPBounds lower_image=
                    evaluate(scalar_procedure,lower_arguments);
                FloatDPBounds upper_image=
                    evaluate(scalar_procedure,upper_arguments);
                monotone_stopwatch.click();
                counts.monotone_seconds+=monotone_stopwatch.elapsed_seconds();
                counts.monotone_evaluations+=2u;

                monotone_image=make_interval(
                    FloatDPBounds(lower_image.lower(),upper_image.upper()));
            }
            counts.monotone_width_sum+=
                monotone_image.width().raw().get_d();
            if(not crosses_zero(monotone_image)) {
                ++counts.monotone_sign_definite;
            }

            UpperIntervalType combined_image=
                intersection(natural_image,centered_image);
            combined_image=intersection(combined_image,monotone_image);
            counts.combined_width_sum+=
                combined_image.width().raw().get_d();
            if(not crosses_zero(combined_image)) {
                ++counts.combined_sign_definite;
            }
        }

        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }
    stopwatch.click();

    auto average=[](double sum,SizeType count) {
        return count==0u ? 0.0 : sum/static_cast<double>(count);
    };

    std::cout << "[lie-gradient-shared-procedure-profile]"
              << " time=" << stopwatch.elapsed_seconds()
              << " derivative-build-time="
              << derivative_build_stopwatch.elapsed_seconds()
              << " procedure-build-time="
              << procedure_build_stopwatch.elapsed_seconds()
              << " separate-instructions=" << separate_instruction_count
              << " shared-instructions=" << shared_instruction_count
              << " processed=" << counts.processed
              << " pruned=" << counts.pruned
              << " split=" << counts.split
              << " natural-cross-zero=" << counts.natural_cross_zero
              << " centered-sign-definite="
              << counts.centered_sign_definite
              << " monotone-sign-definite="
              << counts.monotone_sign_definite
              << " combined-sign-definite="
              << counts.combined_sign_definite
              << " monotone-boxes=" << counts.monotone_boxes
              << " monotone-coordinates=" << counts.monotone_coordinates
              << " avg-natural-width="
              << average(counts.natural_width_sum,counts.natural_cross_zero)
              << " avg-centered-width="
              << average(counts.centered_width_sum,counts.natural_cross_zero)
              << " avg-monotone-width="
              << average(counts.monotone_width_sum,counts.natural_cross_zero)
              << " avg-combined-width="
              << average(counts.combined_width_sum,counts.natural_cross_zero)
              << " shared-evals=" << counts.shared_evaluations
              << " shared-time=" << counts.shared_seconds
              << " midpoint-evals=" << counts.midpoint_evaluations
              << " midpoint-time=" << counts.midpoint_seconds
              << " monotone-evals=" << counts.monotone_evaluations
              << " monotone-time=" << counts.monotone_seconds
              << std::endl;
}

struct GradientSharedProcedureCheckCounts {
    SizeType processed = 0u;
    SizeType pruned = 0u;
    SizeType split = 0u;
    SizeType checked = 0u;
    SizeType natural_mismatch = 0u;
    SizeType dx_mismatch = 0u;
    SizeType dy_mismatch = 0u;
    SizeType negative_raw_widths = 0u;
    double separate_dx_width_sum = 0.0;
    double shared_dx_width_sum = 0.0;
    double separate_dy_width_sum = 0.0;
    double shared_dy_width_sum = 0.0;
    double maximum_endpoint_difference = 0.0;
    Bool printed_first_mismatch = false;
};

double raw_interval_width(UpperIntervalType const& interval)
{
    return interval.upper_bound().raw().get_d()
        - interval.lower_bound().raw().get_d();
}

double endpoint_difference(
    UpperIntervalType const& first,
    UpperIntervalType const& second)
{
    double const lower_difference=std::abs(
        first.lower_bound().raw().get_d()
        - second.lower_bound().raw().get_d());
    double const upper_difference=std::abs(
        first.upper_bound().raw().get_d()
        - second.upper_bound().raw().get_d());
    return std::max(lower_difference,upper_difference);
}

Bool same_endpoints(
    UpperIntervalType const& first,
    UpperIntervalType const& second)
{
    return first.lower_bound().raw().get_d()
            == second.lower_bound().raw().get_d()
        && first.upper_bound().raw().get_d()
            == second.upper_bound().raw().get_d();
}

Void profile_lie_gradient_shared_procedure_check(
    SizeType box_limit,
    UpperBoxType const& domain,
    ValidatedScalarMultivariateFunction const& db_dy_function,
    ValidatedScalarMultivariateFunction const& lie_function)
{
    ValidatedScalarMultivariateFunction derivative_x=
        db_dy_function.derivative(0u);
    ValidatedScalarMultivariateFunction derivative_y=
        db_dy_function.derivative(1u);

    ValidatedProcedure scalar_procedure(db_dy_function);
    ValidatedProcedure derivative_x_procedure(derivative_x);
    ValidatedProcedure derivative_y_procedure(derivative_y);
    ValidatedVectorMultivariateFunction shared_function=
        join(db_dy_function,join(derivative_x,derivative_y));
    Vector<ValidatedProcedure> shared_procedure(shared_function);

    std::vector<UpperBoxType> pending;
    pending.push_back(domain);
    GradientSharedProcedureCheckCounts counts;

    Stopwatch<Milliseconds> stopwatch;
    while(not pending.empty() && counts.processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++counts.processed;

        UpperIntervalType lie_image=apply(lie_function,box);
        if(definitely(lie_image.lower_bound()>=0)) {
            ++counts.pruned;
            continue;
        }

        auto children=box.split();
        if(definitely(children.first==children.second)) {
            continue;
        }
        ++counts.split;

        Vector<FloatDPBounds> arguments=bounds_arguments(box);
        UpperIntervalType direct_image=
            make_interval(evaluate(scalar_procedure,arguments));
        if(crosses_zero(direct_image)) {
            ++counts.checked;

            Vector<FloatDPBounds> shared_images=
                evaluate(shared_procedure,arguments);
            UpperIntervalType shared_natural=make_interval(shared_images[0u]);
            UpperIntervalType shared_dx=make_interval(shared_images[1u]);
            UpperIntervalType shared_dy=make_interval(shared_images[2u]);

            UpperIntervalType separate_dx=
                make_interval(evaluate(derivative_x_procedure,arguments));
            UpperIntervalType separate_dy=
                make_interval(evaluate(derivative_y_procedure,arguments));

            Bool const natural_mismatch=
                not same_endpoints(direct_image,shared_natural);
            Bool const dx_mismatch=
                not same_endpoints(separate_dx,shared_dx);
            Bool const dy_mismatch=
                not same_endpoints(separate_dy,shared_dy);

            if(natural_mismatch) { ++counts.natural_mismatch; }
            if(dx_mismatch) { ++counts.dx_mismatch; }
            if(dy_mismatch) { ++counts.dy_mismatch; }

            double const separate_dx_width=raw_interval_width(separate_dx);
            double const shared_dx_width=raw_interval_width(shared_dx);
            double const separate_dy_width=raw_interval_width(separate_dy);
            double const shared_dy_width=raw_interval_width(shared_dy);

            if(separate_dx_width<0.0 || shared_dx_width<0.0
                    || separate_dy_width<0.0 || shared_dy_width<0.0) {
                ++counts.negative_raw_widths;
            }

            counts.separate_dx_width_sum+=separate_dx_width;
            counts.shared_dx_width_sum+=shared_dx_width;
            counts.separate_dy_width_sum+=separate_dy_width;
            counts.shared_dy_width_sum+=shared_dy_width;

            if((natural_mismatch || dx_mismatch || dy_mismatch)
                    && not counts.printed_first_mismatch) {
                counts.printed_first_mismatch=true;
                std::cout << "[lie-gradient-shared-procedure-check-first-mismatch]"
                          << " box=" << box
                          << " natural-separate=" << direct_image
                          << " natural-shared=" << shared_natural
                          << " dx-separate=" << separate_dx
                          << " dx-shared=" << shared_dx
                          << " dy-separate=" << separate_dy
                          << " dy-shared=" << shared_dy
                          << std::endl;
            }

            counts.maximum_endpoint_difference=std::max(
                counts.maximum_endpoint_difference,
                endpoint_difference(direct_image,shared_natural));
            counts.maximum_endpoint_difference=std::max(
                counts.maximum_endpoint_difference,
                endpoint_difference(separate_dx,shared_dx));
            counts.maximum_endpoint_difference=std::max(
                counts.maximum_endpoint_difference,
                endpoint_difference(separate_dy,shared_dy));
        }

        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }
    stopwatch.click();

    auto average=[](double sum,SizeType count) {
        return count==0u ? 0.0 : sum/static_cast<double>(count);
    };

    std::cout << "[lie-gradient-shared-procedure-check]"
              << " time=" << stopwatch.elapsed_seconds()
              << " processed=" << counts.processed
              << " pruned=" << counts.pruned
              << " split=" << counts.split
              << " checked=" << counts.checked
              << " natural-mismatch=" << counts.natural_mismatch
              << " dx-mismatch=" << counts.dx_mismatch
              << " dy-mismatch=" << counts.dy_mismatch
              << " negative-raw-widths=" << counts.negative_raw_widths
              << " avg-separate-dx-width="
              << average(counts.separate_dx_width_sum,counts.checked)
              << " avg-shared-dx-width="
              << average(counts.shared_dx_width_sum,counts.checked)
              << " avg-separate-dy-width="
              << average(counts.separate_dy_width_sum,counts.checked)
              << " avg-shared-dy-width="
              << average(counts.shared_dy_width_sum,counts.checked)
              << " max-endpoint-difference="
              << counts.maximum_endpoint_difference
              << std::endl;
}

Void profile_lie_gradient_expression_cse(
    RealExpression const& db_dy_expression,
    RealVariable const& x,
    RealVariable const& y,
    RealSpace const& space)
{
    Stopwatch<Milliseconds> derivative_stopwatch;
    Vector<RealExpression> raw_expressions({
        db_dy_expression,
        simplify(derivative(db_dy_expression,x)),
        simplify(derivative(db_dy_expression,y))
    });
    derivative_stopwatch.click();

    Stopwatch<Milliseconds> raw_compile_stopwatch;
    Vector<Formula<EffectiveNumber>> raw_formulae=
        make_formula(raw_expressions,space);
    Vector<EffectiveProcedure> raw_procedure(
        space.dimension(),raw_formulae);
    raw_compile_stopwatch.click();

    Vector<RealExpression> cse_expressions(raw_expressions);
    Stopwatch<Milliseconds> cse_stopwatch;
    eliminate_common_subexpressions(cse_expressions);
    cse_stopwatch.click();

    Stopwatch<Milliseconds> cse_compile_stopwatch;
    Vector<Formula<EffectiveNumber>> cse_formulae=
        make_formula(cse_expressions,space);
    Vector<EffectiveProcedure> cse_procedure(
        space.dimension(),cse_formulae);
    cse_compile_stopwatch.click();

    SizeType raw_node_sum=0u;
    SizeType cse_node_sum=0u;
    SizeType raw_pointer_sum=0u;
    SizeType cse_pointer_sum=0u;
    for(SizeType i=0u;i!=raw_expressions.size();++i) {
        raw_node_sum+=count_nodes(raw_expressions[i]);
        cse_node_sum+=count_nodes(cse_expressions[i]);
        raw_pointer_sum+=count_distinct_node_pointers(raw_expressions[i]);
        cse_pointer_sum+=count_distinct_node_pointers(cse_expressions[i]);
    }

    std::cout << "[lie-gradient-expression-cse-profile]"
              << " derivative-build-time="
              << derivative_stopwatch.elapsed_seconds()
              << " cse-time=" << cse_stopwatch.elapsed_seconds()
              << " raw-compile-time="
              << raw_compile_stopwatch.elapsed_seconds()
              << " cse-compile-time="
              << cse_compile_stopwatch.elapsed_seconds()
              << " raw-node-sum=" << raw_node_sum
              << " cse-node-sum=" << cse_node_sum
              << " raw-pointer-sum=" << raw_pointer_sum
              << " cse-pointer-sum=" << cse_pointer_sum
              << " raw-instructions="
              << raw_procedure.temporaries_size()
              << " cse-instructions="
              << cse_procedure.temporaries_size()
              << std::endl;
}

struct GradientExpressionCseFrontierCounts {
    SizeType processed = 0u;
    SizeType pruned = 0u;
    SizeType split = 0u;
    SizeType natural_cross_zero = 0u;
    SizeType centered_sign_definite = 0u;
    SizeType monotone_sign_definite = 0u;
    SizeType combined_sign_definite = 0u;
    SizeType monotone_boxes = 0u;
    SizeType monotone_coordinates = 0u;
    SizeType natural_evaluations = 0u;
    SizeType cse_evaluations = 0u;
    SizeType midpoint_evaluations = 0u;
    SizeType monotone_evaluations = 0u;
    double natural_width_sum = 0.0;
    double centered_width_sum = 0.0;
    double monotone_width_sum = 0.0;
    double combined_width_sum = 0.0;
    double natural_seconds = 0.0;
    double cse_seconds = 0.0;
    double midpoint_seconds = 0.0;
    double monotone_seconds = 0.0;
};

Void profile_lie_gradient_expression_cse_frontier(
    SizeType box_limit,
    UpperBoxType const& domain,
    RealExpression const& db_dy_expression,
    RealVariable const& x,
    RealVariable const& y,
    RealSpace const& space,
    ValidatedScalarMultivariateFunction const& db_dy_function,
    ValidatedScalarMultivariateFunction const& lie_function)
{
    Stopwatch<Milliseconds> derivative_stopwatch;
    Vector<RealExpression> expressions({
        db_dy_expression,
        simplify(derivative(db_dy_expression,x)),
        simplify(derivative(db_dy_expression,y))
    });
    derivative_stopwatch.click();

    Stopwatch<Milliseconds> cse_stopwatch;
    eliminate_common_subexpressions(expressions);
    cse_stopwatch.click();

    Stopwatch<Milliseconds> procedure_build_stopwatch;
    Vector<Formula<EffectiveNumber>> formulae=
        make_formula(expressions,space);
    Vector<EffectiveProcedure> cse_procedure(
        space.dimension(),formulae);
    ValidatedProcedure natural_procedure(db_dy_function);
    EffectiveProcedure scalar_cse_procedure(
        space.dimension(),formulae[0u]);
    procedure_build_stopwatch.click();

    std::vector<UpperBoxType> pending;
    pending.push_back(domain);
    GradientExpressionCseFrontierCounts counts;

    Stopwatch<Milliseconds> stopwatch;
    while(not pending.empty() && counts.processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++counts.processed;

        UpperIntervalType lie_image=apply(lie_function,box);
        if(definitely(lie_image.lower_bound()>=0)) {
            ++counts.pruned;
            continue;
        }

        auto children=box.split();
        if(definitely(children.first==children.second)) {
            continue;
        }
        ++counts.split;

        Vector<FloatDPBounds> arguments=bounds_arguments(box);

        Stopwatch<Milliseconds> natural_stopwatch;
        UpperIntervalType natural_image=
            make_interval(evaluate(natural_procedure,arguments));
        natural_stopwatch.click();
        counts.natural_seconds+=natural_stopwatch.elapsed_seconds();
        ++counts.natural_evaluations;

        if(crosses_zero(natural_image)) {
            ++counts.natural_cross_zero;
            counts.natural_width_sum+=natural_image.width().raw().get_d();

            Stopwatch<Milliseconds> cse_evaluation_stopwatch;
            Vector<FloatDPBounds> images=
                evaluate(cse_procedure,arguments);
            cse_evaluation_stopwatch.click();
            counts.cse_seconds+=cse_evaluation_stopwatch.elapsed_seconds();
            ++counts.cse_evaluations;

            Vector<FloatDPBounds> midpoint_arguments(arguments);
            for(SizeType i=0u;i!=box.dimension();++i) {
                auto midpoint=box[i].midpoint();
                midpoint_arguments[i]=FloatDPBounds(
                    midpoint.raw(),midpoint.raw());
            }

            Stopwatch<Milliseconds> midpoint_stopwatch;
            FloatDPBounds centered_bounds=
                evaluate(scalar_cse_procedure,midpoint_arguments);
            midpoint_stopwatch.click();
            counts.midpoint_seconds+=midpoint_stopwatch.elapsed_seconds();
            ++counts.midpoint_evaluations;

            for(SizeType i=0u;i!=box.dimension();++i) {
                centered_bounds+=images[i+1u]
                    *(arguments[i]-midpoint_arguments[i]);
            }
            UpperIntervalType centered_image=make_interval(centered_bounds);
            counts.centered_width_sum+=
                centered_image.width().raw().get_d();
            if(not crosses_zero(centered_image)) {
                ++counts.centered_sign_definite;
            }

            Vector<FloatDPBounds> lower_arguments(arguments);
            Vector<FloatDPBounds> upper_arguments(arguments);
            SizeType monotone_coordinates=0u;
            for(SizeType i=0u;i!=box.dimension();++i) {
                UpperIntervalType derivative_image=
                    make_interval(images[i+1u]);
                if(definitely(derivative_image.lower_bound()>0)) {
                    auto lower=box[i].lower_bound().raw();
                    auto upper=box[i].upper_bound().raw();
                    lower_arguments[i]=FloatDPBounds(lower,lower);
                    upper_arguments[i]=FloatDPBounds(upper,upper);
                    ++monotone_coordinates;
                } else if(definitely(derivative_image.upper_bound()<0)) {
                    auto lower=box[i].lower_bound().raw();
                    auto upper=box[i].upper_bound().raw();
                    lower_arguments[i]=FloatDPBounds(upper,upper);
                    upper_arguments[i]=FloatDPBounds(lower,lower);
                    ++monotone_coordinates;
                }
            }

            UpperIntervalType monotone_image=natural_image;
            if(monotone_coordinates!=0u) {
                ++counts.monotone_boxes;
                counts.monotone_coordinates+=monotone_coordinates;

                Stopwatch<Milliseconds> monotone_stopwatch;
                FloatDPBounds lower_image=
                    evaluate(scalar_cse_procedure,lower_arguments);
                FloatDPBounds upper_image=
                    evaluate(scalar_cse_procedure,upper_arguments);
                monotone_stopwatch.click();
                counts.monotone_seconds+=monotone_stopwatch.elapsed_seconds();
                counts.monotone_evaluations+=2u;

                monotone_image=make_interval(
                    FloatDPBounds(lower_image.lower(),upper_image.upper()));
            }
            counts.monotone_width_sum+=
                monotone_image.width().raw().get_d();
            if(not crosses_zero(monotone_image)) {
                ++counts.monotone_sign_definite;
            }

            UpperIntervalType combined_image=
                intersection(natural_image,centered_image);
            combined_image=intersection(combined_image,monotone_image);
            counts.combined_width_sum+=
                combined_image.width().raw().get_d();
            if(not crosses_zero(combined_image)) {
                ++counts.combined_sign_definite;
            }
        }

        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }
    stopwatch.click();

    auto average=[](double sum,SizeType count) {
        return count==0u ? 0.0 : sum/static_cast<double>(count);
    };

    std::cout << "[lie-gradient-expression-cse-frontier]"
              << " time=" << stopwatch.elapsed_seconds()
              << " derivative-build-time="
              << derivative_stopwatch.elapsed_seconds()
              << " cse-time=" << cse_stopwatch.elapsed_seconds()
              << " procedure-build-time="
              << procedure_build_stopwatch.elapsed_seconds()
              << " cse-instructions="
              << cse_procedure.temporaries_size()
              << " processed=" << counts.processed
              << " pruned=" << counts.pruned
              << " split=" << counts.split
              << " natural-cross-zero=" << counts.natural_cross_zero
              << " centered-sign-definite="
              << counts.centered_sign_definite
              << " monotone-sign-definite="
              << counts.monotone_sign_definite
              << " combined-sign-definite="
              << counts.combined_sign_definite
              << " monotone-boxes=" << counts.monotone_boxes
              << " monotone-coordinates=" << counts.monotone_coordinates
              << " avg-natural-width="
              << average(counts.natural_width_sum,counts.natural_cross_zero)
              << " avg-centered-width="
              << average(counts.centered_width_sum,counts.natural_cross_zero)
              << " avg-monotone-width="
              << average(counts.monotone_width_sum,counts.natural_cross_zero)
              << " avg-combined-width="
              << average(counts.combined_width_sum,counts.natural_cross_zero)
              << " natural-evals=" << counts.natural_evaluations
              << " natural-time=" << counts.natural_seconds
              << " cse-evals=" << counts.cse_evaluations
              << " cse-eval-time=" << counts.cse_seconds
              << " midpoint-evals=" << counts.midpoint_evaluations
              << " midpoint-time=" << counts.midpoint_seconds
              << " monotone-evals=" << counts.monotone_evaluations
              << " monotone-time=" << counts.monotone_seconds
              << std::endl;
}

struct GradientCsePruneProfileCounts {
    SizeType processed = 0u;
    SizeType natural_pruned = 0u;
    SizeType split = 0u;
    SizeType db_dy_cross_zero = 0u;
    SizeType recomposed_natural_pruned = 0u;
    SizeType improved_pruned = 0u;
    SizeType improved_intersection_pruned = 0u;
    SizeType centered_sign_definite = 0u;
    SizeType monotone_coordinates = 0u;
    SizeType cse_evaluations = 0u;
    SizeType component_evaluations = 0u;
    double natural_lie_width_sum = 0.0;
    double recomposed_natural_width_sum = 0.0;
    double improved_lie_width_sum = 0.0;
    double intersected_lie_width_sum = 0.0;
    double cse_seconds = 0.0;
    double component_seconds = 0.0;
};

Void profile_lie_gradient_cse_pruning(
    SizeType box_limit,
    UpperBoxType const& domain,
    RealExpression const& db_dy_expression,
    RealVariable const& x,
    RealVariable const& y,
    RealSpace const& space,
    ValidatedScalarMultivariateFunction const& db_dx_function,
    ValidatedScalarMultivariateFunction const& db_dy_function,
    ValidatedScalarMultivariateFunction const& barrier_function,
    ValidatedScalarMultivariateFunction const& dy_function,
    ValidatedScalarMultivariateFunction const& lie_function)
{
    Stopwatch<Milliseconds> derivative_stopwatch;
    Vector<RealExpression> expressions({
        db_dy_expression,
        simplify(derivative(db_dy_expression,x)),
        simplify(derivative(db_dy_expression,y))
    });
    derivative_stopwatch.click();

    Stopwatch<Milliseconds> cse_stopwatch;
    eliminate_common_subexpressions(expressions);
    cse_stopwatch.click();

    Stopwatch<Milliseconds> procedure_build_stopwatch;
    Vector<Formula<EffectiveNumber>> formulae=
        make_formula(expressions,space);
    Vector<EffectiveProcedure> cse_procedure(
        space.dimension(),formulae);
    EffectiveProcedure scalar_cse_procedure(
        space.dimension(),formulae[0u]);
    ValidatedProcedure natural_db_dy_procedure(db_dy_function);
    procedure_build_stopwatch.click();

    std::vector<UpperBoxType> pending;
    pending.push_back(domain);
    GradientCsePruneProfileCounts counts;

    Stopwatch<Milliseconds> stopwatch;
    while(not pending.empty() && counts.processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++counts.processed;

        UpperIntervalType natural_lie=apply(lie_function,box);
        if(definitely(natural_lie.lower_bound()>=0)) {
            ++counts.natural_pruned;
            continue;
        }

        auto children=box.split();
        if(definitely(children.first==children.second)) {
            continue;
        }
        ++counts.split;

        Vector<FloatDPBounds> arguments=bounds_arguments(box);
        UpperIntervalType natural_db_dy=
            make_interval(evaluate(natural_db_dy_procedure,arguments));
        if(crosses_zero(natural_db_dy)) {
            ++counts.db_dy_cross_zero;

            Stopwatch<Milliseconds> cse_eval_stopwatch;
            Vector<FloatDPBounds> images=
                evaluate(cse_procedure,arguments);

            Vector<FloatDPBounds> midpoint_arguments(arguments);
            for(SizeType i=0u;i!=box.dimension();++i) {
                auto midpoint=box[i].midpoint();
                midpoint_arguments[i]=FloatDPBounds(
                    midpoint.raw(),midpoint.raw());
            }
            FloatDPBounds centered_bounds=
                evaluate(scalar_cse_procedure,midpoint_arguments);
            for(SizeType i=0u;i!=box.dimension();++i) {
                centered_bounds+=images[i+1u]
                    *(arguments[i]-midpoint_arguments[i]);
            }
            UpperIntervalType centered_image=make_interval(centered_bounds);
            if(not crosses_zero(centered_image)) {
                ++counts.centered_sign_definite;
            }

            UpperIntervalType improved_db_dy=
                intersection(natural_db_dy,centered_image);

            Vector<FloatDPBounds> lower_arguments(arguments);
            Vector<FloatDPBounds> upper_arguments(arguments);
            SizeType monotone_coordinates=0u;
            for(SizeType i=0u;i!=box.dimension();++i) {
                UpperIntervalType derivative_image=
                    make_interval(images[i+1u]);
                if(definitely(derivative_image.lower_bound()>0)) {
                    auto lower=box[i].lower_bound().raw();
                    auto upper=box[i].upper_bound().raw();
                    lower_arguments[i]=FloatDPBounds(lower,lower);
                    upper_arguments[i]=FloatDPBounds(upper,upper);
                    ++monotone_coordinates;
                } else if(definitely(derivative_image.upper_bound()<0)) {
                    auto lower=box[i].lower_bound().raw();
                    auto upper=box[i].upper_bound().raw();
                    lower_arguments[i]=FloatDPBounds(upper,upper);
                    upper_arguments[i]=FloatDPBounds(lower,lower);
                    ++monotone_coordinates;
                }
            }
            counts.monotone_coordinates+=monotone_coordinates;
            if(monotone_coordinates!=0u) {
                FloatDPBounds lower_image=
                    evaluate(scalar_cse_procedure,lower_arguments);
                FloatDPBounds upper_image=
                    evaluate(scalar_cse_procedure,upper_arguments);
                UpperIntervalType monotone_image=make_interval(
                    FloatDPBounds(lower_image.lower(),upper_image.upper()));
                improved_db_dy=
                    intersection(improved_db_dy,monotone_image);
            }
            cse_eval_stopwatch.click();
            counts.cse_seconds+=cse_eval_stopwatch.elapsed_seconds();
            ++counts.cse_evaluations;

            Stopwatch<Milliseconds> component_stopwatch;
            UpperIntervalType db_dx_image=apply(db_dx_function,box);
            UpperIntervalType barrier_image=apply(barrier_function,box);
            UpperIntervalType dy_image=apply(dy_function,box);
            UpperIntervalType dx_image=box[1u];

            UpperIntervalType recomposed_natural=
                db_dx_image*dx_image
                + natural_db_dy*dy_image
                + barrier_image;
            UpperIntervalType improved_lie=
                db_dx_image*dx_image
                + improved_db_dy*dy_image
                + barrier_image;
            UpperIntervalType intersected_lie=
                intersection(natural_lie,improved_lie);
            component_stopwatch.click();
            counts.component_seconds+=component_stopwatch.elapsed_seconds();
            counts.component_evaluations+=3u;

            if(definitely(recomposed_natural.lower_bound()>=0)) {
                ++counts.recomposed_natural_pruned;
            }
            if(definitely(improved_lie.lower_bound()>=0)) {
                ++counts.improved_pruned;
            }
            if(definitely(intersected_lie.lower_bound()>=0)) {
                ++counts.improved_intersection_pruned;
            }

            counts.natural_lie_width_sum+=
                natural_lie.width().raw().get_d();
            counts.recomposed_natural_width_sum+=
                recomposed_natural.width().raw().get_d();
            counts.improved_lie_width_sum+=
                improved_lie.width().raw().get_d();
            counts.intersected_lie_width_sum+=
                intersected_lie.width().raw().get_d();
        }

        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }
    stopwatch.click();

    auto average=[](double sum,SizeType count) {
        return count==0u ? 0.0 : sum/static_cast<double>(count);
    };

    std::cout << "[lie-gradient-cse-prune-profile|lie-root-cse-range-profile|lie-root-cse-search|lie-algebraic-form-profile|lie-directional-propagation-profile|lie-layer2-directional-profile|lie-gradient-block-reassociation-profile|lie-gradient-width-attribution-profile|lie-gradient-correlation-attribution-profile]"
              << " time=" << stopwatch.elapsed_seconds()
              << " derivative-build-time="
              << derivative_stopwatch.elapsed_seconds()
              << " cse-time=" << cse_stopwatch.elapsed_seconds()
              << " procedure-build-time="
              << procedure_build_stopwatch.elapsed_seconds()
              << " cse-instructions="
              << cse_procedure.temporaries_size()
              << " processed=" << counts.processed
              << " natural-pruned=" << counts.natural_pruned
              << " split=" << counts.split
              << " db/dy-cross-zero=" << counts.db_dy_cross_zero
              << " centered-sign-definite="
              << counts.centered_sign_definite
              << " monotone-coordinates="
              << counts.monotone_coordinates
              << " recomposed-natural-pruned="
              << counts.recomposed_natural_pruned
              << " improved-pruned=" << counts.improved_pruned
              << " improved-intersection-pruned="
              << counts.improved_intersection_pruned
              << " avg-natural-lie-width="
              << average(counts.natural_lie_width_sum,counts.db_dy_cross_zero)
              << " avg-recomposed-natural-width="
              << average(counts.recomposed_natural_width_sum,counts.db_dy_cross_zero)
              << " avg-improved-lie-width="
              << average(counts.improved_lie_width_sum,counts.db_dy_cross_zero)
              << " avg-intersected-lie-width="
              << average(counts.intersected_lie_width_sum,counts.db_dy_cross_zero)
              << " cse-evals=" << counts.cse_evaluations
              << " cse-eval-time=" << counts.cse_seconds
              << " component-evals=" << counts.component_evaluations
              << " component-time=" << counts.component_seconds
              << std::endl;
}

struct RootCseRangeProfileCounts {
    SizeType processed = 0u;
    SizeType natural_pruned = 0u;
    SizeType split = 0u;
    SizeType centered_pruned = 0u;
    SizeType monotone_pruned = 0u;
    SizeType combined_pruned = 0u;
    SizeType monotone_boxes = 0u;
    SizeType monotone_coordinates = 0u;
    SizeType cse_evaluations = 0u;
    SizeType midpoint_evaluations = 0u;
    SizeType monotone_evaluations = 0u;
    double natural_width_sum = 0.0;
    double centered_width_sum = 0.0;
    double monotone_width_sum = 0.0;
    double combined_width_sum = 0.0;
    double cse_seconds = 0.0;
    double midpoint_seconds = 0.0;
    double monotone_seconds = 0.0;
};

Void profile_lie_root_cse_range(
    SizeType box_limit,
    UpperBoxType const& domain,
    RealExpression const& expression,
    RealVariable const& x,
    RealVariable const& y,
    RealSpace const& space,
    ValidatedScalarMultivariateFunction const& function)
{
    Stopwatch<Milliseconds> derivative_stopwatch;
    Vector<RealExpression> expressions({
        expression,
        simplify(derivative(expression,x)),
        simplify(derivative(expression,y))
    });
    derivative_stopwatch.click();

    Stopwatch<Milliseconds> cse_stopwatch;
    eliminate_common_subexpressions(expressions);
    cse_stopwatch.click();

    Stopwatch<Milliseconds> procedure_build_stopwatch;
    Vector<Formula<EffectiveNumber>> formulae=
        make_formula(expressions,space);
    Vector<EffectiveProcedure> cse_procedure(
        space.dimension(),formulae);
    EffectiveProcedure scalar_procedure(
        space.dimension(),formulae[0u]);
    procedure_build_stopwatch.click();

    std::vector<UpperBoxType> pending;
    pending.push_back(domain);
    RootCseRangeProfileCounts counts;

    Stopwatch<Milliseconds> stopwatch;
    while(not pending.empty() && counts.processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++counts.processed;

        UpperIntervalType natural_image=apply(function,box);
        if(definitely(natural_image.lower_bound()>=0)) {
            ++counts.natural_pruned;
            continue;
        }

        auto children=box.split();
        if(definitely(children.first==children.second)) {
            continue;
        }
        ++counts.split;
        counts.natural_width_sum+=natural_image.width().raw().get_d();

        Vector<FloatDPBounds> arguments=bounds_arguments(box);
        Stopwatch<Milliseconds> cse_eval_stopwatch;
        Vector<FloatDPBounds> images=
            evaluate(cse_procedure,arguments);
        cse_eval_stopwatch.click();
        counts.cse_seconds+=cse_eval_stopwatch.elapsed_seconds();
        ++counts.cse_evaluations;

        Vector<FloatDPBounds> midpoint_arguments(arguments);
        for(SizeType i=0u;i!=box.dimension();++i) {
            auto midpoint=box[i].midpoint();
            midpoint_arguments[i]=FloatDPBounds(
                midpoint.raw(),midpoint.raw());
        }

        Stopwatch<Milliseconds> midpoint_stopwatch;
        FloatDPBounds centered_bounds=
            evaluate(scalar_procedure,midpoint_arguments);
        midpoint_stopwatch.click();
        counts.midpoint_seconds+=midpoint_stopwatch.elapsed_seconds();
        ++counts.midpoint_evaluations;

        for(SizeType i=0u;i!=box.dimension();++i) {
            centered_bounds+=images[i+1u]
                *(arguments[i]-midpoint_arguments[i]);
        }
        UpperIntervalType centered_image=make_interval(centered_bounds);
        counts.centered_width_sum+=centered_image.width().raw().get_d();
        if(definitely(centered_image.lower_bound()>=0)) {
            ++counts.centered_pruned;
        }

        Vector<FloatDPBounds> lower_arguments(arguments);
        Vector<FloatDPBounds> upper_arguments(arguments);
        SizeType monotone_coordinates=0u;
        for(SizeType i=0u;i!=box.dimension();++i) {
            UpperIntervalType derivative_image=
                make_interval(images[i+1u]);
            if(definitely(derivative_image.lower_bound()>0)) {
                auto lower=box[i].lower_bound().raw();
                auto upper=box[i].upper_bound().raw();
                lower_arguments[i]=FloatDPBounds(lower,lower);
                upper_arguments[i]=FloatDPBounds(upper,upper);
                ++monotone_coordinates;
            } else if(definitely(derivative_image.upper_bound()<0)) {
                auto lower=box[i].lower_bound().raw();
                auto upper=box[i].upper_bound().raw();
                lower_arguments[i]=FloatDPBounds(upper,upper);
                upper_arguments[i]=FloatDPBounds(lower,lower);
                ++monotone_coordinates;
            }
        }

        UpperIntervalType monotone_image=natural_image;
        if(monotone_coordinates!=0u) {
            ++counts.monotone_boxes;
            counts.monotone_coordinates+=monotone_coordinates;

            Stopwatch<Milliseconds> monotone_stopwatch;
            FloatDPBounds lower_image=
                evaluate(scalar_procedure,lower_arguments);
            FloatDPBounds upper_image=
                evaluate(scalar_procedure,upper_arguments);
            monotone_stopwatch.click();
            counts.monotone_seconds+=monotone_stopwatch.elapsed_seconds();
            counts.monotone_evaluations+=2u;

            monotone_image=make_interval(
                FloatDPBounds(lower_image.lower(),upper_image.upper()));
        }
        counts.monotone_width_sum+=monotone_image.width().raw().get_d();
        if(definitely(monotone_image.lower_bound()>=0)) {
            ++counts.monotone_pruned;
        }

        UpperIntervalType combined_image=
            intersection(natural_image,centered_image);
        combined_image=intersection(combined_image,monotone_image);
        counts.combined_width_sum+=combined_image.width().raw().get_d();
        if(definitely(combined_image.lower_bound()>=0)) {
            ++counts.combined_pruned;
        }

        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }
    stopwatch.click();

    auto average=[](double sum,SizeType count) {
        return count==0u ? 0.0 : sum/static_cast<double>(count);
    };

    std::cout << "[lie-root-cse-range-profile|lie-root-cse-search|lie-algebraic-form-profile|lie-directional-propagation-profile|lie-layer2-directional-profile|lie-gradient-block-reassociation-profile|lie-gradient-width-attribution-profile|lie-gradient-correlation-attribution-profile]"
              << " time=" << stopwatch.elapsed_seconds()
              << " derivative-build-time="
              << derivative_stopwatch.elapsed_seconds()
              << " cse-time=" << cse_stopwatch.elapsed_seconds()
              << " procedure-build-time="
              << procedure_build_stopwatch.elapsed_seconds()
              << " cse-instructions="
              << cse_procedure.temporaries_size()
              << " processed=" << counts.processed
              << " natural-pruned=" << counts.natural_pruned
              << " split=" << counts.split
              << " centered-pruned=" << counts.centered_pruned
              << " monotone-pruned=" << counts.monotone_pruned
              << " combined-pruned=" << counts.combined_pruned
              << " monotone-boxes=" << counts.monotone_boxes
              << " monotone-coordinates=" << counts.monotone_coordinates
              << " avg-natural-width="
              << average(counts.natural_width_sum,counts.split)
              << " avg-centered-width="
              << average(counts.centered_width_sum,counts.split)
              << " avg-monotone-width="
              << average(counts.monotone_width_sum,counts.split)
              << " avg-combined-width="
              << average(counts.combined_width_sum,counts.split)
              << " cse-evals=" << counts.cse_evaluations
              << " cse-eval-time=" << counts.cse_seconds
              << " midpoint-evals=" << counts.midpoint_evaluations
              << " midpoint-time=" << counts.midpoint_seconds
              << " monotone-evals=" << counts.monotone_evaluations
              << " monotone-time=" << counts.monotone_seconds
              << std::endl;
}

struct RootCseSearchCounts {
    SizeType processed = 0u;
    SizeType natural_lie_pruned = 0u;
    SizeType cse_lie_pruned = 0u;
    SizeType barrier_pruned = 0u;
    SizeType split = 0u;
    SizeType epsilon_certified = 0u;
    SizeType natural_lie_evaluations = 0u;
    SizeType cse_evaluations = 0u;
    SizeType midpoint_evaluations = 0u;
    SizeType barrier_evaluations = 0u;
    double natural_lie_seconds = 0.0;
    double cse_seconds = 0.0;
    double midpoint_seconds = 0.0;
    double barrier_seconds = 0.0;
};

Void profile_lie_root_cse_search(
    SizeType box_limit,
    UpperBoxType const& domain,
    RealExpression const& lie_expression,
    RealVariable const& x,
    RealVariable const& y,
    RealSpace const& space,
    ValidatedScalarMultivariateFunction const& lie_function,
    ValidatedScalarMultivariateFunction const& barrier_function,
    String const& literal_order,
    Bool upper_child_first)
{
    Stopwatch<Milliseconds> derivative_stopwatch;
    Vector<RealExpression> expressions({
        lie_expression,
        simplify(derivative(lie_expression,x)),
        simplify(derivative(lie_expression,y))
    });
    derivative_stopwatch.click();

    Stopwatch<Milliseconds> cse_stopwatch;
    eliminate_common_subexpressions(expressions);
    cse_stopwatch.click();

    Stopwatch<Milliseconds> procedure_build_stopwatch;
    Vector<Formula<EffectiveNumber>> formulae=
        make_formula(expressions,space);
    Vector<EffectiveProcedure> cse_procedure(
        space.dimension(),formulae);
    EffectiveProcedure scalar_procedure(
        space.dimension(),formulae[0u]);
    procedure_build_stopwatch.click();

    FloatDP epsilon(0.00001_x,dp);
    std::vector<UpperBoxType> pending;
    pending.push_back(domain);
    RootCseSearchCounts counts;

    Stopwatch<Milliseconds> search_stopwatch;
    while(not pending.empty() && counts.processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++counts.processed;

        Bool pruned=false;
        Bool lie_epsilon_satisfied=false;
        Bool barrier_epsilon_satisfied=false;

        auto classify_lie=[&]() {
            Stopwatch<Milliseconds> natural_stopwatch;
            UpperIntervalType natural_image=apply(lie_function,box);
            natural_stopwatch.click();
            counts.natural_lie_seconds+=natural_stopwatch.elapsed_seconds();
            ++counts.natural_lie_evaluations;

            if(definitely(natural_image.lower_bound()>=0)) {
                ++counts.natural_lie_pruned;
                pruned=true;
                return;
            }

            Vector<FloatDPBounds> arguments=bounds_arguments(box);
            Stopwatch<Milliseconds> cse_eval_stopwatch;
            Vector<FloatDPBounds> images=
                evaluate(cse_procedure,arguments);
            cse_eval_stopwatch.click();
            counts.cse_seconds+=cse_eval_stopwatch.elapsed_seconds();
            ++counts.cse_evaluations;

            Vector<FloatDPBounds> midpoint_arguments(arguments);
            for(SizeType i=0u;i!=box.dimension();++i) {
                auto midpoint=box[i].midpoint();
                midpoint_arguments[i]=FloatDPBounds(
                    midpoint.raw(),midpoint.raw());
            }
            Stopwatch<Milliseconds> midpoint_stopwatch;
            FloatDPBounds centered_bounds=
                evaluate(scalar_procedure,midpoint_arguments);
            midpoint_stopwatch.click();
            counts.midpoint_seconds+=midpoint_stopwatch.elapsed_seconds();
            ++counts.midpoint_evaluations;

            for(SizeType i=0u;i!=box.dimension();++i) {
                centered_bounds+=images[i+1u]
                    *(arguments[i]-midpoint_arguments[i]);
            }
            UpperIntervalType centered_image=make_interval(centered_bounds);
            UpperIntervalType improved_image=
                intersection(natural_image,centered_image);

            if(definitely(improved_image.lower_bound()>=0)) {
                ++counts.cse_lie_pruned;
                pruned=true;
                return;
            }
            lie_epsilon_satisfied=
                definitely(improved_image.upper_bound()<epsilon);
        };

        auto classify_barrier=[&]() {
            Stopwatch<Milliseconds> barrier_stopwatch;
            UpperIntervalType image=apply(barrier_function,box);
            barrier_stopwatch.click();
            counts.barrier_seconds+=barrier_stopwatch.elapsed_seconds();
            ++counts.barrier_evaluations;
            if(definitely(image.upper_bound()<0)) {
                ++counts.barrier_pruned;
                pruned=true;
                return;
            }
            auto relaxed_lower=sub(down,FloatDP(0,dp),epsilon);
            barrier_epsilon_satisfied=
                definitely(image.lower_bound()>=relaxed_lower);
        };

        if(literal_order=="lie-first") {
            classify_lie();
            if(not pruned) { classify_barrier(); }
        } else {
            classify_barrier();
            if(not pruned) { classify_lie(); }
        }

        if(pruned) {
            continue;
        }
        if(lie_epsilon_satisfied && barrier_epsilon_satisfied) {
            ++counts.epsilon_certified;
            continue;
        }

        auto children=box.split();
        if(definitely(children.first==children.second)) {
            continue;
        }
        ++counts.split;
        if(upper_child_first) {
            pending.push_back(std::move(children.first));
            pending.push_back(std::move(children.second));
        } else {
            pending.push_back(std::move(children.second));
            pending.push_back(std::move(children.first));
        }
    }
    search_stopwatch.click();

    std::cout << "[lie-root-cse-search|lie-algebraic-form-profile|lie-directional-propagation-profile|lie-layer2-directional-profile|lie-gradient-block-reassociation-profile|lie-gradient-width-attribution-profile|lie-gradient-correlation-attribution-profile]"
              << " search-time=" << search_stopwatch.elapsed_seconds()
              << " derivative-build-time="
              << derivative_stopwatch.elapsed_seconds()
              << " cse-time=" << cse_stopwatch.elapsed_seconds()
              << " procedure-build-time="
              << procedure_build_stopwatch.elapsed_seconds()
              << " cse-instructions="
              << cse_procedure.temporaries_size()
              << " processed=" << counts.processed
              << " natural-lie-pruned=" << counts.natural_lie_pruned
              << " cse-lie-pruned=" << counts.cse_lie_pruned
              << " barrier-pruned=" << counts.barrier_pruned
              << " total-pruned="
              << counts.natural_lie_pruned+counts.cse_lie_pruned+counts.barrier_pruned
              << " split=" << counts.split
              << " epsilon-certified=" << counts.epsilon_certified
              << " natural-lie-evals=" << counts.natural_lie_evaluations
              << " natural-lie-time=" << counts.natural_lie_seconds
              << " cse-evals=" << counts.cse_evaluations
              << " cse-eval-time=" << counts.cse_seconds
              << " midpoint-evals=" << counts.midpoint_evaluations
              << " midpoint-time=" << counts.midpoint_seconds
              << " barrier-evals=" << counts.barrier_evaluations
              << " barrier-time=" << counts.barrier_seconds
              << std::endl;
}

struct LieAlgebraicFormCounts {
    SizeType processed = 0u;
    SizeType baseline_pruned = 0u;
    SizeType expanded_pruned = 0u;
    SizeType grouped_y_pruned = 0u;
    SizeType split = 0u;
    double baseline_width_sum = 0.0;
    double expanded_width_sum = 0.0;
    double grouped_y_width_sum = 0.0;
    double baseline_seconds = 0.0;
    double expanded_seconds = 0.0;
    double grouped_y_seconds = 0.0;
};

Void profile_lie_algebraic_forms(
    SizeType box_limit,
    UpperBoxType const& domain,
    ValidatedScalarMultivariateFunction const& baseline_function,
    ValidatedScalarMultivariateFunction const& expanded_function,
    ValidatedScalarMultivariateFunction const& grouped_y_function)
{
    UpperIntervalType baseline_root=apply(baseline_function,domain);
    UpperIntervalType expanded_root=apply(expanded_function,domain);
    UpperIntervalType grouped_y_root=apply(grouped_y_function,domain);

    std::cout << "[lie-algebraic-form-root]"
              << " baseline=" << baseline_root
              << " baseline-width=" << baseline_root.width().raw().get_d()
              << " expanded=" << expanded_root
              << " expanded-width=" << expanded_root.width().raw().get_d()
              << " grouped-y=" << grouped_y_root
              << " grouped-y-width=" << grouped_y_root.width().raw().get_d()
              << std::endl;

    std::vector<UpperBoxType> pending;
    pending.push_back(domain);
    LieAlgebraicFormCounts counts;

    Stopwatch<Milliseconds> stopwatch;
    while(not pending.empty() && counts.processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++counts.processed;

        Stopwatch<Milliseconds> baseline_stopwatch;
        UpperIntervalType baseline_image=apply(baseline_function,box);
        baseline_stopwatch.click();
        counts.baseline_seconds+=baseline_stopwatch.elapsed_seconds();

        Stopwatch<Milliseconds> expanded_stopwatch;
        UpperIntervalType expanded_image=apply(expanded_function,box);
        expanded_stopwatch.click();
        counts.expanded_seconds+=expanded_stopwatch.elapsed_seconds();

        Stopwatch<Milliseconds> grouped_y_stopwatch;
        UpperIntervalType grouped_y_image=apply(grouped_y_function,box);
        grouped_y_stopwatch.click();
        counts.grouped_y_seconds+=grouped_y_stopwatch.elapsed_seconds();

        if(definitely(baseline_image.lower_bound()>=0)) {
            ++counts.baseline_pruned;
        }
        if(definitely(expanded_image.lower_bound()>=0)) {
            ++counts.expanded_pruned;
        }
        if(definitely(grouped_y_image.lower_bound()>=0)) {
            ++counts.grouped_y_pruned;
        }

        counts.baseline_width_sum+=baseline_image.width().raw().get_d();
        counts.expanded_width_sum+=expanded_image.width().raw().get_d();
        counts.grouped_y_width_sum+=grouped_y_image.width().raw().get_d();

        if(definitely(baseline_image.lower_bound()>=0)) {
            continue;
        }

        auto children=box.split();
        if(definitely(children.first==children.second)) {
            continue;
        }
        ++counts.split;
        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }
    stopwatch.click();

    auto average=[](double sum,SizeType count) {
        return count==0u ? 0.0 : sum/static_cast<double>(count);
    };

    std::cout << "[lie-algebraic-form-profile]"
              << " time=" << stopwatch.elapsed_seconds()
              << " processed=" << counts.processed
              << " baseline-pruned=" << counts.baseline_pruned
              << " expanded-pruned=" << counts.expanded_pruned
              << " grouped-y-pruned=" << counts.grouped_y_pruned
              << " split=" << counts.split
              << " avg-baseline-width="
              << average(counts.baseline_width_sum,counts.processed)
              << " avg-expanded-width="
              << average(counts.expanded_width_sum,counts.processed)
              << " avg-grouped-y-width="
              << average(counts.grouped_y_width_sum,counts.processed)
              << " baseline-time=" << counts.baseline_seconds
              << " expanded-time=" << counts.expanded_seconds
              << " grouped-y-time=" << counts.grouped_y_seconds
              << std::endl;
}





Vector<RealExpression> gradient_width_attribution_expressions(
    RealExpression const& x0,
    RealExpression const& x1)
{
    auto const& p=TestBarr3Full64::parameters();
    constexpr SizeType width=TestBarr3Full64::width;
    auto c=[&](SizeType i) { return TestBarr3Full64::constant(p[i]); };

    std::array<RealExpression,width> h1;
    std::array<RealExpression,width> a1;
    std::array<RealExpression,width> dh1_dy;
    for(SizeType j=0u;j!=width;++j) {
        RealExpression z=c(TestBarr3Full64::b1_offset+j)
            + c(TestBarr3Full64::w1_offset+2u*j)*x0
            + c(TestBarr3Full64::w1_offset+2u*j+1u)*x1;
        h1[j]=tanh(z);
        a1[j]=1-sqr(h1[j]);
        dh1_dy[j]=
            a1[j]*c(TestBarr3Full64::w1_offset+2u*j+1u);
    }

    Vector<RealExpression> expressions(3u*width+1u);
    RealExpression db_dy=RealExpression(0);
    for(SizeType i=0u;i!=width;++i) {
        RealExpression z=c(TestBarr3Full64::b2_offset+i);
        RealExpression dz_dy=RealExpression(0);
        for(SizeType j=0u;j!=width;++j) {
            RealExpression w=c(TestBarr3Full64::w2_offset+width*i+j);
            z=z+w*h1[j];
            dz_dy=dz_dy+w*dh1_dy[j];
        }
        RealExpression h2=tanh(z);
        RealExpression a2=1-sqr(h2);
        RealExpression term=
            c(TestBarr3Full64::w3_offset+i)*a2*dz_dy;

        expressions[i]=a2;
        expressions[width+i]=dz_dy;
        expressions[2u*width+i]=term;
        db_dy=db_dy+term;
    }
    expressions[3u*width]=db_dy;
    return expressions;
}

Void profile_lie_gradient_width_attribution(
    SizeType box_limit,
    UpperBoxType const& domain,
    RealSpace const& space,
    RealExpression const& lie_expression,
    Vector<RealExpression> const& expressions)
{
    constexpr SizeType width=TestBarr3Full64::width;
    auto const& p=TestBarr3Full64::parameters();

    Stopwatch<Milliseconds> build_stopwatch;
    Vector<Formula<EffectiveNumber>> formulae=
        make_formula(expressions,space);
    Vector<EffectiveProcedure> procedure(
        space.dimension(),formulae);
    ValidatedProcedure lie_procedure(
        make_function(space,lie_expression));
    build_stopwatch.click();

    std::array<double,width> a2_width_sum{};
    std::array<double,width> dz_width_sum{};
    std::array<double,width> term_width_sum{};
    std::array<SizeType,width> a2_touches_zero{};
    std::array<SizeType,width> dz_cross_zero{};
    std::array<SizeType,width> term_cross_zero{};

    SizeType processed=0u;
    SizeType pruned=0u;
    SizeType inspected=0u;
    SizeType split=0u;
    double db_dy_width_sum=0.0;
    double term_width_total=0.0;
    double evaluation_seconds=0.0;

    std::vector<UpperBoxType> pending;
    pending.push_back(domain);

    Stopwatch<Milliseconds> total_stopwatch;
    while(not pending.empty() && processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++processed;

        auto arguments=bounds_arguments(box);
        UpperIntervalType lie_image=
            make_interval(evaluate(lie_procedure,arguments));
        if(definitely(lie_image.lower_bound()>=0)) {
            ++pruned;
            continue;
        }

        auto children=box.split();
        if(definitely(children.first==children.second)) {
            continue;
        }
        ++split;
        ++inspected;

        Stopwatch<Milliseconds> evaluation_stopwatch;
        Vector<FloatDPBounds> images=evaluate(procedure,arguments);
        evaluation_stopwatch.click();
        evaluation_seconds+=evaluation_stopwatch.elapsed_seconds();

        UpperIntervalType db_dy_image=
            make_interval(images[3u*width]);
        db_dy_width_sum+=db_dy_image.width().raw().get_d();

        for(SizeType i=0u;i!=width;++i) {
            UpperIntervalType a2=make_interval(images[i]);
            UpperIntervalType dz=make_interval(images[width+i]);
            UpperIntervalType term=make_interval(images[2u*width+i]);

            double const a2_width=a2.width().raw().get_d();
            double const dz_width=dz.width().raw().get_d();
            double const term_width=term.width().raw().get_d();

            a2_width_sum[i]+=a2_width;
            dz_width_sum[i]+=dz_width;
            term_width_sum[i]+=term_width;
            term_width_total+=term_width;

            if(definitely(a2.lower_bound()<=0)
                    && definitely(a2.upper_bound()>=0)) {
                ++a2_touches_zero[i];
            }
            if(crosses_zero(dz)) {
                ++dz_cross_zero[i];
            }
            if(crosses_zero(term)) {
                ++term_cross_zero[i];
            }
        }

        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }
    total_stopwatch.click();

    std::array<SizeType,width> order{};
    for(SizeType i=0u;i!=width;++i) { order[i]=i; }
    std::sort(
        order.begin(),order.end(),
        [&](SizeType lhs,SizeType rhs) {
            return term_width_sum[lhs]>term_width_sum[rhs];
        });

    auto average=[&](double sum) {
        return inspected==0u ? 0.0 : sum/static_cast<double>(inspected);
    };

    std::cout << "[lie-gradient-width-attribution-profile]"
              << " time=" << total_stopwatch.elapsed_seconds()
              << " build-time=" << build_stopwatch.elapsed_seconds()
              << " evaluation-time=" << evaluation_seconds
              << " procedure-instructions=" << procedure.temporaries_size()
              << " processed=" << processed
              << " pruned=" << pruned
              << " split=" << split
              << " inspected=" << inspected
              << " avg-db-dy-width=" << average(db_dy_width_sum)
              << " avg-sum-term-width=" << average(term_width_total)
              << std::endl;

    double cumulative=0.0;
    SizeType const report_count=std::min(SizeType(12u),width);
    for(SizeType rank=0u;rank!=report_count;++rank) {
        SizeType const i=order[rank];
        cumulative+=term_width_sum[i];
        double const share=
            term_width_total==0.0 ? 0.0 : term_width_sum[i]/term_width_total;
        double const cumulative_share=
            term_width_total==0.0 ? 0.0 : cumulative/term_width_total;
        std::cout << "[lie-gradient-width-attribution-neuron]"
                  << " rank=" << rank+1u
                  << " neuron=" << i
                  << " w3=" << p[TestBarr3Full64::w3_offset+i]
                  << " avg-term-width=" << average(term_width_sum[i])
                  << " width-share=" << share
                  << " cumulative-share=" << cumulative_share
                  << " avg-dz-width=" << average(dz_width_sum[i])
                  << " dz-cross-zero=" << dz_cross_zero[i]
                  << " avg-a2-width=" << average(a2_width_sum[i])
                  << " a2-touches-zero=" << a2_touches_zero[i]
                  << " term-cross-zero=" << term_cross_zero[i]
                  << std::endl;
    }
}


Void profile_lie_gradient_correlation_attribution(
    SizeType box_limit,
    UpperBoxType const& domain,
    RealSpace const& space,
    RealExpression const& lie_expression,
    Vector<RealExpression> const& expressions)
{
    constexpr SizeType width=TestBarr3Full64::width;
    Vector<Formula<EffectiveNumber>> formulae=make_formula(expressions,space);
    Vector<EffectiveProcedure> procedure(space.dimension(),formulae);
    ValidatedProcedure lie_procedure(make_function(space,lie_expression));

    SizeType processed=0u;
    SizeType pruned=0u;
    SizeType split=0u;
    SizeType inspected=0u;
    double direct_db_dy_width_sum=0.0;
    double term_quadrant_hull_width_sum=0.0;
    double db_dy_quadrant_hull_width_sum=0.0;
    double evaluation_seconds=0.0;

    std::vector<UpperBoxType> pending;
    pending.push_back(domain);

    Stopwatch<Milliseconds> total_stopwatch;
    while(not pending.empty() && processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++processed;

        Vector<FloatDPBounds> arguments=bounds_arguments(box);
        UpperIntervalType lie_image=
            make_interval(evaluate(lie_procedure,arguments));
        if(definitely(lie_image.lower_bound()>=0)) {
            ++pruned;
            continue;
        }

        auto children_x=box.split(0u);
        if(definitely(children_x.first==children_x.second)) {
            continue;
        }
        auto low_y=children_x.first.split(1u);
        auto high_y=children_x.second.split(1u);
        std::array<UpperBoxType,4u> quadrants={
            low_y.first,low_y.second,high_y.first,high_y.second
        };

        ++split;
        ++inspected;

        Stopwatch<Milliseconds> evaluation_stopwatch;
        Vector<FloatDPBounds> direct_images=evaluate(procedure,arguments);
        UpperIntervalType direct_db_dy=
            make_interval(direct_images[3u*width]);

        std::array<UpperIntervalType,width> term_hulls;
        Bool first_quadrant=true;
        UpperIntervalType db_dy_hull;

        for(auto const& quadrant:quadrants) {
            Vector<FloatDPBounds> qargs=bounds_arguments(quadrant);
            Vector<FloatDPBounds> qimages=evaluate(procedure,qargs);
            UpperIntervalType q_db_dy=make_interval(qimages[3u*width]);

            if(first_quadrant) {
                db_dy_hull=q_db_dy;
                for(SizeType i=0u;i!=width;++i) {
                    term_hulls[i]=make_interval(qimages[2u*width+i]);
                }
                first_quadrant=false;
            } else {
                db_dy_hull=hull(db_dy_hull,q_db_dy);
                for(SizeType i=0u;i!=width;++i) {
                    term_hulls[i]=hull(
                        term_hulls[i],
                        make_interval(qimages[2u*width+i]));
                }
            }
        }
        evaluation_stopwatch.click();
        evaluation_seconds+=evaluation_stopwatch.elapsed_seconds();

        double term_hull_width=0.0;
        for(SizeType i=0u;i!=width;++i) {
            term_hull_width+=term_hulls[i].width().raw().get_d();
        }

        direct_db_dy_width_sum+=direct_db_dy.width().raw().get_d();
        term_quadrant_hull_width_sum+=term_hull_width;
        db_dy_quadrant_hull_width_sum+=db_dy_hull.width().raw().get_d();

        pending.push_back(std::move(children_x.second));
        pending.push_back(std::move(children_x.first));
    }
    total_stopwatch.click();

    auto average=[&](double sum) {
        return inspected==0u ? 0.0 : sum/static_cast<double>(inspected);
    };
    double const direct=average(direct_db_dy_width_sum);
    double const local=average(term_quadrant_hull_width_sum);
    double const joint=average(db_dy_quadrant_hull_width_sum);

    std::cout << "[lie-gradient-correlation-attribution-profile]"
              << " time=" << total_stopwatch.elapsed_seconds()
              << " evaluation-time=" << evaluation_seconds
              << " processed=" << processed
              << " pruned=" << pruned
              << " split=" << split
              << " inspected=" << inspected
              << " avg-direct-db-dy-width=" << direct
              << " avg-sum-term-quadrant-hull-width=" << local
              << " avg-db-dy-quadrant-hull-width=" << joint
              << " local-recovery="
              << (direct==0.0 ? 0.0 : (direct-local)/direct)
              << " cross-term-recovery="
              << (direct==0.0 ? 0.0 : (local-joint)/direct)
              << " total-recovery="
              << (direct==0.0 ? 0.0 : (direct-joint)/direct)
              << std::endl;
}

struct GradientBlockForms {
    RealExpression barrier;
    RealExpression db_dx;
    RealExpression db_dy_forward;
    std::array<RealExpression,5u> db_dy_blocked;
};

GradientBlockForms gradient_block_forms(
    RealExpression const& x0,
    RealExpression const& x1)
{
    auto const& p=TestBarr3Full64::parameters();
    constexpr SizeType width=TestBarr3Full64::width;
    constexpr std::array<SizeType,5u> block_sizes={2u,4u,8u,16u,32u};
    auto c=[&](SizeType i) { return TestBarr3Full64::constant(p[i]); };

    std::array<RealExpression,width> h1;
    std::array<RealExpression,width> a1;
    std::array<RealExpression,width> dh1_dx;
    std::array<RealExpression,width> dh1_dy;
    for(SizeType j=0u;j!=width;++j) {
        RealExpression wx=c(TestBarr3Full64::w1_offset+2u*j);
        RealExpression wy=c(TestBarr3Full64::w1_offset+2u*j+1u);
        RealExpression z=c(TestBarr3Full64::b1_offset+j)+wx*x0+wy*x1;
        h1[j]=tanh(z);
        a1[j]=1-sqr(h1[j]);
        dh1_dx[j]=a1[j]*wx;
        dh1_dy[j]=a1[j]*wy;
    }

    std::array<RealExpression,width> h2;
    std::array<RealExpression,width> a2;
    std::array<RealExpression,width> dz_dx;
    std::array<RealExpression,width> dz_dy;
    for(SizeType i=0u;i!=width;++i) {
        RealExpression z=c(TestBarr3Full64::b2_offset+i);
        RealExpression local_dx=RealExpression(0);
        RealExpression local_dy=RealExpression(0);
        for(SizeType j=0u;j!=width;++j) {
            RealExpression w=c(TestBarr3Full64::w2_offset+width*i+j);
            z=z+w*h1[j];
            local_dx=local_dx+w*dh1_dx[j];
            local_dy=local_dy+w*dh1_dy[j];
        }
        h2[i]=tanh(z);
        a2[i]=1-sqr(h2[i]);
        dz_dx[i]=local_dx;
        dz_dy[i]=local_dy;
    }

    RealExpression barrier=c(TestBarr3Full64::b3_offset);
    RealExpression db_dx=RealExpression(0);
    RealExpression db_dy_forward=RealExpression(0);
    for(SizeType i=0u;i!=width;++i) {
        RealExpression w3=c(TestBarr3Full64::w3_offset+i);
        barrier=barrier+w3*h2[i];
        db_dx=db_dx+w3*a2[i]*dz_dx[i];
        db_dy_forward=db_dy_forward+w3*a2[i]*dz_dy[i];
    }

    std::array<RealExpression,5u> blocked;
    for(SizeType k=0u;k!=block_sizes.size();++k) {
        RealExpression value=RealExpression(0);
        SizeType const bs=block_sizes[k];
        for(SizeType begin=0u;begin<width;begin+=bs) {
            SizeType const end=std::min(width,begin+bs);
            RealExpression block=RealExpression(0);
            for(SizeType j=0u;j!=width;++j) {
                RealExpression inner=RealExpression(0);
                for(SizeType i=begin;i!=end;++i) {
                    inner=inner
                        + c(TestBarr3Full64::w3_offset+i)
                        * a2[i]
                        * c(TestBarr3Full64::w2_offset+width*i+j);
                }
                block=block
                    + a1[j]
                    * c(TestBarr3Full64::w1_offset+2u*j+1u)
                    * inner;
            }
            value=value+block;
        }
        blocked[k]=value;
    }
    return {barrier,db_dx,db_dy_forward,blocked};
}

Void profile_lie_gradient_block_reassociation(
    SizeType box_limit,
    UpperBoxType const& domain,
    RealSpace const& space,
    RealExpression const& baseline_expression,
    std::array<RealExpression,5u> const& candidates)
{
    constexpr std::array<SizeType,5u> block_sizes={2u,4u,8u,16u,32u};

    ValidatedProcedure baseline_procedure(
        make_function(space,baseline_expression));
    std::vector<ValidatedProcedure> procedures;
    procedures.reserve(candidates.size());
    for(auto const& expression:candidates) {
        procedures.emplace_back(make_function(space,expression));
    }

    std::array<SizeType,5u> pruned{};
    std::array<SizeType,5u> candidate_only{};
    std::array<SizeType,5u> baseline_only{};
    std::array<double,5u> width_sum{};
    std::array<double,5u> eval_seconds{};
    SizeType processed=0u;
    SizeType baseline_pruned=0u;
    SizeType split=0u;
    double baseline_width_sum=0.0;
    double baseline_seconds=0.0;

    auto root_arguments=bounds_arguments(domain);
    UpperIntervalType baseline_root=
        make_interval(evaluate(baseline_procedure,root_arguments));
    std::cout << "[lie-gradient-block-reassociation-root]"
              << " baseline=" << baseline_root
              << " baseline-width=" << baseline_root.width().raw().get_d()
              << " baseline-instructions="
              << baseline_procedure._instructions.size()
              << std::endl;
    for(SizeType k=0u;k!=block_sizes.size();++k) {
        UpperIntervalType image=
            make_interval(evaluate(procedures[k],root_arguments));
        std::cout << "[lie-gradient-block-reassociation-root-candidate]"
                  << " block=" << block_sizes[k]
                  << " image=" << image
                  << " width=" << image.width().raw().get_d()
                  << " instructions=" << procedures[k]._instructions.size()
                  << std::endl;
    }

    std::vector<UpperBoxType> pending;
    pending.push_back(domain);
    Stopwatch<Milliseconds> total;
    while(not pending.empty() && processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++processed;
        auto args=bounds_arguments(box);

        Stopwatch<Milliseconds> bsw;
        UpperIntervalType baseline=
            make_interval(evaluate(baseline_procedure,args));
        bsw.click();
        baseline_seconds+=bsw.elapsed_seconds();
        Bool const bprune=definitely(baseline.lower_bound()>=0);
        if(bprune) { ++baseline_pruned; }
        baseline_width_sum+=baseline.width().raw().get_d();

        for(SizeType k=0u;k!=block_sizes.size();++k) {
            Stopwatch<Milliseconds> csw;
            UpperIntervalType image=
                make_interval(evaluate(procedures[k],args));
            csw.click();
            eval_seconds[k]+=csw.elapsed_seconds();
            Bool const cprune=definitely(image.lower_bound()>=0);
            if(cprune) { ++pruned[k]; }
            if(not bprune && cprune) { ++candidate_only[k]; }
            if(bprune && not cprune) { ++baseline_only[k]; }
            width_sum[k]+=image.width().raw().get_d();
        }

        if(bprune) { continue; }
        auto children=box.split();
        if(definitely(children.first==children.second)) { continue; }
        ++split;
        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }
    total.click();

    auto avg=[&](double sum) {
        return processed==0u ? 0.0 : sum/static_cast<double>(processed);
    };
    std::cout << "[lie-gradient-block-reassociation-profile]"
              << " time=" << total.elapsed_seconds()
              << " processed=" << processed
              << " baseline-pruned=" << baseline_pruned
              << " split=" << split
              << " avg-baseline-width=" << avg(baseline_width_sum)
              << " baseline-time=" << baseline_seconds
              << std::endl;
    for(SizeType k=0u;k!=block_sizes.size();++k) {
        std::cout << "[lie-gradient-block-reassociation-candidate]"
                  << " block=" << block_sizes[k]
                  << " pruned=" << pruned[k]
                  << " candidate-only-pruned=" << candidate_only[k]
                  << " baseline-only-pruned=" << baseline_only[k]
                  << " avg-width=" << avg(width_sum[k])
                  << " eval-time=" << eval_seconds[k]
                  << " instructions=" << procedures[k]._instructions.size()
                  << std::endl;
    }
}

struct Layer2DirectionalNetworkAndLie {
    RealExpression barrier;
    RealExpression lie_factored;
    RealExpression lie_grouped;
};

Layer2DirectionalNetworkAndLie layer2_directional_network_and_lie(
    RealExpression const& x0,
    RealExpression const& x1)
{
    auto const& p=TestBarr3Full64::parameters();
    constexpr SizeType width=TestBarr3Full64::width;
    auto c=[&](SizeType i) {
        return TestBarr3Full64::constant(p[i]);
    };

    RealExpression polynomial=x0*(sqr(x0)/3-1);
    RealExpression dx=x1;
    RealExpression dy=polynomial-x1;

    std::array<RealExpression,width> h1;
    std::array<RealExpression,width> dh1_dx;
    std::array<RealExpression,width> dh1_dy;
    for(SizeType j=0u;j!=width;++j) {
        RealExpression wx=c(TestBarr3Full64::w1_offset+2u*j);
        RealExpression wy=c(TestBarr3Full64::w1_offset+2u*j+1u);
        RealExpression z=c(TestBarr3Full64::b1_offset+j)+wx*x0+wy*x1;
        h1[j]=tanh(z);
        RealExpression factor=1-sqr(h1[j]);
        dh1_dx[j]=factor*wx;
        dh1_dy[j]=factor*wy;
    }

    RealExpression barrier=c(TestBarr3Full64::b3_offset);
    RealExpression lie_factored=RealExpression(0);
    RealExpression lie_grouped=RealExpression(0);

    for(SizeType i=0u;i!=width;++i) {
        RealExpression z=c(TestBarr3Full64::b2_offset+i);
        RealExpression dz_dx=RealExpression(0);
        RealExpression dz_dy=RealExpression(0);
        for(SizeType j=0u;j!=width;++j) {
            RealExpression weight=c(TestBarr3Full64::w2_offset+width*i+j);
            z=z+weight*h1[j];
            dz_dx=dz_dx+weight*dh1_dx[j];
            dz_dy=dz_dy+weight*dh1_dy[j];
        }

        RealExpression h2=tanh(z);
        RealExpression factor=1-sqr(h2);
        RealExpression output_weight=c(TestBarr3Full64::w3_offset+i);
        barrier=barrier+output_weight*h2;

        // Delay the directional combination until each second-layer neuron.
        // This retains the baseline's separate first-layer derivative sums,
        // while exposing the shared second-layer activation factor before the
        // final output aggregation.
        RealExpression local_factored=
            factor*(dz_dx*dx+dz_dy*dy);
        RealExpression local_grouped=
            factor*(x1*(dz_dx-dz_dy)+dz_dy*polynomial);

        lie_factored=lie_factored+output_weight*local_factored;
        lie_grouped=lie_grouped+output_weight*local_grouped;
    }

    return {barrier,lie_factored,lie_grouped};
}

struct LieLayer2DirectionalCounts {
    SizeType processed = 0u;
    SizeType baseline_pruned = 0u;
    SizeType factored_pruned = 0u;
    SizeType grouped_pruned = 0u;
    SizeType factored_only_pruned = 0u;
    SizeType grouped_only_pruned = 0u;
    SizeType baseline_only_vs_factored = 0u;
    SizeType baseline_only_vs_grouped = 0u;
    SizeType split = 0u;
    double baseline_width_sum = 0.0;
    double factored_width_sum = 0.0;
    double grouped_width_sum = 0.0;
    double baseline_seconds = 0.0;
    double factored_seconds = 0.0;
    double grouped_seconds = 0.0;
};

Void profile_lie_layer2_directional(
    SizeType box_limit,
    UpperBoxType const& domain,
    ValidatedScalarMultivariateFunction const& baseline_function,
    ValidatedScalarMultivariateFunction const& factored_function,
    ValidatedScalarMultivariateFunction const& grouped_function)
{
    UpperIntervalType baseline_root=apply(baseline_function,domain);
    UpperIntervalType factored_root=apply(factored_function,domain);
    UpperIntervalType grouped_root=apply(grouped_function,domain);

    std::cout << "[lie-layer2-directional-root]"
              << " baseline=" << baseline_root
              << " baseline-width=" << baseline_root.width().raw().get_d()
              << " factored=" << factored_root
              << " factored-width=" << factored_root.width().raw().get_d()
              << " grouped=" << grouped_root
              << " grouped-width=" << grouped_root.width().raw().get_d()
              << std::endl;

    std::vector<UpperBoxType> pending;
    pending.push_back(domain);
    LieLayer2DirectionalCounts counts;

    Stopwatch<Milliseconds> stopwatch;
    while(not pending.empty() && counts.processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++counts.processed;

        Stopwatch<Milliseconds> baseline_stopwatch;
        UpperIntervalType baseline_image=apply(baseline_function,box);
        baseline_stopwatch.click();
        counts.baseline_seconds+=baseline_stopwatch.elapsed_seconds();

        Stopwatch<Milliseconds> factored_stopwatch;
        UpperIntervalType factored_image=apply(factored_function,box);
        factored_stopwatch.click();
        counts.factored_seconds+=factored_stopwatch.elapsed_seconds();

        Stopwatch<Milliseconds> grouped_stopwatch;
        UpperIntervalType grouped_image=apply(grouped_function,box);
        grouped_stopwatch.click();
        counts.grouped_seconds+=grouped_stopwatch.elapsed_seconds();

        Bool baseline_rejects=definitely(baseline_image.lower_bound()>=0);
        Bool factored_rejects=definitely(factored_image.lower_bound()>=0);
        Bool grouped_rejects=definitely(grouped_image.lower_bound()>=0);

        if(baseline_rejects) { ++counts.baseline_pruned; }
        if(factored_rejects) { ++counts.factored_pruned; }
        if(grouped_rejects) { ++counts.grouped_pruned; }
        if(not baseline_rejects && factored_rejects) {
            ++counts.factored_only_pruned;
        }
        if(not baseline_rejects && grouped_rejects) {
            ++counts.grouped_only_pruned;
        }
        if(baseline_rejects && not factored_rejects) {
            ++counts.baseline_only_vs_factored;
        }
        if(baseline_rejects && not grouped_rejects) {
            ++counts.baseline_only_vs_grouped;
        }

        counts.baseline_width_sum+=baseline_image.width().raw().get_d();
        counts.factored_width_sum+=factored_image.width().raw().get_d();
        counts.grouped_width_sum+=grouped_image.width().raw().get_d();

        if(baseline_rejects) {
            continue;
        }

        auto children=box.split();
        if(definitely(children.first==children.second)) {
            continue;
        }
        ++counts.split;
        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }
    stopwatch.click();

    auto average=[](double sum,SizeType count) {
        return count==0u ? 0.0 : sum/static_cast<double>(count);
    };

    std::cout << "[lie-layer2-directional-profile]"
              << " time=" << stopwatch.elapsed_seconds()
              << " processed=" << counts.processed
              << " baseline-pruned=" << counts.baseline_pruned
              << " factored-pruned=" << counts.factored_pruned
              << " grouped-pruned=" << counts.grouped_pruned
              << " factored-only-pruned=" << counts.factored_only_pruned
              << " grouped-only-pruned=" << counts.grouped_only_pruned
              << " baseline-only-vs-factored=" << counts.baseline_only_vs_factored
              << " baseline-only-vs-grouped=" << counts.baseline_only_vs_grouped
              << " split=" << counts.split
              << " avg-baseline-width="
              << average(counts.baseline_width_sum,counts.processed)
              << " avg-factored-width="
              << average(counts.factored_width_sum,counts.processed)
              << " avg-grouped-width="
              << average(counts.grouped_width_sum,counts.processed)
              << " baseline-time=" << counts.baseline_seconds
              << " factored-time=" << counts.factored_seconds
              << " grouped-time=" << counts.grouped_seconds
              << std::endl;
}

struct LieDirectionalPropagationCounts {
    SizeType processed = 0u;
    SizeType baseline_pruned = 0u;
    SizeType factored_pruned = 0u;
    SizeType grouped_pruned = 0u;
    SizeType factored_only_pruned = 0u;
    SizeType grouped_only_pruned = 0u;
    SizeType baseline_only_vs_factored = 0u;
    SizeType baseline_only_vs_grouped = 0u;
    SizeType split = 0u;
    double baseline_width_sum = 0.0;
    double factored_width_sum = 0.0;
    double grouped_width_sum = 0.0;
    double baseline_seconds = 0.0;
    double factored_seconds = 0.0;
    double grouped_seconds = 0.0;
};

struct DirectionalNetworkAndLie {
    RealExpression barrier;
    RealExpression lie_factored_seed;
    RealExpression lie_grouped_seed;
};

DirectionalNetworkAndLie directional_network_and_lie(
    RealExpression const& x0,
    RealExpression const& x1)
{
    auto const& p=TestBarr3Full64::parameters();
    constexpr SizeType width=TestBarr3Full64::width;
    auto c=[&](SizeType i) {
        return TestBarr3Full64::constant(p[i]);
    };

    RealExpression polynomial=x0*(sqr(x0)/3-1);
    RealExpression dx=x1;
    RealExpression dy=polynomial-x1;

    std::array<RealExpression,width> h1;
    std::array<RealExpression,width> dh1_factored;
    std::array<RealExpression,width> dh1_grouped;
    for(SizeType i=0u;i!=width;++i) {
        RealExpression wx=c(TestBarr3Full64::w1_offset+2u*i);
        RealExpression wy=c(TestBarr3Full64::w1_offset+2u*i+1u);
        RealExpression z=c(TestBarr3Full64::b1_offset+i)+wx*x0+wy*x1;
        h1[i]=tanh(z);
        RealExpression factor=1-sqr(h1[i]);

        // Same directional derivative, with two algebraically equivalent
        // seed forms. The first preserves the already successful factored
        // dynamics. The second groups y only while its coefficients are
        // constants, before any network-dependent factors are introduced.
        dh1_factored[i]=factor*(wx*dx+wy*dy);
        dh1_grouped[i]=factor*(x1*(wx-wy)+wy*polynomial);
    }

    std::array<RealExpression,width> h2;
    std::array<RealExpression,width> dh2_factored;
    std::array<RealExpression,width> dh2_grouped;
    for(SizeType i=0u;i!=width;++i) {
        RealExpression z=c(TestBarr3Full64::b2_offset+i);
        RealExpression dz_factored=RealExpression(0);
        RealExpression dz_grouped=RealExpression(0);
        for(SizeType j=0u;j!=width;++j) {
            RealExpression weight=c(TestBarr3Full64::w2_offset+width*i+j);
            z=z+weight*h1[j];
            dz_factored=dz_factored+weight*dh1_factored[j];
            dz_grouped=dz_grouped+weight*dh1_grouped[j];
        }
        h2[i]=tanh(z);
        RealExpression factor=1-sqr(h2[i]);
        dh2_factored[i]=factor*dz_factored;
        dh2_grouped[i]=factor*dz_grouped;
    }

    RealExpression barrier=c(TestBarr3Full64::b3_offset);
    RealExpression lie_factored=RealExpression(0);
    RealExpression lie_grouped=RealExpression(0);
    for(SizeType i=0u;i!=width;++i) {
        RealExpression weight=c(TestBarr3Full64::w3_offset+i);
        barrier=barrier+weight*h2[i];
        lie_factored=lie_factored+weight*dh2_factored[i];
        lie_grouped=lie_grouped+weight*dh2_grouped[i];
    }

    return {barrier,lie_factored,lie_grouped};
}

Void profile_lie_directional_propagation(
    SizeType box_limit,
    UpperBoxType const& domain,
    ValidatedScalarMultivariateFunction const& baseline_function,
    ValidatedScalarMultivariateFunction const& factored_function,
    ValidatedScalarMultivariateFunction const& grouped_function)
{
    UpperIntervalType baseline_root=apply(baseline_function,domain);
    UpperIntervalType factored_root=apply(factored_function,domain);
    UpperIntervalType grouped_root=apply(grouped_function,domain);

    std::cout << "[lie-directional-propagation-root]"
              << " baseline=" << baseline_root
              << " baseline-width=" << baseline_root.width().raw().get_d()
              << " factored-seed=" << factored_root
              << " factored-seed-width=" << factored_root.width().raw().get_d()
              << " grouped-seed=" << grouped_root
              << " grouped-seed-width=" << grouped_root.width().raw().get_d()
              << std::endl;

    std::vector<UpperBoxType> pending;
    pending.push_back(domain);
    LieDirectionalPropagationCounts counts;

    Stopwatch<Milliseconds> stopwatch;
    while(not pending.empty() && counts.processed<box_limit) {
        UpperBoxType box=std::move(pending.back());
        pending.pop_back();
        ++counts.processed;

        Stopwatch<Milliseconds> baseline_stopwatch;
        UpperIntervalType baseline_image=apply(baseline_function,box);
        baseline_stopwatch.click();
        counts.baseline_seconds+=baseline_stopwatch.elapsed_seconds();

        Stopwatch<Milliseconds> factored_stopwatch;
        UpperIntervalType factored_image=apply(factored_function,box);
        factored_stopwatch.click();
        counts.factored_seconds+=factored_stopwatch.elapsed_seconds();

        Stopwatch<Milliseconds> grouped_stopwatch;
        UpperIntervalType grouped_image=apply(grouped_function,box);
        grouped_stopwatch.click();
        counts.grouped_seconds+=grouped_stopwatch.elapsed_seconds();

        Bool baseline_rejects=definitely(baseline_image.lower_bound()>=0);
        Bool factored_rejects=definitely(factored_image.lower_bound()>=0);
        Bool grouped_rejects=definitely(grouped_image.lower_bound()>=0);

        if(baseline_rejects) { ++counts.baseline_pruned; }
        if(factored_rejects) { ++counts.factored_pruned; }
        if(grouped_rejects) { ++counts.grouped_pruned; }
        if(not baseline_rejects && factored_rejects) {
            ++counts.factored_only_pruned;
        }
        if(not baseline_rejects && grouped_rejects) {
            ++counts.grouped_only_pruned;
        }
        if(baseline_rejects && not factored_rejects) {
            ++counts.baseline_only_vs_factored;
        }
        if(baseline_rejects && not grouped_rejects) {
            ++counts.baseline_only_vs_grouped;
        }

        counts.baseline_width_sum+=baseline_image.width().raw().get_d();
        counts.factored_width_sum+=factored_image.width().raw().get_d();
        counts.grouped_width_sum+=grouped_image.width().raw().get_d();

        // Preserve exactly the established baseline frontier. Alternative
        // forms are measured but never allowed to alter traversal here.
        if(baseline_rejects) {
            continue;
        }

        auto children=box.split();
        if(definitely(children.first==children.second)) {
            continue;
        }
        ++counts.split;
        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }
    stopwatch.click();

    auto average=[](double sum,SizeType count) {
        return count==0u ? 0.0 : sum/static_cast<double>(count);
    };

    std::cout << "[lie-directional-propagation-profile]"
              << " time=" << stopwatch.elapsed_seconds()
              << " processed=" << counts.processed
              << " baseline-pruned=" << counts.baseline_pruned
              << " factored-seed-pruned=" << counts.factored_pruned
              << " grouped-seed-pruned=" << counts.grouped_pruned
              << " factored-seed-only-pruned=" << counts.factored_only_pruned
              << " grouped-seed-only-pruned=" << counts.grouped_only_pruned
              << " baseline-only-vs-factored=" << counts.baseline_only_vs_factored
              << " baseline-only-vs-grouped=" << counts.baseline_only_vs_grouped
              << " split=" << counts.split
              << " avg-baseline-width="
              << average(counts.baseline_width_sum,counts.processed)
              << " avg-factored-seed-width="
              << average(counts.factored_width_sum,counts.processed)
              << " avg-grouped-seed-width="
              << average(counts.grouped_width_sum,counts.processed)
              << " baseline-time=" << counts.baseline_seconds
              << " factored-seed-time=" << counts.factored_seconds
              << " grouped-seed-time=" << counts.grouped_seconds
              << std::endl;
}

TestBarr3Full64::NetworkAndLie reassociated_network_and_lie(
    RealExpression const& x0,
    RealExpression const& x1)
{
    auto const& p=TestBarr3Full64::parameters();
    constexpr SizeType width=TestBarr3Full64::width;
    auto c=[&](SizeType i) {
        return TestBarr3Full64::constant(p[i]);
    };

    std::array<RealExpression,width> h1;
    std::array<RealExpression,width> a1;
    for(SizeType j=0u;j!=width;++j) {
        RealExpression z=c(TestBarr3Full64::b1_offset+j)
            + c(TestBarr3Full64::w1_offset+2u*j)*x0
            + c(TestBarr3Full64::w1_offset+2u*j+1u)*x1;
        h1[j]=tanh(z);
        a1[j]=1-sqr(h1[j]);
    }

    std::array<RealExpression,width> h2;
    std::array<RealExpression,width> a2;
    for(SizeType i=0u;i!=width;++i) {
        RealExpression z=c(TestBarr3Full64::b2_offset+i);
        for(SizeType j=0u;j!=width;++j) {
            z=z+c(TestBarr3Full64::w2_offset+width*i+j)*h1[j];
        }
        h2[i]=tanh(z);
        a2[i]=1-sqr(h2[i]);
    }

    RealExpression barrier=c(TestBarr3Full64::b3_offset);
    for(SizeType i=0u;i!=width;++i) {
        barrier=barrier+c(TestBarr3Full64::w3_offset+i)*h2[i];
    }

    RealExpression db_dx=RealExpression(0);
    RealExpression db_dy=RealExpression(0);
    for(SizeType j=0u;j!=width;++j) {
        RealExpression accumulated=RealExpression(0);
        for(SizeType i=0u;i!=width;++i) {
            accumulated=accumulated
                + c(TestBarr3Full64::w3_offset+i)
                    * a2[i]
                    * c(TestBarr3Full64::w2_offset+width*i+j);
        }
        RealExpression local=accumulated*a1[j];
        db_dx=db_dx+local*c(TestBarr3Full64::w1_offset+2u*j);
        db_dy=db_dy+local*c(TestBarr3Full64::w1_offset+2u*j+1u);
    }

    RealExpression dx=x1;
    RealExpression dy=x0*(sqr(x0)/3-1)-x1;
    return {barrier,db_dx*dx+db_dy*dy,db_dx,db_dy};
}

Void profile_mean_value_range(
    String const& name,
    ValidatedScalarMultivariateFunction const& function,
    ExactBoxType const& domain)
{
    UpperBoxType upper_domain(domain);

    Stopwatch<Milliseconds> interval_stopwatch;
    UpperIntervalType interval_image=apply(function,upper_domain);
    interval_stopwatch.click();

    UpperBoxType midpoint_box(upper_domain);
    for(SizeType i=0u; i!=midpoint_box.dimension(); ++i) {
        auto midpoint=upper_domain[i].midpoint();
        midpoint_box[i]=UpperIntervalType(midpoint.raw(),midpoint.raw());
    }

    Stopwatch<Milliseconds> midpoint_stopwatch;
    UpperIntervalType mean_value_image=apply(function,midpoint_box);
    midpoint_stopwatch.click();

    std::vector<ValidatedScalarMultivariateFunction> derivatives;
    derivatives.reserve(upper_domain.dimension());
    Stopwatch<Milliseconds> derivative_build_stopwatch;
    for(SizeType i=0u; i!=upper_domain.dimension(); ++i) {
        derivatives.push_back(function.derivative(i));
    }
    derivative_build_stopwatch.click();

    Stopwatch<Milliseconds> derivative_eval_stopwatch;
    for(SizeType i=0u; i!=upper_domain.dimension(); ++i) {
        UpperIntervalType derivative_image=apply(derivatives[i],upper_domain);
        auto midpoint=upper_domain[i].midpoint();
        UpperIntervalType midpoint_interval(midpoint.raw(),midpoint.raw());
        mean_value_image+=derivative_image*(upper_domain[i]-midpoint_interval);
    }
    derivative_eval_stopwatch.click();

    std::cout << "[mean-value-profile] " << name
              << " interval-time=" << interval_stopwatch.elapsed_seconds()
              << " interval-image=" << interval_image
              << " midpoint-eval-time=" << midpoint_stopwatch.elapsed_seconds()
              << " derivative-build-time="
              << derivative_build_stopwatch.elapsed_seconds()
              << " derivative-eval-time="
              << derivative_eval_stopwatch.elapsed_seconds()
              << " mean-value-image=" << mean_value_image
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

    if(query=="mean-value") {
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,network.lie+network.barrier);
        profile_mean_value_range(
            "lie+barrier",lie_function,domain);
        return 0;
    }

    if(query=="lie-components") {
        RealExpression dx=ey;
        RealExpression dy=-ex-ey+(ex*ex*ex)/3;
        RealExpression term_x=network.db_dx*dx;
        RealExpression term_y=network.db_dy*dy;
        profile_interval_expression("db/dx",network.db_dx,space,UpperBoxType(domain));
        profile_interval_expression("db/dy",network.db_dy,space,UpperBoxType(domain));
        profile_interval_expression("dx=y",dx,space,UpperBoxType(domain));
        profile_interval_expression("dy=-x-y+x^3/3",dy,space,UpperBoxType(domain));
        profile_interval_expression("db/dx*dx",term_x,space,UpperBoxType(domain));
        profile_interval_expression("db/dy*dy",term_y,space,UpperBoxType(domain));
        profile_interval_expression("lie",term_x+term_y,space,UpperBoxType(domain));
        profile_interval_expression("barrier",network.barrier,space,UpperBoxType(domain));
        profile_interval_expression(
            "lie+barrier",term_x+term_y+network.barrier,
            space,UpperBoxType(domain));
        return 0;
    }

    if(query=="lie-split-profile") {
        RealExpression dy=-ex-ey+(ex*ex*ex)/3;
        RealExpression dominant_term=network.db_dy*dy;
        ValidatedScalarMultivariateFunction dominant_function=
            make_function(space,dominant_term);
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,network.lie+network.barrier);
        UpperBoxType upper_domain(domain);
        for(SizeType coordinate=0u; coordinate!=upper_domain.dimension(); ++coordinate) {
            profile_split_candidate(
                "db/dy*dy",dominant_function,upper_domain,coordinate);
            profile_split_candidate(
                "lie+barrier",lie_function,upper_domain,coordinate);
        }
        return 0;
    }

    if(query=="lie-dynamics-rewrite") {
        RealExpression dx=ey;
        RealExpression dy_original=-ex-ey+(ex*ex*ex)/3;
        RealExpression dy_factored=ex*(sqr(ex)/3-1)-ey;
        RealExpression term_x=network.db_dx*dx;
        RealExpression term_y_original=network.db_dy*dy_original;
        RealExpression term_y_factored=network.db_dy*dy_factored;
        RealExpression lie_original=term_x+term_y_original;
        RealExpression lie_factored=term_x+term_y_factored;
        UpperBoxType upper_domain(domain);

        profile_interval_expression(
            "dy-original",dy_original,space,upper_domain);
        profile_interval_expression(
            "dy-factored",dy_factored,space,upper_domain);
        profile_interval_expression(
            "db/dy*dy-original",term_y_original,space,upper_domain);
        profile_interval_expression(
            "db/dy*dy-factored",term_y_factored,space,upper_domain);
        profile_interval_expression(
            "lie+barrier-original",
            lie_original+network.barrier,space,upper_domain);
        profile_interval_expression(
            "lie+barrier-factored",
            lie_factored+network.barrier,space,upper_domain);
        return 0;
    }

    if(query=="lie-correlation-profile") {
        RealExpression dy=ex*(sqr(ex)/3-1)-ey;
        ValidatedScalarMultivariateFunction db_dy_function=
            make_function(space,network.db_dy);
        ValidatedScalarMultivariateFunction dy_function=
            make_function(space,dy);
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,network.lie+network.barrier);
        profile_lie_correlation_frontier(
            box_limit,UpperBoxType(domain),
            db_dy_function,dy_function,lie_function);
        return 0;
    }

    if(query=="lie-gradient-reassociation") {
        Stopwatch<Milliseconds> reassociation_stopwatch;
        auto reassociated=reassociated_network_and_lie(ex,ey);
        reassociation_stopwatch.click();
        std::cout << "[gradient-reassociation] construction="
                  << reassociation_stopwatch.elapsed_seconds()
                  << " s" << std::endl;
        UpperBoxType upper_domain(domain);
        profile_interval_expression(
            "db/dx-forward",network.db_dx,space,upper_domain);
        profile_interval_expression(
            "db/dx-reassociated",reassociated.db_dx,space,upper_domain);
        profile_interval_expression(
            "db/dy-forward",network.db_dy,space,upper_domain);
        profile_interval_expression(
            "db/dy-reassociated",reassociated.db_dy,space,upper_domain);
        profile_interval_expression(
            "lie+barrier-forward",
            network.lie+network.barrier,space,upper_domain);
        profile_interval_expression(
            "lie+barrier-reassociated",
            reassociated.lie+reassociated.barrier,space,upper_domain);
        return 0;
    }

    if(query=="lie-gradient-sign-profile") {
        ValidatedScalarMultivariateFunction db_dy_function=
            make_function(space,network.db_dy);
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,network.lie+network.barrier);
        profile_lie_gradient_sign_frontier(
            box_limit,UpperBoxType(domain),db_dy_function,lie_function);
        return 0;
    }

    if(query=="lie-gradient-split-profile") {
        ValidatedScalarMultivariateFunction db_dy_function=
            make_function(space,network.db_dy);
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,network.lie+network.barrier);
        profile_lie_gradient_split_frontier(
            box_limit,UpperBoxType(domain),db_dy_function,lie_function);
        return 0;
    }

    if(query=="lie-gradient-mean-value-profile") {
        ValidatedScalarMultivariateFunction db_dy_function=
            make_function(space,network.db_dy);
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,network.lie+network.barrier);
        profile_lie_gradient_mean_value_frontier(
            box_limit,UpperBoxType(domain),db_dy_function,lie_function);
        return 0;
    }

    if(query=="lie-gradient-quadrant-profile") {
        ValidatedScalarMultivariateFunction db_dy_function=
            make_function(space,network.db_dy);
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,network.lie+network.barrier);
        profile_lie_gradient_quadrant_frontier(
            box_limit,UpperBoxType(domain),db_dy_function,lie_function);
        return 0;
    }

    if(query=="lie-gradient-composite-profile") {
        ValidatedScalarMultivariateFunction db_dy_function=
            make_function(space,network.db_dy);
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,network.lie+network.barrier);
        profile_lie_gradient_composite_frontier(
            box_limit,UpperBoxType(domain),db_dy_function,lie_function);
        return 0;
    }

    if(query=="lie-gradient-symbolic-procedure-profile") {
        ValidatedScalarMultivariateFunction db_dy_function=
            make_function(space,network.db_dy);
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,network.lie+network.barrier);
        profile_lie_gradient_symbolic_procedure_frontier(
            box_limit,UpperBoxType(domain),db_dy_function,lie_function);
        return 0;
    }

    if(query=="lie-gradient-shared-procedure-profile") {
        ValidatedScalarMultivariateFunction db_dy_function=
            make_function(space,network.db_dy);
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,network.lie+network.barrier);
        profile_lie_gradient_shared_procedure_frontier(
            box_limit,UpperBoxType(domain),db_dy_function,lie_function);
        return 0;
    }

    if(query=="lie-gradient-shared-procedure-check") {
        ValidatedScalarMultivariateFunction db_dy_function=
            make_function(space,network.db_dy);
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,network.lie+network.barrier);
        profile_lie_gradient_shared_procedure_check(
            box_limit,UpperBoxType(domain),db_dy_function,lie_function);
        return 0;
    }

    if(query=="lie-gradient-expression-cse-profile") {
        profile_lie_gradient_expression_cse(
            network.db_dy,x,y,space);
        return 0;
    }

    if(query=="lie-gradient-expression-cse-frontier") {
        ValidatedScalarMultivariateFunction db_dy_function=
            make_function(space,network.db_dy);
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,network.lie+network.barrier);
        profile_lie_gradient_expression_cse_frontier(
            box_limit,UpperBoxType(domain),
            network.db_dy,x,y,space,
            db_dy_function,lie_function);
        return 0;
    }

    if(query=="lie-gradient-cse-prune-profile") {
        RealExpression dy=x*(sqr(x)/3-1)-y;
        ValidatedScalarMultivariateFunction db_dx_function=
            make_function(space,network.db_dx);
        ValidatedScalarMultivariateFunction db_dy_function=
            make_function(space,network.db_dy);
        ValidatedScalarMultivariateFunction barrier_function=
            make_function(space,network.barrier);
        ValidatedScalarMultivariateFunction dy_function=
            make_function(space,dy);
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,network.lie+network.barrier);
        profile_lie_gradient_cse_pruning(
            box_limit,UpperBoxType(domain),
            network.db_dy,x,y,space,
            db_dx_function,db_dy_function,
            barrier_function,dy_function,lie_function);
        return 0;
    }

    if(query=="lie-root-cse-range-profile") {
        RealExpression lie_barrier=network.lie+network.barrier;
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,lie_barrier);
        profile_lie_root_cse_range(
            box_limit,UpperBoxType(domain),
            lie_barrier,x,y,space,lie_function);
        return 0;
    }

    if(query=="lie-root-cse-search") {
        RealExpression lie_barrier=network.lie+network.barrier;
        ValidatedScalarMultivariateFunction lie_function=
            make_function(space,lie_barrier);
        ValidatedScalarMultivariateFunction barrier_function=
            make_function(space,network.barrier);
        profile_lie_root_cse_search(
            box_limit,UpperBoxType(domain),
            lie_barrier,x,y,space,
            lie_function,barrier_function,
            lie_literal_order,upper_child_first);
        return 0;
    }

    if(query=="lie-gradient-correlation-attribution-profile") {
        Vector<RealExpression> expressions=
            gradient_width_attribution_expressions(ex,ey);
        profile_lie_gradient_correlation_attribution(
            box_limit,UpperBoxType(domain),space,
            network.lie+network.barrier,expressions);
        return 0;
    }

    if(query=="lie-gradient-width-attribution-profile") {
        Vector<RealExpression> expressions=
            gradient_width_attribution_expressions(ex,ey);
        profile_lie_gradient_width_attribution(
            box_limit,UpperBoxType(domain),space,
            network.lie+network.barrier,expressions);
        return 0;
    }

    if(query=="lie-gradient-block-reassociation-profile") {
        auto forms=gradient_block_forms(ex,ey);
        RealExpression dy=ex*(sqr(ex)/3-1)-ey;
        RealExpression baseline=
            forms.db_dx*ey+forms.db_dy_forward*dy+forms.barrier;
        std::array<RealExpression,5u> candidates;
        for(SizeType i=0u;i!=candidates.size();++i) {
            candidates[i]=
                forms.db_dx*ey+forms.db_dy_blocked[i]*dy+forms.barrier;
        }
        profile_lie_gradient_block_reassociation(
            box_limit,UpperBoxType(domain),space,baseline,candidates);
        return 0;
    }

    if(query=="lie-layer2-directional-profile") {
        Stopwatch<Milliseconds> layer2_construction_stopwatch;
        auto layer2=layer2_directional_network_and_lie(ex,ey);
        layer2_construction_stopwatch.click();
        std::cout << "[layer2-directional]"
                  << " construction="
                  << layer2_construction_stopwatch.elapsed_seconds()
                  << " s" << std::endl;

        RealExpression baseline_lie_barrier=
            network.lie+network.barrier;
        RealExpression factored_lie_barrier=
            layer2.lie_factored+layer2.barrier;
        RealExpression grouped_lie_barrier=
            layer2.lie_grouped+layer2.barrier;

        ValidatedScalarMultivariateFunction baseline_function=
            make_function(space,baseline_lie_barrier);
        ValidatedScalarMultivariateFunction factored_function=
            make_function(space,factored_lie_barrier);
        ValidatedScalarMultivariateFunction grouped_function=
            make_function(space,grouped_lie_barrier);

        profile_lie_layer2_directional(
            box_limit,UpperBoxType(domain),
            baseline_function,factored_function,grouped_function);
        return 0;
    }

    if(query=="lie-directional-propagation-profile") {
        Stopwatch<Milliseconds> directional_construction_stopwatch;
        auto directional=directional_network_and_lie(ex,ey);
        directional_construction_stopwatch.click();
        std::cout << "[directional-propagation]"
                  << " construction="
                  << directional_construction_stopwatch.elapsed_seconds()
                  << " s" << std::endl;

        RealExpression baseline_lie_barrier=
            network.lie+network.barrier;
        RealExpression factored_lie_barrier=
            directional.lie_factored_seed+directional.barrier;
        RealExpression grouped_lie_barrier=
            directional.lie_grouped_seed+directional.barrier;

        ValidatedScalarMultivariateFunction baseline_function=
            make_function(space,baseline_lie_barrier);
        ValidatedScalarMultivariateFunction factored_function=
            make_function(space,factored_lie_barrier);
        ValidatedScalarMultivariateFunction grouped_function=
            make_function(space,grouped_lie_barrier);

        profile_lie_directional_propagation(
            box_limit,UpperBoxType(domain),
            baseline_function,factored_function,grouped_function);
        return 0;
    }

    if(query=="lie-algebraic-form-profile") {
        RealExpression polynomial=x*(sqr(x)/3-1);
        RealExpression baseline=
            network.db_dx*y
            + network.db_dy*(polynomial-y)
            + network.barrier;
        RealExpression expanded=
            network.db_dx*y
            + network.db_dy*polynomial
            - network.db_dy*y
            + network.barrier;
        RealExpression grouped_y=
            y*(network.db_dx-network.db_dy)
            + network.db_dy*polynomial
            + network.barrier;

        ValidatedScalarMultivariateFunction baseline_function=
            make_function(space,baseline);
        ValidatedScalarMultivariateFunction expanded_function=
            make_function(space,expanded);
        ValidatedScalarMultivariateFunction grouped_y_function=
            make_function(space,grouped_y);

        profile_lie_algebraic_forms(
            box_limit,UpperBoxType(domain),
            baseline_function,expanded_function,grouped_y_function);
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
