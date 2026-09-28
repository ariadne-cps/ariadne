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
    if(argument=="all" || argument=="lie" || argument=="lie-only" || argument=="eval" || argument=="taylor" || argument=="affine" || argument=="mean-value" || argument=="lie-components" || argument=="lie-split-profile" || argument=="lie-dynamics-rewrite" || argument=="lie-correlation-profile" || argument=="lie-gradient-reassociation" || argument=="lie-gradient-sign-profile" || argument=="lie-gradient-split-profile" || argument=="lie-gradient-mean-value-profile") { return argument; }
    throw std::runtime_error(
        "Usage: benchmark_smt_barr3_verification "
        "[positive-box-limit|full] [sensitivity|geometric|lookahead] "
        "[witness|no-witness] [shaving|no-shaving] [hull|no-hull] "
        "[all|lie|lie-only|eval|taylor|affine|mean-value|lie-components|lie-split-profile|lie-dynamics-rewrite|lie-correlation-profile|lie-gradient-reassociation|lie-gradient-sign-profile|lie-gradient-split-profile|lie-gradient-mean-value-profile]");
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
        "[all|lie|lie-only|eval|taylor|affine|mean-value|lie-components|lie-split-profile|lie-dynamics-rewrite|lie-correlation-profile|lie-gradient-reassociation|lie-gradient-sign-profile|lie-gradient-split-profile|lie-gradient-mean-value-profile] [monotone|no-monotone]");
}

String lie_literal_order_from_argument(Int argc,const char* argv[]) {
    if(argc<=8) { return "barrier-first"; }
    String argument(argv[8]);
    if(argument=="barrier-first" || argument=="lie-first") { return argument; }
    throw std::runtime_error(
        "Usage: benchmark_smt_barr3_verification "
        "[positive-box-limit|full] [sensitivity|geometric|lookahead] "
        "[witness|no-witness] [shaving|no-shaving] [hull|no-hull] "
        "[all|lie|lie-only|eval|taylor|affine|mean-value|lie-components|lie-split-profile|lie-dynamics-rewrite|lie-correlation-profile|lie-gradient-reassociation|lie-gradient-sign-profile|lie-gradient-split-profile|lie-gradient-mean-value-profile] [monotone|no-monotone] "
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
        "[all|lie|lie-only|eval|taylor|affine|mean-value|lie-components|lie-split-profile|lie-dynamics-rewrite|lie-correlation-profile|lie-gradient-reassociation|lie-gradient-sign-profile|lie-gradient-split-profile|lie-gradient-mean-value-profile] [monotone|no-monotone] "
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
