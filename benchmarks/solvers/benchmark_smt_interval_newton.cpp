/***************************************************************************
 *            benchmark_smt_interval_newton.cpp
 *
 *  Copyright  2026  Luca Geretti
 ****************************************************************************/

#include <chrono>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "function/function.hpp"
#include "geometry/box.hpp"
#include "solvers/solver.hpp"
#include "solvers/smt_solver.hpp"

using namespace Ariadne;

namespace {

String mode_from_argument(Int argc,const char* argv[])
{
    if(argc<=1) { return "contract"; }
    String const mode(argv[1]);
    if(mode=="contract" || mode=="infeasible" || mode=="singular"
       || mode=="smt-contract" || mode=="smt-infeasible"
       || mode=="smt-singular" || mode=="smt-mixed-contract"
       || mode=="smt-overdetermined-infeasible"
       || mode=="smt-underdetermined"
       || mode=="smt-first-subsystem-singular"
       || mode=="smt-boolean-contract"
       || mode=="smt-search-compare"
       || mode=="smt-search-compare-wide"
       || mode=="smt-search-compare-4d") {
        return mode;
    }
    throw std::runtime_error(
        "Usage: benchmark_smt_interval_newton "
        "[contract|infeasible|singular|smt-contract|smt-infeasible|smt-singular|"
        "smt-mixed-contract|smt-overdetermined-infeasible|smt-underdetermined|"
        "smt-first-subsystem-singular|smt-boolean-contract|smt-search-compare|"
        "smt-search-compare-wide|smt-search-compare-4d]");
}

double elapsed_seconds(std::chrono::steady_clock::time_point const& start)
{
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now()-start).count();
}

double width_sum(Vector<SolverInterface::ValidatedNumericType> const& box)
{
    double result=0.0;
    for(SizeType i=0u;i!=box.size();++i) {
        result+=box[i].upper().raw().get_d()-box[i].lower().raw().get_d();
    }
    return result;
}

} // namespace

Int main(Int argc,const char* argv[])
{
    String const mode=mode_from_argument(argc,argv);

    if(mode=="smt-search-compare-4d") {
        RealVariable x0("newton_4d_x0"), x1("newton_4d_x1");
        RealVariable x2("newton_4d_x2"), x3("newton_4d_x3");
        RealExpression e0=x0, e1=x1, e2=x2, e3=x3;
        RealSpace space({x0,x1,x2,x3});
        List<SmtTheoryPrimitiveLiteral> literals({
            SmtTheoryPrimitiveLiteral(
                sqr(e0)+sqr(e1)+sqr(e2)+sqr(e3)-2,
                SmtTheoryPrimitiveRelation::EQ_ZERO),
            SmtTheoryPrimitiveLiteral(
                e0-e1,SmtTheoryPrimitiveRelation::EQ_ZERO),
            SmtTheoryPrimitiveLiteral(
                e1-e2,SmtTheoryPrimitiveRelation::EQ_ZERO),
            SmtTheoryPrimitiveLiteral(
                e2-e3,SmtTheoryPrimitiveRelation::EQ_ZERO)
        });
        ExactBoxType domain({
            ExactIntervalType(0.5_x,1.0_x),
            ExactIntervalType(0.5_x,1.0_x),
            ExactIntervalType(0.5_x,1.0_x),
            ExactIntervalType(0.5_x,1.0_x)
        });

        auto run=[&](Bool interval_newton_enabled,String const& label) {
            SmtSolver solver(SmtSolverConfiguration(
                1e-5_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                65536u,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                interval_newton_enabled));

            auto const start=std::chrono::steady_clock::now();
            SmtResult result=solver.solve(space,domain,literals);
            double const seconds=elapsed_seconds(start);
            auto const& statistics=result.statistics();

            std::cout << "[smt-interval-newton-search]"
                      << " variant=" << label
                      << " dimension=4"
                      << " status=" << result.status();
            if(result.is_unknown()) {
                std::cout << " reason=" << result.unknown_reason();
            }
            std::cout << " time=" << seconds
                      << " boxes=" << statistics.boxes_processed
                      << " pruned=" << statistics.boxes_pruned
                      << " split=" << statistics.boxes_split
                      << " epsilon-certified="
                      << statistics.epsilon_box_certifications
                      << " newton-attempts="
                      << statistics.interval_newton_attempts
                      << " newton-effective="
                      << statistics.interval_newton_effective_reductions
                      << " newton-infeasible="
                      << statistics.interval_newton_infeasible
                      << " newton-singular="
                      << statistics.interval_newton_singular
                      << " newton-time="
                      << statistics.interval_newton_seconds
                      << std::endl;
        };

        run(false,"baseline");
        run(true,"newton");
        return 0;
    }

    if(mode=="smt-search-compare" || mode=="smt-search-compare-wide") {
        RealVariable rx("newton_search_x"), ry("newton_search_y");
        RealExpression ex=rx;
        RealExpression ey=ry;
        RealSpace space({rx,ry});
        List<SmtTheoryPrimitiveLiteral> literals({
            SmtTheoryPrimitiveLiteral(
                sqr(ex)+sqr(ey)-1,
                SmtTheoryPrimitiveRelation::EQ_ZERO),
            SmtTheoryPrimitiveLiteral(
                ex-ey,
                SmtTheoryPrimitiveRelation::EQ_ZERO)
        });
        ExactBoxType domain=
            mode=="smt-search-compare-wide"
                ? ExactBoxType({
                    ExactIntervalType(-2.0_x,2.0_x),
                    ExactIntervalType(-2.0_x,2.0_x)})
                : ExactBoxType({
                    ExactIntervalType(0.5_x,1.0_x),
                    ExactIntervalType(0.5_x,1.0_x)});


        auto run=[&](Bool interval_newton_enabled,String const& label) {
            SmtSolver solver(SmtSolverConfiguration(
                1e-5_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                4096u,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                interval_newton_enabled));

            auto const start=std::chrono::steady_clock::now();
            SmtResult result=solver.solve(space,domain,literals);
            double const seconds=elapsed_seconds(start);
            auto const& statistics=result.statistics();

            std::cout << "[smt-interval-newton-search]"
                      << " variant=" << label
                      << " status=" << result.status();
            if(result.is_unknown()) {
                std::cout << " reason=" << result.unknown_reason();
            }
            std::cout << " time=" << seconds
                      << " boxes=" << statistics.boxes_processed
                      << " pruned=" << statistics.boxes_pruned
                      << " split=" << statistics.boxes_split
                      << " epsilon-certified="
                      << statistics.epsilon_box_certifications
                      << " newton-attempts="
                      << statistics.interval_newton_attempts
                      << " newton-effective="
                      << statistics.interval_newton_effective_reductions
                      << " newton-infeasible="
                      << statistics.interval_newton_infeasible
                      << " newton-singular="
                      << statistics.interval_newton_singular
                      << " newton-time="
                      << statistics.interval_newton_seconds
                      << std::endl;
        };

        run(false,"baseline");
        run(true,"newton");
        return 0;
    }

    if(mode=="smt-boolean-contract") {
        RealVariable rx("newton_boolean_x"), ry("newton_boolean_y");
        RealExpression ex=rx;
        RealExpression ey=ry;
        RealSpace space({rx,ry});
        ContinuousPredicate formula=
            ((sqr(ex)+sqr(ey)-1)==0) && ((ex-ey)==0);
        ExactBoxType domain({
            ExactIntervalType(0.5_x,1.0_x),
            ExactIntervalType(0.5_x,1.0_x)
        });

        SmtSolver solver(SmtSolverConfiguration(
            1e-5_x,
            std::numeric_limits<SizeType>::max(),
            std::numeric_limits<SizeType>::max(),
            1u,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            true));

        auto const start=std::chrono::steady_clock::now();
        SmtResult result=solver.solve(space,domain,formula);
        double const seconds=elapsed_seconds(start);
        auto const& statistics=result.statistics();

        std::cout << "[smt-interval-newton-integration]"
                  << " mode=" << mode
                  << " status=" << result.status();
        if(result.is_unknown()) {
            std::cout << " reason=" << result.unknown_reason();
        }
        std::cout << " time=" << seconds
                  << " theory-checks=" << statistics.theory_checks
                  << " boxes=" << statistics.boxes_processed
                  << " pruned=" << statistics.boxes_pruned
                  << " split=" << statistics.boxes_split
                  << " newton-attempts=" << statistics.interval_newton_attempts
                  << " newton-effective="
                  << statistics.interval_newton_effective_reductions
                  << " newton-infeasible="
                  << statistics.interval_newton_infeasible
                  << " newton-singular=" << statistics.interval_newton_singular
                  << " newton-time=" << statistics.interval_newton_seconds
                  << std::endl;
        return 0;
    }

    if(mode=="smt-contract" || mode=="smt-infeasible" || mode=="smt-singular"
       || mode=="smt-mixed-contract"
       || mode=="smt-overdetermined-infeasible"
       || mode=="smt-underdetermined"
       || mode=="smt-first-subsystem-singular") {
        RealVariable rx("newton_x"), ry("newton_y");
        RealExpression ex=rx;
        RealExpression ey=ry;
        RealSpace space({rx,ry});

        SmtTheoryPrimitiveLiteral circle(
            sqr(ex)+sqr(ey)-1,
            SmtTheoryPrimitiveRelation::EQ_ZERO);
        SmtTheoryPrimitiveLiteral diagonal(
            ex-ey,
            SmtTheoryPrimitiveRelation::EQ_ZERO);
        List<SmtTheoryPrimitiveLiteral> literals;
        if(mode=="smt-mixed-contract") {
            // Put a non-equality first to verify that subsystem extraction
            // scans for EQ_ZERO literals rather than taking the first n.
            literals.append(SmtTheoryPrimitiveLiteral(
                ex,SmtTheoryPrimitiveRelation::GEQ_ZERO));
            literals.append(circle);
            literals.append(diagonal);
        } else if(mode=="smt-overdetermined-infeasible") {
            literals.append(circle);
            literals.append(diagonal);
            literals.append(SmtTheoryPrimitiveLiteral(
                2*(ex-ey),SmtTheoryPrimitiveRelation::EQ_ZERO));
        } else if(mode=="smt-underdetermined") {
            literals.append(circle);
            literals.append(SmtTheoryPrimitiveLiteral(
                ex,SmtTheoryPrimitiveRelation::GEQ_ZERO));
        } else if(mode=="smt-first-subsystem-singular") {
            // The first two equalities are linearly dependent, while replacing
            // the second by the circle equation yields a regular subsystem on
            // [0.5,1]^2. This exposes the limitation of first-n selection.
            literals.append(diagonal);
            literals.append(SmtTheoryPrimitiveLiteral(
                2*(ex-ey),SmtTheoryPrimitiveRelation::EQ_ZERO));
            literals.append(circle);
        } else {
            literals.append(circle);
            literals.append(diagonal);
        }

        ExactBoxType domain=
            (mode=="smt-contract" || mode=="smt-mixed-contract"
             || mode=="smt-underdetermined"
             || mode=="smt-first-subsystem-singular")
                ? ExactBoxType({
                    ExactIntervalType(0.5_x,1.0_x),
                    ExactIntervalType(0.5_x,1.0_x)})
                : (mode=="smt-infeasible"
                   || mode=="smt-overdetermined-infeasible")
                    ? ExactBoxType({
                        ExactIntervalType(0.8_x,1.0_x),
                        ExactIntervalType(0.8_x,1.0_x)})
                    : ExactBoxType({
                        ExactIntervalType(-1.0_x,1.0_x),
                        ExactIntervalType(-1.0_x,1.0_x)});

        SmtSolver solver(SmtSolverConfiguration(
            1e-5_x,
            std::numeric_limits<SizeType>::max(),
            std::numeric_limits<SizeType>::max(),
            1u,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            true));

        auto const start=std::chrono::steady_clock::now();
        SmtResult result=solver.solve(space,domain,literals);
        double const seconds=elapsed_seconds(start);
        auto const& statistics=result.statistics();

        std::cout << "[smt-interval-newton-integration]"
                  << " mode=" << mode
                  << " status=" << result.status();
        if(result.is_unknown()) {
            std::cout << " reason=" << result.unknown_reason();
        }
        std::cout << " time=" << seconds
                  << " boxes=" << statistics.boxes_processed
                  << " pruned=" << statistics.boxes_pruned
                  << " split=" << statistics.boxes_split
                  << " newton-attempts=" << statistics.interval_newton_attempts
                  << " newton-effective="
                  << statistics.interval_newton_effective_reductions
                  << " newton-infeasible="
                  << statistics.interval_newton_infeasible
                  << " newton-singular=" << statistics.interval_newton_singular
                  << " newton-time=" << statistics.interval_newton_seconds
                  << std::endl;
        return 0;
    }

    auto x=ValidatedScalarMultivariateFunction::coordinates(2u);
    ValidatedVectorMultivariateFunction function(
        List<ValidatedScalarMultivariateFunction>({
            sqr(x[0u])+sqr(x[1u])-1,
            x[0u]-x[1u]
        }));

    ExactBoxType exact_domain=
        mode=="contract"
            ? ExactBoxType({
                ExactIntervalType(0.5_x,1.0_x),
                ExactIntervalType(0.5_x,1.0_x)})
            : mode=="infeasible"
                ? ExactBoxType({
                    ExactIntervalType(0.8_x,1.0_x),
                    ExactIntervalType(0.8_x,1.0_x)})
                : ExactBoxType({
                    ExactIntervalType(-1.0_x,1.0_x),
                    ExactIntervalType(-1.0_x,1.0_x)});

    Vector<SolverInterface::ValidatedNumericType> current=
        cast_singleton(exact_domain);
    IntervalNewtonSolver solver(1e-12,1u);

    auto const start=std::chrono::steady_clock::now();
    try {
        Vector<SolverInterface::ValidatedNumericType> image=
            solver.step(function,current);
        double const seconds=elapsed_seconds(start);
        Bool const disjoint=not consistent(image,current);

        Vector<SolverInterface::ValidatedNumericType> contracted=current;
        if(not disjoint) {
            contracted=refinement(image,current);
        }

        std::cout << "[smt-interval-newton]"
                  << " mode=" << mode
                  << " time=" << seconds
                  << " disjoint=" << disjoint
                  << " input-width-sum=" << width_sum(current)
                  << " newton-width-sum=" << width_sum(image)
                  << " contracted-width-sum="
                  << (disjoint ? 0.0 : width_sum(contracted))
                  << " input=" << current
                  << " newton=" << image;
        if(not disjoint) {
            std::cout << " contracted=" << contracted;
        }
        std::cout << std::endl;
        return 0;
    }
    catch(const SingularMatrixException& exception) {
        double const seconds=elapsed_seconds(start);
        std::cout << "[smt-interval-newton]"
                  << " mode=" << mode
                  << " time=" << seconds
                  << " singular=true"
                  << " input-width-sum=" << width_sum(current)
                  << " message=" << exception.what()
                  << std::endl;
        return mode=="singular" ? 0 : 1;
    }
}
