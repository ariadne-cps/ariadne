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
       || mode=="smt-singular") {
        return mode;
    }
    throw std::runtime_error(
        "Usage: benchmark_smt_interval_newton "
        "[contract|infeasible|singular|smt-contract|smt-infeasible|smt-singular]");
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

    if(mode=="smt-contract" || mode=="smt-infeasible" || mode=="smt-singular") {
        RealVariable rx("newton_x"), ry("newton_y");
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
            mode=="smt-contract"
                ? ExactBoxType({
                    ExactIntervalType(0.5_x,1.0_x),
                    ExactIntervalType(0.5_x,1.0_x)})
                : mode=="smt-infeasible"
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
