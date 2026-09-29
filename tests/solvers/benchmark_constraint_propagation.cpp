/***************************************************************************
 *            benchmark_constraint_propagation.cpp
 *
 *  Copyright  2026  Luca Geretti
 ****************************************************************************/

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>

#include "function/function.hpp"
#include "function/procedure.hpp"
#include "geometry/box.hpp"
#include "solvers/constraint_solver.hpp"

using namespace Ariadne;

namespace {

SizeType chain_size_from_argument(Int argc,const char* argv[])
{
    if(argc<=1) { return 64u; }
    char* end=nullptr;
    unsigned long long parsed=std::strtoull(argv[1],&end,10);
    if(end==argv[1] || *end!='\0' || parsed<2u) {
        throw std::runtime_error(
            "Usage: benchmark_constraint_propagation [chain-size>=2]");
    }
    return static_cast<SizeType>(parsed);
}

double elapsed_seconds(std::chrono::steady_clock::time_point const& start)
{
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now()-start).count();
}

} // namespace

Int main(Int argc,const char* argv[])
{
    SizeType const n=chain_size_from_argument(argc,argv);
    auto x=ValidatedScalarMultivariateFunction::coordinates(n);

    ExactBoxType exact_domain(n,ExactIntervalType(0.0_x,1.0_x));
    UpperBoxType domain(exact_domain);

    std::vector<ConstraintPropagationConstraint> constraints;
    constraints.reserve(n);

    // Deliberately reverse the sparse chain. A sequential full-scan fixpoint
    // can propagate only one new coordinate per round:
    // x0=0, x1=x0, x2=x1, ... .
    for(SizeType i=n;i-- > 1u;) {
        auto function=x[i]-x[i-1u];
        constraints.push_back({
            function,
            ExactIntervalType(0.0_x,0.0_x),
            {},
            false,
            false,
            std::make_shared<ValidatedProcedure>(function)
        });
    }
    {
        auto function=x[0u];
        constraints.push_back({
            function,
            ExactIntervalType(0.0_x,0.0_x),
            {},
            false,
            false,
            std::make_shared<ValidatedProcedure>(function)
        });
    }

    ConstraintSolver solver;
    ConstraintPropagationStatistics statistics;
    auto const start=std::chrono::steady_clock::now();
    Bool const empty=solver.propagate(
        domain,constraints,false,statistics,false,true);
    double const seconds=elapsed_seconds(start);

    double final_width_sum=0.0;
    for(SizeType i=0u;i!=n;++i) {
        final_width_sum+=domain[i].width().raw().get_d();
    }

    std::cout << "[constraint-propagation-chain]"
              << " size=" << n
              << " time=" << seconds
              << " empty=" << empty
              << " hull-rounds=" << statistics.hull_rounds
              << " hull-effective=" << statistics.hull_effective
              << " hull-contractor-calls=" << statistics.hull_contractor_calls
              << " agenda-pushes=" << statistics.hull_agenda_pushes
              << " agenda-pops=" << statistics.hull_agenda_pops
              << " agenda-effective=" << statistics.hull_agenda_effective_calls
              << " hull-procedure-builds=" << statistics.hull_procedure_builds
              << " hull-contract-time=" << statistics.hull_contraction_seconds
              << " hull-forward-time=" << statistics.hull_forward_execution_seconds
              << " hull-backward-time=" << statistics.hull_backward_propagation_seconds
              << " final-width-sum=" << final_width_sum
              << std::endl;

    return empty ? 1 : 0;
}
