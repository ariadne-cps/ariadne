/***************************************************************************
 * benchmark_constraint_shaving.cpp
 *
 * Sparse shaving baseline for adaptive ACID-like scheduling.
 ***************************************************************************/

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "function/function.hpp"
#include "function/procedure.hpp"
#include "geometry/box.hpp"
#include "solving/constraint_solver.hpp"

using namespace Ariadne;

namespace {

SizeType size_from_argument(Int argc,const char* argv[])
{
    if(argc<=1) { return 32u; }
    char* end=nullptr;
    unsigned long long parsed=std::strtoull(argv[1],&end,10);
    if(end==argv[1] || *end!='\0' || parsed==0u) {
        throw std::runtime_error(
            "Usage: benchmark_constraint_shaving [size>=1] [diagonal|selective]");
    }
    return static_cast<SizeType>(parsed);
}

String mode_from_argument(Int argc,const char* argv[])
{
    if(argc<=2) { return "diagonal"; }
    String mode(argv[2]);
    if(mode=="diagonal" || mode=="selective") {
        return mode;
    }
    throw std::runtime_error(
        "Usage: benchmark_constraint_shaving [size>=1] [diagonal|selective]");
}

double elapsed_seconds(std::chrono::steady_clock::time_point const& start)
{
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now()-start).count();
}

} // namespace

Int main(Int argc,const char* argv[])
{
    SizeType const n=size_from_argument(argc,argv);
    String const mode=mode_from_argument(argc,argv);
    SizeType const nuisance_count=(mode=="selective") ? 3u : 0u;
    SizeType const dimension=n+nuisance_count;
    auto x=ValidatedScalarMultivariateFunction::coordinates(dimension);

    UpperBoxType domain(
        ExactBoxType(dimension,ExactIntervalType(-1.0_x,1.0_x)));

    std::vector<ConstraintPropagationConstraint> constraints;
    constraints.reserve(n);
    for(SizeType i=0u;i!=n;++i) {
        ValidatedScalarMultivariateFunction function=x[i];
        if(mode=="selective") {
            function=function
                +(x[n]+x[n+1u]+x[n+2u])/64;
        }
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
        domain,constraints,false,statistics,true,false);
    double const seconds=elapsed_seconds(start);

    double final_width_sum=0.0;
    for(SizeType i=0u;i!=dimension;++i) {
        final_width_sum+=domain[i].width().raw().get_d();
    }

    std::cout << "[constraint-shaving-sparse]"
              << " mode=" << mode
              << " constraints=" << n
              << " dimension=" << dimension
              << " time=" << seconds
              << " empty=" << empty
              << " shaving-rounds=" << statistics.shaving_rounds
              << " shaving-effective-rounds=" << statistics.shaving_effective
              << " shaving-attempts=" << statistics.shaving_coordinate_attempts
              << " shaving-effective-attempts=" << statistics.shaving_coordinate_effective
              << " shaving-dependency-skipped=" << statistics.shaving_dependency_skipped
              << " shaving-adaptive-skipped=" << statistics.shaving_adaptive_skipped
              << " shaving-refresh-rounds=" << statistics.shaving_refresh_rounds
              << " shaving-active-rounds=" << statistics.shaving_active_rounds
              << " shaving-evals=" << statistics.shaving_function_evaluations
              << " shaving-time=" << statistics.shaving_seconds
              << " final-width-sum=" << final_width_sum
              << std::endl;

    return empty ? 1 : 0;
}
