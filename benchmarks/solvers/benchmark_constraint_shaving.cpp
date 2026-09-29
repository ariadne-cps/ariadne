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
#include "solvers/constraint_solver.hpp"

using namespace Ariadne;

namespace {

SizeType size_from_argument(Int argc,const char* argv[])
{
    if(argc<=1) { return 32u; }
    char* end=nullptr;
    unsigned long long parsed=std::strtoull(argv[1],&end,10);
    if(end==argv[1] || *end!='\0' || parsed==0u) {
        throw std::runtime_error(
            "Usage: benchmark_constraint_shaving [size>=1]");
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
    SizeType const n=size_from_argument(argc,argv);
    auto x=ValidatedScalarMultivariateFunction::coordinates(n);

    UpperBoxType domain(
        ExactBoxType(n,ExactIntervalType(-1.0_x,1.0_x)));

    std::vector<ConstraintPropagationConstraint> constraints;
    constraints.reserve(n);
    for(SizeType i=0u;i!=n;++i) {
        auto function=x[i];
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
    for(SizeType i=0u;i!=n;++i) {
        final_width_sum+=domain[i].width().raw().get_d();
    }

    std::cout << "[constraint-shaving-sparse]"
              << " size=" << n
              << " time=" << seconds
              << " empty=" << empty
              << " shaving-rounds=" << statistics.shaving_rounds
              << " shaving-effective-rounds=" << statistics.shaving_effective
              << " shaving-attempts=" << statistics.shaving_coordinate_attempts
              << " shaving-effective-attempts=" << statistics.shaving_coordinate_effective
              << " shaving-dependency-skipped=" << statistics.shaving_dependency_skipped
              << " shaving-evals=" << statistics.shaving_function_evaluations
              << " shaving-time=" << statistics.shaving_seconds
              << " final-width-sum=" << final_width_sum
              << std::endl;

    return empty ? 1 : 0;
}
