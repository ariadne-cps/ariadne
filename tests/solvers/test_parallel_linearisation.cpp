/***************************************************************************
 *            test_parallel_linearisation.cpp
 *
 *  Copyright  2023  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "config.hpp"

#include "solvers/solver.hpp"
#include "solvers/linear_programming.hpp"
#include "function/function.hpp"
#include "function/function_patch.hpp"
#include "function/taylor_function.hpp"
#include "algebra/vector.hpp"
#include "algebra/algebra.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/space.hpp"
#include "function/formula.hpp"
#include "io/command_line_interface.hpp"

#include "../test.hpp"

using namespace Ariadne;
using namespace std;

typedef Vector<DyadicInterval> DyadicIntervalVector;

Tuple<Matrix<FloatDPBounds>,Vector<FloatDPBounds>,Vector<FloatDPBounds>,Vector<FloatDPBounds>> construct_problem(ValidatedVectorMultivariateFunction const& f, ExactBoxType const& d) {

    auto x0 = midpoint(d);
    auto Jx0 = f.jacobian(x0);
    ARIADNE_TEST_PRINT(Jx0)

    auto J_rng = jacobian_range(f,cast_vector(d));
    auto b_rng = f(x0) - Jx0 * Vector<FloatDPUpperInterval>(x0) + (J_rng-Jx0) * (d-x0);
    ARIADNE_TEST_PRINT(b_rng)

    auto n = f.result_size();
    auto m = f.argument_size();
    auto nv = m+2*n;

    Matrix<FloatDPBounds> A(2*n,nv,FloatDP(0,DoublePrecision()));
    for (SizeType i=0; i<n; ++i) {
        for (SizeType j=0; j<m; ++j)
            A.at(i,j) = -Jx0.at(i,j);
        A.at(i,m+i) = 1;
    }
    for (SizeType i=0; i<n; ++i) {
        for (SizeType j=0; j<m; ++j)
            A.at(n+i,j) = Jx0.at(i,j);
        A.at(n+i,m+n+i) = 1;
    }
    ARIADNE_TEST_PRINT(A)

    Vector<FloatDPBounds> b(2*n,FloatDP(0,DoublePrecision()));
    for (SizeType i=0; i<n; ++i) {
        b.at(i) = -b_rng.at(i).lower_bound().raw();
    }
    for (SizeType i=0; i<n; ++i) {
        b.at(n+i) = b_rng.at(i).upper_bound().raw();
    }
    ARIADNE_TEST_PRINT(b)

    ARIADNE_TEST_EQUAL(A.row_size(),b.size())

    Vector<FloatDPBounds> xl(nv,FloatDP(0,DoublePrecision()));
    for (SizeType i=0; i<m; ++i) {
        xl.at(i) = -1;
    }
    ARIADNE_TEST_PRINT(xl)

    Vector<FloatDPBounds> xu(nv,FloatDP(0,DoublePrecision()));
    for (SizeType i=0; i<m; ++i) {
        xu.at(i) = 1;
    }
    for (SizeType i=m; i<nv; ++i) {
        xu.at(i) = inf;
    }
    ARIADNE_TEST_PRINT(xu)

    return std::make_tuple(A,b,xl,xu);
}

ExactBoxType intersection_domain(ValidatedVectorMultivariateFunction const& f, ExactBoxType const& d) {
    auto problem = construct_problem(f,d);

    auto const& A = get<0>(problem);
    auto const& b = get<1>(problem);
    auto const& xl = get<2>(problem);
    auto const& xu = get<3>(problem);

    SimplexSolver<FloatDPBounds> solver;

    auto n = f.result_size();
    auto m = f.argument_size();
    auto nv = m+2*n;

    ExactBoxType q(n,ExactIntervalType(0,0,DoublePrecision()));
    for (SizeType p=0;p<n;++p) {
        Vector<FloatDPBounds> c(nv,FloatDP(0,DoublePrecision()));
        c[p] = 1;
        auto sol_lower = solver.minimise(c,xl,xu,A,b);
        ARIADNE_TEST_PRINT(sol_lower)
        q[p].set_lower_bound(sol_lower[p].lower_raw());
        c[p] = -1;
        auto sol_upper = solver.minimise(c,xl,xu,A,b);
        ARIADNE_TEST_PRINT(sol_upper)
        q[p].set_upper_bound(sol_upper[p].upper_raw());
    }
    return q;
}

class TestParallelLinearisation
{
  public:
    Int test() {
        ARIADNE_TEST_CALL(test_parallel_linearization());
        return 0;
    }

    Void test_parallel_linearization() {
        EffectiveScalarMultivariateFunction x1=EffectiveScalarMultivariateFunction::coordinate(7,0);
        EffectiveScalarMultivariateFunction x2=EffectiveScalarMultivariateFunction::coordinate(7,1);
        EffectiveScalarMultivariateFunction x3=EffectiveScalarMultivariateFunction::coordinate(7,2);
        EffectiveScalarMultivariateFunction x4=EffectiveScalarMultivariateFunction::coordinate(7,3);
        ExactBoxType d({{-1.0_x,1.0_x},{-1.0_x,1.0_x},{-1.0_x,1.0_x},{-1.0_x,1.0_x},{-1.0_x,1.0_x},{-1.0_x,1.0_x},{-1.0_x,1.0_x}});

        EffectiveVectorMultivariateFunction outer({6.39_dec+1.06_dec*x1+0.5_dec*x2-0.02_dec*sqr(x1)-0.01_dec*x1*x2+0.05_dec*x3,
                                                   5.6_dec+0.08_dec*x1+0.92_dec*x2-0.07_dec*sqr(x1)-0.06_dec*x1*x2+0.04_dec*x4});

        EffectiveScalarMultivariateFunction x1b=EffectiveScalarMultivariateFunction::coordinate(7,4);
        EffectiveScalarMultivariateFunction x3b=EffectiveScalarMultivariateFunction::coordinate(7,5);
        EffectiveScalarMultivariateFunction x4b=EffectiveScalarMultivariateFunction::coordinate(7,6);
        EffectiveVectorMultivariateFunction outer_b_x2_1({6.89_dec+1.05_dec*x1b-0.02_dec*sqr(x1b)+0.05_dec*x3b,
                                                          6.52_dec+0.02_dec*x1b-0.07_dec*sqr(x1b)+0.04_dec*x4b});

        auto f = outer - outer_b_x2_1;

        auto q = intersection_domain(f,d);
        ARIADNE_TEST_PRINT(q)
    }

};

Int main(Int argc, const char **argv) {

    if (not CommandLineInterface::instance().acquire(argc,argv)) return -1;

    TestParallelLinearisation().test();

    return ARIADNE_TEST_FAILURES;
}
