/***************************************************************************
 *            test_inner_approximation.cpp
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

#include <fstream>
#include <iostream>

#include "config.hpp"
#include "utility/tuple.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/constraint.hpp"
#include "function/taylor_function.hpp"
#include "dynamics/enclosure.hpp"
#include "dynamics/orbit.hpp"
#include "dynamics/vector_field.hpp"
#include "dynamics/vector_field_evolver.hpp"
#include "geometry/box.hpp"
#include "geometry/list_set.hpp"
#include "geometry/function_set.hpp"
#include "solvers/integrator.hpp"
#include "solvers/linear_programming.hpp"
#include "symbolic/expression_set.hpp"
#include "io/figure.hpp"
#include "io/drawer.hpp"
#include "io/graphics_manager.hpp"
#include "conclog/logging.hpp"

#include "../test.hpp"

using namespace ConcLog;

using namespace Ariadne;
using namespace std;

double gamma(LabelledEnclosure const& encl, SizeType idx) {
    auto rng = encl.bounding_box().euclidean_set()[idx];
    return rng.width().get_d();
}

double gamma_min(LabelledEnclosure const& inner, LabelledEnclosure const& outer) {
    double result = std::numeric_limits<double>::infinity();
    for (SizeType i=0; i<inner.dimension(); ++i)
        result = min(result,gamma(inner,i)/gamma(outer,i));
    return result;
}

List<LabelledEnclosure> boundary(LabelledEnclosure const& enclosure) {
    List<LabelledEnclosure> result;
    auto nv = enclosure.state_space().size();
    for (SizeType i=0; i<nv; ++i) {
        {
            auto boundary_piece  = enclosure;
            auto bounding_domain = boundary_piece.domain();
            bounding_domain[i].set_lower_bound(bounding_domain[i].upper_bound());
            boundary_piece.restrict(bounding_domain);
            result.push_back(boundary_piece);
        }
        {
            auto boundary_piece  = enclosure;
            auto bounding_domain = boundary_piece.domain();
            bounding_domain[i].set_upper_bound(bounding_domain[i].lower_bound());
            boundary_piece.restrict(bounding_domain);
            result.push_back(boundary_piece);
        }
    }
    return result;
}

Tuple<Matrix<FloatDPBounds>,Vector<FloatDPBounds>,Vector<FloatDPBounds>,Vector<FloatDPBounds>> construct_problem(ValidatedVectorMultivariateFunction const& f, ExactBoxType const& d) {

    auto x0 = midpoint(d);
    auto Jx0 = f.jacobian(x0);

    auto J_rng = jacobian_range(f,cast_vector(d));
    auto b_rng = f(x0) - Jx0 * Vector<FloatDPUpperInterval>(x0) + (J_rng-Jx0) * (d-x0);

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

    Vector<FloatDPBounds> b(2*n,FloatDP(0,DoublePrecision()));
    for (SizeType i=0; i<n; ++i) {
        b.at(i) = -b_rng.at(i).lower_bound().raw();
    }
    for (SizeType i=0; i<n; ++i) {
        b.at(n+i) = b_rng.at(i).upper_bound().raw();
    }

    Vector<FloatDPBounds> xl(nv,FloatDP(0,DoublePrecision()));
    for (SizeType i=0; i<m; ++i) {
        xl.at(i) = -1;
    }

    Vector<FloatDPBounds> xu(nv,FloatDP(0,DoublePrecision()));
    for (SizeType i=0; i<m; ++i) {
        xu.at(i) = 1;
    }
    for (SizeType i=m; i<nv; ++i) {
        xu.at(i) = inf;
    }

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

    ExactBoxType q(m,ExactIntervalType(0,0,DoublePrecision()));
    for (SizeType p=0;p<m;++p) {
        auto c = Vector<FloatDPBounds>::unit(nv,p,DoublePrecision());
        auto sol_lower = solver.minimise(c,xl,xu,A,b);
        q[p].set_lower_bound(sol_lower[p].lower_raw());
        c = -c;
        auto sol_upper = solver.minimise(c,xl,xu,A,b);
        q[p].set_upper_bound(sol_upper[p].upper_raw());
    }
    return q;
}

ExactBoxType inner_difference(ExactBoxType const& bx1, ExactBoxType const& bx2) {
    ARIADNE_PRECONDITION(bx1.dimension() == bx2.dimension())

    auto result = bx1;

    auto n = bx1.dimension();

    SizeType change_idx = 0;
    Interval<FloatDPBounds> change_range(0,0);
    FloatDPBounds vmax = FloatDP(-inf,DoublePrecision());

    for (SizeType i=0; i<n; ++i) {
        auto dl = bx2[i].lower_bound()-bx1[i].lower_bound();
        auto du = bx1[i].upper_bound()-bx2[i].upper_bound();
        auto v = max(dl,du);
        if (definitely(v > vmax)) {
            change_idx = i;
            change_range = (definitely(dl >= du) ? Interval<FloatDPBounds>(bx1[i].lower_bound(),bx2[i].lower_bound()) : Interval<FloatDPBounds>(bx2[i].upper_bound(),bx1[i].upper_bound()));
            vmax = v;
        }
    }

    result[change_idx].set_lower_bound(change_range.lower_bound().lower_raw());
    result[change_idx].set_upper_bound(change_range.upper_bound().upper_raw());

    return result;
}

LabelledEnclosure inner_approximation(LabelledEnclosure const& outer) {
    auto result = outer;
    result.uniform_error_recondition();
    auto const& outer_function = result.state_function();
    auto outer_domain = result.domain();

    List<ValidatedVectorMultivariateFunction> boundaries;
    for (SizeType i=0;i<outer_function.result_size();++i) {
        auto new_domain_lower = outer_domain;
        new_domain_lower[i].set_upper_bound(new_domain_lower[i].lower_bound());
        boundaries.push_back(restriction(outer_function,new_domain_lower));
        auto new_domain_upper = outer_domain;
        new_domain_upper[i].set_lower_bound(new_domain_upper[i].upper_bound());
        boundaries.push_back(restriction(outer_function,new_domain_upper));
    }

    ExactBoxType I = outer_domain;

    for (auto const& boundary : boundaries) {
        auto f = outer_function - boundary;
        auto intersection = intersection_domain(f,I);
        I = inner_difference(I,intersection);
    }

    result.restrict(I);

    return result;
}

class TestInnerApproximation
{

  public:

    void test() const {
        ARIADNE_TEST_CALL(test_inner_approximation());
    }

    void test_inner_approximation() const {
        RealVariable x1("x1"), x2("x2");

        VectorField dynamics({dot(x1)=x2/2+5, dot(x2)= x1/200*(100-x1*(10+x2))+5});

        StepMaximumError max_err=1e-6;
        TaylorPicardIntegrator integrator(max_err);

        VectorFieldEvolver evolver(dynamics,integrator);
        evolver.configuration().set_maximum_enclosure_radius(1.0);
        evolver.configuration().set_maximum_step_size(0.01);

        RealExpressionBoundedConstraintSet initial_set({-1<=x1<=1,-1<=x2<=1});

        Real evolution_time = 1;

        LabelledFigure fig=LabelledFigure({4<=x1<=8,4<=x2<=8});

        auto evolution = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);

        auto outer_final = evolution.final()[0];

        auto inner_final = inner_approximation(outer_final);

        auto gamma = gamma_min(inner_final,outer_final);
        ARIADNE_TEST_PRINT(gamma)

        GraphicsManager::instance().set_drawer(AffineDrawer(7));
        fig << fill_colour(lightgrey) << outer_final << fill_colour(red) << line_colour(red) << line_width(3.0);

        auto outer_final_boundary = boundary(outer_final);
        for (auto const& encl : outer_final_boundary)
            fig << encl;

        fig << line_colour(black) << line_width(1.0) << fill_colour(orange) << inner_final;
        fig.write("test_inner_approximation");
    }
};

Int main()
{
    ARIADNE_TEST_CALL(TestInnerApproximation().test());
    return ARIADNE_TEST_FAILURES;
}
