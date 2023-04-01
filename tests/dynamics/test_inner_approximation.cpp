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

        ARIADNE_TEST_PRINT(outer_final)

        auto outer_final_boundary = boundary(outer_final);
        auto inner_final_domain = outer_final.domain();
        inner_final_domain[0].set_lower_bound(FloatDP(-0.84_x,DoublePrecision()));
        inner_final_domain[0].set_upper_bound(FloatDP(0.84_x,DoublePrecision()));
        inner_final_domain[1].set_lower_bound(FloatDP(-0.84_x,DoublePrecision()));
        inner_final_domain[1].set_upper_bound(FloatDP(0.85_x,DoublePrecision()));
        auto inner_final = outer_final;
        inner_final.restrict(inner_final_domain);

        auto gamma = gamma_min(inner_final,outer_final);

        ARIADNE_TEST_PRINT(gamma)

        GraphicsManager::instance().set_drawer(AffineDrawer(6));
        fig << fill_colour(lightgrey) << outer_final << fill_colour(red) << line_colour(red) << line_width(3.0);
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
