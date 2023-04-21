/***************************************************************************
 *            inner.cpp
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

#include "ariadne.hpp"
#include "ariadne_main.hpp"
#include "dynamics/inner_approximation.hpp"
#include "io/drawer.hpp"
#include "utility/stopwatch.hpp"

using namespace ConcLog;

using namespace Ariadne;
using namespace std;

double gamma(LabelledEnclosure const& inner, LabelledEnclosure const& outer) {
    auto inner_domain = project(inner.domain(),Range(0,outer.domain().dimension()));
    return (inner_domain.volume()/outer.domain().volume()).get_d();
}

LabelledFigure bounded_figure(LabelledEnclosure const& e) {
    auto bx = e.bounding_box().euclidean_set();
    auto vars = List<RealVariable>(e.bounding_box().variables());
    auto x = vars[0];
    auto y = vars[1];
    auto xlb = bx[0].lower_bound().get_d();
    auto xub = bx[0].upper_bound().get_d();
    auto ylb = bx[1].lower_bound().get_d();
    auto yub = bx[1].upper_bound().get_d();

    auto xw = xub-xlb;
    auto yw = yub-ylb;
    xlb = xlb - 0.1*xw;
    xub = xub + 0.1*xw;
    ylb = ylb - 0.1*yw;
    yub = yub + 0.1*yw;

    return LabelledFigure(Axes2d(xlb<=x<=xub,ylb<=y<=yub));
}

LabelledEnclosure brusselator_sample() {
    RealVariable x("x"), y("y");
    VectorField dynamics({dot(x)=-y-1.5_dec*pow(x,2)-0.5_dec*pow(x,3)-0.5_dec,dot(y)=3*x-y});

    Real e1=5/100_q; Real e2=7/100_q;
    RealExpressionBoundedConstraintSet initial_set({1-e1<=x<=1+e1,1-e2<=y<=1+e2});

    StepMaximumError max_err=1e-6;
    TaylorPicardIntegrator integrator(max_err);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.02);

    Real evolution_time = 1.0_dec;

    LabelledFigure fig=LabelledFigure({0.8_dec<=x<=1.1_dec,0.9_dec<=y<=1.2_dec});

    return evolver.orbit(initial_set,evolution_time,Semantics::UPPER).final()[0];
}

LabelledEnclosure article_sample() {
    using VFT = ValidatedVectorMultivariateTaylorFunctionModelDP;
    using SFT = ValidatedScalarMultivariateTaylorFunctionModelDP;

    RealVariable x1("x1"), x2("x2");
    RealSpace spc({x1,x2});

    ExactBoxType domain({{-1,1},{-1,1},{-1,1},{-1,1}});

    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),1e-9);

    auto p0 = SFT::coordinate(domain,0,sweeper);
    auto p1 = SFT::coordinate(domain,1,sweeper);
    auto p2 = SFT::coordinate(domain,2,sweeper);
    auto p3 = SFT::coordinate(domain,3,sweeper);

    auto f = VFT(2,domain,sweeper);
    f[0] = 6.39_x + 1.06_x*p0 + 0.5_x*p1 - 0.02_x*p0*p0 - 0.01_x*p0*p1 + 0.05_x*p2;
    f[1] = 5.6_x + 0.08_x*p0 + 0.92_x*p1 - 0.07_x*p0*p0 - 0.06_x*p0*p1 + 0.04_x*p3;

    auto factory = TaylorFunctionFactory(sweeper);
    EnclosureConfiguration config(factory);

    return {Enclosure(domain,f,config),spc};
}

LabelledEnclosure vanderpol_sample() {
    RealConstant mu("mu",1);
    RealVariable x("x"), y("y");

    VectorField dynamics({dot(x)=y, dot(y)= mu*y*(1-sqr(x))-x});

    StepMaximumError max_err=1e-6;
    GradedTaylorSeriesIntegrator integrator(max_err);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.02);
    evolver.configuration().set_maximum_spacial_error(1e2);
    CONCLOG_PRINTLN(evolver.configuration());

    Real x0 = 1.40_dec;
    Real y0 = 2.40_dec;
    Real eps_x0 = 0.15_dec;
    Real eps_y0 = 0.05_dec;

    RealExpressionBoundedConstraintSet initial_set({x0-eps_x0<=x<=x0+eps_x0,y0-eps_y0<=y<=y0+eps_y0});

    CONCLOG_PRINTLN("Initial set: " << initial_set);
    Real evolution_time = 0.3_dec;

    return evolver.orbit(initial_set,evolution_time,Semantics::UPPER).final()[0];
}

void ariadne_main() {

    //auto linear_solver = NativeSimplex();
    //auto linear_solver = NativeIPM();
    auto linear_solver = GLPKSimplex();
    //auto linear_solver = GLPKIPM();

    auto approximator = NonlinearCandidateValidationInnerApproximator(ParallelLinearisationContractor(linear_solver,2,1));

    auto outer_final = vanderpol_sample();

    CONCLOG_PRINTLN_AT(1,"enclosure function = " << outer_final.state_function())

    auto fig = bounded_figure(outer_final);

    bool inner_found = false;
    auto inner_final = outer_final;
    try {
        Stopwatch<Milliseconds> sw;
        inner_final = approximator.compute_from(outer_final);
        sw.click();
        CONCLOG_PRINTLN("Done in " << sw.elapsed_seconds() << " seconds.");

        inner_found = true;
        auto gamma_value = gamma(inner_final, outer_final);
        CONCLOG_PRINTLN_VAR(gamma_value)

    } catch (std::exception& e) {
        CONCLOG_PRINTLN("Inner approximation could not be found")
    }

    GraphicsManager::instance().set_drawer(AffineDrawer(7));
    fig << fill_colour(lightgrey) << outer_final << fill_colour(red) << line_colour(red) << line_width(3.0);

    auto outer_final_boundary = boundary(outer_final);
    for (auto const &encl: outer_final_boundary)
        fig << encl;

    if (inner_found) fig << line_colour(black) << line_width(1.0) << fill_colour(orange) << inner_final;
    CONCLOG_RUN_AT(2,fig.write("inner"));

}