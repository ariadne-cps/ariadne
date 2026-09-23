/***************************************************************************
 *            vanderpol.cpp
 *
 *  Copyright  2017-20  Luca Geretti
 *
 ****************************************************************************/

#include "utility/stopwatch.hpp"
#include "function/taylor_function.hpp"
#include "function/constraint.hpp"
#include "dynamics/enclosure.hpp"
#include "ariadne_main.hpp"

void ariadne_main()
{
    CONCLOG_PRINTLN("van der Pol oscillator");

    RealConstant mu("mu",1);
    RealVariable x("x"), y("y");
    VectorField dynamics({dot(x)=y, dot(y)=mu*y*(1-sqr(x))-x});

    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),1e-12);

    Real x0=1.40_dec;
    Real y0=2.30_dec;
    Real eps_x0=0.15_dec;
    Real eps_y0=0.05_dec;
    RealExpressionBoundedConstraintSet initial_set({
        x0-eps_x0<=x<=x0+eps_x0,
        y0-eps_y0<=y<=y0+eps_y0
    });

    PreconditionedGradedTaylorSeriesIntegrator integrator(
        StepMaximumError(1e-6),sweeper,lipschitz_tolerance=0.5_x,
        minimum_spacial_order=5,minimum_temporal_order=5,
        maximum_spacial_order=5,maximum_temporal_order=5);
    integrator.set_preconditioning(TaylorSeriesPreconditioning::QR);
    integrator.set_diagnostics(true);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.02);
    evolver.configuration().set_maximum_spacial_error(1e-6);
    evolver.configuration().set_enable_reconditioning(false);

    auto orbit=evolver.orbit(initial_set,Real(0.04_dec),Semantics::UPPER);
    std::cerr << "[PreconditionedSecondStepSummary]"
              << " reach_sets=" << orbit.reach().size()
              << " intermediate_sets=" << orbit.intermediate().size()
              << std::endl;

}
