/***************************************************************************
 *            vanderpol.cpp
 *
 *  Copyright  2017-20  Luca Geretti
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

using namespace Ariadne;

int main(int argc, const char* argv[])
{
    ARIADNE_LOG_SET_VERBOSITY(get_verbosity(argc,argv));
    Logger::instance().configuration().set_theme(TT_THEME_DARK);
    Logger::instance().configuration().set_thread_name_printing_policy(ThreadNamePrintingPolicy::BEFORE);

    ConcurrencyManager::instance().set_concurrency(4);

    ARIADNE_LOG_PRINTLN("van der Pol oscillator");

    RealConstant mu("mu",1);
    RealVariable x("x"), y("y");

    VectorField dynamics({dot(x)=y, dot(y)=mu*y*(1-sqr(x))-x});

    MaximumError max_err=1e-6;
    StartingStepSizeNumRefinements num_ref=2;
    TaylorPicardIntegrator integrator(max_err,
                                      ThresholdSweeper<FloatDP>(DoublePrecision(),max_err.value()/1024),
                                      LipschitzConstant(0.5),
                                      num_ref);

    VectorFieldEvolver evolver(dynamics,integrator);
    ARIADNE_LOG_PRINTLN_VAR_AT(1,evolver.configuration());

    typedef TaskInput<VectorFieldEvolver> I;
    typedef TaskOutput<VectorFieldEvolver> O;
    auto verification_parameter = ScalarAppraisalParameter<VectorFieldEvolver>("y",TaskAppraisalParameterOptimisation::MINIMISE,[y](I const& i,O const& o,DurationType const& d) {
        return o.evolve.bounding_box()[y].upper_bound().get_d(); });
    auto verification_constraint = TaskAppraisalConstraint<VectorFieldEvolver>(verification_parameter,2.75,AppraisalConstraintSeverity::CRITICAL);
    VerificationManager::instance().add_safety_specification(evolver, {verification_constraint});

    Real x0 = 1.4_dec;
    Real y0 = 2.4_dec;
    Real eps_x0 = 0.15_dec;
    Real eps_y0 = 0.05_dec;

    EnclosureConfiguration enclosure_config(evolver.function_factory());
    enclosure_config.set_reconditioning_num_blocks(3);
    auto initial_set = evolver.enclosure({x0-eps_x0<=x<=x0+eps_x0,y0-eps_y0<=y<=y0+eps_y0},enclosure_config);
    ARIADNE_LOG_PRINTLN_VAR_AT(1,initial_set);

    Real evolution_time = 7;

    auto start = std::chrono::high_resolution_clock::now();
    ARIADNE_LOG_PRINTLN("Computing orbit... ");
    try {
        auto orbit = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
        auto end = std::chrono::high_resolution_clock::now();
        ARIADNE_LOG_PRINTLN_AT(1,"Done in " << ((double)std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count())/1000 << " seconds.");

        ARIADNE_LOG_PRINTLN_VAR_AT(1,ConcurrencyManager::instance().last_search_best_points());

        ARIADNE_LOG_PRINTLN("Plotting...");
        LabelledFigure fig({-2.5<=x<=2.5,-3<=y<=3});
        fig << fill_colour(1.0,0.75,0.5);
        fig.draw(orbit.reach());
        fig.write("vanderpol");
    } catch (CriticalAppraisalFailureException& ex) {
        ARIADNE_LOG_PRINTLN("Safety verification failure: " << ex.what());
    }
}
