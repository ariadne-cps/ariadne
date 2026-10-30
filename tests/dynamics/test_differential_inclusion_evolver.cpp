/***************************************************************************
 *            test_differential_inclusion_evolver.cpp
 *
 *  Copyright  2008-21  Luca Geretti, Pieter Collins, Sanja Zivanovic
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

#include "dynamics/differential_inclusion_evolver.hpp"
#include "dynamics/orbit.hpp"
#include "algebra/sweeper.hpp"
#include "integration/integrator.hpp"
#include "geometry/box.hpp"
#include "function/function.hpp"
#include "function/formula.hpp"
#include "function/taylor_model.hpp"
#include "algebra/algebra.hpp"
#include "set/function_set.hpp"
#include "labelling/figure.hpp"
#include "drawing/drawing_command_line_interface.hpp"
#include "expression-set/expression_set.hpp"

#include "macro/testing.hpp"

using namespace Ariadne;

class TestDifferentialInclusionEvolver {
    // One scalar affine inclusion: cheap even under MSVC Debug, but still
    // exercises uncertain inputs and periodic parameter reconditioning.
    void run_case(const List<InputApproximation>& approximations,
                  Nat reduction_period, ExactDouble ratio,
                  Real evolution_time, ExactDouble maximum_step) const {
        RealVariable x("x"), u("u");
        DottedRealAssignments dynamics={dot(x)=u-x};
        RealVariablesBox inputs={0.5_dec<=u<=1.5_dec};
        Real e=1/32_q;
        RealVariablesBox initial={-e<=x<=e};
        DifferentialInclusion inclusion(dynamics,inputs);

        ThresholdSweeperDP sweeper(DoublePrecision(),1e-8_pr);
        TaylorPicardIntegrator integrator(
            step_maximum_error=1e-5_pr,
            sweeper,
            lipschitz_tolerance=0.5_x,
            minimum_temporal_order=4,
            maximum_temporal_order=8);
        LohnerReconditioner reconditioner(
            static_cast<Nat>(initial.variables().size()),
            static_cast<Nat>(inputs.variables().size()),
            reduction_period,ratio);

        for(auto approximation : approximations) {
            DifferentialInclusionEvolver evolver(inclusion,integrator,reconditioner);
            List<InputApproximation> selected={approximation};
            evolver.configuration().set_approximations(selected);
            evolver.configuration().set_maximum_step_size(maximum_step);
            auto orbit=evolver.orbit(initial,evolution_time);
            ARIADNE_TEST_ASSERT(orbit.reach().size()>0u);
            ARIADNE_TEST_ASSERT(orbit.final().size()>0u);
        }
    }

  public:
    void test_uncertain_input() const {
        // One step with each of the two input approximations.
        run_case({ConstantApproximation(),AffineApproximation()},
                 3,3.0_x,1/32_q,1.0_x/pow(two,5));
    }

    void test_periodic_reconditioning() const {
        // Four short steps and parameter reduction every two steps.
        run_case({AffineApproximation()},
                 2,1.25_x,1/8_q,1.0_x/pow(two,5));
    }

    void test() const {
        ARIADNE_TEST_CALL(test_uncertain_input());
        ARIADNE_TEST_CALL(test_periodic_reconditioning());
    }
};

int main(int argc, const char* argv[]) {
    if (not DrawingCommandLineInterface().acquire(argc,argv)) return -1;
    TestDifferentialInclusionEvolver().test();
    return ARIADNE_TEST_FAILURES;
}
