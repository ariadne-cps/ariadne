/***************************************************************************
 *            test_differential_inclusion.cpp
 *
 *  Copyright  2008-18 Luca Geretti, Pieter Collins, Sanja Zivanovic
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

#include "../../source/dynamics/inclusion_evolver.hpp"
#include "algebra/sweeper.hpp"
#include "solvers/integrator_interface.hpp"
#include "geometry/box.hpp"
#include "function/function.hpp"
#include "function/formula.hpp"
#include "function/taylor_model.hpp"
#include "algebra/algebra.hpp"
#include "geometry/function_set.hpp"
#include "output/graphics.hpp"
#include "symbolic/expression_set.hpp"

#include "../test.hpp"

using namespace Ariadne;

class TestInclusionIntegrator {

    Void run_each_approximation(String name, InclusionVectorField const& ivf, BoxDomainType const& initial, Real evolution_time, StepSizeType step, List<InputApproximation> approximations, SweeperDP sweeper, Reconditioner const& reconditioner, unsigned int verbosity) const
    {
        for (auto appro: approximations) {
            List<InputApproximation> singleapproximation = {appro};
            std::cout << appro << std::endl;
            run_single_test(name,ivf,initial,evolution_time,step,singleapproximation,sweeper,reconditioner,verbosity);
        }
    }

    Void run_single_test(String name, InclusionVectorField const& ivf, BoxDomainType const& initial, Real evolution_time, StepSizeType step, List<InputApproximation> approximations, SweeperDP sweeper, Reconditioner const& reconditioner, unsigned int verbosity) const
    {
        auto evolver = InclusionEvolver(approximations,sweeper,step,reconditioner);
        evolver.verbosity = verbosity;
        List<ValidatedVectorMultivariateFunctionModelType> flow_functions = evolver.flow(ivf,initial,evolution_time);
    }

    Void run_test(String name, const DottedRealAssignments& dynamics, const RealVariablesBox& inputs,
                  const RealVariablesBox& initial, Real evolution_time, StepSizeType step) const {

        InclusionVectorField ivf(dynamics,inputs);

        SizeType freq=12;
        ThresholdSweeperDP sweeper(DoublePrecision(),1e-8);
        unsigned int verbosity = 0;

        LohnerReconditioner reconditioner(freq,10000);

        List<InputApproximation> approximations = {ZeroApproximation(),ConstantApproximation(),AffineApproximation(),SinusoidalApproximation(),PiecewiseApproximation()};

        this->run_each_approximation(name,ivf,initial_ranges_to_box(initial),evolution_time,step,approximations,sweeper,reconditioner,verbosity);
    }

  public:
    void test() const;

    void test_higgins_selkov() const;
    void test_jet_engine() const;
    void test_rossler() const;
};



void TestInclusionIntegrator::test_higgins_selkov() const {
    RealVariable S("S"), P("P"), v0("v0"), k1("k1"), k2("k2");
    DottedRealAssignments dynamics={dot(S)=v0-S*k1*pow(P,2),dot(P)=S*k1*pow(P,2)-k2*P};
    RealVariablesBox inputs={0.9998_dec<=v0<=1.0002_dec,0.9998_dec<=k1<=1.0002_dec,0.99981_dec<=k2<=1.00021_dec};

    Real e=1/100_q;
    RealVariablesBox initial={{2-e<=S<=2+e},{1-e<=P<=1+e}};

    Real evolution_time=1/25_q;
    StepSizeType step=0.015625_dy;

    this->run_test("higgins-selkov",dynamics,inputs,initial,evolution_time,step);
}


void TestInclusionIntegrator::test_jet_engine() const {
    RealVariable x("x"), y("y"), u1("u1"), u2("u2");
    DottedRealAssignments dynamics={dot(x)=-y-1.5_dec*pow(x,2)-0.5_dec*pow(x,3)-0.5_dec+u1,dot(y)=3*x-y+u2};
    RealVariablesBox inputs={-5/1000_q<=u1<=5/1000_q,-5/1000_q<=u2<=5/1000_q};

    Real e1=5/100_q; Real e2=7/100_q;
    RealVariablesBox initial={{1-e1<=x<=1+e1},{1-e2<=y<=1+e2}};

    Real evolution_time=1/25_q;
    StepSizeType step=0.015625_dy;

    this->run_test("jet-engine",dynamics,inputs,initial,evolution_time,step);
}

void TestInclusionIntegrator::test_rossler() const {
    RealVariable x("x"), y("y"), z("z"), u("u");
    DottedRealAssignments dynamics={dot(x)=-y-z,dot(y)=x+y*0.1_dec,dot(z)=z*(x-6)+u};
    RealVariablesBox inputs={0.099_dec<=u<=0.101_dec};

    Real e=1/1024_q;
    RealVariablesBox initial={{-9-e<=x<=-9+e},{-e<=y<=e},{0.01_dec-e<=z<=0.01_dec+e}};

    Real evolution_time=1/128_q;
    StepSizeType step=0.00390625_dy;

    this->run_test("rossler",dynamics,inputs,initial,evolution_time,step);
}


void TestInclusionIntegrator::test() const {
    ARIADNE_TEST_CALL(test_higgins_selkov());
    ARIADNE_TEST_CALL(test_jet_engine());
    ARIADNE_TEST_CALL(test_rossler());
}

int main() {

    TestInclusionIntegrator().test();
    return ARIADNE_TEST_FAILURES;
}
