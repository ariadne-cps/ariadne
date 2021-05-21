/***************************************************************************
 *            solver_demonstration.cpp
 *
 *  Copyright  2009-21  Pieter Collins
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

void print() { ARIADNE_LOG_PRINTLN(""); }
template<class T> void print(const char* label, T const& expr) { ARIADNE_LOG_PRINTLN(label << ": " << (expr)) }


void algebraic_solver_demonstration() {
    //! [Algebraic Solver demonstration]

    // Compute the solution h to the vector equation f(x,h(x))=0
    auto slv=IntervalNewtonSolver(1e-8,6);

    auto dom=BoxDomainType({{-1,+1},{-1,+1},{-1,+1}});
    auto x=EffectiveScalarMultivariateFunction::coordinate(3,0);
    auto y0=EffectiveScalarMultivariateFunction::coordinate(3,1);
    auto y1=EffectiveScalarMultivariateFunction::coordinate(3,2);
    auto f=join(x+4*y0+y1,y0+y1);
    print("f:",f);
    auto hf=slv.implicit(f,BoxDomainType({{-1,+1}}),BoxDomainType({{-1,+1},{-1,+1}}));
    print("implicit(f):",hf);

    // Compute the solution h to the scalar equation g(x,hg(x))=0
    // with f(x,y)=4+x-y^2, so y=sqrt(4+x);
    dom=BoxDomainType({{-1,+1},{-1,+1}});
    x=EffectiveScalarMultivariateFunction::coordinate(2,0);
    auto y=EffectiveScalarMultivariateFunction::coordinate(2,1);
    auto g=x-4*y+y*y;
    print("g:",g);
    auto hg=slv.implicit(g, BoxDomainType({{-1,+1}}),IntervalDomainType({-1,+1}));
    print("implicit(g):",hg);
    //! [Algebraic Solver demonstration]
}

void differential_solver_demonstration() {
    //! [Differential Solver demonstration]

    // Compute the flow of the Taylor function f starting in the domain dom for time interval [-h,+h]
    auto integrator=GradedTaylorSeriesIntegrator(1e-8);

    auto dom=BoxDomainType({{-1,+1}});
    auto bbx=BoxDomainType({{-4,+4}});
    auto h=1/two;
    auto f=ValidatedVectorMultivariateFunction::identity(1);

    auto phis=integrator.flow(f,dom,h);
    assert(phis.size()==1);
    auto phi=phis[0];
    print("phi:",phi);

    // Compute two time steps of the flow of the Taylor function f starting in domain D for the interval [h,2h]
    auto phi0=phi;
    print("phi.domain():",phi.domain());
    print("h:",h);
    auto phi0h=partial_evaluate(phi,1,FloatDPBounds(h,dp));
    auto dom1=phi0h.codomain();
    phi=integrator.flow(f,dom1,h)[0];
    print("phi:",phi);
    SweeperDP swp=GradedSweeperDP(dp,6);
    auto tr=ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate({{0,2*h}},0,swp)-h;
    tr=ValidatedScalarMultivariateFunctionModelDP(tr);
    auto phi1=compose(phi,combine(phi0h,tr));
    print("phi0:",phi0);
    print("phi1:",phi1);
    //! [Differential Solver demonstration]
}


int main(int argc, const char* argv[]) {
    if (not CommandLineInterface::instance().acquire(argc,argv)) return -1;

    algebraic_solver_demonstration();
    print();
    differential_solver_demonstration();
}



