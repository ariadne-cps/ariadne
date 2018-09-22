/***************************************************************************
 *            laub-loomis-nonoise.cpp
 *
 *  Copyright  2018  Luca Geretti
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

#include <cstdarg>
#include "ariadne.hpp"

using namespace Ariadne;


int main()
{
    RealVariable x1("x1"),x2("x2"),x3("x3"),x4("x4"),x5("x5"),x6("x6"),x7("x7");

    VectorField dynamics({dot(x1)=Real(1.4)*x3-Real(0.9)*x1,
                          dot(x2)=Real(2.5)*x5-Real(1.5)*x2,
                          dot(x3)=Real(0.6)*x7-Real(0.8)*x2*x3,
                          dot(x4)=Real(2)-Real(1.3)*x3*x4,
                          dot(x5)=Real(0.7)*x1-x4*x5,
                          dot(x6)=Real(0.3)*x1-Real(3.1)*x6,
                          dot(x7)=Real(1.8)*x6-Real(1.5)*x2*x7
                         });

    MaximumError max_err=1e-3;
    TaylorSeriesIntegrator integrator(max_err);
    integrator.set_maximum_step_size(0.05);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().maximum_enclosure_radius(1.0);
    evolver.configuration().maximum_step_size(0.05);
    evolver.configuration().maximum_spacial_error(1e-6);
    evolver.verbosity = 1;
    std::cout <<  evolver.configuration() << std::endl;

    Real x1_0(1.2);
    Real x2_0(1.05);
    Real x3_0(1.5);
    Real x4_0(2.4);
    Real x5_0(1.0);
    Real x6_0(0.1);
    Real x7_0(0.45);
    Real eps = 1/10_q;

    Box<RealInterval> initial_set({{x1_0-eps,x1_0+eps},{x2_0-eps,x2_0+eps},{x3_0-eps,x3_0+eps},{x4_0-eps,x4_0+eps},{x5_0-eps,x5_0+eps},{x6_0-eps,x6_0+eps},{x7_0-eps,x7_0+eps}});

    std::cout << "Initial set: " << initial_set << std::endl;
    Real evolution_time(20.0);

    std::cout << "Computing orbit... " << std::flush;
    auto orbit = evolver.orbit(evolver.enclosure(initial_set),evolution_time,Semantics::UPPER);
    std::cout << "done." << std::endl;
}
