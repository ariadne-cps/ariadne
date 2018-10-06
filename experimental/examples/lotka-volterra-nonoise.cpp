/***************************************************************************
 *            lotka-volterra-nonoise.cpp
 *
 *  Copyright  2017  Luca Geretti
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

    RealConstant u1("u1",3.0_dec);
    RealConstant u2("u2",1.0_dec);
    RealVariable x1("x1"), x2("x2");

    VectorField dynamics({dot(x1)=u1*x1*(1-x2), dot(x2)= u2*x2*(x1-1)});

    MaximumError max_err=0.01;
    TaylorSeriesIntegrator integrator(max_err);
    std::cout << integrator << std::endl;
    TaylorPicardIntegrator integrator2(max_err);
    std::cout << integrator2 << std::endl;

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().maximum_enclosure_radius(1.0);
    evolver.configuration().maximum_step_size(1.0/50);
    evolver.configuration().maximum_spacial_error(1e-3);
    evolver.verbosity = 1;
    std::cout <<  evolver.configuration() << std::endl;

    Real x1_0(1.2);
    Real x2_0(1.1);
    Real eps = 1/100000000_q;

    Box<RealInterval> initial_set({{x1_0-eps,x1_0+eps},{x2_0-eps,x2_0+eps}});

    std::cout << "Initial set: " << initial_set << std::endl;
    Real evolution_time(10.0);

    std::cout << "Computing orbit... " << std::flush;
    auto orbit = evolver.orbit(evolver.enclosure(initial_set),evolution_time,Semantics::UPPER);
    std::cout << "done." << std::endl;

    plot("lotka-volterra",ApproximateBoxType({{0.5,1.5}, {0.5,1.5}}), Colour(1.0,0.75,0.5), orbit);
}
