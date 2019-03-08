/***************************************************************************
 *            lotka-volterra.cpp
 *
 *  Copyright  2008-18 Luca Geretti
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

#include "lotka-volterra.hpp"
#include "noisy-utilities.hpp"

using namespace Ariadne;


int main()
{
    // run_noisy_system(LV());
    auto tmp = LV();

    auto dynamics = std::get<1>(tmp);

    MaximumError max_err=1e-6;
    TaylorSeriesIntegrator integrator(max_err);
    //integrator.set_maximum_step_size(0.02);

    VectorFieldEvolver evolver(dynamics,integrator);

    evolver.configuration().maximum_enclosure_radius(1.0);
    evolver.configuration().maximum_step_size(0.02);
    evolver.configuration().maximum_spacial_error(1e-6);
    evolver.verbosity = 1;
    std::cout <<  evolver.configuration() << std::endl;

    // RealVariablesBox initial={{x==1.2_dec},{y==1.1_dec}};
    Real x0(1.2_dec);
    Real y0(1.1_dec);
    Box<RealInterval> initial_set({{x0,x0},{y0,y0}});
    std::cout << "Initial set: " << initial_set << std::endl;

    auto evolution_time = std::get<4>(tmp);

    std::cout << "Computing orbit... " << std::flush;
    auto orbit = evolver.orbit(evolver.enclosure(initial_set),evolution_time,Semantics::UPPER);
    std::cout << "done." << std::endl;
}
