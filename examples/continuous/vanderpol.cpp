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

#include <cstdarg>
#include "ariadne.hpp"

using namespace Ariadne;


int main()
{

    RealConstant mu("mu",1.0_dec);
    RealVariable x("x"), y("y");

    VectorField dynamics({dot(x)=y, dot(y)= mu*y*(1-sqr(x))-x});

    MaximumError max_err=1e-6;

    GradedTaylorSeriesIntegrator integrator(max_err);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.02);
    evolver.configuration().set_maximum_spacial_error(1e-6);
    evolver.verbosity = 1;
    std::cout <<  evolver.configuration() << std::endl;

    Real x0(1.40);
    Real y0(2.40);
    Real eps_x0 = 15/100_q;
    Real eps_y0 = 5/100_q;

    Box<RealInterval> initial_set({{x0-eps_x0,x0+eps_x0},{y0-eps_y0,y0+eps_y0}});

    std::cout << "Initial set: " << initial_set << std::endl;
    Real evolution_time(7.0);

    std::cout << "Computing orbit... " << std::flush;
    auto orbit = evolver.orbit(evolver.enclosure(initial_set),evolution_time,Semantics::UPPER);
    std::cout << "done." << std::endl;

    std::cout << "plotting..." << std::endl;
    Box<FloatDPUpperInterval> graphics_box(2);
    graphics_box[0] = FloatDPUpperInterval(-2.5,2.5);
    graphics_box[1] = FloatDPUpperInterval(-3.0,3.0);
    Figure fig=Figure();
    fig.set_bounding_box(graphics_box);
    fig.set_line_colour(0.0,0.0,0.0);
    fig.set_line_style(false);
    fig.set_fill_colour(0.5,0.5,0.5);
    fig.set_fill_colour(1.0,0.75,0.5);
    fig.draw(orbit.reach());
    fig.write("vanderpol");
}
