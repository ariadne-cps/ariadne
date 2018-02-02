/***************************************************************************
 *            vanderpol.cpp
 *
 *  Copyright  2017  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
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
    integrator.set_maximum_step_size(0.0625);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().maximum_enclosure_radius(1.0);
    evolver.configuration().maximum_step_size(0.0625);
    evolver.configuration().maximum_spacial_error(0.0001);
    evolver.verbosity = 1;
    std::cout <<  evolver.configuration() << std::endl;

    Real x1_0(1.2);
    Real x2_0(1.1);
    Real eps = 1/1024_q;

    Box<RealInterval> initial_set({{x1_0-eps,x1_0+eps},{x2_0-eps,x2_0+eps}});

    std::cout << "Initial set: " << initial_set << std::endl;
    Real evolution_time(3.6);

    std::cout << "Computing orbit... " << std::flush;
    auto orbit = evolver.orbit(evolver.enclosure(initial_set),evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    plot("lotka-volterra",ApproximateBoxType({{0.5,1.5}, {0.5,1.5}}), Colour(1.0,0.75,0.5), orbit);
}
