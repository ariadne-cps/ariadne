/***************************************************************************
 *            harmonic-oscillator-nonoise.cpp
 *
 *  Copyright  2018  Luca Geretti
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

    RealVariable x1("x1"), x2("x2");

    VectorField dynamics({dot(x1)=x2, dot(x2)= -x1});

    MaximumError max_err=0.01;
    TaylorSeriesIntegrator integrator(max_err);
    std::cout << integrator << std::endl;
    TaylorPicardIntegrator integrator2(max_err);
    std::cout << integrator2 << std::endl;

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().maximum_enclosure_radius(1.0);
    evolver.configuration().maximum_step_size(1.0/64);
    evolver.configuration().maximum_spacial_error(1e-3);
    evolver.verbosity = 1;
    std::cout <<  evolver.configuration() << std::endl;

    Real x1_0(1.0);
    Real x2_0(0.0);
    Real eps = 1/1024_q;

    Box<RealInterval> initial_set({{x1_0-eps,x1_0+eps},{x2_0-eps,x2_0+eps}});

    std::cout << "Initial set: " << initial_set << std::endl;
    Real evolution_time(2.0*3.141592);

    std::cout << "Computing orbit... " << std::flush;
    auto orbit = evolver.orbit(evolver.enclosure(initial_set),evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    plot("harmonic-oscillator",ApproximateBoxType({{-1.5,1.5}, {-1.5,1.5}}), Colour(1.0,0.75,0.5), orbit);
}
