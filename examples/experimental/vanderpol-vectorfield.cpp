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

    RealConstant mu("mu",1.0_dec);
    RealVariable x("x"), y("y");

    VectorField dynamics({dot(x)=y, dot(y)= mu*y*(1-sqr(x))-x});

    MaximumError max_err=0.01;
    TaylorSeriesIntegrator integrator(max_err);
    integrator.set_maximum_step_size(0.0625);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().maximum_enclosure_radius(1.0);
    evolver.configuration().maximum_step_size(0.0625);
    evolver.configuration().maximum_spacial_error(0.0001);
    evolver.verbosity = 0;
    std::cout <<  evolver.configuration() << std::endl;

    Real x0(2.01);
    Real eps = 1/1024_q;

    Box<RealInterval> initial_set({{x0-eps,x0+eps},{-eps,eps}});

    std::cout << "Initial set: " << initial_set << std::endl;
    Real evolution_time(6.68);

    std::cout << "Computing orbit... " << std::flush;
    auto orbit = evolver.orbit(evolver.enclosure(initial_set),evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;
/*
    plot("vanderpol-vf",ApproximateBoxType({{-2.1,2.1}, {-3.0,3.0}}), Colour(0.0,0.5,1.0), orbit);

    std::cout << "Discretising orbit" << std::flush;
    Grid grid(2);
    GridTreeSet gts(grid);

    for (ListSet<Enclosure>::ConstIterator it = orbit.reach().begin(); it != orbit.reach().end(); it++)
    {
        std::cout<<"."<<std::flush;
        it->state_auxiliary_set().adjoin_outer_approximation_to(gts,4);
    }
    std::cout << "done." << std::endl;

    // The following currently fails since auxiliary variables are not tracked
    plot("vanderpol-vf-reach",ApproximateBoxType({{-2.1,2.1}, {-3.0,3.0}}), Colour(0.0,0.5,1.0), gts);
*/
}
