/***************************************************************************
 *            vanderpol-vectorfield.cpp
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

    RealConstant mu("mu",1.0_dec);
    RealVariable x("x"), y("y");

    VectorField dynamics({dot(x)=y, dot(y)= mu*y*(1-sqr(x))-x});

    MaximumError max_err=1e-5;
    TaylorSeriesIntegrator integrator(max_err, Order(5u));
    //integrator.set_maximum_step_size(0.02);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().maximum_enclosure_radius(1.0);
    evolver.configuration().maximum_step_size(0.02);
    evolver.configuration().maximum_spacial_error(2e-4);
    evolver.verbosity = 1;
    std::cout <<  evolver.configuration() << std::endl;

    Box<RealInterval> initial_set({{1.25_dec,1.55_dec},{2.35_dec,2.45_dec}});

    std::cout << "Initial set: " << initial_set << std::endl;
    Real evolution_time(7.0);

    std::cout << "Computing orbit... " << std::flush;
    auto orbit = evolver.orbit(evolver.enclosure(initial_set),evolution_time,Semantics::UPPER);
    std::cout << "done." << std::endl;

    for (auto set : orbit.reach()) {
        if (definitely(set.bounding_box()[1].upper().raw() >= 2.75))
            std::cout << "set with upper value " << set.bounding_box()[1].upper().raw() << " is outside the safe set." << std::endl;
    }

    std::cout << "plotting..." << std::endl;
    Box<FloatDPUpperInterval> graphics_box(2);
    graphics_box[0] = FloatDPUpperInterval(-2.5,2.5);
    graphics_box[1] = FloatDPUpperInterval(-3.0,3.0);
    Figure fig=Figure();
    fig.set_bounding_box(graphics_box);
    fig.set_line_colour(0.0,0.0,0.0);
    fig.set_line_style(true);
    fig.set_fill_colour(0.5,0.5,0.5);
    fig.set_fill_colour(1.0,0.75,0.5);
    fig.draw(orbit.reach());
    fig.write("vanderpol-vectorfield");
/*
    plot("vanderpol-vf",Axes2d(-2.5,x,2.5, -3.0,y,3.0), Colour(0.0,0.5,1.0), orbit);

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
