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

inline char activity_symbol(SizeType step) {
    switch (step % 4) {
    case 0: return '\\';
    case 1: return '|';
    case 2: return '/';
    default: return '-';
    }
}

void discretize(GridTreePaving& hgts, Ariadne::Orbit<Ariadne::Enclosure>& orbit, unsigned precision)
{
  int oSize=orbit.reach().size();
  std::cerr<<"\n";
  int index=1;
  for (ListSet<Enclosure>::ConstIterator it = orbit.reach().begin(); it != orbit.reach().end(); it++,index++)
  {
      std::cerr << "\r[" << activity_symbol(static_cast<unsigned>(index)) << "] " << static_cast<int>((index*100)/oSize) << "% " << std::flush;
      it->state_auxiliary_set().adjoin_outer_approximation_to(hgts,precision);
  }
  fprintf(stderr,"\n");
}

int main()
{

    RealConstant mu("mu",1.0_dec);
    RealVariable x("x"), y("y");

    VectorField dynamics({dot(x)=y, dot(y)= mu*y*(1-sqr(x))-x});

    MaximumError max_err=1e-6;
    GradedTaylorSeriesIntegrator integrator(max_err);
    //integrator.set_maximum_step_size(0.02);

    VectorFieldEvolver evolver(dynamics,integrator);
    evolver.configuration().maximum_enclosure_radius(1.0);
    evolver.configuration().maximum_step_size(0.02);
    evolver.configuration().maximum_spacial_error(1e-6);
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

    // std::cout << "plotting..." << std::endl;
    // Box<FloatDPUpperInterval> graphics_box(2);
    // graphics_box[0] = FloatDPUpperInterval(-2.5,2.5);
    // graphics_box[1] = FloatDPUpperInterval(-3.0,3.0);
    // Figure fig=Figure();
    // fig.set_bounding_box(graphics_box);
    // fig.set_line_colour(0.0,0.0,0.0);
    // fig.set_line_style(false);
    // fig.set_fill_colour(0.5,0.5,0.5);
    // fig.set_fill_colour(1.0,0.75,0.5);
    // fig.draw(orbit.reach());
    // fig.write("vanderpol-vectorfield");

    // plot("vanderpol-vf",Axes2d(-2.5,x,2.5, -3.0,y,3.0), Colour(0.0,0.5,1.0), orbit);

    std::cout << "Discretising orbit" << std::flush;
    Grid grid(2);
    GridTreePaving gts(grid);

    for(unsigned i=2;i<=4;++i)
    {
      auto h = gts;
      clock_t s_time = clock();
      // run code
      discretize(h,orbit,i);
      // End time
      clock_t e_time = clock();
      float elapsed_time =static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
      std::cout << "instance "<<i<<" in "<<elapsed_time<<" sec" << std::endl;
      char title[32];
      sprintf(title,"%d",i);
      plot(title,ApproximateBoxType({{-2.1,2.1}, {-3.0,3.0}}), Colour(0.0,0.6,1.0), h);
    }
    std::cout << "done." << std::endl;

    // The following currently fails since auxiliary variables are not tracked
}
