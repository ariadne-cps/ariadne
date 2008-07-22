/***************************************************************************
 *            test_discretiser.cc
 *
 *  Copyright  2006-8  Pieter Collins
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

#include <fstream>

#include "test_float.h"

#include "ariadne.h"
#include "base/tuple.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"
#include "geometry/empty_set.h"
#include "evaluation/standard_approximator.h"
#include "system/affine_vector_field.h"
#include "system/numerical_system.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/vector_field_evolver.h"
#include "evaluation/discrete_evolver.h"
#include "evaluation/affine_integrator.h"
#include "evaluation/standard_integrator.h"
#include "evaluation/standard_subdivider.h"
#include "evaluation/cascade_reducer.h"
#include "models/vanderpol.h"
#include "output/epsstream.h"
#include "output/logging.h"

#include "test.h"

using namespace Ariadne;
using namespace std;
using Models::VanDerPolEquation;

template<class R> 
class TestDiscreteEvolver
{
 public:
  void test() const;
};

int main() 
{
  TestDiscreteEvolver<Flt>().test();
  return ARIADNE_TEST_FAILURES;
}

template<class R> 
void TestDiscreteEvolver<R>::test() const
{
  cout << __PRETTY_FUNCTION__ << endl;
  typedef Interval<R> I;

  // Set up the evolution parameters and grid
  time_type time(6.0);
  time_type step_size(0.0625);
  R grid_size(0.125);
  R enclosure_radius(0.25);
    
  EvolutionParameters<R> parameters;
  parameters.set_maximum_enclosure_radius(enclosure_radius);
  parameters.set_maximum_step_size(step_size);
  Grid<R> grid(2,grid_size);

  // Set up the evaluators
  StandardIntegrator< Zonotope<R> > integrator;
  StandardSubdivider< Zonotope<R> > subdivider;
  CascadeReducer< Zonotope<R> > reducer(3);
  Evolver< VectorField<R>, Zonotope<R> > evolver(parameters,integrator,subdivider,reducer);
  StandardApproximator< Zonotope<R> > approximator;
  DiscreteEvolver< VectorField<R>, GridApproximationScheme<R>, Zonotope<R> > discretiser(evolver,approximator);

  // Define the initial box
  Box<R> initial_box=Box<R>("[1.01,1.02]x[0.51,0.52]");
  cout << "initial_box=" << initial_box << endl;

  // Set up the vector field
  R mu=0.865;
  VanDerPolEquation<R> vdp=VanDerPolEquation<R>(Point<R>(1,&mu));
  cout << "vdp=" << vdp << endl;

  //Function evaluation sanity check
  cout << "vdp.evaluate(" << initial_box << ") = " << vdp.evaluate(initial_box) << endl;
  cout << "vdp.jacobian(" << initial_box << ") = " << vdp.jacobian(initial_box) << endl;
  cout << endl;
  
  // Define a bounding box for the evolution
  Box<R> bounding_box=Box<R>("[-4,4]x[-4,4]") ;
  Box<R> eps_bounding_box=bounding_box.neighbourhood(0.1);
 
  // Over-approximate the initial set by a grid cell
  GridCellListSet<R> initial_cells(grid);
  initial_cells.adjoin(over_approximation(initial_box,grid));
  cout << "initial_cells=" << initial_cells << endl << endl;

  // Compute the reachable sets
  GridCellListSet<R> evolve_set = discretiser.upper_evolve(vdp,initial_cells,time,grid);
  cout << "evolve_set=" << evolve_set << endl;
  GridCellListSet<R> reach_set = discretiser.upper_reach(vdp,initial_cells,time,grid);
  cout << "reach_set=" << reach_set << endl;
  
  // Print the intial, evolve and reach sets
  epsfstream eps;
  eps.open("test_discrete_evolver-vdp.eps",eps_bounding_box);
  eps << line_style(true);
  eps << fill_colour(cyan) << reach_set;
  eps << fill_colour(yellow) << evolve_set;
  eps << fill_colour(blue) << initial_cells;
  eps.close();



}
