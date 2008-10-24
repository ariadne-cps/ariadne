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

#include "tuple.h"
#include "vector.h"
#include "matrix.h"
#include "zonotope.h"
#include "list_set.h"
#include "grid_set.h"
#include "hybrid_set.h"
#include "hybrid_automaton.h"
#include "hybrid_evolver.h"
#include "orbit.h"
#include "discretiser.h"
#include "graphics.h"

#include "models.h"

#include "logging.h"

#include "test.h"

using namespace Ariadne;
using namespace std;
using Models::VanDerPol;

class TestDiscretiser
{
 public:
  void test() const;
};

int main() 
{
  TestDiscretiser().test();
  return ARIADNE_TEST_FAILURES;
}

void TestDiscretiser::test() const
{
  cout << __PRETTY_FUNCTION__ << endl;

  // Set up the evolution parameters and grid
  Float time(6.0);
  Float maximum_step_size(0.0625);
  Float maximum_enclosure_radius(0.25);
    
  EvolutionParameters parameters;
  parameters.maximum_enclosure_radius=maximum_enclosure_radius;
  parameters.maximum_step_size=maximum_step_size;
  Grid grid(2);

  // Set up the evaluators
  HybridEvolver evolver(parameters);
  HybridDiscretiser< ApproximateTaylorModel > discretiser(evolver);

  // Define the initial box
  Box initial_box=make_box("[1.01,1.02]x[0.51,0.52]");
  cout << "initial_box=" << initial_box << endl;

  // Set up the vector field
  Float mu=0.865;
  Function<VanDerPol> vdp(Vector<Float>(1,&mu));
  cout << "vdp=" << vdp << endl;
  HybridAutomaton hvdp("Van der Pol");
  hvdp.new_mode(DiscreteState(1),vdp);

  //Function evaluation sanity check
  cout << "vdp.evaluate(" << initial_box << ") = " << vdp.evaluate(initial_box) << endl;
  cout << "vdp.jacobian(" << initial_box << ") = " << vdp.jacobian(initial_box) << endl;
  cout << endl;
  
  // Define a bounding box for the evolution
  Box bounding_box=make_box("[-4,4]x[-4,4]") ;
  //Box eps_bounding_box=bounding_box.neighbourhood(0.1);
 
  // Over-approximate the initial set by a grid cell
  GridCell initial_cell(over_approximation(initial_box,grid));
  HybridGridCell hybrid_initial_cell(DiscreteState(1),initial_cell);
  cout << "hybrid_initial_cell=" << hybrid_initial_cell << endl << endl;

  HybridTime htime(1,time);

  // Compute the reachable sets
  Orbit<HybridGridCell> evolve_orbit
    = discretiser.upper_evolve(hvdp,hybrid_initial_cell,htime);
  cout << "evolve_orbit=" << evolve_orbit << endl;
  
  // Print the intial, evolve and reach sets
  Graphic fig;
  fig << line_style(true);
  fig << fill_colour(blue) << initial_cell;
  fig.write("test_discretiser-vdp");



}
