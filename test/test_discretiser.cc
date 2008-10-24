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
  typedef DefaultEnclosureType ApproximateTaylorModel;
  typedef std::pair<DiscreteState,DefaultEnclosureType> DefaultHybridEnclosureType;

  cout << __PRETTY_FUNCTION__ << endl;

  // Set up the evolution parameters and grid
  Float time(6.0);
  Float maximum_step_size(0.125);
  Float maximum_enclosure_radius(0.25);
  uint depth=10;
  DiscreteState location(1);

  EvolutionParameters parameters;
  parameters.maximum_enclosure_radius=maximum_enclosure_radius;
  parameters.maximum_step_size=maximum_step_size;
  Grid grid(2);

  // Set up the evaluators
  HybridEvolver evolver(parameters);
  HybridDiscretiser< DefaultEnclosureType > discretiser(evolver);

  
  // Set up the vector field
  Float mu=0.865;
  Function<VanDerPol> vdp(Vector<Float>(1,&mu));
  cout << "vdp=" << vdp << endl;
  HybridAutomaton hvdp("Van der Pol");
  hvdp.new_mode(DiscreteState(1),vdp);

  // Define a bounding box for the evolution
  Box bounding_box=make_box("[-4,4]x[-4,4]") ;
  //Box eps_bounding_box=bounding_box.neighbourhood(0.1);
 
  // Define the initial cell
  Box box=make_box("[1.01,1.02]x[0.51,0.52]");
  cout << "box=" << box << endl;
  GridTreeSet approx_tree_set=outer_approximation(box,grid,depth);
  GridCell initial_cell=*approx_tree_set.begin();
  HybridGridCell hybrid_initial_cell(location,initial_cell);
  cout << "hybrid_initial_cell=" << hybrid_initial_cell << endl << endl;
  HybridBox hybrid_initial_set=hybrid_initial_cell.box();
  cout << "hybrid_initial_set=" << hybrid_initial_set << endl << endl;

  HybridTime htime(time,1);
  cout << "hybrid_time=" << htime << endl << endl;

  // Compute the reachable sets
  Orbit<DefaultHybridEnclosureType> evolve_orbit
    = evolver.orbit(hvdp,hybrid_initial_set,htime,UPPER_SEMANTICS);
  DefaultEnclosureType const& initial_set=evolve_orbit.initial().second;
  ListSet<DefaultEnclosureType> const& reach_set=evolve_orbit.reach()[location];
  ListSet<DefaultEnclosureType> const& intermediate_set=evolve_orbit.intermediate()[location];
  ListSet<DefaultEnclosureType> const& final_set=evolve_orbit.final()[location];

  // Compute the reachable sets
  Orbit<HybridGridCell> discrete_orbit
    = discretiser.upper_evolve(hvdp,hybrid_initial_cell,htime);
  GridCellListSet const& reach_cells=discrete_orbit.reach()[location];
  GridCellListSet const& intermediate_cells=discrete_orbit.intermediate()[location];
  GridCellListSet const& final_cells=discrete_orbit.final()[location];
  
  cout << "initial_set=" << initial_set.range() << endl << endl;
  cout << "initial_cell=" << initial_cell.box() << endl << endl;
  cout << "final_set=" << final_set << endl << endl;
  cout << "final_cells=" << final_cells << endl << endl;

  // Print the intial, evolve and reach sets
  {
    Graphic fig;
    fig.set_bounding_box(Box(2,Interval(-3,3)));
    fig << line_style(true);
    fig << line_style(true);
    fig << fill_colour(cyan) << reach_cells;
    fig << fill_colour(magenta) << intermediate_cells;
    fig << fill_colour(blue) << initial_cell;
    fig << fill_colour(blue) << final_cells;
    fig.write("test_discretiser-vdp");
  }

  // Print the intial, evolve and reach sets
  {
    Graphic fig;
    fig.set_bounding_box(Box(2,Interval(-3,3)));
    fig << line_style(true);
    fig << fill_colour(cyan) << reach_set;
    fig << fill_colour(magenta) << intermediate_set;
    fig << fill_colour(blue) << initial_set;
    fig << fill_colour(blue) << final_set;
    fig.write("test_discretiser-orbit-vdp");
  }


}
